#include <unordered_map>
//#include <iostream>
#include <typeinfo>
#include <random>
#include <algorithm>
#include <cstring>
#include "helper.hpp"
#include "kernel.hpp"
#define TOT_THREADS 1
#define KMER_SZ 4



//TODO: DO it such that contigs with now left or righ reads are offloaded to kernels, then try to make separate left and right kernels so that contigs only right reads are launched in right kernel
// and contigs with only left are launched in left kernels.
int main (int argc, char* argv[]){
    std::string in_file = argv[1];
    std::vector<CtgWithReads> data_in;
    int32_t max_ctg_size, total_r_reads, total_l_reads, max_read_size, max_r_count, max_l_count;
    read_locassm_data(&data_in, in_file, max_ctg_size, total_r_reads, total_l_reads, max_read_size,max_r_count, max_l_count);
    int32_t vec_size = data_in.size();
    if(DEBUG_PRINT_CPU){
        print_vals("max_l_count:",max_l_count,"max_r_count:", max_r_count);
        print_loc_data(&data_in);
        }
    int32_t max_read_count = max_r_count>max_l_count ? max_r_count : max_l_count;

    //host allocations for converting loc_assm_data to prim types
    int32_t *cid_h = new int32_t[vec_size];
    char *ctg_seqs_h = new char[max_ctg_size * vec_size];
    int32_t *ctg_seq_offsets_h = new int32_t[vec_size];
    double *depth_h = new double[vec_size];
    char *reads_left_h = new char[total_l_reads * max_read_size];
    char *reads_right_h = new char[total_r_reads * max_read_size];
    char *quals_right_h = new char[total_r_reads * max_read_size];
    char *quals_left_h = new char[total_l_reads * max_read_size];
    int32_t *reads_l_offset_h = new int32_t[total_l_reads];
    int32_t *reads_r_offset_h = new int32_t[total_r_reads];
    int32_t *rds_l_cnt_offset_h = new int32_t[vec_size];
    int32_t *rds_r_cnt_offset_h = new int32_t[vec_size];
    int32_t *quals_l_offset_h = new int32_t[total_l_reads];
    int32_t *quals_r_offset_h = new int32_t[total_r_reads];
    int64_t *term_counts_h = new int64_t[3];

    int max_mer_len = 22;

    //compute total device memory required:
    size_t total_dev_mem = sizeof(int32_t) * vec_size * 4 + sizeof(int32_t) * total_l_reads
                           + sizeof(int32_t) * total_r_reads + sizeof(char) * max_ctg_size * vec_size
                           + sizeof(char) * total_l_reads * max_read_size + sizeof(char) * total_r_reads * max_read_size
                           + sizeof(double) * vec_size + sizeof(char) * total_r_reads * max_read_size 
                           + sizeof(char) * total_l_reads * max_read_size + sizeof(int64_t)*3
                           + sizeof(loc_ht)*(max_read_size*max_read_count)*vec_size
                           + sizeof(char)*vec_size * MAX_WALK_LEN
                           + (max_mer_len + MAX_WALK_LEN) * sizeof(char) * vec_size;
    print_vals("Total GPU Mem. requested:", total_dev_mem);

    //device allocations for loc_assm_data
    int32_t *cid_d, *ctg_seq_offsets_d, *reads_l_offset_d, *reads_r_offset_d; 
    int32_t *rds_l_cnt_offset_d, *rds_r_cnt_offset_d;
    char *ctg_seqs_d, *reads_left_d, *reads_right_d, *quals_left_d, *quals_right_d;
    char *longest_walks_d, *mer_walk_temp_d;
    double *depth_d;
    int64_t *term_counts_d;
    loc_ht *d_ht;

    CUDA_CHECK(cudaMalloc(&cid_d, sizeof(int32_t) * vec_size));
    CUDA_CHECK(cudaMalloc(&ctg_seq_offsets_d, sizeof(int32_t) * vec_size));
    CUDA_CHECK(cudaMalloc(&reads_l_offset_d, sizeof(int32_t) * total_l_reads));
    CUDA_CHECK(cudaMalloc(&reads_r_offset_d, sizeof(int32_t) * total_r_reads));
    CUDA_CHECK(cudaMalloc(&rds_l_cnt_offset_d, sizeof(int32_t) * vec_size));
    CUDA_CHECK(cudaMalloc(&rds_r_cnt_offset_d, sizeof(int32_t) * vec_size));
    CUDA_CHECK(cudaMalloc(&ctg_seqs_d, sizeof(char) * max_ctg_size * vec_size));
    CUDA_CHECK(cudaMalloc(&reads_left_d, sizeof(char) * total_l_reads * max_read_size));
    CUDA_CHECK(cudaMalloc(&reads_right_d, sizeof(char) * total_r_reads * max_read_size));
    CUDA_CHECK(cudaMalloc(&depth_d, sizeof(double) * vec_size));
    CUDA_CHECK(cudaMalloc(&quals_right_d, sizeof(char) * total_r_reads * max_read_size));
    CUDA_CHECK(cudaMalloc(&quals_left_d, sizeof(char) * total_l_reads * max_read_size));
    CUDA_CHECK(cudaMalloc(&term_counts_d, sizeof(int64_t)*3));
    // if we separate out kernels for right and left walks then we can use r_count/l_count separately but for now use the max of two
    // also subtract the appropriate kmer length from max_read_size to reduce memory footprint of global ht_loc.
    // one local hashtable for each thread, so total hash_tables equal to vec_size i.e. total contigs
    // TODO: to account for overfilling of the hashtable, consider assuming load factor of 0.8 and add a cushion of memory in hashtable
    CUDA_CHECK(cudaMalloc(&d_ht, sizeof(loc_ht)*(max_read_size*max_read_count)*vec_size)); 
    //TODO: come back to this and see if we can find a better approximation of longest walk size
    CUDA_CHECK(cudaMalloc(&longest_walks_d, sizeof(char)*vec_size * MAX_WALK_LEN));
    CUDA_CHECK(cudaMalloc(&mer_walk_temp_d, (max_mer_len + MAX_WALK_LEN) * sizeof(char) * vec_size));


    

    //convert the loc_assem data to primitive structures for device
    int32_t ctgs_offset_sum = 0;
    int32_t reads_r_offset_sum = 0;
    int32_t reads_l_offset_sum = 0;
    int read_l_index = 0, read_r_index = 0;
    for(int i = 0; i < data_in.size(); i++){
        CtgWithReads temp_data = data_in[i];
        cid_h[i] = temp_data.cid;
        depth_h[i] = temp_data.depth;
        //convert string to c-string
        char *ctgs_ptr = ctg_seqs_h + ctgs_offset_sum;
        memcpy(ctgs_ptr, temp_data.seq.c_str(), temp_data.seq.size());
        ctgs_offset_sum += temp_data.seq.size();
        ctg_seq_offsets_h[i] = ctgs_offset_sum;

        for(int j = 0; j < temp_data.reads_left.size(); j++){
            char *reads_l_ptr = reads_left_h + reads_l_offset_sum;
            char *quals_l_ptr = quals_left_h + reads_l_offset_sum;
            memcpy(reads_l_ptr, temp_data.reads_left[j].seq.c_str(), temp_data.reads_left[j].seq.size());
            //quals offsets will be same as reads offset because quals and reads have same length
            memcpy(quals_l_ptr, temp_data.reads_left[j].quals.c_str(), temp_data.reads_left[j].quals.size());
            reads_l_offset_sum += temp_data.reads_left[j].seq.size();
            reads_l_offset_h[read_l_index] = reads_l_offset_sum;
            read_l_index++;
        }
        rds_l_cnt_offset_h[i] = read_l_index; // running sum of left reads count

        for(int j = 0; j < temp_data.reads_right.size(); j++){
            char *reads_r_ptr = reads_right_h + reads_r_offset_sum;
            char *quals_r_ptr = quals_right_h + reads_r_offset_sum;
            memcpy(reads_r_ptr, temp_data.reads_right[j].seq.c_str(), temp_data.reads_right[j].seq.size());
            //quals offsets will be same as reads offset because quals and reads have same length
            memcpy(quals_r_ptr, temp_data.reads_right[j].quals.c_str(), temp_data.reads_right[j].quals.size());
            reads_r_offset_sum += temp_data.reads_right[j].seq.size();
            reads_r_offset_h[read_r_index] = reads_r_offset_sum;
            read_r_index++;
        }
        rds_r_cnt_offset_h[i] = read_r_index; // running sum of right reads count
    }// data conversion for loop ends
    if(DEBUG_PRINT_CPU)
        print_vals("Completed reading in file and converting data to primitive types...");
    for(int i = 0; i < 3; i++){
        term_counts_h[i] = 0;
    }

    //move the data to device
    CUDA_CHECK(cudaMemcpy(cid_d, cid_h, sizeof(int32_t) * vec_size, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(ctg_seq_offsets_d, ctg_seq_offsets_h, sizeof(int32_t) * vec_size, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(reads_l_offset_d, reads_l_offset_h, sizeof(int32_t) * total_l_reads, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(reads_r_offset_d, reads_r_offset_h, sizeof(int32_t) * total_r_reads, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(rds_l_cnt_offset_d, rds_l_cnt_offset_h, sizeof(int32_t) * vec_size, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(rds_r_cnt_offset_d, rds_r_cnt_offset_h, sizeof(int32_t) * vec_size, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(ctg_seqs_d, ctg_seqs_h, sizeof(char) * max_ctg_size * vec_size, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(reads_left_d, reads_left_h, sizeof(char) * total_l_reads * max_read_size, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(reads_right_d, reads_right_h, sizeof(char) * total_r_reads * max_read_size, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(depth_d, depth_h, sizeof(double) * vec_size, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(quals_right_d, quals_right_h, sizeof(char) * total_r_reads * max_read_size, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(quals_left_d, quals_left_h, sizeof(char) * total_l_reads * max_read_size, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(term_counts_d, term_counts_h, sizeof(int64_t)*3, cudaMemcpyHostToDevice));

    //call kernel here, one thread per contig
    unsigned total_threads = vec_size;
    
    if(DEBUG_PRINT_CPU)
        print_vals("Calling Kernel...");
    iterative_walks_kernel<<<1,2>>>(cid_d, ctg_seq_offsets_d, ctg_seqs_d, reads_left_d, reads_right_d, quals_right_d, quals_left_d, reads_l_offset_d, reads_r_offset_d, rds_l_cnt_offset_d, rds_r_cnt_offset_d, 
    depth_d, d_ht, max_mer_len, 22, 0, term_counts_d, 0, 0, 0, max_read_size, max_read_count, longest_walks_d, mer_walk_temp_d);

    CUDA_CHECK(cudaFree(term_counts_d));
    if(DEBUG_PRINT_CPU)
        print_vals("Exiting...");
    //free up the memory

    
    



    return 0;
}


    // //allocate memory for host and device char arrays
    // h_contigs = new char[max_ctg_size*contigs.size()];
    // CUDA_CHECK(cudaMalloc(&d_contigs, sizeof(char)*max_ctg_size*contigs.size()));
    // //memory for offset array
    // h_offset_arr  = new int[contigs.size()];
    // CUDA_CHECK(cudaMalloc(&d_offset_arr, sizeof(int)*contigs.size()));
    // //memory allocation on device for local hashtables
    // CUDA_CHECK(cudaMalloc(&d_ht, sizeof(loc_ht)*HT_SIZE*TOT_THREADS));


    // //convert strings to char* and prepare offset array
    // for(int i = 0; i < contigs.size(); i++)
    // {
    //     char* seq_ptr = h_contigs + offset_sum;
    //     memcpy(seq_ptr, contigs[i].c_str(), contigs[i].size());
    //     offset_sum += contigs[i].size();
    //     h_offset_arr[i] = offset_sum;
    // }
    // //moving data to device
    // CUDA_CHECK(cudaMemcpy(d_contigs, h_contigs, sizeof(char)*max_ctg_size*contigs.size(), cudaMemcpyHostToDevice));
    // CUDA_CHECK(cudaMemcpy(d_offset_arr, h_offset_arr, sizeof(int)*contigs.size(), cudaMemcpyHostToDevice));
    // //print_vals("launch kernel", "another print");
    // ht_kernel<<<1,1>>>(d_ht, d_contigs, d_offset_arr, KMER_SZ);

    // CUDA_CHECK(cudaFree(d_ht));

    // std::unordered_map<std::string, int> my_map;

    // my_map.insert({"this_map",10});
    // my_map.insert({"map",100});
    // my_map.insert({"this_map", 15});
    // print_vals("total_insertions:", "count");
   // print_vals(my_map["this_map"], my_map["map"]);
