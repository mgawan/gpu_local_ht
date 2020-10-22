#include <unordered_map>
//#include <iostream>
#include <typeinfo>
#include <random>
#include <algorithm>
#include <cstring>
#include <fstream>
#include "helper.hpp"
#include "kernel.hpp"
#define TOT_THREADS 1
#define KMER_SZ 4

double total_data_time = 0;
double total_kernel_time = 0;

std::ofstream ofile("contig-test-1.dat");
void call_kernel(std::vector<CtgWithReads>& data_in, int32_t max_ctg_size, int32_t total_r_reads, int32_t total_l_reads, int32_t max_read_size, int32_t max_r_count, int32_t max_l_count, int insert_avg, int insert_stddev);
//TODO: DO it such that contigs with now left or righ reads are offloaded to kernels, then try to make separate left and right kernels so that contigs only right reads are launched in right kernel
// and contigs with only left are launched in left kernels.
int main (int argc, char* argv[]){

    std::string in_file = argv[1];
    std::vector<CtgWithReads> data_in;
    int32_t max_ctg_size, total_r_reads, total_l_reads, max_read_size, max_r_count, max_l_count;
    read_locassm_data(&data_in, in_file, max_ctg_size, total_r_reads, total_l_reads, max_read_size,max_r_count, max_l_count);
    print_vals("total exts:",data_in.size());
        timer final_time;
    final_time.timer_start();
    int slice_size = 10000;
    int iterations = (data_in.size() + slice_size)/slice_size;
    // print_vals("Total Contigs:", data_in.size());
    // print_vals("Slices:", iterations);
    
    for(int i = 0; i < iterations; i++){
        int left_over;
        if(iterations - 1 == i)
            left_over = data_in.size() % slice_size;
        else
            left_over = slice_size;
        std::vector<CtgWithReads> slice_data (&data_in[i*slice_size], &data_in[i*slice_size + left_over]);

        print_vals("*** NEW SLICE LAUNCHING ***","Current slice size:", slice_data.size(), "Slice ID:",i);
        print_vals("max:",max_r_count);
        call_kernel(slice_data, max_ctg_size, total_r_reads, total_l_reads, max_read_size, max_r_count, max_l_count, 121, 251);
    }
    

  //TODO: free up memory
  final_time.timer_end();
  print_vals("Total Overall Time:", final_time.get_total_time());
  print_vals("Total Data Transfer Time:", total_data_time);
  print_vals("Total Kernel Time:", total_kernel_time);
    return 0;
}

void call_kernel(std::vector<CtgWithReads>& data_in, int32_t max_ctg_size, int32_t total_r_reads, int32_t total_l_reads, int32_t max_read_size, int32_t max_r_count, int32_t max_l_count, int insert_avg, int insert_stddev){
    int32_t vec_size = data_in.size();
    int32_t max_read_count = max_r_count>max_l_count ? max_r_count : max_l_count;
    int max_walk_len = insert_avg + 2 * insert_stddev;

    //host allocations for converting loc_assm_data to prim types
    int32_t *cid_h = new int32_t[vec_size];
    char *ctg_seqs_h = new char[max_ctg_size * vec_size];
    char * ctgs_seqs_rc_h = new char[max_ctg_size * vec_size];// not requried on GPU, ctg space will be re-used
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
    int32_t *term_counts_h = new int32_t[3];
    char* longest_walks_r_h = new char[vec_size * max_walk_len];
    char* longest_walks_l_h = new char[vec_size * max_walk_len]; // not needed on device, will re-use right walk memory
    int* final_walk_lens_r_h = new int[vec_size];
    int* final_walk_lens_l_h = new int[vec_size]; // not needed on device, will re use right walk memory

    int max_mer_len = 33;

    //compute total device memory required:
    size_t total_dev_mem = sizeof(int32_t) * vec_size * 4 + sizeof(int32_t) * total_l_reads
                           + sizeof(int32_t) * total_r_reads + sizeof(char) * max_ctg_size * vec_size
                           + sizeof(char) * total_l_reads * max_read_size + sizeof(char) * total_r_reads * max_read_size
                           + sizeof(double) * vec_size + sizeof(char) * total_r_reads * max_read_size 
                           + sizeof(char) * total_l_reads * max_read_size + sizeof(int64_t)*3
                           + sizeof(loc_ht)*(max_read_size*max_read_count)*vec_size
                           + sizeof(char)*vec_size * max_walk_len
                           + (max_mer_len + max_walk_len) * sizeof(char) * vec_size
                           + sizeof(int) * vec_size;
    print_vals("Total GPU Mem. (GB) requested:", (double)total_dev_mem/(1024*1024*1024));
    print_vals("max_l_count:",max_l_count,"max_r_count:", max_r_count);
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    print_vals("GPU Mem Avail (GB):", (double)prop.totalGlobalMem/(1024*1024*1024));

    timer gpu_total;
    gpu_total.timer_start();
    //device allocations for loc_assm_data
    int32_t *cid_d, *ctg_seq_offsets_d, *reads_l_offset_d, *reads_r_offset_d; 
    int32_t *rds_l_cnt_offset_d, *rds_r_cnt_offset_d;
    char *ctg_seqs_d, *reads_left_d, *reads_right_d, *quals_left_d, *quals_right_d;
    char *longest_walks_d, *mer_walk_temp_d;
    double *depth_d;
    int32_t *term_counts_d;
    loc_ht *d_ht;
    loc_ht_bool *d_ht_bool;
    int* final_walk_lens_d;

    timer cuda_alloc_time;
    cuda_alloc_time.timer_start();
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
    CUDA_CHECK(cudaMalloc(&longest_walks_d, sizeof(char)*vec_size * max_walk_len));
    CUDA_CHECK(cudaMalloc(&mer_walk_temp_d, (max_mer_len + max_walk_len) * sizeof(char) * vec_size));
    CUDA_CHECK(cudaMalloc(&d_ht_bool, sizeof(loc_ht_bool) * vec_size * max_walk_len));
    CUDA_CHECK(cudaMalloc(&final_walk_lens_d, sizeof(int) * vec_size));
    cuda_alloc_time.timer_end();
    //convert the loc_assem data to primitive structures for device
    int32_t ctgs_offset_sum = 0;
    int32_t reads_r_offset_sum = 0;
    int32_t reads_l_offset_sum = 0;
    int read_l_index = 0, read_r_index = 0;
    timer data_packing;
    data_packing.timer_start();

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
    #ifdef DEBUG_PRINT_CPU
        print_vals("Completed reading in file and converting data to primitive types...");
    #endif
    for(int i = 0; i < 3; i++){
        term_counts_h[i] = 0;
    }
    data_packing.timer_end();
    timer gpu_transfer;

     print_vals("Host to Device Transfer...");
    //move the data to device
    gpu_transfer.timer_start();
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
    CUDA_CHECK(cudaMemcpy(term_counts_d, term_counts_h, sizeof(int32_t)*3, cudaMemcpyHostToDevice));

    gpu_transfer.timer_end();
    double transfer_time = gpu_transfer.get_total_time();
    

    //call kernel here, one thread per contig
    unsigned total_threads = vec_size;
    unsigned thread_per_blk = 512;
    unsigned blocks = (vec_size + thread_per_blk)/thread_per_blk;
    

    print_vals("Calling Kernel with blocks:", blocks, "Threads:", thread_per_blk);
    int64_t sum_ext, num_walks;
    timer kernel_time;

    kernel_time.timer_start();
    iterative_walks_kernel<<<blocks,thread_per_blk>>>(cid_d, ctg_seq_offsets_d, ctg_seqs_d, reads_left_d, reads_right_d, quals_right_d, quals_left_d, reads_l_offset_d, reads_r_offset_d, rds_l_cnt_offset_d, rds_r_cnt_offset_d, 
    depth_d, d_ht, d_ht_bool, max_mer_len, term_counts_d, num_walks, max_walk_len, sum_ext, max_read_size, max_read_count, longest_walks_d, mer_walk_temp_d, final_walk_lens_d, vec_size);
    CUDA_CHECK(cudaDeviceSynchronize());
    kernel_time.timer_end();
    double right_kernel_time = kernel_time.get_total_time();

    print_vals("Device to Host Transfer...", "Copying back right walks");
    gpu_transfer.timer_start();
    CUDA_CHECK(cudaMemcpy(longest_walks_r_h, longest_walks_d, sizeof(char) * vec_size * max_walk_len, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(final_walk_lens_r_h, final_walk_lens_d, sizeof(int32_t) * vec_size, cudaMemcpyDeviceToHost)); 
    gpu_transfer.timer_end();
    transfer_time += gpu_transfer.get_total_time();

    #ifdef DEBUG_PRINT_CPU
        print_vals("printing walks for right extension:");
        print_vals("longest walk for contig 0 has length:", final_walk_lens_r_h[0]);
        for(int i = 0; i < final_walk_lens_r_h[0]; i++){
            std::cout<<longest_walks_r_h[i];
        }
        std::cout<<std::endl;
        print_vals("longest walk for contig 1 has length:", final_walk_lens_r_h[1]);
        for(int i = 0; i < final_walk_lens_r_h[1]; i++){
            std::cout<<longest_walks_r_h[max_walk_len*1 + i];
        }
        std::cout<<std::endl;

        print_vals("longest walk for contig 2 has length:", final_walk_lens_r_h[2]);
        for(int i = 0; i < final_walk_lens_r_h[2]; i++){
            std::cout<<longest_walks_r_h[max_walk_len*2 + i];
        }
        std::cout<<std::endl;
        print_vals("longest walk for contig 3 has length:", final_walk_lens_r_h[3]);
        for(int i = 0; i < final_walk_lens_r_h[3]; i++){
            std::cout<<longest_walks_r_h[max_walk_len*3 + i];
        }
        std::cout<<std::endl;
                std::cout<<std::endl;
                print_vals("longest walk for contig 4 has length:", final_walk_lens_r_h[4]);
        for(int i = 0; i < final_walk_lens_r_h[4]; i++){
            std::cout<<longest_walks_r_h[max_walk_len*4 + i];
        }
        std::cout<<std::endl;

        print_vals("longest walk for contig 5 has length:", final_walk_lens_r_h[5]);
        for(int i = 0; i < final_walk_lens_r_h[5]; i++){
            std::cout<<longest_walks_r_h[max_walk_len*5 + i];
        }
        std::cout<<std::endl;
        print_vals("longest walk for contig 6 has length:", final_walk_lens_r_h[6]);
        for(int i = 0; i < final_walk_lens_r_h[6]; i++){
            std::cout<<longest_walks_r_h[max_walk_len*6 + i];
        }
        std::cout<<std::endl;
    #endif
    // perform revcomp of contig sequences and launch kernel with left reads, TODO: move right and left reads data separately
    print_vals("revcomp-ing the contigs for next kernel");
    timer rev_comp_;
    rev_comp_.timer_start();
    for(int j = 0; j < vec_size; j++){
        int size_lst;
        char* curr_seq;
        char* curr_seq_rc;
        if(j == 0){
            size_lst = ctg_seq_offsets_h[j];
            curr_seq = ctg_seqs_h;
            curr_seq_rc = ctgs_seqs_rc_h;
        }
        else{
            size_lst = ctg_seq_offsets_h[j] - ctg_seq_offsets_h[j-1];
            curr_seq = ctg_seqs_h + ctg_seq_offsets_h[j - 1];
            curr_seq_rc = ctgs_seqs_rc_h + ctg_seq_offsets_h[j - 1];
        }
        revcomp(curr_seq, curr_seq_rc, size_lst);
        #ifdef DEBUG_PRINT_CPU   
        print_vals("orig seq:");
        for(int h = 0; h < size_lst; h++)
            std::cout<<curr_seq[h];
        std::cout << std::endl;
        print_vals("recvomp seq:");
        for(int h = 0; h < size_lst; h++)
            std::cout<<curr_seq_rc[h];
        std::cout << std::endl;    
        #endif   

    }

    rev_comp_.timer_end();

    //cpying rev comped ctgs to device on same memory as previous ctgs
    gpu_transfer.timer_start();
    CUDA_CHECK(cudaMemcpy(ctg_seqs_d, ctgs_seqs_rc_h, sizeof(char) * max_ctg_size * vec_size, cudaMemcpyHostToDevice));
    gpu_transfer.timer_end();
    transfer_time += gpu_transfer.get_total_time();
    // launching kernel by swapping right and left reads, TODO: make this correct
    kernel_time.timer_start();
    iterative_walks_kernel<<<blocks,thread_per_blk>>>(cid_d, ctg_seq_offsets_d, ctg_seqs_d, reads_right_d, reads_left_d, quals_left_d, quals_right_d, reads_r_offset_d, reads_l_offset_d, rds_r_cnt_offset_d, rds_l_cnt_offset_d, 
    depth_d, d_ht, d_ht_bool, max_mer_len, term_counts_d, num_walks, max_walk_len, sum_ext, max_read_size, max_read_count, longest_walks_d, mer_walk_temp_d, final_walk_lens_d, vec_size);
    CUDA_CHECK(cudaDeviceSynchronize());
    kernel_time.timer_end();
    double left_kernel_time = kernel_time.get_total_time();
    print_vals("Device to Host Transfer...", "Copying back left walks");
    
    gpu_transfer.timer_start();
    CUDA_CHECK(cudaMemcpy(longest_walks_l_h, longest_walks_d, sizeof(char) * vec_size * max_walk_len, cudaMemcpyDeviceToHost)); // copy back left walks
    CUDA_CHECK(cudaMemcpy(final_walk_lens_l_h, final_walk_lens_d, sizeof(int32_t) * vec_size, cudaMemcpyDeviceToHost)); 
    gpu_transfer.timer_end();
    transfer_time += gpu_transfer.get_total_time();
    #ifdef DEBUG_PRINT_CPU
        print_vals("printing walks for left extension:");
        print_vals("longest walk for contig 0 has length:", final_walk_lens_l_h[0]);
        for(int i = 0; i < final_walk_lens_l_h[0]; i++){
            std::cout<<longest_walks_l_h[i];
        }
        std::cout<<std::endl;
        print_vals("longest walk for contig 1 has length:", final_walk_lens_l_h[1]);
        for(int i = 0; i < final_walk_lens_l_h[1]; i++){
            std::cout<<longest_walks_l_h[max_walk_len*1 + i];
        }
        std::cout<<std::endl;

        print_vals("longest walk for contig 2 has length:", final_walk_lens_l_h[2]);
        for(int i = 0; i < final_walk_lens_l_h[2]; i++){
            std::cout<<longest_walks_l_h[max_walk_len*2 + i];
        }
        std::cout<<std::endl;
        print_vals("longest walk for contig 3 has length:", final_walk_lens_l_h[3]);
        for(int i = 0; i < final_walk_lens_l_h[3]; i++){
            std::cout<<longest_walks_l_h[max_walk_len*3 + i];
        }
        std::cout<<std::endl;
                print_vals("longest walk for contig 4 has length:", final_walk_lens_l_h[4]);
        for(int i = 0; i < final_walk_lens_l_h[4]; i++){
            std::cout<<longest_walks_l_h[max_walk_len*4 + i];
        }
        std::cout<<std::endl;

        print_vals("longest walk for contig 5 has length:", final_walk_lens_l_h[5]);
        for(int i = 0; i < final_walk_lens_l_h[5]; i++){
            std::cout<<longest_walks_l_h[max_walk_len*5 + i];
        }
        std::cout<<std::endl;
        print_vals("longest walk for contig 6 has length:", final_walk_lens_l_h[6]);
        for(int i = 0; i < final_walk_lens_l_h[6]; i++){
            std::cout<<longest_walks_l_h[max_walk_len*6 + i];
        }
        std::cout<<std::endl;
    #endif
    
    CUDA_CHECK(cudaFree(term_counts_d));
    CUDA_CHECK(cudaFree(cid_d));
    CUDA_CHECK(cudaFree(ctg_seq_offsets_d));
    CUDA_CHECK(cudaFree(reads_l_offset_d));
    CUDA_CHECK(cudaFree(reads_r_offset_d));
    CUDA_CHECK(cudaFree(rds_l_cnt_offset_d));
    CUDA_CHECK(cudaFree(rds_r_cnt_offset_d));
    CUDA_CHECK(cudaFree(ctg_seqs_d));
    CUDA_CHECK(cudaFree(reads_left_d));
    CUDA_CHECK(cudaFree(reads_right_d));
    CUDA_CHECK(cudaFree(depth_d));
    CUDA_CHECK(cudaFree(quals_right_d));
    CUDA_CHECK(cudaFree(quals_left_d));
    CUDA_CHECK(cudaFree(d_ht)); 
    CUDA_CHECK(cudaFree(longest_walks_d));
    CUDA_CHECK(cudaFree(mer_walk_temp_d));
    CUDA_CHECK(cudaFree(d_ht_bool));
    CUDA_CHECK(cudaFree(final_walk_lens_d));
    //adding pre and post walks in strings and printing the results out in a file.
    
    for(int i = 0; i< vec_size; i++){
        if(final_walk_lens_l_h[i] != 0){
            std::string left(&longest_walks_l_h[max_walk_len*i],final_walk_lens_l_h[i]);
            std::string left_rc = revcomp(left);
            data_in[i].seq.insert(0,left_rc);  
        }
        if(final_walk_lens_r_h[i] != 0){
            std::string right(&longest_walks_r_h[max_walk_len*i],final_walk_lens_r_h[i]);
            data_in[i].seq += right;
        }
        ofile << data_in[i].cid<<" "<<data_in[i].seq<<std::endl;
    }

    gpu_total.timer_end();
    print_vals("Total Time:", gpu_total.get_total_time() );
    //print_vals("Total Data Transfer Time:", transfer_time);
    total_data_time += transfer_time;
   // print_vals("Total Right Kernel Time:", right_kernel_time);
    total_kernel_time += right_kernel_time;
   // print_vals("Total Left Kernel Time:", left_kernel_time);
    total_kernel_time += left_kernel_time;
   // print_vals("Total rec comp Time:", rev_comp_.get_total_time());
   // print_vals("Total Data packing time:", data_packing.get_total_time());
    print_vals("Total Cuda malloc time:", cuda_alloc_time.get_total_time());

}
