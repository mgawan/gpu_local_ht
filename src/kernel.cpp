#include "kernel.hpp"


//TODO: all the hashtable entries need to be set to empty, figure that out
__device__ void print_mer(cstr_type& mer){
    for(int i = 0; i < mer.length; i++){
        printf("%c",mer.start_ptr[i]);
    }
    printf("\n");
}

__device__ void cstr_copy(cstr_type& str1, cstr_type& str2){

    for(int i = 0; i < str1.length; i++){
        str1.start_ptr[i] = str2.start_ptr[i];
    }
    str1.length = str2.length;
}
// __global__ void ht_kernel(loc_ht* ht, char* contigs, int* offset_sum, int kmer_size){
//     int idx = blockIdx.x*gridDim.x + threadIdx.x;
//     cstr_type loc_contig;
//     loc_ht* thread_ht;
//     thread_ht = ht + idx*HT_SIZE;

//     for(int i = 0; i < HT_SIZE; i++){
//         thread_ht[i].key.length = EMPTY;
//     }

//     if(idx == 0){
//         loc_contig.start_ptr = contigs;
//         loc_contig.length = offset_sum[idx];
//     }else{
//         loc_contig.start_ptr = contigs + offset_sum[idx];
//         loc_contig.length = offset_sum[idx] - offset_sum[idx-1];
//     }
    
//       cstr_type kmer;
//       kmer.start_ptr = loc_contig.start_ptr;
//       kmer.length = kmer_size;
//       printf("kmer size:%d\n",kmer.length);
//       printf("contig size:%d\n",loc_contig.length);
//     for(int i = 0; i < (loc_contig.length - (kmer_size-1)); i++){
//         char *mystr = kmer.start_ptr;
//         for(int j = 0; j < kmer.length; j++){
//             printf("%c",mystr[j]);
//         }
//         printf("\n");
//         cstr_type temp = ht_get(thread_ht, kmer);
//         if(temp.length == -1)
//             ht_insert(thread_ht, kmer, loc_contig);
//         else{
//             printf("key exists:\n");
//             for(int j = 0; j < temp.length; j++){
//                 printf("%c",temp.start_ptr[j]);
//         }
//         printf("\n");
//         }

//         kmer.start_ptr =  kmer.start_ptr + 1;
//     }


// }

//TODO: make sure that the returned hash value is wthin the range of HT size
//TODO: find a better hash func for strings
__device__ unsigned hash_func(cstr_type key, uint32_t max_size){
    unsigned hash, i;
    for(hash = i = 0; i < key.length; ++i)
    {
        hash += key.start_ptr[i];
        hash += (hash << 10);
        hash ^= (hash >> 6);
    }
    hash += (hash << 3);
    hash ^= (hash >> 11);
    hash += (hash << 15);
    //TODO: this way of limiting hash value is not good, try to find a better way
    return hash%max_size;//(hash & (HT_SIZE - 1));
}

__device__ void ht_insert(loc_ht* thread_ht, cstr_type kmer_key, MerFreqs mer_val, uint32_t max_size){
    unsigned hash_val = hash_func(kmer_key, max_size);
    unsigned orig_hash = hash_val;
    //int count = 0; // for debugging
    while(true){
        int if_empty = thread_ht[hash_val].key.length; // length is set to some unimaginable number to indicate if its empty
        if(if_empty == EMPTY){ //the case where there is a key but no val, will not happen

           // printf("hash_val:%d, orig_hash:%d, attemp:%d\n",hash_val, orig_hash, count); // for debugging
            thread_ht[hash_val].key = kmer_key;
            thread_ht[hash_val].val = mer_val;
            return;
        }
        hash_val = (hash_val +1 ) % max_size;//(hash_val + 1) & (HT_SIZE-1);
        //count++; //for debugging

    }
}
//overload for bool vals
__device__ void ht_insert(loc_ht_bool* thread_ht, cstr_type kmer_key, bool bool_val, uint32_t max_size){
    unsigned hash_val = hash_func(kmer_key, max_size);
    unsigned orig_hash = hash_val;
    //int count = 0; // for debugging
    while(true){
        int if_empty = thread_ht[hash_val].key.length; // length is set to some unimaginable number to indicate if its empty
        if(if_empty == EMPTY){ //the case where there is a key but no val, will not happen

           // printf("hash_val:%d, orig_hash:%d, attemp:%d\n",hash_val, orig_hash, count); // for debugging
            thread_ht[hash_val].key = kmer_key;
            thread_ht[hash_val].val = bool_val;
            return;
        }
        hash_val = (hash_val +1 ) % max_size;//(hash_val + 1) & (HT_SIZE-1);
        //count++; //for debugging

    }
}

// __device__ void ht_delete(loc_ht* thread_ht, cstr_type kmer_key){
//     int hash_val = hash_func(kmer_key);
//     while(true){
//         if(thread_ht[hash_val].key == kmer_key){
//             thread_ht[hash_val].key.length = EMPTY;
//             return;
//         }
//         if(thread_ht[hash_val].key.length == EMPTY){
//             return;
//         }
//         hash_val = (hash_val + 1) & (HT_SIZE -1);
//     }
// }

__device__ 
loc_ht& ht_get(loc_ht* thread_ht, cstr_type kmer_key, uint32_t max_size){
    unsigned hash_val = hash_func(kmer_key, max_size);
    unsigned orig_hash = hash_val;
    
    while(true){
        if(thread_ht[hash_val].key.length == EMPTY){
            return thread_ht[hash_val];
        }
        else if(thread_ht[hash_val].key == kmer_key){
            //printf("key found, returning\n");// keep this for debugging
            return thread_ht[hash_val];
        }
        hash_val = (hash_val +1 ) %max_size;//hash_val = (hash_val + 1) & (HT_SIZE -1);
        if(hash_val == orig_hash){ // loop till you reach the same starting positions and then return error
            printf("*****end reached, hashtable full*****\n"); // for debugging
            printf("*****end reached, hashtable full*****\n");
            printf("*****end reached, hashtable full*****\n");
           // return loc_ht(cstr_type(NULL,-1), MerFreqs());
        }
    }

}

//TODO:use some OOP technique for implementing the bool table, may be inheritence?
//overload for bool vals
__device__ 
loc_ht_bool& ht_get(loc_ht_bool* thread_ht, cstr_type kmer_key, uint32_t max_size){
    unsigned hash_val = hash_func(kmer_key, max_size);
    unsigned orig_hash = hash_val;
    
    while(true){
        if(thread_ht[hash_val].key.length == EMPTY){
            return thread_ht[hash_val];
        }
        else if(thread_ht[hash_val].key == kmer_key){
            //printf("key found, returning\n");// keep this for debugging
            return thread_ht[hash_val];
        }
        hash_val = (hash_val +1 ) %max_size;//hash_val = (hash_val + 1) & (HT_SIZE -1);
        if(hash_val == orig_hash){ // loop till you reach the same starting positions and then return error
            printf("*****end reached, hashtable full*****\n"); // for debugging
            printf("*****end reached, hashtable full*****\n");
            printf("*****end reached, hashtable full*****\n");
           // return loc_ht(cstr_type(NULL,-1), MerFreqs());
        }
    }

}

//TODO: intialize the bool table in kernel main
//TODO: check if we need longest walk in this function
__device__ char walk_mers(loc_ht* thrd_loc_ht, loc_ht_bool* thrd_ht_bool, uint32_t max_ht_size, uint32_t& mer_len, char* mer_walk_temp, char* longest_walk, const int idx, int max_walk_len){
    char walk_result = 'X';
    int walk_length = 0;
    cstr_type mer(mer_walk_temp, mer_len);
    cstr_type walk(mer_walk_temp + mer_len, walk_length); // walk pointer starts at the end of initial mer pointer

    for( int nsteps = 0; nsteps < max_walk_len; nsteps++){
        //check if there is a cycle in graph
        loc_ht_bool &temp_mer_loop = ht_get(thrd_ht_bool, mer, max_walk_len);
        if(temp_mer_loop.key.length == EMPTY){ // if the mer has not been visited, add it to the table and mark visited
            temp_mer_loop.key = mer;
            temp_mer_loop.val = true;
        }else{
            walk_result = 'R'; // if the table already contains this mer then a cycle exits, return the walk with repeat.
            break;
        }

        loc_ht &temp_mer = ht_get(thrd_loc_ht, mer, max_ht_size);
        if(temp_mer.key.length == EMPTY){//if mer is not found then dead end reached, terminate the walk
            walk_result = 'X';
            break;
        }
        char ext = temp_mer.val.ext;
        if(ext == 'F' || ext == 'X'){ // if the table points that ext is fork or dead end the terminate the walk
            walk_result = ext;
            break;
        }
        mer.start_ptr = mer.start_ptr + 1; // increment the mer pointer and append the ext
        mer.start_ptr[mer.length-1] = ext; // walk pointer points at the end of initial mer point.
        walk.length++;
        
    }

    return walk_result;
}
// // return the result of the walk (f, r or x)
// static char walk_mers(MerMap &mers_ht, string &mer, string &walk, int mer_len, int walk_len_limit) {
//   HASH_TABLE<string, bool> loop_check_ht;
//   char walk_result = 'X';
//   for (int nsteps = 0; nsteps < walk_len_limit; nsteps++) {
//     if (!(nsteps % 10)) progress();
//     // check for a cycle in the graph
//     if (loop_check_ht.find(mer) != loop_check_ht.end()) {
//       walk_result = 'R';
//       break;
//     } else {
//       loop_check_ht.insert({mer, true});
//     }
//     auto it = mers_ht.find(mer);
//     if (it == mers_ht.end()) {
//       walk_result = 'X';
//       break;
//     }
//     char ext = it->second.ext;
//     if (ext == 'F' || ext == 'X') {
//       walk_result = ext;
//       break;
//     }
//     mer.erase(0, 1);
//     mer += ext;
//     walk += ext;
//   }
//   return walk_result;
// }

__device__ 
void count_mers(loc_ht* thrd_loc_ht, char* loc_r_reads, uint32_t max_ht_size, char* loc_r_quals, int32_t* reads_r_offset, int32_t& r_rds_cnt, 
int32_t* rds_count_r_sum, double& loc_ctg_depth, uint32_t& mer_len, uint32_t& qual_offset, int64_t& excess_reads, const int idx){
    cstr_type read;
    cstr_type qual;
    uint32_t running_sum_len = 0;
    int test = 1;
    if(DEBUG_PRINT_GPU && idx == test)
        printf("inside_count_mers\n");
    for(int i = 0; i < r_rds_cnt; i++){
        if(DEBUG_PRINT_GPU)
            printf("read loop iter:%d\n",i);
        //TODO: pass idx here
        read.start_ptr = loc_r_reads + running_sum_len;
        qual.start_ptr = loc_r_quals + running_sum_len;
        if(i == 0){
            if(idx == 0){
                read.length = reads_r_offset[(rds_count_r_sum[idx] - r_rds_cnt) + i];
                qual.length = reads_r_offset[(rds_count_r_sum[idx] - r_rds_cnt) + i];
                if(DEBUG_PRINT_GPU && idx == test)
                    printf("rds_count_r_sum[idx]:%d, rds_cnt:%d, reads_offset_0:%d\n",rds_count_r_sum[idx], r_rds_cnt, read.length);
                }
            else{   
                read.length = reads_r_offset[(rds_count_r_sum[idx] - r_rds_cnt) + i] - reads_r_offset[(rds_count_r_sum[idx - 1] -1)];
                qual.length = reads_r_offset[(rds_count_r_sum[idx] - r_rds_cnt) + i] - reads_r_offset[(rds_count_r_sum[idx - 1] -1)];
                if(DEBUG_PRINT_GPU && idx == test)
                    printf("rds_count_r_sum[idx]:%d, rds_cnt:%d, reads_offset_0:%d\n",rds_count_r_sum[idx], r_rds_cnt, read.length);
                
                }
            }
        else{
            read.length = reads_r_offset[(rds_count_r_sum[idx] - r_rds_cnt) + i] - reads_r_offset[(rds_count_r_sum[idx] - r_rds_cnt) + (i-1)];
            qual.length = reads_r_offset[(rds_count_r_sum[idx] - r_rds_cnt) + i] - reads_r_offset[(rds_count_r_sum[idx] - r_rds_cnt) + (i-1)];
            if(DEBUG_PRINT_GPU && idx == test)
                 printf("rds_count_r_sum[idx]:%d, rds_cnt:%d, reads_offset_0:%d\n",rds_count_r_sum[idx], r_rds_cnt, read.length);
                
            }
        if(DEBUG_PRINT_GPU && idx == test){
            printf("mer_len:%d, read_len:%d\n",mer_len, read.length);
            printf("read from idx:%d\n", idx);
            print_mer(read);
          }
        if (mer_len > read.length) // skip the read which is smaller than merlen
            continue;
        int num_mers = read.length - mer_len;
        cstr_type mer(read.start_ptr, mer_len);
        for( int start = 0; start < num_mers; start++){
            // if(DEBUG_PRINT_GPU)
            //     print_mer(mer);
            //TODO: on cpu side add a check that if a certain read contains 'N', that is not included, check this with steve, 
            // because searching a single mer for an N is going to be too slow
            loc_ht &temp_Mer = ht_get(thrd_loc_ht, mer, max_ht_size);
            if(temp_Mer.key.length == EMPTY){
                temp_Mer.key = mer;
                temp_Mer.val = {.hi_q_exts = {0}, .low_q_exts = {0}, .ext = 0, .count = 0};
            }
            int ext_pos = start + mer_len;
          //  assert(ext_pos < (int)read.length); // TODO: verify that assert works on gpu, for now commenting it out and replacing with printf
          if(ext_pos >= (int) read.length)
            printf("*********ASSERTION FAILUER IN COUNT_MERS****");
            char ext = read.start_ptr[ext_pos];
            if (ext == 'N') continue; // TODO: why the redundant check?
            int qual_diff = qual.start_ptr[ext_pos] - qual_offset;
            if (qual_diff >= LASSM_MIN_QUAL) temp_Mer.val.low_q_exts.inc(ext, 1);
            if (qual_diff >= LASSM_MIN_HI_QUAL) temp_Mer.val.hi_q_exts.inc(ext, 1);
            
            mer.start_ptr = mer.start_ptr + 1;
        }
       running_sum_len += read.length; // right before the for loop ends, update the prev_len to offset next read correctly
    }

    //setting extension by traversing the completed table
    // TODO: think of a better way to do this
    for (int k = 0; k < max_ht_size; k++) {
        if( thrd_loc_ht[k].key.length != EMPTY){
            thrd_loc_ht[k].val.set_ext(loc_ctg_depth);
            if(DEBUG_PRINT_GPU && idx == test){
                printf("from ht:\n");
                print_mer(thrd_loc_ht[k].key);
                printf("MerFreq.ext:%c, MerFreq.count:%d\n",thrd_loc_ht[k].val.ext,thrd_loc_ht[k].val.count);
                thrd_loc_ht[k].val.hi_q_exts.print();
                thrd_loc_ht[k].val.low_q_exts.print();
            }
        }
    }
}

//same kernel will be used for right and left walks
__global__ void iterative_walks_kernel(int32_t* cid, int32_t* ctg_offsets, char* contigs, 
char* reads_l, char* reads_r, char* quals_r, char* quals_l, int32_t* reads_l_offset, int32_t* reads_r_offset, int32_t* rds_count_l_sum, int32_t* rds_count_r_sum, double* ctg_depth, loc_ht* global_ht, loc_ht_bool* global_ht_bool,
int max_mer_len, int kmer_len, int walk_len_limit, int64_t *term_counts, int64_t num_walks, int64_t max_walk_len, int64_t sum_ext, int32_t max_read_size, int32_t max_read_count, char* longest_walks, char* mer_walk_temp){
    const int idx = threadIdx.x + blockIdx.x * gridDim.x;
    cstr_type loc_ctg;
    char *loc_r_reads, *loc_l_reads, *loc_r_quals, *loc_l_quals;
    int32_t r_rds_cnt, l_rds_cnt, loc_rds_r_offset, loc_rds_l_offset;
    loc_ht* loc_mer_map = global_ht + idx * max_read_size * max_read_count;
    loc_ht_bool* loc_bool_map = global_ht_bool + idx * max_walk_len;
    double loc_ctg_depth = ctg_depth[idx];
    int64_t excess_reads;
    uint32_t qual_offset = 0, max_ht_size = max_read_size * max_read_count;
    char* loc_mer_walk_temp = mer_walk_temp + idx*MAX_WALK_LEN;

    for(uint32_t k = 0; k < max_ht_size; k++){
        loc_mer_map[k].key.length = EMPTY;
    }
    //TODO: initalize hash table find a faster way of doing this

    if(idx == 0){
        loc_ctg.start_ptr = contigs;
        loc_ctg.length = ctg_offsets[idx];
        r_rds_cnt = rds_count_r_sum[idx];
        l_rds_cnt = rds_count_l_sum[idx];
        loc_r_reads = reads_r;
        loc_l_reads = reads_l;
        loc_r_quals = quals_r;
        loc_l_quals = quals_l;
    }else{
        loc_ctg.start_ptr = contigs + ctg_offsets[idx-1];
        loc_ctg.length = ctg_offsets[idx] - ctg_offsets[idx - 1];
        r_rds_cnt = rds_count_r_sum[idx] - rds_count_r_sum[idx - 1];
        l_rds_cnt = rds_count_l_sum[idx] - rds_count_l_sum[idx - 1];
        if (rds_count_r_sum[idx - 1] == 0)
            loc_r_reads = reads_r;
        else
            loc_r_reads = reads_r + reads_r_offset[rds_count_r_sum[idx - 1] - 1]; // you want to start from where previous contigs, last read ends.

        if (rds_count_l_sum[idx - 1] == 0)
            loc_l_reads = reads_l;
        else
            loc_l_reads = reads_l + reads_l_offset[rds_count_l_sum[idx - 1] - 1]; // you want to start from where previous contigs, last read ends.
        
        if (rds_count_r_sum[idx - 1] == 0)
            loc_r_quals = quals_r;
        else
            loc_r_quals = quals_r + reads_r_offset[rds_count_r_sum[idx - 1] - 1]; // you want to start from where previous contigs, last read ends.

        if (rds_count_l_sum[idx - 1] == 0)
            loc_l_quals = quals_l;
        else
            loc_l_quals = quals_l + reads_l_offset[rds_count_l_sum[idx - 1] - 1]; // you want to start from where previous contigs, last read ends. 
    }

    cstr_copy(loc_mer_walk_temp, loc_ctg);

    //main for loop
    //TODO: commenting out the main for loop for testing count_mers
    //for(int mer_len = kmer_len; mer_len >= min_mer_len && mer_len <= max_mer_len; mer_len += shift){
          //TODO: add a check if total number of reads exceeds a certain number/too large, skip that one, may be do this on cpu 
          // to preserve memory on GPU
          //TODO: need to reinitialize the hashtable after each kmer size is done
        uint32_t mer_len = 21;
        if(DEBUG_PRINT_GPU)
            printf("read_count:%d, idx:%d\n",r_rds_cnt, idx);
        if(r_rds_cnt != 0){    //if count is zero, no need to count
            if(DEBUG_PRINT_GPU)
                printf("calling count_mers\n");
            count_mers(loc_mer_map, loc_r_reads, max_ht_size, loc_r_quals, reads_r_offset, r_rds_cnt, rds_count_r_sum, loc_ctg_depth, mer_len, qual_offset, excess_reads, idx);
            __device__ char walk_mers(loc_ht* thrd_loc_ht, loc_ht_bool* thrd_ht_bool, uint32_t max_ht_size, uint32_t& mer_len, char* mer_walk_temp, char* longest_walk, const int idx, int max_walk_len);
            }
   // }

}