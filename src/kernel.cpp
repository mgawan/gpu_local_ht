#include "kernel.hpp"


//TODO: all the hashtable entries need to be set to empty, figure that out
__device__ void print_mer(cstr_type& mer){
    for(int i = 0; i < mer.length; i++){
        printf("%c",mer.start_ptr[i]);
    }
    printf("\n");
}

__device__ void cstr_copy(cstr_type& str1, cstr_type& str2){

    for(int i = 0; i < str2.length; i++){
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
__device__ char walk_mers(loc_ht* thrd_loc_ht, loc_ht_bool* thrd_ht_bool, uint32_t max_ht_size, int& mer_len, cstr_type& mer_walk_temp, cstr_type& longest_walk, cstr_type& walk, const int idx, int max_walk_len){
    char walk_result = 'X';
    int test = 629;
    int walk_length = 0;
    //cstr_type mer(mer_walk_temp, mer_len);
    //cstr_type walk(mer_walk_temp + mer_len, walk_length); // walk pointer starts at the end of initial mer pointer

    for( int nsteps = 0; nsteps < max_walk_len; nsteps++){
        //check if there is a cycle in graph
        loc_ht_bool &temp_mer_loop = ht_get(thrd_ht_bool, mer_walk_temp, max_walk_len);
        if(temp_mer_loop.key.length == EMPTY){ // if the mer has not been visited, add it to the table and mark visited
            temp_mer_loop.key = mer_walk_temp;
            temp_mer_loop.val = true;
        }else{
            walk_result = 'R'; // if the table already contains this mer then a cycle exits, return the walk with repeat.
            #ifdef DEBUG_PRINT_GPU
            if(idx == test)
                printf("breaking at cycle found, res: %c\n", walk_result);
            #endif
            break;
        }

        loc_ht &temp_mer = ht_get(thrd_loc_ht, mer_walk_temp, max_ht_size);
        if(temp_mer.key.length == EMPTY){//if mer is not found then dead end reached, terminate the walk
            walk_result = 'X';
            #ifdef DEBUG_PRINT_GPU
            if(idx == test)
               printf("breaking at mer not found,res: %c\n", walk_result);
            #endif
            break;
        }
        char ext = temp_mer.val.ext;
        if(ext == 'F' || ext == 'X'){ // if the table points that ext is fork or dead end the terminate the walk
            walk_result = ext;
            #ifdef DEBUG_PRINT_GPU
                if(idx == test){
                    printf("breaking at dead end, res: %c\n", walk_result);
                    printf("Mer Looked up:\n");
                    print_mer(mer_walk_temp);
                    printf("ext:%c\n",temp_mer.val.ext);
                    printf("walk with mer_len:%d\n", mer_len);
                    print_mer(walk);
                }
            #endif
            break;
        }

        #ifdef DEBUG_PRINT_GPU
        if(test == idx){
            printf("Mer Looked up:\n");
            print_mer(mer_walk_temp);
            printf("ext:%c\n",temp_mer.val.ext);
            printf("walk with mer_len:%d\n", mer_len);
            print_mer(walk);
        }
        #endif
        mer_walk_temp.start_ptr = mer_walk_temp.start_ptr + 1; // increment the mer pointer and append the ext
        mer_walk_temp.start_ptr[mer_walk_temp.length-1] = ext; // walk pointer points at the end of initial mer point.
        walk.length++;

        
    }
    
    #ifdef DEBUG_PRINT_GPU
    if(idx == test)
        for (int k = 0; k < max_walk_len; k++) {
        if(thrd_ht_bool[k].key.length != EMPTY){
            printf("from bool ht:\n");
            print_mer(thrd_ht_bool[k].key);
            printf("Bool:%d\n",thrd_ht_bool[k].val);
        }
        }
    #endif
    return walk_result;
}

__device__ 
void count_mers(loc_ht* thrd_loc_ht, char* loc_r_reads, uint32_t max_ht_size, char* loc_r_quals, uint32_t* reads_r_offset, uint32_t& r_rds_cnt, 
uint32_t* rds_count_r_sum, double& loc_ctg_depth, int& mer_len, uint32_t& qual_offset, int64_t& excess_reads, const int idx){
    cstr_type read;
    cstr_type qual;
    uint32_t running_sum_len = 0;
    int test = 629;
    #ifdef DEBUG_PRINT_GPU
    if(DEBUG_PRINT_GPU && idx == test)
        printf("inside_count_mers\n");
    #endif
    for(int i = 0; i < r_rds_cnt; i++){
        #ifdef DEBUG_PRINT_GPU
        if(DEBUG_PRINT_GPU && idx == test)
            printf("read loop iter:%d\n",i);
        #endif
        //TODO: pass idx here
        read.start_ptr = loc_r_reads + running_sum_len;
        qual.start_ptr = loc_r_quals + running_sum_len;
        if(i == 0){
            if(idx == 0){
                read.length = reads_r_offset[(rds_count_r_sum[idx] - r_rds_cnt) + i];
                qual.length = reads_r_offset[(rds_count_r_sum[idx] - r_rds_cnt) + i];
                #ifdef DEBUG_PRINT_GPU
                if(DEBUG_PRINT_GPU && idx == test)
                    printf("rds_count_r_sum[idx]:%d, rds_cnt:%d, read_length:%d\n",rds_count_r_sum[idx], r_rds_cnt, read.length);
                #endif
                }
            else{   
                read.length = reads_r_offset[(rds_count_r_sum[idx] - r_rds_cnt) + i] - reads_r_offset[(rds_count_r_sum[idx - 1] -1)];
                qual.length = reads_r_offset[(rds_count_r_sum[idx] - r_rds_cnt) + i] - reads_r_offset[(rds_count_r_sum[idx - 1] -1)];
                #ifdef DEBUG_PRINT_GPU
                if(DEBUG_PRINT_GPU && idx == test)
                    printf("rds_count_r_sum[idx]:%d, rds_cnt:%d, reads_offset_0:%d\n",rds_count_r_sum[idx], r_rds_cnt, read.length);
                #endif
                
                }
            }
        else{
            read.length = reads_r_offset[(rds_count_r_sum[idx] - r_rds_cnt) + i] - reads_r_offset[(rds_count_r_sum[idx] - r_rds_cnt) + (i-1)];
            qual.length = reads_r_offset[(rds_count_r_sum[idx] - r_rds_cnt) + i] - reads_r_offset[(rds_count_r_sum[idx] - r_rds_cnt) + (i-1)];
            #ifdef DEBUG_PRINT_GPU
            if(DEBUG_PRINT_GPU && idx == test)
                 printf("rds_count_r_sum[idx]:%d, rds_cnt:%d, reads_offset_0:%d\n",rds_count_r_sum[idx], r_rds_cnt, read.length);
            #endif
                
            }
        #ifdef DEBUG_PRINT_GPU
        if(DEBUG_PRINT_GPU && idx == test){
            printf("mer_len:%d, read_len:%d\n",mer_len, read.length);
            printf("read from idx:%d\n", idx);
            print_mer(read);
          }
          #endif
        if (mer_len > read.length) // skip the read which is smaller than merlen
            continue;
        int num_mers = read.length - mer_len;
        cstr_type mer(read.start_ptr, mer_len);
        for( int start = 0; start < num_mers; start++){
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
            #ifdef DEBUG_PRINT_GPU
            if(DEBUG_PRINT_GPU && idx == test){
                printf("from ht:\n");
                print_mer(thrd_loc_ht[k].key);
                printf("MerFreq.ext:%c, MerFreq.count:%d\n",thrd_loc_ht[k].val.ext,thrd_loc_ht[k].val.count);
                thrd_loc_ht[k].val.hi_q_exts.print();
                thrd_loc_ht[k].val.low_q_exts.print();
            }
            #endif
        }
    }
}

//same kernel will be used for right and left walks
__global__ void iterative_walks_kernel(uint32_t* cid, uint32_t* ctg_offsets, char* contigs, 
char* reads_l, char* reads_r, char* quals_r, char* quals_l, uint32_t* reads_l_offset, uint32_t* reads_r_offset, uint32_t* rds_count_l_sum, uint32_t* rds_count_r_sum, 
double* ctg_depth, loc_ht* global_ht, loc_ht_bool* global_ht_bool, int kmer_len, uint32_t max_mer_len_off, uint32_t *term_counts, int64_t num_walks, int64_t max_walk_len, 
int64_t sum_ext, int32_t max_read_size, int32_t max_read_count, uint32_t qual_offset, char* longest_walks, char* mer_walk_temp, uint32_t* final_walk_lens, int tot_ctgs)
{
    const int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if(idx < tot_ctgs){
    cstr_type loc_ctg;
    char *loc_r_reads, *loc_l_reads, *loc_r_quals, *loc_l_quals;
    uint32_t r_rds_cnt, l_rds_cnt, loc_rds_r_offset, loc_rds_l_offset;
    loc_ht* loc_mer_map = global_ht + idx * max_read_size * max_read_count;
    loc_ht_bool* loc_bool_map = global_ht_bool + idx * max_walk_len;
    double loc_ctg_depth = ctg_depth[idx];
    int64_t excess_reads;
    uint32_t max_ht_size = max_read_size * max_read_count;
    char* longest_walk_loc = longest_walks + idx * max_walk_len;
    int test = 629;



      int min_mer_len = LASSM_MIN_KMER_LEN;
      int max_mer_len = LASSM_MAX_KMER_LEN;
      

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
    max_mer_len = min(max_mer_len, loc_ctg.length);
    char* loc_mer_walk_temp = mer_walk_temp + idx * (max_walk_len + max_mer_len_off);

   // uint32_t mer_len = 21;
    // cstr_type ctg_mer(loc_ctg.start_ptr + (loc_ctg.length - mer_len), mer_len);
    // cstr_type loc_mer_walk(loc_mer_walk_temp, 0);
    // cstr_copy(loc_mer_walk, ctg_mer);
    // cstr_type longest_walk_thread(longest_walk_loc,0);
    // cstr_type walk(loc_mer_walk.start_ptr + mer_len, 0);

    cstr_type longest_walk_thread(longest_walk_loc,0);

    //main for loop
    int shift = 0;
    for(int mer_len = kmer_len; mer_len >= min_mer_len && mer_len <= max_mer_len; mer_len += shift){
            #ifdef DEBUG_PRINT_GPU
            if(idx == test){
               printf("GPU: shift:%d, mer_len:%d, min_mer_len:%d, idx:%d, max_mer_len:%d\n", shift, mer_len, min_mer_len, idx, max_mer_len);
               printf("contig:\n");
               print_mer(loc_ctg);
               }
            #endif
        
          //TODO: add a check if total number of reads exceeds a certain number/too large, skip that one, may be do this on cpu 
          // to preserve memory on GPU
          //TODO: need to reinitialize the hashtable after each kmer size is done

        if(r_rds_cnt != 0){    //if read count is zero, no need to count mers

            for(uint32_t k = 0; k < max_ht_size; k++){ // resetting hash table for next go
                loc_mer_map[k].key.length = EMPTY;
            }
            count_mers(loc_mer_map, loc_r_reads, max_ht_size, loc_r_quals, reads_r_offset, r_rds_cnt, rds_count_r_sum, loc_ctg_depth, mer_len, qual_offset, excess_reads, idx);

            cstr_type ctg_mer(loc_ctg.start_ptr + (loc_ctg.length - mer_len), mer_len);
            cstr_type loc_mer_walk(loc_mer_walk_temp, 0);
            cstr_copy(loc_mer_walk, ctg_mer);
            cstr_type walk(loc_mer_walk.start_ptr + mer_len, 0);

            #ifdef DEBUG_PRINT_GPU
            if(idx == test){
                printf("read_count:%d, idx:%d\n",r_rds_cnt, idx);
                printf("mer ctg len:%d mer_walk before:\n",loc_mer_walk.length);
                print_mer(loc_mer_walk);

               printf("ctg mer:\n");
               print_mer(ctg_mer);
            }
            #endif

            for(uint32_t k = 0; k < max_walk_len; k++){ // resetting bool map for next go
                loc_bool_map[k].key.length = EMPTY;
            }
            //TODO: initalize hash table find a faster way of doing this
            char walk_res = walk_mers(loc_mer_map, loc_bool_map, max_ht_size, mer_len, loc_mer_walk, longest_walk_thread, walk, idx, max_walk_len);
            #ifdef DEBUG_PRINT_GPU
            if(idx == test){
                printf("walk_res:%c, idx:%d\n",walk_res, idx);
                printf("GPU: walk.len:%d, longest.len:%d, idx:%d\n", walk.length, longest_walk_thread.length, idx);
            }
            #endif
            //int walk_len = walk.length
            if (walk.length > longest_walk_thread.length){ // this walk is longer than longest then copy it to longest walk
                cstr_copy(longest_walk_thread, walk);
            }
            if (walk_res == 'X') {
                atomicAdd(&term_counts[0], 1);
                // walk reaches a dead-end, downshift, unless we were upshifting
                if (shift == LASSM_SHIFT_SIZE) 
                    break;
                shift = -LASSM_SHIFT_SIZE;
            }else {
                if (walk_res == 'F') 
                    atomicAdd(&term_counts[1], 1);
                else 
                    atomicAdd(&term_counts[2], 1);
                // otherwise walk must end with a fork or repeat, so upshift
                if (shift == -LASSM_SHIFT_SIZE){
                    #ifdef DEBUG_PRINT_GPU
                    printf("breaking at shift neg:%d\n", shift);
                    #endif
                    break;
                    }
                if (mer_len > loc_ctg.length){
                    #ifdef DEBUG_PRINT_GPU
                    printf("breaking at mer_len too large:%d\n", mer_len);
                    #endif
                    break;
                }
                //if(walk_res == 'R') // if the walk ends at repeat, then the walk with repeat may have been longest so reset longest walk to zero
                  //  longest_walk_thread.length = 0;// TODO: DOUBLE CHECK WITH STEVE IF THIS IS CORRECT
                shift = LASSM_SHIFT_SIZE;
            }

        }else{
            break;
        }
    }
    if(longest_walk_thread.length > 0){
        final_walk_lens[idx] = longest_walk_thread.length;
        // printf("final longest walk len:%d/n", longest_walk_thread.length);
        // print_mer(longest_walk_thread);
       // atomicAdd(num_walks, 1);
     //   atomicAdd(sum_ext, longest_walk_thread.length);
    }else{
        final_walk_lens[idx] = 0;
    }

    #ifdef DEBUG_PRINT_GPU
    if(idx == test){
       // printf("walk:\n");
       // print_mer(walk);
       // printf("walk len:%d\n", walk.length);
       // printf("mer_walk after:\n");
       // print_mer(loc_mer_walk);
       // printf("mer_walk after, len:%d\n", loc_mer_walk.length);
        }
        //printf("walk result:%c\n", walk_res);
    #endif
}//end if to check if idx exceeds contigs
}