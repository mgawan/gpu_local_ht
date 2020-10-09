#include "kernel.hpp"


//TODO: all the hashtable entries need to be set to empty, figure that out
__device__ print_mer(cstr_type& mer){
    for(int i = 0; i < mer.length; i++){
        printf("%c",mer.start_ptr[i]);
    }
    printf("\n");
}
__global__ void ht_kernel(loc_ht* ht, char* contigs, int* offset_sum, int kmer_size){
    int idx = blockIdx.x*gridDim.x + threadIdx.x;
    cstr_type loc_contig;
    loc_ht* thread_ht;
    thread_ht = ht + idx*HT_SIZE;

    for(int i = 0; i < HT_SIZE; i++){
        thread_ht[i].key.length = EMPTY;
    }

    if(idx == 0){
        loc_contig.start_ptr = contigs;
        loc_contig.length = offset_sum[idx];
    }else{
        loc_contig.start_ptr = contigs + offset_sum[idx];
        loc_contig.length = offset_sum[idx] - offset_sum[idx-1];
    }
    
      cstr_type kmer;
      kmer.start_ptr = loc_contig.start_ptr;
      kmer.length = kmer_size;
      printf("kmer size:%d\n",kmer.length);
      printf("contig size:%d\n",loc_contig.length);
    for(int i = 0; i < (loc_contig.length - (kmer_size-1)); i++){
        char *mystr = kmer.start_ptr;
        for(int j = 0; j < kmer.length; j++){
            printf("%c",mystr[j]);
        }
        printf("\n");
        cstr_type temp = ht_get(thread_ht, kmer);
        if(temp.length == -1)
            ht_insert(thread_ht, kmer, loc_contig);
        else{
            printf("key exists:\n");
            for(int j = 0; j < temp.length; j++){
                printf("%c",temp.start_ptr[j]);
        }
        printf("\n");
        }

        kmer.start_ptr =  kmer.start_ptr + 1;
    }


}

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

__device__ void ht_insert(loc_ht* thread_ht, cstr_type kmer_key, MerFreq mer_val, uint32_t max_size){
    unsigned hash_val = hash_func(kmer_key);
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
ht_loc& ht_get(loc_ht* thread_ht, cstr_type kmer_key, uint32_t max_size){
    unsigned hash_val = hash_func(kmer_key);
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
            return ht_loc(cstr_type(NULL,-1), MerFreqs());
        }
    }

}



__device__ 
void count_mers(ht_loc* thrd_loc_ht, char* loc_r_reads, uint32_t max_ht_size, char* loc_r_quals, uint32_t* reads_r_offset, uint32_t& r_rds_cnt, 
uint32_t* rds_count_r_sum, uint32_t& loc_ctg_depth, uint32_t& mer_len, uint32_t& qual_offset, uint32_t& excess_reads){
    cstr_type read;
    cstr_type qual;
    uint32_t running_sum_len = 0;
    for(int i = 0; i < r_rds_cnt; ++){
        read.start_ptr = loc_r_reads + running_sum_len;
        qual.start_ptr = loc_r_quals + running_sum_len;
        if(i == 0){
            if(threadIdx.x == 0){
                read.length = reads_r_offset[(rds_count_r_sum[threadIdx.x] - r_rds_cnt) + i];
                qual.length = reads_r_offset[(rds_count_r_sum[threadIdx.x] - r_rds_cnt) + i];
                }
            else{   
                read.length = reads_r_offset[(rds_count_r_sum[threadIdx.x] - r_rds_cnt) + i] - reads_r_offset[(rds_count_r_sum[threadIdx.x - 1] -1)];
                qual.length = reads_r_offset[(rds_count_r_sum[threadIdx.x] - r_rds_cnt) + i] - reads_r_offset[(rds_count_r_sum[threadIdx.x - 1] -1)];
                }
            }
        else{
            read.length = reads_r_offset[(rds_count_r_sum[threadIdx.x] - r_rds_cnt) + i] - reads_r_offset[(rds_count_r_sum[threadIdx.x] - r_rds_cnt) + (i-1)];
            qual.length = reads_r_offset[(rds_count_r_sum[threadIdx.x] - r_rds_cnt) + i] - reads_r_offset[(rds_count_r_sum[threadIdx.x] - r_rds_cnt) + (i-1)];
            }

        if (mer_len > read.length) // skip the read which is smaller than merlen
            continue;
        int num_mers = read.length - mer_len;
        cstr_type mer(read.cstr_type, mer_len)
        for( int start = 0; start < num_mers; start++){
            print_mer(mer);
            //TODO: on cpu side add a check that if a certain read contains 'N', that is not included, check this with steve, 
            // because searching a single mer for an N is going to be too slow
            ht_loc &temp_Mer = ht_get(thrd_loc_ht, mer, max_ht_size);
            if(temp_Mer.key.length == EMPTY){
                temp_Mer.key = mer;
                temp_Mer.val = {.hi_q_exts = {0}, .low_q_exts = {0}, .ext = 0, .count = 0}; // TODO: verify that this constructor works on GPU
            }
            int ext_pos = start + mer_len;
            assert(ext_pos < (int)read.length); // TODO: verify that assert works on gpu
            char ext = read[ext_pos];
            if (ext == 'N') continue; // TODO: why the redundant check?
            int qual_diff = qual[ext_pos] - qual_offset;
            if (qual_diff >= LASSM_MIN_QUAL) temp_Mer.val.low_q_exts.inc(ext, 1);
            if (qual_diff >= LASSM_MIN_HI_QUAL) temp_Mer.val.hi_q_exts.inc(ext, 1);
            
            mer.start_ptr = mer.start_ptr + 1;
        }
       running_sum_len += read.length; // right before the for loop ends, update the prev_len to offset next read correctly
    }

    //setting extension by traversing the completed table
    // TODO: think of a better way to do this
    for (int k = 0; k < max_table_size; k++) {
        if( thrd_loc_ht[k].key.length != EMPTY)
            thrd_loc_ht[k].val.set_ext(loc_ctg_depth);
    }
}

//same kernel will be used for right and left walks
__global__ void iterative_walks_kernel(uint32_t* cid, uint32_t* ctg_offsets, char* contigs, 
char* reads_l, char* reads_r, char* quals_r, char* quals_l, uint32_t* reads_l_offset, uint32_t* reads_r_offset, uint32_t* rds_count_l_sum, uint32_t* rds_count_r_sum, uint32_t* ctg_depth, ht_loc* global_ht,
int max_mer_len, int kmer_len, int walk_len_limit, int64_t *term_counts, int64_t num_walks, int64_t max_walk_len, int64_t sum_ext, int32_t max_read_size, int32_t max_read_count){
    unsigned idx = threadIdx.x + blockIdx.x * blockDimx.x;
    cstr_type loc_ctg;
    char *loc_r_reads, *loc_l_reads, *loc_r_quals, *loc_l_quals;
    uint32_t r_rds_cnt, l_rds_cnt, loc_rds_r_offset, loc_rds_l_offset;
    loc_ht* loc_mer_map = global_ht + idx * max_read_size * max_read_count;
    uint32_t loc_ctg_depth = ctg_depth[idx];
    int64_t excess_reads;
    uint32_t qual_offset, max_ht_size = max_read_size * max_read_count;

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
        if (reads_count_r_sum[idx - 1] == 0)
            loc_r_reads = reads_r;
        else
            loc_r_reads = reads_r + reads_r_offset[reads_count_r_sum[idx - 1] - 1]; // you want to start from where previous contigs, last read ends.

        if (reads_count_l_sum[idx - 1] == 0)
            loc_l_reads = reads_l;
        else
            loc_l_reads = reads_l + reads_l_offset[reads_count_l_sum[idx - 1] - 1]; // you want to start from where previous contigs, last read ends.
        
        if (reads_count_r_sum[idx - 1] == 0)
            loc_r_quals = quals_r;
        else
            loc_r_quals = quals_r + reads_r_offset[reads_count_r_sum[idx - 1] - 1]; // you want to start from where previous contigs, last read ends.

        if (reads_count_l_sum[idx - 1] == 0)
            loc_l_quals = quals_l;
        else
            loc_l_quals = quals_l + reads_l_offset[reads_count_l_sum[idx - 1] - 1]; // you want to start from where previous contigs, last read ends. 
    }

    //main for loop
    //TODO: commenting out the main for loop for testing count_mers
    //for(int mer_len = kmer_len; mer_len >= min_mer_len && mer_len <= max_mer_len; mer_len += shift){
          //TODO: add a check if total number of reads exceeds a certain number/too large, skip that one, may be do this on cpu 
          // to preserve memory on GPU
          int mer_len = 21;
        if(r_rds_cnt != 0)    //if count is zero, no need to count
            count_mers(loc_mer_map, loc_r_reads, max_ht_size, loc_r_quals, reads_r_offset, r_rds_cnt, rds_count_r_sum, loc_ctg_depth, 
            mer_len, qual_offset, excess_reads);
   // }

}