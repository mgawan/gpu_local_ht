#include<stdio.h>
#define EMPTY 0xFFFFFFFF
#define HT_SIZE 300
struct cstr_type{
    char* start_ptr;
    int length;
    __device__ cstr_type(){}
    __device__ cstr_type(char* ptr, int len){
        start_ptr = ptr;
        length = len;
    }

    __device__ bool operator==(const cstr_type& in2){
        bool str_eq = true;
        if(length != EMPTY && in2.length != EMPTY)
            for(int i = 0; i < in2.length; i++){
                if(start_ptr[i] != in2.start_ptr[i]){
                    str_eq = false;
                    break;
                }
            }
        return (str_eq && (length == in2.length));
    }
};

struct loc_ht{
    cstr_type key;
    cstr_type val;
};
__global__ void ht_kernel(loc_ht* ht, char* contigs, int* offset_sum, int kmer_size);
__device__ void ht_insert(loc_ht* thread_ht, cstr_type kmer_key, cstr_type ctg_val);
__device__ void ht_delete(loc_ht* thread_ht, cstr_type kmer_key);
__device__ cstr_type ht_get(loc_ht* thread_ht, cstr_type kmer_key);
__device__ unsigned hash_func(cstr_type key);
__global__ void iterative_walks_kernel(uint32_t* cid, char *contigs, uint32_t* ctg_depth, char* reads_seqs, int max_mer_len, int kmer_len,
                                      int qual_offset, int walk_len_limit, int64_t *term_counts,
                                      int64_t num_walks, int64_t max_walk_len, int64_t sum_ext, int64_t excess_reads);