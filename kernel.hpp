#define EMPTY 0xFFFFFFFF
#define HT_SIZE 300
struct cstr_type{
    char* start_ptr;
    int length;

    __device__ bool operator==(const cstr_type& in2){
        return ((start_ptr == in2.start_ptr) && (length == in2.length));
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
__device__ int hash_func(cstr_type key);