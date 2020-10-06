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

struct cu_pair{
    char x;
    int y;
}

struct ExtCounts {
  uint16_t count_A;
  uint16_t count_C;
  uint16_t count_G;
  uint16_t count_T;

  __device__  
  cu_pr_comp(cu_pair& elem1, cu_pair& elem2){
    if(elem1.y == elem2.y)
        return elem1.x > elem2.x;
    else
        return elem1.y > elem1.y;
  }
  __device__
  sort(cu_pair (&counts)[4]){
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            if(cu_pr_comp(&counts[i], &counts[j])){
                
            }
        }
    }
  }

  get_sorted(cu_pair (&counts)[4]) {
    cu_pair temp_pair;
    temp_pair.x = 'A';
    temp_pair.y = (int)count_A;
    counts[0] = temp_pair;

    temp_pair.x = 'C';
    temp_pair.y = (int)count_C;
    counts[1] = temp_pair;

    temp_pair.x = 'G';
    temp_pair.y = (int)count_G;
    counts[2] = temp_pair;

    temp_pair.x = 'T';
    temp_pair.y = (int)count_T;
    counts[3] = temp_pair;

    sort(std::begin(counts), std::end(counts),
         [](const auto &elem1, const auto &elem2) {
           if (elem1.second == elem2.second) return elem1.first > elem2.first;
           else return elem1.second > elem2.second;
         });
    return counts;
  }

  bool is_zero() {
    if (count_A + count_C + count_G + count_T == 0) return true;
    return false;
  }

  void inc(char ext, int count) {
    switch (ext) {
      case 'A':
        count += count_A;
        count_A = (count < numeric_limits<ext_count_t>::max()) ? count : numeric_limits<ext_count_t>::max();
        break;
      case 'C':
        count += count_C;
        count_C = (count < numeric_limits<ext_count_t>::max()) ? count : numeric_limits<ext_count_t>::max();
        break;
      case 'G':
        count += count_G;
        count_G = (count < numeric_limits<ext_count_t>::max()) ? count : numeric_limits<ext_count_t>::max();
        break;
      case 'T':
        count += count_T;
        count_T = (count < numeric_limits<ext_count_t>::max()) ? count : numeric_limits<ext_count_t>::max();
        break;
    }
  }

  char get_ext(uint16_t count) {
    auto sorted_counts = get_sorted();
    int top_count = sorted_counts[0].second;
    int runner_up_count = sorted_counts[1].second;
    // set dynamic_min_depth to 1.0 for single depth data (non-metagenomes)
    int dmin_dyn = max((int)((1.0 - _dynamic_min_depth) * count), _dmin_thres);
    if (top_count < dmin_dyn) return 'X';
    if (runner_up_count >= dmin_dyn) return 'F';
    return sorted_counts[0].first;
    /*
    // FIXME: this is not very helpful. With qual_cutoff = 20) it increases ctgy & coverage a little bit, but at a cost
    // of increased msa. We really need to try both low q (qual cutoff 10) and hi q, as we do with localassm.
    double dmin_dyn = max(2.0, LASSM_MIN_EXPECTED_DEPTH * count);
    if ((top_count < dmin_dyn && runner_up_count > 0) || (top_count >= dmin_dyn && runner_up_count >= dmin_dyn)) return 'F';
    return sorted_counts[0].first;
    */
  }

};

__global__ void ht_kernel(loc_ht* ht, char* contigs, int* offset_sum, int kmer_size);
__device__ void ht_insert(loc_ht* thread_ht, cstr_type kmer_key, cstr_type ctg_val);
__device__ void ht_delete(loc_ht* thread_ht, cstr_type kmer_key);
__device__ cstr_type ht_get(loc_ht* thread_ht, cstr_type kmer_key);
__device__ unsigned hash_func(cstr_type key);
__global__ void iterative_walks_kernel(uint32_t* cid, char *contigs, uint32_t* ctg_depth, char* reads_seqs, int max_mer_len, int kmer_len,
 int qual_offset, int walk_len_limit, int64_t *term_counts, int64_t num_walks, int64_t max_walk_len, int64_t sum_ext, int64_t excess_reads);