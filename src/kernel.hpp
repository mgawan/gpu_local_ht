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

struct ExtCounts {
  uint16_t count_A;
  uint16_t count_C;
  uint16_t count_G;
  uint16_t count_T;

//   __device__  
//   void  cu_pr_comp(cu_pair& elem1, cu_pair& elem2){
//     if(elem1.y == elem2.y)
//         return elem1.x > elem2.x;
//     else
//         return elem1.y > elem1.y;
//   }
//   __device__
//   void sort(cu_pair (&counts)[4]){
//     for(int i = 0; i < 4; i++){
//         for(int j = 0; j < 4; j++){
//             if(cu_pr_comp(counts[i], counts[j])){
//                 cu_pr_comp temp = counts[i];
//                 counts[i] = counts[j];
//                 counts[j] = temp;
//             }
//         }
//     }
//   }

//   __device__ 
//   void get_sorted(cu_pair (&counts)[4]) {
//     cu_pair temp_pair;
//     temp_pair.x = 'A';
//     temp_pair.y = (int)count_A;
//     counts[0] = temp_pair;

//     temp_pair.x = 'C';
//     temp_pair.y = (int)count_C;
//     counts[1] = temp_pair;

//     temp_pair.x = 'G';
//     temp_pair.y = (int)count_G;
//     counts[2] = temp_pair;

//     temp_pair.x = 'T';
//     temp_pair.y = (int)count_T;
//     counts[3] = temp_pair;

//     sort(counts);// sorts in place on the array, passed by reference
//   }

//   __device__ 
//   bool is_zero() {
//     if (count_A + count_C + count_G + count_T == 0) return true;
//     return false;
//   }

//TODO: replace numeric_limits by either a suitable constant/predefined value or find a device alternate
  __device__
  void inc(char ext, int count) {
    switch (ext) {
      case 'A':
        count += count_A;
        count_A = (count < numeric_limits<uint16_t>::max()) ? count : numeric_limits<uint16_t>::max();
        break;
      case 'C':
        count += count_C;
        count_C = (count < numeric_limits<uint16_t>::max()) ? count : numeric_limits<uint16_t>::max();
        break;
      case 'G':
        count += count_G;
        count_G = (count < numeric_limits<uint16_t>::max()) ? count : numeric_limits<uint16_t>::max();
        break;
      case 'T':
        count += count_T;
        count_T = (count < numeric_limits<uint16_t>::max()) ? count : numeric_limits<uint16_t>::max();
        break;
    }
  }
};

struct MerFreqs {
  // how many times this kmer has occurred: don't need to count beyond 65536
  // count of high quality extensions and low quality extensions - structure comes from kmer_dht.hpp
  ExtCounts hi_q_exts, low_q_exts;
  // the final extensions chosen - A,C,G,T, or F,X
  char ext;
  // the count of the final extension
  int count;

  struct MerBase {
    char base;
    uint16_t nvotes_hi_q, nvotes, rating;

    __device__
    uint16_t get_base_rating(int depth) {
      double min_viable = max(LASSM_MIN_VIABLE_DEPTH * depth, 2.0);
      double min_expected_depth = max(LASSM_MIN_EXPECTED_DEPTH * depth, 2.0);
      if (nvotes == 0) return 0;
      if (nvotes == 1) return 1;
      if (nvotes < min_viable) return 2;
      if (min_expected_depth > nvotes && nvotes >= min_viable && nvotes_hi_q < min_viable) return 3;
      if (min_expected_depth > nvotes && nvotes >= min_viable && nvotes_hi_q >= min_viable) return 4;
      if (nvotes >= min_expected_depth && nvotes_hi_q < min_viable) return 5;
      if (nvotes >= min_expected_depth && min_viable < nvotes_hi_q && nvotes_hi_q < min_expected_depth) return 6;
      return 7;
    }
  };

  void set_ext(int seq_depth) {
    // set extension similarly to how it is done with localassm in mhm
    MerBase mer_bases[4] = {{.base = 'A', .nvotes_hi_q = hi_q_exts.count_A, .nvotes = low_q_exts.count_A},
                            {.base = 'C', .nvotes_hi_q = hi_q_exts.count_C, .nvotes = low_q_exts.count_C},
                            {.base = 'G', .nvotes_hi_q = hi_q_exts.count_G, .nvotes = low_q_exts.count_G},
                            {.base = 'T', .nvotes_hi_q = hi_q_exts.count_T, .nvotes = low_q_exts.count_T}};
    for (int i = 0; i < 4; i++) {
      mer_bases[i].rating = mer_bases[i].get_base_rating(seq_depth);
    }
    // sort bases in descending order of quality
    sort(mer_bases, mer_bases + sizeof(mer_bases) / sizeof(mer_bases[0]),
         [](const auto &elem1, const auto &elem2) -> bool {
           if (elem1.rating != elem2.rating) return elem1.rating > elem2.rating;
           if (elem1.nvotes_hi_q != elem2.nvotes_hi_q) return elem1.nvotes_hi_q > elem2.nvotes_hi_q;
           if (elem1.nvotes != elem2.nvotes) return elem1.nvotes > elem2.nvotes;
           return true;
         });
    int top_rating = mer_bases[0].rating;
    int runner_up_rating = mer_bases[1].rating;
    if (top_rating < runner_up_rating) DIE("top_rating ", top_rating, " < ", runner_up_rating, "\n");
    assert(top_rating >= runner_up_rating);
    int top_rated_base = mer_bases[0].base;
    ext = 'X';
    count = 0;
    // no extension (base = 0) if the runner up is close to the top rating
    // except, if rating is 7 (best quality), then all bases of rating 7 are forks
    if (top_rating > LASSM_RATING_THRES) {         // must have at least minViable bases
      if (top_rating <= 3) {    // must be uncontested
        if (runner_up_rating == 0) ext = top_rated_base;
      } else if (top_rating < 6) {
        if (runner_up_rating < 3) ext = top_rated_base;
      } else if (top_rating == 6) {  // viable and fair hiQ support
        if (runner_up_rating < 4) ext = top_rated_base;
      } else {                     // strongest rating trumps
        if (runner_up_rating < 7) {
          ext = top_rated_base;
        } else {
          if (mer_bases[2].rating == 7 || mer_bases[0].nvotes == mer_bases[1].nvotes) ext = 'F';
          else if (mer_bases[0].nvotes > mer_bases[1].nvotes) ext = mer_bases[0].base;
          else if (mer_bases[1].nvotes > mer_bases[0].nvotes) ext = mer_bases[1].base;
        }
      }
    }
    for (int i = 0; i < 4; i++) {
      if (mer_bases[i].base == ext) {
        count = mer_bases[i].nvotes;
        break;
      }
    }
  }

};

__global__ void ht_kernel(loc_ht* ht, char* contigs, int* offset_sum, int kmer_size);
__device__ void ht_insert(loc_ht* thread_ht, cstr_type kmer_key, cstr_type ctg_val);
__device__ void ht_delete(loc_ht* thread_ht, cstr_type kmer_key);
__device__ cstr_type ht_get(loc_ht* thread_ht, cstr_type kmer_key);
__device__ unsigned hash_func(cstr_type key);
__global__ void iterative_walks_kernel(uint32_t* cid, char *contigs, uint32_t* ctg_depth, char* reads_seqs, int max_mer_len, int kmer_len,
 int qual_offset, int walk_len_limit, int64_t *term_counts, int64_t num_walks, int64_t max_walk_len, int64_t sum_ext, int64_t excess_reads);