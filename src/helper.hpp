
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>

#define CUDA_CHECK(ans)                                                                  \
    {                                                                                    \
        gpuAssert((ans), __FILE__, __LINE__);                                            \
    }
inline void
gpuAssert(cudaError_t code, const char* file, int line, bool abort = true)
{
    if(code != cudaSuccess)
    {
        fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if(abort)
            exit(code);
    }
}

struct ReadSeq {
  std::string read_id;
  std::string seq;
  std::string quals;
};

struct CtgWithReads {
  int32_t cid;
  std::string seq;
  double depth;
  std::vector<ReadSeq> reads_left;
  std::vector<ReadSeq> reads_right;
};




std::vector<std::string> read_fasta(std::string in_file, int &largest);
void read_locassm_data(std::vector<CtgWithReads> *data_in, std::string fname, 
int32_t& max_ctg_size, int32_t& total_r_reads, int32_t& total_l_reads, int32_t& max_read_size);
//templated functions needs to be in the same translation unit
template<typename T>
void print_log(T _log){
    std::cout<<_log<<std::endl;
}

template<typename T>
void print_vals(T val){
    print_log(val);
}

template<typename T, typename... Types>
void print_vals(T val, Types... val_){
    if(sizeof...(val_) == 0){
        print_vals(val);
    }else{
        print_vals(val);
        print_vals(val_...);
        }
}

int32_t print_loc_data(std::vector<CtgWithReads> *data_in);

