#include <iostream>
#include <unordered_map>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <typeinfo>
#include <random>
#include <algorithm>
#include <cstring>

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


std::vector<std::string> read_fasta(std::string in_file, int &largest);
template<typename T>
void print_log(T _log);

template<typename T>
void print_vals(T val);

template<typename T, typename... Types>
void print_vals(T val, Types... val_);

__global__ void ht_kernel(int x, int y){
    return;
}

int main (int argc, char* argv[]){
    std::string in_file = argv[1];
    int max_ctg_size;
    std::vector<std::string> contigs = read_fasta(in_file, max_ctg_size);
    char *d_contigs, *h_contigs;
    int *h_offset_arr, *d_offset_arr;
    int offset_sum = 0;
    //allocate memory for host and device char arrays
    h_contigs = new char[max_ctg_size*contigs.size()];
    CUDA_CHECK(cudaMalloc(&d_contigs, sizeof(char)*max_ctg_size*contigs.size()));
    //memory for offset array
    h_offset_arr  = new int[contigs.size()];
    CUDA_CHECK(cudaMalloc(&d_offset_arr, sizeof(char)*contigs.size()));

    //convert strings to char* and prepare offset array
    for(int i = 0; i < contigs.size(); i++)
    {
        char* seq_ptr = h_contigs + offset_sum;
        memcpy(seq_ptr, contigs[i].c_str(), contigs[i].size());
        offset_sum += contigs[i].size();
        h_offset_arr[i] = offset_sum;
    }
    //moving data to device
    CUDA_CHECK(cudaMemcpy(d_contigs, h_contigs, sizeof(char)*offset_sum, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_offset_arr, h_offset_arr, sizeof(int)*contigs.size(), cudaMemcpyHostToDevice));

    ht_kernel<<<1,1>>>(1,1);

    return 0;
}

std::vector<std::string> read_fasta(std::string in_file, int &largest){
    std::ifstream in_stream(in_file);
    std::string in_line;
    std::vector<std::string> in_sequences;
    int max_size = 0;

    if(in_stream.is_open())
    {
      while(getline(in_stream, in_line))
      {
          if(in_line[0] == '>'){
            continue;
          }else{
            std::string seq = in_line;
            in_sequences.push_back(seq);
            if(max_size < seq.length())
                max_size = seq.length();
          }
      }
      in_stream.close();
    }
    largest = max_size;
    return in_sequences;
}

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
