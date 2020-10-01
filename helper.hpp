
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

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

