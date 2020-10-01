#include "kernel.hpp"

__global__ void ht_kernel(loc_ht* ht, char* contigs, int* offset_sum, int kmer_size){
    int idx = blockIdx.x*gridDim.x + threadIdx.x;
    int total_threads = blockDim.x * gridDim.x;
    char* loc_contig;
    loc_ht* thread_ht;
    thread_ht = ht+ total_threads*idx;
    if(idx == 0){
        loc_contig = contigs;
    }else{
        loc_contig = contigs + offset_sum[idx];
    }

    for(int i = 0; i < )
}