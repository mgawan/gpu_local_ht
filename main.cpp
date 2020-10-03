#include <unordered_map>
#include <iostream>
#include <typeinfo>
#include <random>
#include <algorithm>
#include <cstring>
#include "helper.hpp"
#include "kernel.hpp"
#define TOT_THREADS 1
#define KMER_SZ 4


int main (int argc, char* argv[]){
    std::string in_file = argv[1];
    int max_ctg_size;
    std::vector<std::string> contigs = read_fasta(in_file, max_ctg_size);
    char *d_contigs, *h_contigs;
    int *h_offset_arr, *d_offset_arr;
    loc_ht *d_ht;
    int offset_sum = 0;
    //allocate memory for host and device char arrays
    h_contigs = new char[max_ctg_size*contigs.size()];
    CUDA_CHECK(cudaMalloc(&d_contigs, sizeof(char)*max_ctg_size*contigs.size()));
    //memory for offset array
    h_offset_arr  = new int[contigs.size()];
    CUDA_CHECK(cudaMalloc(&d_offset_arr, sizeof(int)*contigs.size()));
    //memory allocation on device for local hashtables
    CUDA_CHECK(cudaMalloc(&d_ht, sizeof(loc_ht)*HT_SIZE*TOT_THREADS));


    //convert strings to char* and prepare offset array
    for(int i = 0; i < contigs.size(); i++)
    {
        char* seq_ptr = h_contigs + offset_sum;
        memcpy(seq_ptr, contigs[i].c_str(), contigs[i].size());
        offset_sum += contigs[i].size();
        h_offset_arr[i] = offset_sum;
    }
    //moving data to device
    CUDA_CHECK(cudaMemcpy(d_contigs, h_contigs, sizeof(char)*max_ctg_size*contigs.size(), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_offset_arr, h_offset_arr, sizeof(int)*contigs.size(), cudaMemcpyHostToDevice));
    //print_vals("launch kernel", "another print");
    ht_kernel<<<1,1>>>(d_ht, d_contigs, d_offset_arr, KMER_SZ);

    CUDA_CHECK(cudaFree(d_ht));

    std::unordered_map<std::string, int> my_map;

    my_map.insert({"this_map",10});
    my_map.insert({"map",100});
    my_map.insert({"this_map", 15});
    print_vals("total_insertions:", "count");
   // print_vals(my_map["this_map"], my_map["map"]);

    return 0;
}
