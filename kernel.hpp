
struct loc_ht{
    int key;
    int val;
};
__global__ void ht_kernel(loc_ht* ht, char* contigs, int* offset_sum);