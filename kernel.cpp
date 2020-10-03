#include "kernel.hpp"


//TODO: all the hashtable entries need to be set to empty, figure that out
__global__ void ht_kernel(loc_ht* ht, char* contigs, int* offset_sum, int kmer_size){
    int idx = blockIdx.x*gridDim.x + threadIdx.x;
    cstr_type loc_contig;
    loc_ht* thread_ht;
    thread_ht = ht + idx*HT_SIZE;

    for(int i = 0; i < HT_SIZE; i++){
        thread_ht[i].key.length = EMPTY;
    }

    if(idx == 0){
        loc_contig.start_ptr = contigs;
        loc_contig.length = offset_sum[idx];
    }else{
        loc_contig.start_ptr = contigs + offset_sum[idx];
        loc_contig.length = offset_sum[idx] - offset_sum[idx-1];
    }
    
      cstr_type kmer;
      kmer.start_ptr = loc_contig.start_ptr;
      kmer.length = kmer_size;
      printf("kmer size:%d\n",kmer.length);
      printf("contig size:%d\n",loc_contig.length);
    for(int i = 0; i < (loc_contig.length - (kmer_size-1)); i++){
        char *mystr = kmer.start_ptr;
        for(int j = 0; j < kmer.length; j++){
            printf("%c",mystr[j]);
        }
        printf("\n");
        cstr_type temp = ht_get(thread_ht, kmer);
        if(temp.length == -1)
            ht_insert(thread_ht, kmer, loc_contig);
        else{
            printf("key exists:\n");
            for(int j = 0; j < temp.length; j++){
                printf("%c",temp.start_ptr[j]);
        }
        printf("\n");
        }

        kmer.start_ptr =  kmer.start_ptr + 1;
    }


}

//TODO: make sure that the returned hash value is wthin the range of HT size
__device__ unsigned hash_func(cstr_type key){
    unsigned hash, i;
    for(hash = i = 0; i < key.length; ++i)
    {
        hash += key.start_ptr[i];
        hash += (hash << 10);
        hash ^= (hash >> 6);
    }
    hash += (hash << 3);
    hash ^= (hash >> 11);
    hash += (hash << 15);
    return hash%HT_SIZE;//(hash & (HT_SIZE - 1));
}

__device__ void ht_insert(loc_ht* thread_ht, cstr_type kmer_key, cstr_type ctg_val){
    unsigned hash_val = hash_func(kmer_key);
    unsigned orig_hash = hash_val;
    int count = 0;
    while(true){
        int if_empty = thread_ht[hash_val].key.length; // length is set to some unimaginable number to indicate if its empty
        if(if_empty == EMPTY || (thread_ht[hash_val].key == kmer_key)){

            printf("hash_val:%d, orig_hash:%d, attemp:%d\n",hash_val, orig_hash, count);
            thread_ht[hash_val].key = kmer_key;
            thread_ht[hash_val].val = ctg_val;

            return;
        }
        hash_val = (hash_val +1 ) %HT_SIZE;//(hash_val + 1) & (HT_SIZE-1);
        count++;

    }
}

__device__ void ht_delete(loc_ht* thread_ht, cstr_type kmer_key){
    int hash_val = hash_func(kmer_key);
    while(true){
        if(thread_ht[hash_val].key == kmer_key){
            thread_ht[hash_val].key.length = EMPTY;
            return;
        }
        if(thread_ht[hash_val].key.length == EMPTY){
            return;
        }
        hash_val = (hash_val + 1) & (HT_SIZE -1);
    }
}

__device__ cstr_type ht_get(loc_ht* thread_ht, cstr_type kmer_key){
    unsigned hash_val = hash_func(kmer_key);
    unsigned orig_hash = hash_val;
    
    while(true){
        if(thread_ht[hash_val].key == kmer_key){
            printf("key found, returning\n");
            return thread_ht[hash_val].val;
        }
        // if(thread_ht[hash_val].key.length == EMPTY){

        //     printf("empty key found, returning\n");
        //     return thread_ht[hash_val].val;
        // }

        hash_val = (hash_val +1 ) %HT_SIZE;//hash_val = (hash_val + 1) & (HT_SIZE -1);
        if(hash_val == orig_hash){
            printf("end reached\n");
            return cstr_type(NULL, -1);
            
            }
    }

}
