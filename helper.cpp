#include "helper.hpp"

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