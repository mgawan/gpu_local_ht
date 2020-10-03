#include "helper.hpp"

struct CtgWithReads {
  int32_t cid;
  std::string seq;
  double depth;
  std::vector<ReadSeq> reads_left;
  std::vector<ReadSeq> reads_right;
};

struct ReadSeq {
  std::string read_id;
  std::string seq;
  std::string quals;
};


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

void read_locassm_data(vector<CtgWithReads> &data_in, string fname){
    std::ifstream f(fname);
    std::string line;
    ctgs_map->clear();
    while(getline(f, line)) {
      std::stringstream ss(line);
      int lsize = 0, rsize = 0;
      //CtgWithReads ctg_with_reads;
      ss >> data_in.cid >> data_in.seq >> data_in.depth >> lsize >> rsize;
      for (int i = 0; i < lsize; i++) {
        ReadSeq read_seq;
        ss >> read_seq.read_id >> read_seq.seq >> read_seq.quals;
        data_in.reads_left.push_back(read_seq);
      }
      for (int i = 0; i < rsize; i++) {
        ReadSeq read_seq;
        ss >> read_seq.read_id >> read_seq.seq >> read_seq.quals;
        data_in.reads_right.push_back(read_seq);
      }
      //ctgs_map->insert({ctg_with_reads.cid, ctg_with_reads});
    }
  }



