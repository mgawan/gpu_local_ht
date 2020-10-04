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

void read_locassm_data(std::vector<CtgWithReads> *data_in, std::string fname){
    std::ifstream f(fname);
    std::string line;
    //ctgs_map->clear();
    while(getline(f, line)) {
      std::stringstream ss(line);
      CtgWithReads temp_in;
      int lsize = 0, rsize = 0;
      //CtgWithReads ctg_with_reads;
      ss >> temp_in.cid >> temp_in.seq >> temp_in.depth >> lsize >> rsize;
      for (int i = 0; i < lsize; i++) {
        ReadSeq read_seq;
        ss >> read_seq.read_id >> read_seq.seq >> read_seq.quals;
        temp_in.reads_left.push_back(read_seq);
      }
      for (int i = 0; i < rsize; i++) {
        ReadSeq read_seq;
        ss >> read_seq.read_id >> read_seq.seq >> read_seq.quals;
        temp_in.reads_right.push_back(read_seq);
      }
      //ctgs_map->insert({ctg_with_reads.cid, ctg_with_reads});
      data_in->push_back(temp_in);
    }
  }

  void print_loc_data(std::vector<CtgWithReads> *data_in){
        for(int i = 0; i < 7; i++){
        print_vals("contig_id: ", (*data_in)[i].cid, "\n seq: ", (*data_in)[i].seq, "\n depth: ", (*data_in)[i].depth, "\n right: ", (*data_in)[i].reads_left.size(), "\n left: ",(*data_in)[i].reads_right.size());
        print_vals("**READS**");
        for(int j = 0; j< (*data_in)[i].reads_left.size(); j++){
            ReadSeq read = (*data_in)[i].reads_left[j];
            print_vals(read.read_id, " ", read.seq, " ", read.quals);
        }
        for(int j = 0; j< (*data_in)[i].reads_right.size(); j++){
            ReadSeq read = (*data_in)[i].reads_right[j];
            print_vals(read.read_id, " ", read.seq, " ", read.quals);
        }

    }
  }



