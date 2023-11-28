#ifndef _MM_BUILD_
#define _MM_BUILD_

#include "loader.h"
#include "util_func.h"
#include "minimizer.h"

#include <boost/filesystem.hpp>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <map>
#include <vector>
#include <tuple>
#include <string>
#include <list>
#include <sys/stat.h>

#define MAX_BLOCK 128



class MMBuild {
 public:
  MMBuild(void);
  MMBuild(BioAlphabet& alpha, bool& reduced_alpha);

  ~MMBuild(void);
  void LoadSequences(std::string& seq_file, const bool &rev_comp);
  void BuildIndex(int& kmer, int& window, std::string &dir, std::string &file_stem);
  void MapQuery(int& kmer, int& window, int& id);
  void MapQuery(MMBuild& query_seq, int& kmer, int& window, int& id);
  void PrintAllSeqsInfo(void);
  void DumpIndextoFile(std::string& dir, std::string& file_stem);
  void LoadIndex(std::string& dir, std::string& file_stem);
  void LoadMultiIndex(std::string& dir, std::string& file_stem);
  void SetBlockConfig(const int& num_blocks, std::string &dir, std::string &file_stem);
  void SplitSequence(std::vector<MMBuild>& block_mm);
  void PrintOverlaps();
  void WritePAF(std::ofstream& outfile);
  
  int GetNumSeqs()  {
    assert(is_sequence_loaded_);
    return num_seqs_;
  }

  int GetAlphabet() {
    assert(is_sequence_loaded_);
    return alpha_;
  }

 protected:
  int num_seqs_;
  double db_size_MB_;
  BioAlphabet alphabet_;
  int alpha_;
  char** header_;
  char** sequence_;
  int* seq_len_;
  uint32_t* seq_id_;
  std::vector<FastaInfo> seq_entry;
  MMSketch* mm_index_{};
  std::vector<MMSketch*> block_mm_;
  std::vector<int> block_size_;
  bool is_header_loaded_, is_sequence_loaded_, is_block_loaded_;
  bool is_mm_built_, is_k_array_built_;
  bool is_mm_loaded_;
  bool is_size_counted_;
  bool is_multi_;
  bool is_rev_comp_;
  int num_blocks_;
  bool is_contained_init_;
  bool is_alphabet_set_;
  std::vector<bool> is_contained_;
  
};

#endif
