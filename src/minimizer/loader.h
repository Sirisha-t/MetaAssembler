#ifndef _LOADER_H_
#define _LOADER_H_

#include <iostream>
#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstring>
#include <cassert>
#include <list>
#include <unordered_map>
#include <ctype.h>

#include "bio_alphabet.h"

typedef uint64_t RIDTYPE;
typedef uint64_t POSTYPE;

#ifndef MAXSEQLEN
#define MAXSEQLEN 1000000
#endif

struct FastaInfo{
    std::string header_;
    std::string sequence_;
    int seq_len_;
    uint32_t rid_;
    FastaInfo(const std::string& header , const std::string& sequence, const int& seq_len, const uint32_t& rid){
        header_ = header;
        sequence_ = sequence;
        seq_len_ = seq_len;
        rid_ = rid;
    };
};

class Loader  {
 public:
  Loader();
  ~Loader();
  // return the number of sequences in the FASTA file
  // file_name: the name of the file to be read
  int CountFastaNumSeqs(const char *file_name);
  /* loading the FASTA file
     @return the number of sequences
     @param alphabet: the alphabet for the file
     @param file_name: the name of the file to be read
     @param header: the two-dimensional array that stores the header information
     @param seq: the two-dimensional array that stores the sequence information*/
  void LoadFasta(BioAlphabet &alphabet, const char* file_name, char** header, char** seq, int* seq_len, uint32_t* seq_id);
  int LoadFasta(BioAlphabet &alphabet, const char* file_name, char** seq, int* seq_len, uint32_t* seq_id);
  void ReverseComplement(
    char *target, char *source, const int& len
  );

 private:
  // checking for special characters in the string
  // alphabet: the alphabet for the string
  // sseq: the string to be checked
  // freq_cutoff: the portion of non-standard char cutoff to drop the string (return FALSE)
  // otherwise modify the string and returns TRUE
  bool CheckSpecialChar(BioAlphabet &alphabet, std::string &sseq, float freq_cutoff = 0.9);
  bool IsLowComplexity(std::string &seq);
  void RecordSequence(
    char **header, char **seq,
    std::string &single_header, std::string &single_seq,
    const int index
  );
  void RecordSequence(
    char** seq, std::string& single_seq, const int& index
  );

};

#endif
