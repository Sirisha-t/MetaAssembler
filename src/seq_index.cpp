#include "seq_index.h"

using namespace std;

MMBuild::MMBuild(void)  {
  is_header_loaded_ = is_sequence_loaded_ = false;
  is_mm_built_ = is_k_array_built_ = false;
  is_size_counted_ = false;
  is_alphabet_set_ = false;
  is_block_loaded_ = false;
  return;
}

MMBuild::MMBuild(BioAlphabet& alpha, bool &reduced_alpha)  {
  if(alpha.GetSeqType()== DNA) alpha_ = 1;
  else if(!reduced_alpha) alpha_ = 0;
  else  alpha_ = 2;
  alphabet_ = alpha;
  is_alphabet_set_ = true;
  is_block_loaded_ = false;
  is_header_loaded_ = is_sequence_loaded_ = is_mm_built_ = false;
  is_size_counted_ = false;
  return;
}


MMBuild::~MMBuild(void) {
  return;
}
void MMBuild::PrintAllSeqsInfo(void)  {
  cout << "Num  Header  Sequences  SeqLen RID" << endl;
  for(int i = 0; i < num_seqs_; ++ i) {
    cout << i << "  " << sequence_[i]<<"  "<<header_[i]<<"  "<<seq_len_[i]<<"  "<<seq_id_[i]<<" "<<endl;
  }
  return;
}


// Accessing the contents in the object
void MMBuild::LoadSequences(std::string& seq_file, const bool &rev_comp)  {
  assert(is_alphabet_set_);
  if(rev_comp && alphabet_.GetSeqType() != DNA)  {
    cout << "MMBuild::LoadSequence: reverse complementary is only available for DNA (nucl) sequences!!!" << endl;
    exit(0);
  }
  Loader seq_loader;
  num_seqs_ = seq_loader.CountFastaNumSeqs(seq_file.c_str());
  cout << "# of seqs:  " << num_seqs_ << endl;
  // double the number of sequences if we need to consider reverse complementary
  uint32_t ns = num_seqs_;   // the original number of sequences
  if(rev_comp) {
    ns = num_seqs_;
    num_seqs_ *= 2;
  }
  sequence_ = new char* [num_seqs_];
  header_ = new char* [num_seqs_];
  seq_len_ = new int [num_seqs_];
  seq_id_ = new uint32_t [num_seqs_];

  for(int i = 0; i < num_seqs_; ++ i) {
    sequence_[i] = header_[i] = NULL;
    seq_len_[i] = 0;
    seq_id_[i] = 0;
  }
  seq_loader.LoadFasta(alphabet_, seq_file.c_str(), header_, sequence_, seq_len_, seq_id_);

  if(rev_comp)  {
    for(int i = 0; i < ns; ++ i) {
      sequence_[ns + i] = new char[seq_len_[i] + 1];
      seq_loader.ReverseComplement(sequence_[ns + i], sequence_[i], seq_len_[i]);
      seq_len_[ns + i] = seq_len_[i];
      seq_id_[ns + i] = seq_id_[i];
      header_[ns + i] = new char[strlen(header_[i]) + 1];
      strcpy(header_[ns + i], header_[i]);
    }
    //cout << "sequence reverse complemented." << endl;
  }
  is_sequence_loaded_ = true;
  is_header_loaded_ = true;
  is_rev_comp_ = rev_comp;
  return;

}


void MMBuild::BuildIndex(int& kmer, int& window, std::string& dir, std::string& file_stem)  {
  if(is_mm_built_) { delete mm_index_;}
  //check if index has already been built and just load it
  string idx_file = dir + "/" + file_stem + ".idx";
  std::ifstream in_fh(idx_file);
  if (in_fh.good()) {
    std::cout<<"Index has already been built. Loading it now."<<endl;
    mm_index_ = new MMSketch((char**) sequence_, (char**) header_, seq_len_, seq_id_, num_seqs_, kmer, window, alpha_, 0);
    this->mm_index_->LoadMMIndexFile(idx_file.c_str());
    is_mm_loaded_ = true;
  }
  else{
    mm_index_ = new MMSketch((char**) sequence_, (char**) header_, seq_len_, seq_id_, num_seqs_, kmer, window, alpha_, 1);
  }
  is_mm_built_ = true;
  return;
}

void MMBuild::DumpIndextoFile(std::string& dir, std::string& file_stem)  {
  if(!is_mm_built_){
    std::cout<<"Cannot load index to file. Please build the index first\n\n";
    exit(0);
  }
  std::string idx_file = dir + "/" + file_stem + ".idx";
  this->mm_index_->DumpAllMMIndex(idx_file.c_str());
  return;
}

void MMBuild::LoadIndex(std::string& dir, std::string& file_stem){
  string idx_file = dir + "/" + file_stem + ".idx";
  this->mm_index_->LoadMMIndexFile(idx_file.c_str());
  is_mm_loaded_ = true;
}

void MMBuild::LoadMultiIndex(std::string& dir, std::string& file_stem) {
  assert(is_multi_ && is_block_loaded_);
  block_mm_.resize(num_blocks_);
  for(int i = 0; i < num_blocks_; ++ i) {
    string idx_file = dir + "/" + file_stem + "." + std::to_string(i) + ".midx";

    if(!boost::filesystem::exists(idx_file))  {
      cerr << "Error: LoadMultiIndex: Minimizers indexing file does not exist. Abort." << endl;
      exit(0);
    }
    // load the existing multi index
    block_mm_[i] = new MMSketch();
    block_mm_[i]->setReadCount(block_size_[i + 1] - block_size_[i] + 1);
    block_mm_[i]->setSequences(sequence_ + block_size_[i]);
    block_mm_[i]->LoadMMIndexFile(idx_file.c_str());
    
    //block_mm_[i]->PrintMinimizers();
  }
  
  is_mm_built_ = true;
  return;
}

void MMBuild::MapQuery(int& kmer, int& window, int& id)  {
  if(!is_mm_built_) { std::cout<<"Index needs to be built before mapping. \n"; exit(0);}
  
  this->mm_index_->Map(kmer, window, 33, id);
  
  return;
}

void MMBuild::MapQuery(MMBuild& query_seq, int& kmer, int& window, int& id )  {
  if(!is_mm_built_) { std::cout<<"Index needs to be built before mapping. \n"; exit(0);}
  
  //Map each query to reference index
  this->mm_index_->Map(query_seq.sequence_[id], query_seq.seq_len_[id], query_seq.header_[id], query_seq.seq_id_[id], kmer, window, 33);
  
  return;
}

void MMBuild::PrintOverlaps(){
  mm_index_->PrintOverlapInfo();
  return;
}



void MMBuild::WritePAF(std::ofstream& outfile){
  mm_index_->WriteOverlaps(outfile);
  return;
}


void MMBuild::SetBlockConfig(const int& num_blocks, std::string &dir, std::string &file_stem) {
  uint16_t begin = 0;
  block_size_.push_back(0);
  
  //std::cout<<"num seqs : "<<num_seqs_<<"\n";
  uint16_t npb = num_seqs_ / num_blocks;
  //std::cout<<"num seqs per block : "<<npb<<"\n";
  if(npb < 1) npb = 1;  // set at least one sequence per block
  for(int i = 1; i < num_blocks; ++ i) {
    uint16_t acc_ns = npb * i;
    if(acc_ns < num_seqs_)  {
      block_size_.push_back(npb * i);
    } else  {
      break;  // break the loop if we have more block than sequences
    }
  }

  num_blocks_ = num_blocks;
  // write the block information
  std::string out_file = dir + "/" + file_stem + ".bsz";
  std::ofstream out_fh(out_file.c_str(), ios_base::out | ios_base::binary);
  if(!out_fh.good())  {
    cout << "Error:: Cannot write block-size index file: " << out_file << endl;
    exit(1);
  }
  out_fh.write((char*) &num_blocks_, sizeof(uint16_t));
  for(int i = 0; i < block_size_.size(); ++ i) {
    out_fh.write((char*) &block_size_[i], sizeof(uint16_t));
  }
  out_fh.close();

  is_multi_ = is_block_loaded_ = true;
  return;
}

void MMBuild::SplitSequence(std::vector<MMBuild>& block_mm)  {
  assert(this->is_sequence_loaded_ && this->is_header_loaded_);
  assert(this->is_alphabet_set_);
  assert(this->is_block_loaded_);
  assert(this->is_multi_);
  assert(block_mm.size() == this->num_blocks_);

  //for(int i = 0; i < this->block_size_.size(); ++ i) {
   // cout << "bs:  " << i << " " << block_size_[i] << endl;
  //}

  for(int i = 0; i < this->num_blocks_; ++ i) {
    block_mm[i].alphabet_ = this->alphabet_;
    block_mm[i].is_alphabet_set_ = true;
    block_mm[i].is_block_loaded_ = false;
    block_mm[i].is_multi_ = false;
    block_mm[i].is_mm_built_ =  false;
    block_mm[i].is_size_counted_ = false;
    if(i < this->num_blocks_ - 1)  {
      block_mm[i].num_seqs_ = this->block_size_[i + 1] - this->block_size_[i];
    } else{
      block_mm[i].num_seqs_ = this->num_seqs_ - this->block_size_[i];
    }
    block_mm[i].sequence_ = new char* [block_mm[i].num_seqs_];
    block_mm[i].header_ = new char* [block_mm[i].num_seqs_];
    block_mm[i].seq_len_ = new int [block_mm[i].num_seqs_];
    block_mm[i].seq_id_ = new uint32_t [block_mm[i].num_seqs_];
    for(int j = 0; j < block_mm[i].num_seqs_; ++ j) {
      block_mm[i].sequence_[j] = new char[strlen(this->sequence_[this->block_size_[i] + j]) + 1];
      strcpy(block_mm[i].sequence_[j], this->sequence_[this->block_size_[i] + j]);
      block_mm[i].seq_len_[j] = strlen(block_mm[i].sequence_[j]);
      block_mm[i].seq_id_[j] = this->seq_id_[this->block_size_[i] + j];
      block_mm[i].header_[j] = new char[strlen(this->header_[this->block_size_[i] + j]) + 1];
      strcpy(block_mm[i].header_[j], this->header_[this->block_size_[i] + j]);
    }
    block_mm[i].is_sequence_loaded_ = block_mm[i].is_header_loaded_ = true;

    //cout << "done copying sequences" << endl;
  }
  return;
}





