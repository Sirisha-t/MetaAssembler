#include "loader.h"

using namespace std;

Loader::Loader()  {
  return;
}

Loader::~Loader() {
  return;
}


int Loader::CountFastaNumSeqs(const char *file_name)  {
  std::ifstream ifstrm(file_name, std::ios_base::in);
	if (!ifstrm) {
    std::cerr << "Cannot open file: " << file_name << std::endl;
    exit (1);
  }
  int count = 0;
  std::string line;
  while (std::getline(ifstrm, line)) {
    if(line[0] == '>') count++;
  }

  ifstrm.close();
  return count;
}

void Loader::RecordSequence(char **header, char **seq, std::string &single_header, std::string &single_seq, const int index) {
  assert(single_header.length() > 0);
  assert(single_seq.length() > 0);
  //std::cout<<index<<"\n";
  header[index] = new char[single_header.length() + 1];
  strcpy(header[index], single_header.c_str());
  seq[index] = new char[single_seq.length() + 1];
  strcpy(seq[index], single_seq.c_str());
  return;
}

bool Loader::CheckSpecialChar(BioAlphabet &alphabet, std::string &sseq, float freq_cutoff) {
  int num_special = 0;
  for(int i = 0; i < sseq.length(); ++ i)  {
    if(!alphabet.IsValid(sseq[i]))  {
      if(alphabet.IsValid(toupper(sseq[i])))  {
        sseq[i] = toupper(sseq[i]);
      } else  {
        ++ num_special;
        sseq[i] = alphabet.RandomChar();
      }
    }
  }
  if(num_special / sseq.length() > 1 - freq_cutoff) return false;
  else return true;
}

void Loader::ReverseComplement(char *target, char *source, const int& len)  {
  int i;
  for(i = 0; i < len; ++ i) {
    switch(source[len - i - 1]) {
      case 'A':
      case 'a':
        target[i] = 'T';
        break;
      case 'C':
      case 'c':
        target[i] = 'G';
        break;
      case 'G':
      case 'g':
        target[i] = 'C';
        break;
      case 'T':
      case 't':
        target[i] = 'A';
        break;
      case 'N':
        target[i] = 'N';
        break;
      default:
        cout << "Loader::ReverseComplement: Invalid letter: " << source[len - i - 1] << endl;
        cout << i << "  " << len << " " << strlen(source) << endl;
        exit(0);
    }
  }
  target[i] = '\0';
  return;
}

// with header loaded
void Loader::LoadFasta(BioAlphabet &alphabet, const char* file_name, char** header, char** seq, int* seq_len, uint32_t* seq_id) {
  // opens the file and reads line-by-line
  std::ifstream ifstrm(file_name, std::ios_base::in);
  std::string line, fasta_tag, fasta_seq;
  int count = 0;
  while (std::getline(ifstrm, line)) {
    if (line[0] == '>') {
      if (fasta_tag != "" && fasta_seq != "") {
        CheckSpecialChar(alphabet, fasta_seq);
        RecordSequence(header, seq, fasta_tag, fasta_seq, count);
        seq_len[count] = fasta_seq.length();
        //std::cout<<count<<"\n";
        seq_id[count] = count;
        ++count;
      }
      fasta_tag = line.substr(1, line.find(' ')); fasta_seq = "";
    } else fasta_seq += line;
  }
  ifstrm.close();
  // handle the last sequence
  if (fasta_tag != "" && fasta_seq != "") {
    CheckSpecialChar(alphabet, fasta_seq);
    RecordSequence(header, seq, fasta_tag, fasta_seq, count);
    seq_len[count] = fasta_seq.length();
    seq_id[count] = count;
    ++ count;
  }

  
  //for(int i = 0; i < count; ++ i) {
    //cout << i << "  " << seq[i]<<"  "<<header[i]<<"  "<<seq_len[i]<<"  "<<seq_id[i]<<" "<<endl;
  //}
  return ;
}
