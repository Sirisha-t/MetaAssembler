#ifndef _MM_SKETCH_
#define _MM_SKETCH_

#include "loader.h"
#include "util_func.h"
#include "murmur3.hpp"
#include "xxhash64.h"
#include "ksort.h"


#include <boost/filesystem.hpp>
#include <algorithm>
#include <functional>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <vector>
#include <tuple>
#include <string>
#include <list>
#include <time.h>
#include <sys/stat.h>
#include<set>
#include <cstring>
#include <sstream>
#include <locale>
#include <cstdint>
#include <iostream>
#include <algorithm>
#include <assert.h>
#include <cmath>
#include <limits>
#include <omp.h>
#include <bitset>


#define MAX_BLOCK 128

// emulate 128-bit integers and arrays
typedef struct { uint64_t x, y; } mm128_t;
typedef struct {uint32_t x, y;} mm64_t;
typedef struct { size_t n, m; mm128_t *a; } mm128_v;
typedef struct { size_t n, m; mm64_t *a; } mm64_v;


struct MMIndex{
    uint64_t info_;

    MMIndex(const uint64_t& info ){
      info_ = info;
    }
};
 
//NEW MMHASH impl -- same as minimap2's mm128_t struct
struct MMHash{  
    uint64_t x, y;
};


//Structure to store the query and target matches
struct MMMatch{
    uint32_t id;
    int strand;
    int start;
    int end;

    MMMatch(const uint32_t& id_, const int& strand_, const int& start_, const int& end_){
        id = id_;
        strand = strand_;
        start = start_;
        end = end_;
    }
};

// Structure to store the overlap output information
struct MMOverlapInfo{
    std::string qname;
    int qstart;
    int qend;
    std::string tname;
    int tstart;
    int tend;
    std::string strand;

    MMOverlapInfo(const std::string& qname_, const int& qstart_, const int& qend_, const std::string& tname_, const int& tstart_, const int& tend_, const std::string& strand_){
        qname = qname_;
        qstart = qstart_;
        qend = qend_;
        tname = tname_;
        tstart = tstart_;
        tend = tend_;
        strand = strand_;
    }
};

struct hash_comparer{
    inline bool operator()(const MMHash& one, const MMHash& two){

            return (one.x < two.x || one.x==two.x && one.y < two.y);

    }
};

struct comparer{
    inline bool operator()(const MMIndex& one, const MMIndex& two){

            return (one.info_ <= two.info_);

    }
};

struct qcomparer{
    inline bool operator()(const MMMatch& one, const MMMatch& two){

            return (one.id < two.id || one.id==two.id && one.strand < two.strand ||
                one.id==two.id && one.strand == two.strand && one.start < two.start ||
               one.id==two.id && one.strand == two.strand && one.start == two.start && one.end < two.end);

    }
};

struct Cell {
    int score;
    int i;
    int j;
};



class MMSketch {
 public:
    //Default constructor
    MMSketch();
    MMSketch(char **seq, char** header, int* seqlen, uint32_t* seqid, int n, int kmer, int window,int& alpha, bool idx_flag);
    ~MMSketch(void);
    void PrintMinimizers(void);
    void FindSeedMatches(int& kmer, int& window, int epsilon, int& id);
    void Map(char* seq, int seqlen, char* header, uint32_t& seqid, int& kmer, int& window, int epsilon);
    void DumpAllMMIndex(const char *idx_file);
    void LoadMMIndexFile(const std::string& idxname);
    void PrintOverlapInfo();
    void WriteOverlaps(std::ofstream& outfile);

    /** Set sequence reads */
	void setSequences( char **s ) { seqs   = s; }

	/** Set read count */
	void setReadCount( RIDTYPE n )    { nreads = n; }

 
    std::unordered_map<uint32_t, std::vector<MMIndex>>& getMMIndex() {
            return MMHashTable;
        }
    

    struct PairHash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1, T2>& p) const {
        // Hash both elements and combine the hashes
        std::size_t h1 = std::hash<T1>{}(p.first);
        std::size_t h2 = std::hash<T2>{}(p.second);
        return h1 ^ h2;
    }
};
   

 private:
    void init(char** seq, char** header, int* seqlen, uint32_t* seqid, int n, int kmer, int window, int alpha);
    void buildSketch();
    std::vector<std::pair<char, int>> get_cigar(const std::string& s1, const std::string& s2, const std::vector<Cell>& traceback);
    std::vector<MMMatch> longest_increasing_subset(std::vector<MMMatch>& matches, std::string orientation);
    std::vector<MMMatch> slice(std::vector<MMMatch>const & matches, int start, int end);
    void computeMinimizerSketch(char* seq, int seqlen, char* header, uint32_t& id, int&alphabet, int& kmer, int &window,std::vector<MMHash>& hashes);   
    void xDropAlignmentExtend(const std::string& seq1, const std::string& seq2, int q_start, int q_end, int t_start, int t_end, std::string& aligned_seq1, std::string& aligned_seq2);
    inline uint32_t invertableHash(uint32_t hash_val, uint32_t mask);
    inline std::string reverseComplement(const std::string& seq);
    inline int hashValue(const char& n);
    uint64_t hash(const std::string & kmer);
    



 protected:
    int k, w;
    //std::vector<FastaInfo> seqs;
    char** seqs;
    char** headers;
    int* seqlens;
    uint32_t* seqids;
    bool* qmap;
    int nreads;
    int alphabet;
    std::vector<MMHash> mm_hashes;
    std::unordered_map<uint32_t, std::vector<MMIndex>> MMHashTable;
    //std::unordered_set<std::pair<uint16_t, uint16_t>, PairHash> overlap_pair;
    std::vector<MMOverlapInfo> MMOverlaps;
    bool is_kmer_set_;
    bool is_mm_built_;
    bool is_window_set_;
 


};

#endif
