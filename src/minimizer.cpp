#include "minimizer.h"
#include "kvec.h"


using namespace std;

#ifdef __GNUC__
#define FORCE_INLINE __attribute__((always_inline)) inline
#else
#define FORCE_INLINE inline
#endif


// Scoring parameters
const int match_score = 1;
const int mismatch_score = -1;
const int gap_score = -1;

// X-drop parameter
const int GAP_PENALTY = 1; // gap penalty
const int X_DROP_THRESHOLD = 50; // X-drop threshold


//DNA overlap parameters=
int MAPLEN = 25;
int MINHANG = 10;

unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

static FORCE_INLINE uint64_t rotl64 ( uint64_t x, int8_t r )
{
  return (x << r) | (x >> (64 - r));
}

#define ROTL64(x,y)	rotl64(x,y)

#define BIG_CONSTANT(x) (x##LLU)
//-----------------------------------------------------------------------------
// Block read - if your platform needs to do endian-swapping or can only
// handle aligned reads, do the conversion here

#define getblock(p, i) (p[i])

static FORCE_INLINE uint64_t fmix64 ( uint64_t k )
{
  k ^= k >> 33;
  k *= BIG_CONSTANT(0xff51afd7ed558ccd);
  k ^= k >> 33;
  k *= BIG_CONSTANT(0xc4ceb9fe1a85ec53);
  k ^= k >> 33;

  return k;
}


MMSketch::MMSketch()
{
    init(NULL, NULL, 0, 0, 0, 0, 0, 0);
}

MMSketch::MMSketch(char **seq, char **header, int* seqlen, uint32_t* seqid, int n, int kmer, int window, int& alpha, bool idx_flag)  {
 init(seq, header, seqlen, seqid, n, kmer, window, alpha);
 std::cout<<"MMSketch initialized successfully. Building sketch now for k="<<kmer<<",w="<<window<<"\n";
 if(idx_flag)
    MMSketch::buildSketch();
 return;
}

void MMSketch::init( char **seq, char **header, int* seqlen, uint32_t* seqid, int n,int kmer, int window, int alpha){
  seqs = seq;
  headers = header;
  seqlens = seqlen;
  seqids = seqid;
  alphabet = alpha;
  k = kmer;
  w = window;
  nreads = n;
  is_mm_built_ = false;
  is_kmer_set_ = true;
  is_window_set_ = true;

  if(alphabet == 0 || alphabet == 2){
    MAPLEN = 10;
    MINHANG = 5;
  }
  
  qmap = new bool [nreads];   
  for(int i = 0; i < nreads; ++ i) {
    qmap[i] = false;
  }
}

MMSketch::~MMSketch(void) {
  return;
}


/** start of murmur hash implementation --unused **/

inline string str_canonical(const std::string& kmer) {
    auto kmer_rev = kmer;
    std::reverse(kmer_rev.begin(), kmer_rev.end());
    for (size_t j = 0; j < kmer_rev.length(); ++j) {
        if (kmer_rev[j] == 'A') kmer_rev[j] = 'T';
        else if (kmer_rev[j] == 'T') kmer_rev[j] = 'A';
        else if (kmer_rev[j] == 'C') kmer_rev[j] = 'G';
        else if (kmer_rev[j] == 'G') kmer_rev[j] = 'C';
    }
    return kmer < kmer_rev ? kmer : kmer_rev;
}

uint64_t MurmurHash3_x64_128 ( const void * key, int len,
                           unsigned int seed )
{
  const uint8_t * data = (const uint8_t*)key;
  const int nblocks = len / 4;
  int i;

  uint64_t h1 = seed;

  uint64_t c1 = BIG_CONSTANT(0x87c37b91114253d5);
  uint64_t c2 = BIG_CONSTANT(0x4cf5ad432745937f); 


  const uint64_t * blocks = (const uint64_t *)(data);

  for(i = 0; i < nblocks; i++)
  {
    uint64_t k1 = getblock(blocks,i);

    k1 *= c1; k1  = ROTL64(k1,31); k1 *= c2; h1 ^= k1;

    h1 = ROTL64(h1,27); h1 = h1*5+0x52dce729;

  }

  const uint8_t * tail = (const uint8_t*)(data + nblocks*16);

  uint64_t k1 = 0;

  switch(len & 8)
  {

  case  8: k1 ^= (uint64_t)(tail[ 7]) << 56;
  case  7: k1 ^= (uint64_t)(tail[ 6]) << 48;
  case  6: k1 ^= (uint64_t)(tail[ 5]) << 40;
  case  5: k1 ^= (uint64_t)(tail[ 4]) << 32;
  case  4: k1 ^= (uint64_t)(tail[ 3]) << 24;
  case  3: k1 ^= (uint64_t)(tail[ 2]) << 16;
  case  2: k1 ^= (uint64_t)(tail[ 1]) << 8;
  case  1: k1 ^= (uint64_t)(tail[ 0]) << 0;
           k1 *= c1; k1  = ROTL64(k1,31); k1 *= c2; h1 ^= k1;
  };

  //----------
  // finalization
    h1 ^= len;
    h1 = fmix64(h1);

    return h1;
}

uint64_t MurmurHash64A ( const void * key, int len, unsigned int seed )
{
	const uint64_t m = 0xc6a4a7935bd1e995;
	const int r = 47;

	uint64_t h = seed ^ (len * m);

	const uint64_t * data = (const uint64_t *)key;
	const uint64_t * end = data + (len/8);

	while(data != end)
	{
		uint64_t k = *data++;

		k *= m; 
		k ^= k >> r; 
		k *= m; 
		
		h ^= k;
		h *= m; 
	}

	const unsigned char * data2 = (const unsigned char*)data;

	switch(len & 7)
	{
	case 7: h ^= uint64_t(data2[6]) << 48;
	case 6: h ^= uint64_t(data2[5]) << 40;
	case 5: h ^= uint64_t(data2[4]) << 32;
	case 4: h ^= uint64_t(data2[3]) << 24;
	case 3: h ^= uint64_t(data2[2]) << 16;
	case 2: h ^= uint64_t(data2[1]) << 8;
	case 1: h ^= uint64_t(data2[0]);
	        h *= m;
	};
 
	h ^= h >> r;
	h *= m;
	h ^= h >> r;

	return h;
} 

uint64_t MMSketch::hash(const std::string & kmer) {
    string canonical_kmer = str_canonical(kmer);
    const char *c = canonical_kmer.c_str();
    return MurmurHash3_x64_128(c, canonical_kmer.size(),42);
}

/** end of murmurhash implementation **/

inline uint32_t MMSketch::invertableHash(uint32_t key, uint32_t mask){
  key = (~key) + (key << 21) & mask; // key = (key << 21) - key - 1;
  key = key ^ (key >> 24);
  key = (key + (key << 3)) + (key << 8) & mask; // key * 265
  key = key ^ (key >> 14);
  key = (key + (key << 2)) + (key << 4); // key * 21
  key = key ^ (key >> 28);
  key = key + (key << 31) & mask;
  return key;
}


 inline int MMSketch::hashValue(const char& n){
    switch(n)
    {
        case 'A':
            return '0'; break;
        case 'C':
            return '1' ; break;
        case 'G':
            return '2'; break;
        case 'T':
            return '3'; break;
        default:
            break;

    }
    assert(true) ;return (' ');

}

 inline std::string MMSketch::reverseComplement(const std::string& seq){
    std::string rev_comp = seq;
    int seqlen = seq.length();
    std::reverse(rev_comp.begin(), rev_comp.end());
    for (std::size_t i = 0; i < seqlen; ++i){
        switch (seq[i]){
        case 'A':
            rev_comp[i] = 'T';
            break;    
        case 'C':
            rev_comp[i] = 'G';
            break;
        case 'G':
            rev_comp[i] = 'C';
            break;
        case 'T':
            rev_comp[i] = 'A';
            break;
        }
    }
    return rev_comp;
}


/* Compute the minimizer sketch for each sequence -- Minimap2's impl */
void MMSketch::computeMinimizerSketch(char* seq, int seqlen, char* header, uint32_t& id, int& alphabet, int& k, int& w, std::vector<MMHash>& hashes ){
    std::vector<MMHash> tmp_hashes;
    //std::cout<<header<<"\t"<<id<<"\n";
    uint32_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1;
    uint32_t kmer[2] = {0,0}; //to store the forward and reverse kmers
	int i, j, l, buf_pos, min_pos, kmer_span = 0; // i is used to store current position in input seq, l stores length of current kmer
	MMHash buf[256], min = { UINT64_MAX, UINT64_MAX };
    int c;
    assert(seqlen > 0 && (w > 0 && w < 256) && (k > 0 && k <= 28)); 
	memset(buf, 0xff, w * 16);
   
    for (i = l = buf_pos = min_pos = 0; i < seqlen; ++i) {	
        MMHash mm_info = { UINT64_MAX, UINT64_MAX };
        if (alphabet == 1) 
        {
            int c = seq_nt4_table[(uint8_t)seq[i]];  
            if (c < 4 ) { // not an ambiguous base
                int z;
                kmer_span = l + 1 < k? l + 1 : k;
                kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer -- here mask is used for clearing out any bits outside the given bit range
                kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
                if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" 
                z = kmer[0] < kmer[1]? 0 : 1; // strand
                ++l;

                if (l >= k && kmer_span < 256) {
                    //std::cout<<"entered if\n";
                    mm_info.x = invertableHash(kmer[z], mask) << 8 | kmer_span; //use upto 24 bits for hash value
                    mm_info.y = (uint64_t)id<<32 | (uint32_t)i<<1 | z; //use only right most 16 bits of rid
                }
            }
        }
        
        buf[buf_pos] = mm_info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
        if (l == w + k - 1 && min.x != UINT64_MAX) { // special case for the first window - because identical k-mers are not stored yet
			for (j = buf_pos + 1; j < w; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y) 
                {
                    hashes.push_back(buf[j]);
                    tmp_hashes.push_back(buf[j]);
                }    
			for (j = 0; j < buf_pos; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y)
                {  
                    hashes.push_back(buf[j]);
                    tmp_hashes.push_back(buf[j]);        
                }      
		}
        if (mm_info.x <= min.x) { // a new minimum; then write the old min
			if (l >= w + k && min.x != UINT64_MAX) 
            {
                hashes.push_back(min);
                tmp_hashes.push_back(min);
            }
			min = mm_info, min_pos = buf_pos;
		}
        else if (buf_pos == min_pos) { // old min has moved outside the window
			if (l >= w + k - 1 && min.x != UINT64_MAX) 
            {
                hashes.push_back(min);
                tmp_hashes.push_back(min);
			}
			for (j = buf_pos + 1, min.x = UINT64_MAX; j < w; ++j) // the two loops are necessary when there are identical k-mers
				if (min.x >= buf[j].x) min = buf[j], min_pos = j; // >= is important s.t. min is always the closest k-mer
			for (j = 0; j <= buf_pos; ++j)
				if (min.x >= buf[j].x) min = buf[j], min_pos = j;
			if (l >= w + k - 1 && min.x != UINT64_MAX) { // write identical k-mers
				for (j = buf_pos + 1; j < w; ++j) // these two loops make sure the output is sorted
					if (min.x == buf[j].x && min.y != buf[j].y)
                        {
                            hashes.push_back(buf[j]);
                            tmp_hashes.push_back(buf[j]);
                        }
				for (j = 0; j <= buf_pos; ++j)
					if (min.x == buf[j].x && min.y != buf[j].y) 
                        {
                            hashes.push_back(buf[j]);
                            tmp_hashes.push_back(buf[j]);
                        }
			}
		}
        if (++buf_pos == w) buf_pos = 0;
	}
	if (min.x != UINT64_MAX){
        hashes.push_back(min);
        tmp_hashes.push_back(min);
    }

    
   /*if(std::string(header) == "NC_015678.1_823046_823145_2:0:0_1:0:0_34c/1"){
   // if(id == 213633){
        std::cout<<seq<<"  "<<header<<"  "<<id<<"\n";
        for(auto &hash: tmp_hashes){
            std::cout<<hash.x<<", "<<hash.y<<"\n";
    }
    std::cout<<"\n\n";
    }

        if(std::string(header) == "NC_015678.1_823047_823146_0:0:0_0:0:0_10969/1"){
        //if(id == 195977){
        std::cout<<seq<<"  "<<header<<"  "<<id<<"\n";
        for(auto &hash: tmp_hashes){
            std::cout<<hash.x<<", "<<hash.y<<"\n";
    }
    std::cout<<"\n\n";
    }*/
  

    return ;

}


/* Build  minimizer sketch for each target sequence
  param target_seq: target sequence 
  param int kmer: kmer size
  param int window : window size
  return a HashMap, with key as hash values and values as {targetSeq, position, strand}
 */
void MMSketch::buildSketch(){
    assert(is_kmer_set_);
    assert(is_window_set_);
    int count = 0;

    clock_t t;
    t = clock();
    
    //Compute mm sketch for each sequence 
    for(int i=0; i < nreads; i++){
        //std::cout<<"i: "<<i<<"\n";
        computeMinimizerSketch(seqs[i], seqlens[i], headers[i], seqids[i], alphabet, k, w, mm_hashes);   
        count++;
    }
    t = clock() - t;
    std::cout<<"Collected minimizers in "<<((double)t)/CLOCKS_PER_SEC << " s\n";
    
    t=clock();
    std::sort(mm_hashes.begin(), mm_hashes.end(), hash_comparer());
    t = clock() - t;
    std::cout<<"Sorted minimizers in "<<((double)t)/CLOCKS_PER_SEC << " s\n";
    PrintMinimizers();

    t=clock();

    //Insert hashes into hash table
    //Reserve memory to reduce reallocations
    MMHashTable.reserve(mm_hashes.size());
    for (const MMHash& mm_entry : mm_hashes) {
        auto it = MMHashTable.find(mm_entry.x);
        if (it != MMHashTable.end()) {
            if (it->second.empty() || it->second.back().info_ != mm_entry.y) {
                    it->second.emplace_back(MMIndex(mm_entry.y));
            }
            if (!std::is_sorted(it->second.begin(), it->second.end(), comparer())) {
                std::sort(it->second.begin(), it->second.end(), comparer());
            }
        } else {
            MMHashTable[mm_entry.x] = { MMIndex(mm_entry.y) };
            }
        }

    t = clock() - t;
    std::cout<<"Created hash table index in "<<((double)t)/CLOCKS_PER_SEC << " s\n";

    mm_hashes.clear();

    //PrintMinimizers();
 
}

void MMSketch::PrintMinimizers(void){
    for(auto &elem: MMHashTable ){
        std::cout<<elem.first<<":\t";
        for(MMIndex &item : elem.second){
            std::cout<<"("<<item.info_<<")"<<", ";
        }
        std::cout<<"\n";
    }
    return;
}


std::vector<MMMatch> MMSketch::slice(std::vector<MMMatch>const & matches, int start, int end)
{  
    auto first = matches.cbegin() +start;
    auto last = matches.cbegin() + end + 1;

    std::vector<MMMatch> sub_matches(first, last);

    return sub_matches;
}

std::vector<MMMatch> MMSketch::longest_increasing_subset(std::vector<MMMatch>& matches, std::string orientation){
    std::vector<MMMatch> longest_subset;

    // Iterate through the sorted matches and find the longest increasing subset
    for (const MMMatch& match : matches) {
        //std::vector<MMMatch> current_subset;

        if (longest_subset.empty() || match.end > longest_subset.back().end) {
            longest_subset.push_back(match);
        }
    }
    
    return longest_subset;
}

/*std::vector<MMMatch> MMSketch::longest_increasing_subset(std::vector<MMMatch>& matches, std::string orientation)
{
    std::vector<int> map_len;
    for(MMMatch m : matches)    map_len.push_back(m.end);
    if(orientation == "-")  std::reverse(map_len.begin(), map_len.end());

    std::vector<int> P;
    std::vector<int> P_reversed;

    int n = map_len.size();
    for (int i = 0; i < n; i++) {
        auto it = lower_bound(P.begin(), P.end(), map_len[i]);
        if (it == P.end()) {
            P.push_back(map_len[i]);
        }
        else {
            *it = map_len[i];
        }
    }
    if(orientation == "+" || orientation == "-")  return P;
    
    
}*/

std::vector<std::pair<char, int>> MMSketch :: get_cigar(const std::string& s1, const std::string& s2, const std::vector<Cell>& traceback) {
    std::vector<std::pair<char, int>> cigar;
    int i = s1.size();
    int j = s2.size();
    for (int k = traceback.size() - 1; k >= 0; k--) {
        int i_prev = traceback[k].i;
        int j_prev = traceback[k].j;
        //std::cout<<i_prev<<"\t"<<j_prev<<"\n";
        if (i_prev == i) {
            cigar.emplace_back('D', j_prev - j);
        } else if (j_prev == j) {
            cigar.emplace_back('I', i_prev - i);
        } else {
            if (s1[i_prev-1] == s2[j_prev-1]) {
                cigar.emplace_back('M', 1);
            } else {
                cigar.emplace_back('X', 1);
            }
        }
        i = i_prev;
        j = j_prev;
    }
    reverse(cigar.begin(), cigar.end());
    
    return cigar;
}


// Perform x-drop alignment for sequence extension
void MMSketch::xDropAlignmentExtend(const std::string& seq1, const std::string& seq2,
                          int q_start, int q_end, int t_start, int t_end,
                          std::string& aligned_seq1, std::string& aligned_seq2) {

    // Extend the alignment from the known overlap positions
    int q_len = q_end - q_start + 1;
    int t_len = t_end - t_start + 1;

    // Initialize the scoring matrix
    std::vector<std::vector<int>> score(q_len + 1, std::vector<int>(t_len + 1, 0));
    // Fill in the scoring matrix using x-drop
    for (int i = 1; i <= q_len; ++i) {
        for (int j = 1; j <= t_len; ++j) {
            int q_pos = q_start + i - 1;
            int t_pos = t_start + j - 1;
            int match = score[i - 1][j - 1] + (seq1[q_pos] == seq2[t_pos] ? match_score : mismatch_score);
            int gap1 = score[i - 1][j] + gap_score;
            int gap2 = score[i][j - 1] + gap_score;

            score[i][j] = std::max({0, match, gap1, gap2}); // X-drop threshold check
        }
    }

    int max_score = 0;
    int best_i = 0, best_j = 0;
    for (int i = 1; i <= q_len; ++i) {
        if (score[i][t_len] > max_score) {
            max_score = score[i][t_len];
            best_i = i;
            best_j = t_len;
        }
    }
    for (int j = 1; j <= t_len; ++j) {
        if (score[q_len][j] > max_score) {
            max_score = score[q_len][j];
            best_i = q_len;
            best_j = j;
        }
    }
    // Traceback 
    aligned_seq1 = "";
    aligned_seq2 = "";
    int i = best_i, j = best_j;
    while (i > 0 && j > 0 && score[i][j] > 0) {
        int q_pos = q_start + i - 1;
        int t_pos = t_start + j - 1;

        if (score[i][j] == score[i - 1][j - 1] + (seq1[q_pos] == seq2[t_pos] ? match_score : mismatch_score)) {
            aligned_seq1 = seq1[q_pos] + aligned_seq1;
            aligned_seq2 = seq2[t_pos] + aligned_seq2;
            i--;
            j--;
        } else if (score[i][j] == score[i - 1][j] + gap_score) {
            aligned_seq1 = seq1[q_pos] + aligned_seq1;
            aligned_seq2 = "-" + aligned_seq2;
            i--;
        } else {
            aligned_seq1 = "-" + aligned_seq1;
            aligned_seq2 = seq2[t_pos] + aligned_seq2;
            j--;
        }
    }
}

/*
Mapping each query sequence to minimizer index hash table --align mode (read vs ref)
*/
void MMSketch::Map(char* seq, int seqlen, char* header, uint32_t& seqid, int& kmer, int& window, int epsilon)
{
    std::vector <MMMatch> query_matches;
    std::vector<MMMatch> potential_olap;
    std::vector<MMMatch> indices;
    std::string orientation = "";
    std::string qname = header;

    int tstart, tend, qstart, qend, orient, mapped_bp;
    std::string tname, first_target, last_target;
    int first_start = 0;
    int first_end = 0;
    int last_end = 0;
    int last_start = 0;
    uint32_t qid, tid;
    

    // Get mm sketch for each query sequence
    std::vector<MMHash> qhash ;
    computeMinimizerSketch(seq, seqlen, header, seqid, alphabet, kmer, window, qhash);   
    for(MMHash q : qhash){
        //if(qname == "146504" )
           std::cout<<"Query "<<q.x<<"\t"<<q.y<<"\n";
        try{
            for(MMIndex t: MMHashTable[q.x]){  
                qid = (q.y>>32) & 0xFFFFFFFF;
                int qstrand = q.y & 0x01;
                uint16_t qpos = (q.y & 0xFFFFFFFF) >> 1;

                tid = (t.info_>>32) & 0xFFFFFFFF;
                int tstrand = t.info_ & 0x01;
                uint16_t tpos = (t.info_ & 0xFFFFFFFF) >> 1;
                //std::cout<<tpos<<"\t"<<tid<<"\t"<<qpos<<"\t"<<qid<<"\n";

                //ignore same read pair overlaps
                if(tid == qid)  
                    continue;

                if( tstrand == qstrand){
                    query_matches.emplace_back(MMMatch(tid, 0, qpos - tpos, tpos));
                } 
                else{
                    query_matches.emplace_back(MMMatch(tid, 1, qpos + tpos, tpos));
                }
            }    
             //if(std::string(header) == "NC_015678.1_1459772_1459871_1:0:0_0:0:0_452b/1")    std::cout<<seqid<<"\t"<<header<<"\t"<<tid<<"\n";      
        }
        catch(const std::out_of_range& oor) {};
    }
    std::sort(query_matches.begin(), query_matches.end(), qcomparer());

    int b=0;

    if(query_matches.size()>0){
        for(int e = 0; e < query_matches.size()-1; ++e){
            if(e == query_matches.size()-1 | 
                query_matches[e+1].id != query_matches[e].id | 
                query_matches[e+1].strand != query_matches[e].strand | 
                query_matches[e+1].start - query_matches[e].start >= epsilon)
                {   
                    //extract potential overlaps for each query match
                    potential_olap = slice(query_matches, b, e);
                    
                    /*if(std::string(header) == "146504"){
                        for(auto &p: potential_olap){
                        std::cout<<p.id<<"  "<<p.start<<"  "<<p.end<<"  "<<p.strand<<"\n";
                    }
                     }*/
                    
                    if(query_matches[e].strand == 0)
                        orientation = "+"; 
                    
                    else if (query_matches[e].strand == 1)
                        orientation = "-";                 
                
                    b = e+1;
                
                    if(orientation == "+"){
                        indices = longest_increasing_subset(potential_olap, orientation); 
                    }
                    else{
                        std::reverse(potential_olap.begin(), potential_olap.end());
                         
                        indices = longest_increasing_subset(potential_olap, orientation); 
                    }
                    
                    if(indices.size() <= 3) continue;

                    //int first_mm = indices[0];
                    //int last_mm = indices[indices.size()-1];

                    MMMatch first_mm = indices[0];
                    MMMatch last_mm = indices[indices.size()-1];
                    tname = headers[tid];
                    qname = headers[qid];
                         
                    int mapped_bp = indices.size() * kmer;
                    int qlen = seqlen;
                    int tlen = seqlens[first_mm.id];

                   
                    if(orientation == "+"){
                        tstart = first_mm.end - kmer;
                        tend = last_mm.end ;
                        qstart = first_mm.start ;     
                        qend = last_mm.start + last_mm.end ;
                         
                        
                    }
                    else{
                        tstart = tlen - first_mm.end - kmer;
                        tend = tlen - last_mm.end;
                        qstart = first_mm.start - first_mm.end;    
                        qend = last_mm.start - last_mm.end + kmer;
                    }

                    int overhang = std::min(qstart,tstart) + std::min(qlen-qend, tlen-tend);
                    int maplen = std::max((qend - qstart), (tend - tstart));
                    
                    //Do not consider overlaps lesser than MAPLEN
                    if(maplen < MAPLEN) continue;
                        
                    //overhang is the bases after and before the computed overlaps
                    /*if(overhang > std::min(MINHANG, int(0.8*maplen))){
                       continue;
                    }*/
                        
                    if(qstart <= tstart && (qlen - qend) <= (tlen-tend))
                        continue;
                    else if(qstart >= tstart && (qlen - qend) >= (tlen-tend))
                        continue;
                    
                    //Output the overlaps and store into vector
                    //if(qid == 17707)
                           // std::cout<<seqid<<"\t"<<tname<<"  "<<tid<<"  "<<tstart<<"  "<<tend<<"  "<<qname<<"  "<<qid<<"  "<<qstart<<"  "<<qend<<"\n"; 

                    //outfile<<qname<<"\t"<<qlen<<"\t"<<qstart<<"\t"<<qend<<"\t"<<orientation<<"\t"<<tname<<"\t"<<tlen<<"\t"<<tstart<<"\t"<<tend<<"\n";
                    //MMOverlapInfo olap_info = MMOverlapInfo(qname,qstart,qend,tname,tstart,tend,orientation);
                    //MMOverlaps.push_back({olap_info});
           }
        }
    } 
   return;
}  


/*
Mapping each query sequence to minimizer index hash table --overlap mode (read vs read)
*/
void MMSketch::FindSeedMatches(int& kmer, int& window, int epsilon, int& id)
{
    //std::cout<<"Map for query : "<<id<<"\n";
    std::vector <MMMatch> query_matches;
    std::vector<MMMatch> potential_olap;
    std::vector<MMMatch> indices;
    std::string orientation = "";
    std::string qname = headers[id];

    int tstart, tend, qstart, qend, orient, mapped_bp;
    std::string tname, first_target, last_target;
    uint32_t qid, tid;

    // Get mm sketch for each query sequence
    std::vector<MMHash> qhash ;
    computeMinimizerSketch(seqs[id], seqlens[id], headers[id], seqids[id], alphabet, kmer, window, qhash);   
    //std::cout<<"Computed mm sketch\n";
    for(MMHash q : qhash){
        try{
            for(MMIndex t: MMHashTable[q.x]){  
                qid = (q.y>>32) & 0xFFFFFFFF;
                int qstrand = q.y & 0x01;
                uint16_t qpos = (q.y & 0xFFFFFFFF) >> 1;

                tid = (t.info_>>32) & 0xFFFFFFFF;
                int tstrand = t.info_ & 0x01;
                uint16_t tpos = (t.info_ & 0xFFFFFFFF) >> 1;

               // std::cout<<qname<<"\t"<<qid<<"\t"<<tid<<"\n";
                //ignore same read pair overlaps
                if(tid == qid)  
                    continue;

                //ignore tid mapping if it has previously been computed 
                if(qmap[tid])
                    continue;

                if( tstrand == qstrand){
                    query_matches.emplace_back(MMMatch(tid, 0, qpos - tpos, tpos));
                } 
                else{
                    query_matches.emplace_back(MMMatch(tid, 1, qpos + tpos, tpos));
                }

              //if(qname == "NC_015678.1_823046_823145_2:0:0_1:0:0_34c/1")    std::cout<<id<<"\t"<<qid<<"\t"<<tid<<"\n";

            }          
        }
        catch(const std::out_of_range& oor) {};
    }
    std::sort(query_matches.begin(), query_matches.end(), qcomparer());

    int b=0;
  

    if(query_matches.size()>0){
        for(int e = 0; e < query_matches.size()-1; ++e){
            if(e == query_matches.size()-1 | 
                query_matches[e+1].id != query_matches[e].id | 
                query_matches[e+1].strand != query_matches[e].strand | 
                query_matches[e+1].start - query_matches[e].start >= epsilon)
                {   
                    //extract potential overlaps for each query match
                    potential_olap = slice(query_matches, b, e);

                     /* if(qname == "NC_015678.1_823046_823145_2:0:0_1:0:0_34c/1"){
                        for(auto &x: potential_olap){
                            if(x.id == 67945)
                        std::cout<<x.id<<"\t"<<x.start<<"\t"<<x.end<<"\t"<<x.strand<<"\n";
                        }
                        }*/

                   /* if(seqids[id] == 17707){
                        for(auto &p: potential_olap){
                            if(p.id == 65520){
                                std::cout<<p.start<<"\t"<<p.end<<"\n";
                            }
                        }
                    }*/

                    if(query_matches[e].strand == 0)
                        orientation = "+"; 
                    
                    else if (query_matches[e].strand == 1)
                        orientation = "-";                 
                
                    b = e+1;
                
                    if(orientation == "+"){
                        indices = longest_increasing_subset(potential_olap, orientation);
                    }
                    else{
                        //std::reverse(potential_olap.begin(), potential_olap.end());
                        indices = longest_increasing_subset(potential_olap, orientation); 
                        /*if(qname == "NC_015678.1_823046_823145_2:0:0_1:0:0_34c/1"){
                        for(auto &x: indices){
                            if(x.id == 67945)
                            {
                                    std::cout<<"Indices: "<<x.id<<"\t"<<x.start<<"\t"<<x.end<<"\t"<<x.strand<<"\n";
                            }
                        
                        }
                        }*/

                    }
                    
                    /*if(seqids[id] == 17707){
                        for(auto &p: indices){
                            if(p.id == 65520){
                                std::cout<<headers[id]<<"\t"<<id<<"\t"<<qid<<"\t"<<tid<<"\t"<<p.start<<"\t"<<p.end<<"\n";
                            }
                        }
                    }*/
                    
                    if(indices.size() <= 3) continue;

                    MMMatch first_mm = indices[0];
                    MMMatch last_mm = indices[indices.size()-1];
                    tname = headers[first_mm.id];
                    //tname = "test";

                   /* if(qname == "NC_015678.1_823046_823145_2:0:0_1:0:0_34c/1")
                    {
                        std::cout<<tname<<"\t"<<tid<<"\t"<<first_mm.id<<"\t"<<first_mm.end<<"\t"<<last_mm.end<<"\n";
                    }*/
                         
                    int mapped_bp = indices.size() * kmer;
                    int qlen = seqlens[id];
                    int tlen = seqlens[id];;
                   
                    /*if(orientation == "+"){
                        tstart = first_mm.end;
                        tend = last_mm.end + kmer;
                        qstart = first_mm.start + first_mm.end;     
                        qend = last_mm.start + last_mm.end + kmer;    

                    }*/
                    if(orientation == "+"){
                        tstart = first_mm.end - kmer;
                        tend = last_mm.end;
                        qstart = first_mm.start;     
                        qend = last_mm.start + last_mm.end;    
                        //if(qid == 844 && tid == 67945 )
                           // std::cout<<id<<"\t"<<tname<<"  "<<tid<<"  "<<tstart<<"  "<<tend<<"  "<<qname<<"  "<<qid<<"  "<<qstart<<"  "<<qend<<"\n"; 
                       // if(qname == "NC_015678.1_823046_823145_2:0:0_1:0:0_34c/1" )
                            //std::cout<<id<<"\t"<<tname<<"  "<<first_mm.id<<"  "<<tstart<<"  "<<tend<<"  "<<qname<<"  "<<qid<<"  "<<qstart<<"  "<<qend<<"\n"; 
                    }
                    else{
                        //if(qname == "NC_015678.1_823046_823145_2:0:0_1:0:0_34c/1" && tname == "NC_015678.1_823047_823146_0:0:0_0:0:0_10969/1"  )
                            //std::cout<<first_mm.start<<"\t"<<first_mm.end<<"\t"<<last_mm.start<<"\t"<<last_mm.end<<"\n";
                        tstart = first_mm.end - kmer;
                        tend = last_mm.end;
                        qstart = last_mm.start - last_mm.end ;
                        qend = first_mm.start - first_mm.end - kmer;    
                        
                       //if(qname == "NC_015678.1_823046_823145_2:0:0_1:0:0_34c/1" )
                            //std::cout<<id<<"\t"<<tname<<"  "<<first_mm.id<<"  "<<tstart<<"  "<<tend<<"  "<<qname<<"  "<<qid<<"  "<<qstart<<"  "<<qend<<"\n"; 
                        //std::cout<<tname<<"  "<<tid<<"  "<<tstart<<"  "<<tend<<"  "<<qname<<"  "<<qid<<"  "<<qstart<<"  "<<qend<<"\n"; 
                    }

                    int overhang = std::min(qstart,tstart) + std::min(qlen-qend, tlen-tend);
                    int maplen = std::max(std::abs(qend - qstart), std::abs(tend - tstart));
                    
                    //Do not consider overlaps lesser than MAPLEN
                    if(std::abs(maplen) < MAPLEN){
                         continue;
                    } 
                        
                    //overhang is the bases after and before the computed overlaps
                    /*if(overhang > std::min(MINHANG, int(0.8*maplen))){
                        continue;
                    }*/
                        
                    //if(qstart <= tstart && (qlen - qend) <= (tlen-tend))
                      // continue;
                   // else if(qstart >= tstart && (qlen - qend) >= (tlen-tend))
                      //  continue;
                    
                     //if(qname == "NC_015678.1_823046_823145_2:0:0_1:0:0_34c/1")
                    //{
                     //   std::cout<<"FINAL olap: "<<tname<<"\t"<<first_mm.id<<"\t"<<first_mm.id<<"\t"<<first_mm.end<<"\t"<<last_mm.end<<"\n";
                   // }

                    //Output the overlaps and store into vector
                    //outfile<<qname<<"\t"<<qlen<<"\t"<<qstart<<"\t"<<qend<<"\t"<<orientation<<"\t"<<tname<<"\t"<<tlen<<"\t"<<tstart<<"\t"<<tend<<"\n";
                    
                    MMOverlapInfo olap_info = MMOverlapInfo(qname,qstart,qend,tname,tstart,tend,orientation);
                    MMOverlaps.push_back({olap_info});
                    qmap[qid] = true;
                    
           }
        }
        //std::cout<<"End of qmatches\n";
    }
    
    return;
}  


void MMSketch::DumpAllMMIndex(const char *idx_file){
    std::ofstream out_fh(idx_file, ios_base::out | ios_base::binary);
    if(!out_fh.good()){
        std::cout<<"ERROR: Cannot write index to file "<<idx_file<<endl;
        exit(1);
    }
    size_t hashTableSize = MMHashTable.size();
    out_fh.write(reinterpret_cast<const char*>(&hashTableSize), sizeof(hashTableSize));
    //Write each key-value pair in the hash table
    for (const auto& entry : MMHashTable) {
        // Write the key
        uint32_t key = entry.first;
        out_fh.write(reinterpret_cast<const char*>(&key), sizeof(key));

        // Write the size of the vector
        size_t vectorSize = entry.second.size();
        out_fh.write(reinterpret_cast<const char*>(&vectorSize), sizeof(vectorSize));

        // Write each MMIndex in the vector
        for (const MMIndex& mm_index : entry.second) {
            out_fh.write(reinterpret_cast<const char*>(&mm_index.info_), sizeof(mm_index.info_));
        }
    }

    out_fh.close();
    MMHashTable.clear();
}


void MMSketch::LoadMMIndexFile(const std::string& idxname) {
    std::ifstream in_fh(idxname, std::ios_base::in | std::ios_base::binary);
    if (!in_fh.good()) {
        std::cout << "ERROR: Cannot read index from file " << idxname << std::endl;
        exit(1);
    }
    // Clear the existing hash table before loading
    MMHashTable.clear();

    size_t hashTableSize;
    in_fh.read(reinterpret_cast<char*>(&hashTableSize), sizeof(hashTableSize));

    // read each key-value pair in the hash table
    for (size_t i = 0; i < hashTableSize; i++) {
        uint32_t key;
        in_fh.read(reinterpret_cast<char*>(&key), sizeof(key));

        size_t vectorSize;
        in_fh.read(reinterpret_cast<char*>(&vectorSize), sizeof(vectorSize));

        std::vector<MMIndex> mm_indices;
        for (size_t j = 0; j < vectorSize; j++) {
            uint32_t mm_info;
            in_fh.read(reinterpret_cast<char*>(&mm_info), sizeof(mm_info));
            mm_indices.emplace_back(MMIndex(mm_info));
        }

        MMHashTable[key] = mm_indices;

        mm_indices.clear();
    }

    in_fh.close();
}

void MMSketch::PrintOverlapInfo(){
    if(MMOverlaps.size()>0)
    {
        cout << "qname  qstart  qend  tname tstart  tend  strand" << endl;
        for(auto& olap : MMOverlaps){
        std::cout<<olap.qname<<"\t"<<olap.qstart<<"\t"<<olap.qend<<"\t"<<olap.tname<<"\t"<<olap.tstart<<"\t"<<olap.tend<<"\t"<<olap.strand;
        std::cout<<endl;
        }
    }
    else{
        std::cout<<"ERROR: Cannot print overlaps. Please perform indexing and mapping first\n\n";
        exit(1);
    }
}


void MMSketch::WriteOverlaps(std::ofstream& outfile){
    if(MMOverlaps.size()>0)
    {
        outfile << "#qname\tqstart\tqend\ttname\ttstart\ttend\tstrand" << endl;
        for(auto& olap : MMOverlaps){
        outfile<<olap.qname<<"\t"<<olap.qstart<<"\t"<<olap.qend<<"\t"<<olap.tname<<"\t"<<olap.tstart<<"\t"<<olap.tend<<"\t"<<olap.strand;
        outfile<<endl;
        }
    }
    else{
        std::cout<<"ERROR: Cannot write overlaps. Please perform indexing and mapping first\n\n";
        exit(1);
    }
}





    

