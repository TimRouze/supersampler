#ifndef UTIL
#define UTIL
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <atomic>
#include <mutex>
#include <stdint.h>
#include <unordered_map>
#include <pthread.h>
#include <chrono>
#include <omp.h>
#include <tmmintrin.h>
#include <math.h>
#include <algorithm>
#include <sys/stat.h>
#include <unistd.h>
#include "include/robin_hood.h"
#include "include/unordered_dense.h"
#include "include/zstr.hpp"
#include "include/MurmurHash3.h"
#include "include/xxhash64.h"



#define kmer __uint128_t


template <>
struct ankerl::unordered_dense::hash<__uint128_t> {
    using is_avalanching = void;

    [[nodiscard]] auto operator()(__uint128_t const& x) const noexcept -> uint64_t {
        return detail::wyhash::hash(x>>64)+detail::wyhash::hash((uint64_t)x);
    }
};


kmer hash64shift(kmer key);

// Hash function for kint in robin_hood
namespace robin_hood {
 template <> struct hash<__uint128_t>{
    size_t operator()(const __uint128_t & x) const{
      return ((hash64shift(x)) ^ hash64shift(x>>64));
    }
  };
}



using namespace std;



kmer nuc2int(char c);
kmer nuc2intrc(char c);
char int2nuc(unsigned char n);
void updateK(kmer& min, char nuc, uint64_t& k);
void read_vector_bool(vector<bool>& V, zstr::ifstream* out, uint64_t n_bits);
void dump_vector_bool(const vector<bool>& V, ostream* out);
string intToString(uint64_t n);
bool kmer_in_superkmer(const kmer canon, const vector<kmer>& V);
kmer min_k(const kmer& k1, const kmer& k2);
kmer str2num(const string& str);
uint64_t revhash(uint64_t x);
uint16_t parseCoverage(const string& str);
uint16_t parseCoverage_log2(const string& str);
uint16_t parseCoverage_exact(const string& str);
uint16_t parseCoverage_bool(const string& str);
string color_coverage2str(const vector<uint16_t>& V);
vector<string> split(const string& s, char delim);
bool exists_test(const string& name);
__m128i mm_bitshift_right(__m128i x, unsigned count);
uint64_t rcbc(uint64_t in, uint64_t n);
string revComp(const string& s);
void decompress_file(const string& file, const string& output_file);
vector<bool> str2boolv(const string& str);
string bool2strv(const vector<bool>& v);
void split(const string& s, char delim, vector<string>& res);
zstr::ifstream* openFile(const string& input_file);
string num2str(kmer num,uint k);
void print_bin(kmer n);
string strCompressor(const string& str);
string strDecompressor(const string* str);

vector<string> splitSTR(const string& s, char delim);
void Biogetline(zstr::ifstream* in,string& result,char type,uint K);
void clean_dna(string& str);
string getLineFasta(zstr::ifstream* in);
template<typename T>
inline T xs(const T& x) {
	return unrevhash(x);
}
uint64_t canonize(uint64_t x, uint64_t n);
__uint128_t canonize(__uint128_t x, uint64_t n);



#endif
