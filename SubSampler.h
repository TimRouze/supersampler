#ifndef SUB
#define SUB


#include <stdio.h>
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <atomic>
#include <mutex>
#include <stdint.h>
#include <unordered_map>
#include <pthread.h>
#include <filesystem>



using namespace std;



class Subsampler {
  public:
  	uint64_t k, minimizer_size;
    uint64_t minimizer_number;
    uint64_t actual_minimizer_number;
	uint64_t coreNumber;
    uint64_t total_kmer_number;
    uint64_t total_superkmer_number;
    uint64_t selected_kmer_number;
    uint64_t selected_superkmer_number;
    uint64_t offsetUpdateAnchor,offsetUpdateMinimizer;
    uint64_t subsampling_rate;
    uint64_t max_superkmer_size;
    double estimated_subrate;
    Subsampler(uint64_t ik, uint64_t i_minimizer,uint64_t isubsampling_rate,uint64_t icore){
        k=ik;
        minimizer_size=i_minimizer;
        coreNumber=icore;
        minimizer_number=(uint64_t)1<<(2*minimizer_size);
        offsetUpdateAnchor=(uint64_t)1<<(2*k);
        offsetUpdateMinimizer=minimizer_number;
        subsampling_rate=isubsampling_rate;
        total_kmer_number=selected_kmer_number=selected_superkmer_number=total_superkmer_number=0;
        max_superkmer_size=k-minimizer_size+1;
        estimated_subrate = actual_minimizer_number = 0;
        
    }
    void parse_fasta(const string& input_file);
    void parse_fasta_test(const string& input_stream);
    void updateK(uint64_t & min, char nuc);
	void updateRCK(uint64_t& min, char nuc);
	void updateM(uint64_t& min, char nuc);
	void updateRCM(uint64_t& min, char nuc);
    uint64_t regular_minimizer_pos(uint64_t seq, uint64_t& position);
    uint64_t canonize(uint64_t x, uint64_t n);
    void estimate_sub_rate(const string& input_file);
    void compress_files(const string & file, const string& output_file);
    void store_kmers(const string& input_file);

};

 #endif