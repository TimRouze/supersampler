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
#include "Decycling.h"
#include "utils.h"



using namespace std;



class Subsampler {
  public:
  //CONSTANTS
  	uint64_t k, minimizer_size;
    uint64_t minimizer_number;
	uint64_t coreNumber;
    kmer  offsetUpdateAnchor,offsetUpdateMinimizer;
    double  subsampling_rate;
    uint64_t max_superkmer_size;
    uint64_t selection_threshold;
    uint type;
    //VARIABLES
    uint64_t total_kmer_number;
    uint64_t cursed_kmer_number;
    uint64_t actual_minimizer_number;
    uint64_t total_superkmer_number;
    uint64_t selected_kmer_number;
    uint64_t selected_superkmer_number;
    uint64_t count_maximal_skmer;
    uint64_t nb_mmer_selected;
    uint64_t first1;
    uint64_t second1;
    uint64_t mask;
    string subsampled_file;
    DecyclingSet* velo;
    Subsampler(uint64_t ik, uint64_t i_minimizer,double isubsampling_rate,uint64_t icore, uint itype){
		velo= new DecyclingSet(i_minimizer);
        k=ik;
        minimizer_size=i_minimizer;
        coreNumber=icore;
        first1=(uint64_t)1<<63;
        second1=(uint64_t)1<<62;
        mask=second1-1;
        minimizer_number=(uint64_t)1<<(2*minimizer_size);
        offsetUpdateAnchor=((kmer)1<<(2*k))-1;
        offsetUpdateMinimizer=minimizer_number-1;
        subsampling_rate=isubsampling_rate;
        cursed_kmer_number=count_maximal_skmer=total_kmer_number=selected_kmer_number=selected_superkmer_number=total_superkmer_number=0;
        max_superkmer_size=k-minimizer_size+1;
        type = itype;
        if(subsampling_rate>1){
            selection_threshold=compute_threshold(subsampling_rate);
        }else{
            selection_threshold=-1;
        }
        //~ selection_threshold=-1;
        actual_minimizer_number = 0;
    }
    void parse_fasta(const string& input_file);
    void parse_fasta_test(const string& input_file, const string& prefix);
    void updateK(kmer & min, char nuc);
	void updateRCK(kmer& min, char nuc);
	void updateM(uint64_t& min, char nuc);
	void updateRCM(uint64_t& min, char nuc);
    uint64_t regular_minimizer_pos(kmer seq, uint64_t& position, bool& is_rev, bool& multiple_min);
    void handle_superkmer(string& superkmer,map<uint32_t,pair<vector<bool>,string>>& sketch_max,kmer input_minimizer, bool inputrev, bool inputmultiple);
    void store_kmers(const string& input_file);
    uint64_t compute_threshold(double sampling_rate);
    void print_stat();
    uint64_t unrevhash(uint64_t x);
};



 #endif
