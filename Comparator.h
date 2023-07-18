#include <stdio.h>
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <iomanip>

#include "include/robin_hood.h"
#include "utils.h"



using namespace std;



class Comparator{
    public:
        uint64_t skmer_size, k, m, nb_kmer_tot, nb_files_eof,nb_files,query_size;
        uint64_t nb_kmer_seen,precision, sub_rate;
        double min_threshold;

        vector<uint64_t> minimizers, nb_kmer_seen_infile;
        bool run;
        vector<string> files_names;
        ankerl::unordered_dense::map<uint32_t, uint32_t> score_A;
        //zstr::ofstream* kmers_comp;
        Comparator(uint p,double mt){
            skmer_size = m = nb_kmer_seen = nb_kmer_tot = k =nb_files= nb_files_eof= 0;
            run = true;
            precision=p;
            min_threshold=mt;
        }
        void print_containment(const string& outfile);
        void print_jaccard(const string& outfile);
        void compare_files(const string& fileofile);
        bool findMin(const vector<uint64_t>& minims,vector<uint64_t>& min_vector);
        void increment_files(const vector<istream*>& files, const vector<uint64_t>& indices);
        void count_intersection(const vector<istream*>& files,const vector<uint64_t>& indices,const string& strminimizer);//, zstr::ofstream* out_kmer);
        void compare_sketches(uint size_query);
        void skip_bucket(const vector<istream*>& files,const vector<uint64_t>& indices,const string& strminimizer);//, zstr::ofstream* out_kmer);
        void compute_scores(const ankerl::unordered_dense::map<kmer,  vector<bool>>& color_map,vector<kmer>& interesting_hits);
        void getfilesname(const string& fof, vector<string>& result);
        void get_header_info(const vector<istream*>& files);
        string inject_minimizer(const string* str, const string& minimizerstr);
};


void nuc422nuc(string& str);
string print_color(uint64_t c, int n);
string kmer2str(kmer num, uint l);
