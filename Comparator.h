#include <stdio.h>
#include <fstream>
#include <string>
#include <iostream>
#include <vector>

#include "include/robin_hood.h"



using namespace std;

class Comparator{
    public:
        uint64_t skmerWoM_size, k, m, nb_kmer_tot, nb_files_eof;
        uint32_t nb_kmer_seen;
        vector<uint64_t> minimizers, nb_kmer_seen_infile;
        bool run;
        robin_hood::unordered_map<string, int> score_map;
        robin_hood::unordered_map<uint64_t, uint64_t> color_map;
        Comparator(){
            skmerWoM_size = m = nb_kmer_seen = nb_kmer_tot = k = nb_files_eof= 0;
            run = true;
        }


        void compare_files(const string& fileofile);
        vector<uint64_t> findMin(const vector<uint64_t>& minims);
        void increment_files(const vector<istream*>& files, vector<uint64_t> indices);
        void update_colormap(const vector<istream*>& files, vector<uint64_t> indices);
        void compare_buckets(const string& fileofile);
        void skip_bucket(const vector<istream*>& files, vector<uint64_t> indices);
        //string find(istream* file2, uint32_t minimizer);

};
string print_color(uint64_t c, int n);
