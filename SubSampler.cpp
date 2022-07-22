#include <algorithm>
#include <atomic>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <mutex>
#include <omp.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <sys/stat.h>
#include <tmmintrin.h>
#include <unistd.h>
#include <unordered_map>
#include <vector>
#include <filesystem>
#include "SubSampler.h"
#include "utils.h"

#include "include/zstr.hpp"
#include "include/robin_hood.h"




using namespace std;



void Subsampler::updateK(uint64_t& min, char nuc) {
	min <<= 2;
	min += nuc2int(nuc);
	min %= offsetUpdateAnchor;
}



void Subsampler::updateM(uint64_t& min, char nuc) {
	min <<= 2;
	min += nuc2int(nuc);
	min %= offsetUpdateMinimizer;
}



void Subsampler::updateRCK(uint64_t& min, char nuc) {
	min >>= 2;
	min += (nuc2intrc(nuc) << (2 * k - 2));
}



void Subsampler::updateRCM(uint64_t& min, char nuc) {
	min >>= 2;
	min += (nuc2intrc(nuc) << (2 * minimizer_size - 2));
}


uint64_t Subsampler::canonize(uint64_t x, uint64_t n) {
	return min(x, rcbc(x, n));
}


uint64_t Subsampler::regular_minimizer_pos(uint64_t seq, uint64_t& position) {
	uint64_t mini, mmer;
	mmer = seq % minimizer_number;
	mini = mmer        = canonize(mmer, minimizer_size);
	uint64_t hash_mini = (unrevhash(mmer));
	position           = 0;
	for (uint64_t i(1); i <= k - minimizer_size; i++) {
		seq >>= 2;
		mmer          = seq % minimizer_number;
		mmer          = canonize(mmer, minimizer_size);
		uint64_t hash = (unrevhash(mmer));
		if (hash_mini > hash) {
			position  = k - minimizer_size - i;
			mini      = mmer;
			hash_mini = hash;
		}
	}
    return revhash((uint64_t)mini) % minimizer_number;
}

/*void Subsampler::estimate_sub_rate(const string& input_file){
	istream* input_stream = openFile(input_file);
	uint64_t count_kmer_in_max_skmer(0), nb_kmer(0);
	string ref, useless;
	uint32_t old_minimizer, minimizer;
	uint32_t cpt(0);
	while(not input_stream->eof()){
		ref = useless = "";
		getline(*input_stream, useless);
		getline(*input_stream, ref);
		if (ref.size() < k){
			ref = "";
		}
		if(not ref.empty() and not useless.empty()){
			old_minimizer = minimizer = minimizer_number;
			uint64_t last_position(0);
			uint64_t seq(str2num(ref.substr(0, k)));
			uint64_t position_min;
			uint64_t min_seq = (str2num(ref.substr(k - minimizer_size, minimizer_size))),
			min_rcseq(rcbc(min_seq, minimizer_size)),
			min_canon(min(min_seq, min_rcseq));
			minimizer         = regular_minimizer_pos(seq, position_min);
			old_minimizer     = minimizer;
			uint64_t hash_min = unrevhash(minimizer);
			uint64_t i(0);
			for (; i + k < ref.size(); ++i) {
				if(cpt >= 1000000){
					cout << "total kmer seen: " << intToString(nb_kmer) << " kmers in max skmers: " << intToString(count_kmer_in_max_skmer) << endl;
					estimated_subrate = (double)nb_kmer/count_kmer_in_max_skmer;
					return;
				}
				updateK(seq, ref[i + k]);
				updateM(min_seq, ref[i + k]);
				updateRCM(min_rcseq, ref[i + k]);
				min_canon      = (min(min_seq, min_rcseq));
				uint64_t new_h = unrevhash(min_canon);
				// THE NEW mmer is a MINIMIZor
				if (new_h < hash_min) {
					minimizer    = (min_canon);
					hash_min     = new_h;
					position_min = i + k - minimizer_size + 1;
				} else {
					// the previous minimizer is outdated
					if (i >= position_min) {
						minimizer = regular_minimizer_pos(seq, position_min);
						hash_min  = unrevhash(minimizer);
						position_min += (i + 1);
					}
				}
				// COMPUTE KMER MINIMIZER
				if (revhash(old_minimizer) % minimizer_number != revhash(minimizer) % minimizer_number) {
					old_minimizer = (revhash(old_minimizer) % minimizer_number);
					if((i - last_position + 1)==max_superkmer_size){
						count_kmer_in_max_skmer+= (i - last_position + 1);
					}
					nb_kmer += (i - last_position + 1);
					last_position = i + 1;
					old_minimizer = minimizer;
					cpt = k+nb_kmer-1;
				}
			}
			if (ref.size() - last_position > k - 1) {
				old_minimizer = (revhash(old_minimizer) % minimizer_number);
				if((ref.size() - last_position + 1)==max_superkmer_size){
						count_kmer_in_max_skmer+= (i - last_position + 1);
				}
				nb_kmer += (i - last_position + 1);
			}
		}
	}
	cout << "total kmers seen: " << intToString(nb_kmer) << "\nkmers in max skmers: " << intToString(count_kmer_in_max_skmer) << endl;
	estimated_subrate = (double)nb_kmer/count_kmer_in_max_skmer;
}*/

void Subsampler::parse_fasta_test(const string& input_file, const string& output_file) {
	uint64_t total_nuc_number(0);
    uint64_t read_kmer(0);
	string tmp;
	array<mutex, 1024> mutexArray;
	istream* input_stream = openFile(input_file);

	
	std::unique_ptr< std::ostream > out_file_skmer = std::unique_ptr< std::ostream >(new zstr::ofstream(output_file + ".gz"));
	ofstream out_file_skmer_big = ofstream(output_file);
	//robin_hood::unordered_map<uint32_t, vector<bool>> sketch;
	vector<bool>* sketch =  new vector<bool>[minimizer_number];
	//robin_hood::unordered_map<uint32_t, vector<bool>>::iterator val;
	
	omp_lock_t lock_Array[1024];
	for (uint64_t i = 0; i < 1024; i++) {
		omp_init_lock(&(lock_Array[i]));
	}
	#pragma omp parallel num_threads(coreNumber)
	{
		string ref, useless;
		uint32_t old_minimizer, minimizer;
		mutex mutexRead;
		while (not input_stream->eof()) {
			ref = useless = "";
			#pragma omp critical(dataupdate)
			{
				mutexRead.lock();
				getline(*input_stream, useless);
				getline(*input_stream, ref);
				if (ref.size() < k) {
					ref = "";
				} else {
					read_kmer += ref.size() - k + 1;
				}
				mutexRead.unlock();
			}
			// FOREACH UNITIG
			if (not ref.empty() and not useless.empty()) {
				uint64_t count_maximal_skmer(0);
				old_minimizer = minimizer = minimizer_number;
				uint64_t last_position(0);
				uint64_t seq(str2num(ref.substr(0, k)));
				uint64_t position_min;
				uint64_t min_seq = (str2num(ref.substr(k - minimizer_size, minimizer_size))),
                min_rcseq(rcbc(min_seq, minimizer_size)),
				min_canon(min(min_seq, min_rcseq));
				minimizer         = regular_minimizer_pos(seq, position_min);
				old_minimizer     = minimizer;
				uint64_t hash_min = unrevhash(minimizer);
				uint64_t i(0);
				// FOREACH KMER
				for (; i + k < ref.size(); ++i) {
					updateK(seq, ref[i + k]);
					updateM(min_seq, ref[i + k]);
					updateRCM(min_rcseq, ref[i + k]);
					min_canon      = (min(min_seq, min_rcseq));
					uint64_t new_h = unrevhash(min_canon);
					// THE NEW mmer is a MINIMIZor
					if (new_h < hash_min) {
						minimizer    = (min_canon);
						hash_min     = new_h;
						position_min = i + k - minimizer_size + 1;
					} else {
						// the previous minimizer is outdated
						if (i >= position_min) {
							minimizer = regular_minimizer_pos(seq, position_min);
							hash_min  = unrevhash(minimizer);
							position_min += (i + 1);
						} else {
						}
					}
					// COMPUTE KMER MINIMIZER
					if (revhash(old_minimizer) % minimizer_number != revhash(minimizer) % minimizer_number) {
						old_minimizer = (revhash(old_minimizer) % minimizer_number);
						#pragma omp atomic
						count_maximal_skmer++;
						if(old_minimizer <= (double)minimizer_number/subsampling_rate){
							vector<bool> skmer = str2boolv(ref.substr(last_position, ((2*k-minimizer_size)/2)-minimizer_size/2) + ref.substr(last_position + (((2*k-minimizer_size)/2)+minimizer_size/2), ((2*k-minimizer_size)/2) - minimizer_size/2));
							if(sketch[old_minimizer].empty()){
								omp_set_lock(&(lock_Array[old_minimizer%1024]));
								actual_minimizer_number++;
								sketch[old_minimizer] = skmer;
								omp_unset_lock(&(lock_Array[old_minimizer%1024]));
							}else{
								omp_set_lock(&(lock_Array[old_minimizer%1024]));
								sketch[old_minimizer].insert(sketch[old_minimizer].end(), skmer.begin(), skmer.end());
								omp_unset_lock(&(lock_Array[old_minimizer%1024]));
							}
							#pragma omp atomic
							selected_kmer_number+=(i - last_position + 1);
							#pragma omp atomic
							selected_superkmer_number++;
						}
						#pragma omp atomic
						total_kmer_number += (i - last_position + 1);
						#pragma omp atomic
                        total_superkmer_number++;
						last_position = i + 1;
						old_minimizer = minimizer;
					}
				}
				if (ref.size() - last_position > k - 1) {
					old_minimizer = (revhash(old_minimizer) % minimizer_number);
					#pragma omp atomic
					count_maximal_skmer++;
					if(old_minimizer <= (double)minimizer_number/subsampling_rate){
						vector<bool> skmer = str2boolv(ref.substr(last_position, ((2*k-minimizer_size)/2)-minimizer_size/2) + ref.substr(last_position + (((2*k-minimizer_size)/2)+minimizer_size/2), ((2*k-minimizer_size)/2) - minimizer_size/2));
						if(sketch[old_minimizer].empty()){
							omp_set_lock(&(lock_Array[old_minimizer%1024]));
							actual_minimizer_number++;
							sketch[old_minimizer] = skmer;
							omp_unset_lock(&(lock_Array[old_minimizer%1024]));
						}else{
							omp_set_lock(&(lock_Array[old_minimizer%1024]));
							sketch[old_minimizer].insert(sketch[old_minimizer].end(), skmer.begin(), skmer.end());
							omp_unset_lock(&(lock_Array[old_minimizer%1024]));
						}
						#pragma omp atomic
						selected_kmer_number+=(ref.size() - last_position + 1);
						#pragma omp atomic
						selected_superkmer_number++;
					}
					#pragma omp atomic
                    total_kmer_number += (ref.size() - last_position + 1);
					#pragma omp atomic
                    total_superkmer_number++;
				}
			}
		}
	}

	tmp = to_string(k-1+max_superkmer_size) + " " + to_string(minimizer_size) + " " + to_string(selected_kmer_number) + "\n";
	out_file_skmer->write(tmp.c_str(), tmp.size());
	out_file_skmer_big.write(tmp.c_str(), tmp.size());
	for(uint64_t i = 0; i <= minimizer_number-1; ++i){
		if(!sketch[i].empty()){
				//cout << sizeof(sketch) << endl;
				//cout << "coucou" << endl;
			tmp = to_string(i) + " " + to_string(sketch[i].size()/2) + " " + bool2strv(sketch[i]) + "\n";
			out_file_skmer->write(tmp.c_str(), tmp.size());
			//out_file_skmer_big.write(tmp.c_str(), tmp.size());
		}
	}
}


int main(int argc, char** argv) {
	char ch;
	string input, inputfof, query, output;
	uint k(31);
	uint m1(8);
	uint c(1);
    uint64_t s(10); 

	while ((ch = getopt(argc, argv, "hdag:q:k:m:n:s:t:b:e:f:i:o:")) != -1) {
		switch (ch) {
			case 'i': input = optarg; break;
			case 'f': inputfof = optarg; break;
			case 'k': k = stoi(optarg); break;
			case 'm': m1 = stoi(optarg); break;
			case 't': c = stoi(optarg); break;
            case 's': s = stoi(optarg); break;
			case 'o': output = optarg; break;
		}
	}


	if ((input == "" && inputfof == "")) {
		cout << "Core arguments:" << endl
		     << "	-i input file" << endl
			 << "	-f input file of file" << endl
			 << "	-o output file, if fof give output file of files name (.txt)" << endl
		     << "	-k kmer size used  (31) " << endl
             << "	-s subsampling used  (10) " << endl
             << "	-m minimize size used  (8) " << endl;
		return 0;
	}else{
        cout<<" I use k="<<k<<" m="<<m1<<" s="<<s<<endl;
        cout<<"Maximal super kmer are of length "<<2*k-m1<<" or "<<k-m1+1<<" kmers" <<endl;
        Subsampler ss(k,m1,s,c);
		if(input != ""){
			//ss.estimate_sub_rate(input);
			ss.parse_fasta_test(input, output);
			cout<<"I have seen "<<intToString(ss.total_kmer_number)<<" kmers and I selected "<<intToString(ss.selected_kmer_number)<<" kmers"<<endl;
			cout<<"This means a practical subsampling rate of "<<(double)ss.total_kmer_number/ss.selected_kmer_number<<endl;
			cout<<"I have seen "<<intToString(ss.total_superkmer_number)<<" superkmers and I selected "<<intToString(ss.selected_superkmer_number)<<" superkmers"<<endl;
			cout<<"This means a practical subsampling rate of "<<(double)ss.total_superkmer_number/ss.selected_superkmer_number<<endl;
			cout<<"This means a mean superkmer size of "<<(double)ss.total_kmer_number/ss.total_superkmer_number<<" kmer per superkmer in the input"<<endl;
			cout<<"This means a mean superkmer size of "<<(double)ss.selected_kmer_number/ss.selected_superkmer_number<<" kmer per superkmer in the output"<<endl;
			
			//ofstream out_file_res = ofstream("metrics.txt", ios_base::app);
			cout << "Optimal memory size of the output file should be: " << intToString((((4+4+84*(ss.selected_superkmer_number/ss.actual_minimizer_number))*ss.actual_minimizer_number)/8)/1000) << " kilo-octets." << endl;
			//OUTPUT TO FIX
			cout <<"Actual output file size is " << intToString(std::filesystem::file_size(output+".gz")/1000) << "KB" << endl;
			//out_file_res << intToString((((4+4+84*(ss.selected_superkmer_number/ss.actual_minimizer_number))*ss.actual_minimizer_number)/8)/1000) << "KB " << to_string(std::filesystem::file_size(output + ".gz")/1000) << "KB.\n";
			cout << "Minimizer number: " << ss.actual_minimizer_number << " Skmer/minimizer: " << ss.selected_superkmer_number/ss.actual_minimizer_number << endl;
			cout << "I have stored " << intToString(ss.selected_superkmer_number * (2*k-m1))<< " nucleotides in the output file containing the superkmers" << endl;

		}else{
			string output_filename, curr_filename;
			istream* fof = openFile(inputfof);
			int i(0);
			ofstream out_fof = ofstream(output);
			while(not fof ->eof()){
        		Subsampler ss(k,m1,s,c);
				getline(*fof, curr_filename);
				if(curr_filename != ""){
					output_filename = "subsampled_file_" + to_string(i) + ".fa";
					//ss.estimate_sub_rate(curr_filename);
					ss.parse_fasta_test(curr_filename, output_filename);
					out_fof << output_filename << ".gz\n";
					i++;
					cout<<"I have seen "<<intToString(ss.total_kmer_number)<<" kmers and I selected "<<intToString(ss.selected_kmer_number)<<" kmers"<<endl;
					cout<<"This means a practical subsampling rate of "<<(double)ss.total_kmer_number/ss.selected_kmer_number<<endl;
					cout<<"I have seen "<<intToString(ss.total_superkmer_number)<<" superkmers and I selected "<<intToString(ss.selected_superkmer_number)<<" superkmers"<<endl;
					cout<<"This means a practical subsampling rate of "<<(double)ss.total_superkmer_number/ss.selected_superkmer_number<<endl;
					cout<<"This means a mean superkmer size of "<<(double)ss.total_kmer_number/ss.total_superkmer_number<<" kmer per superkmer in the input"<<endl;
					cout<<"This means a mean superkmer size of "<<(double)ss.selected_kmer_number/ss.selected_superkmer_number<<" kmer per superkmer in the output"<<endl;
					
					//ofstream out_file_res = ofstream("metrics.txt", ios_base::app);
					cout << "Optimal memory size of the output file should be: " << intToString((((4+4+84*(ss.selected_superkmer_number/ss.actual_minimizer_number))*ss.actual_minimizer_number)/8)/1000) << " kilo-octets." << endl;
					//OUTPUT TO FIX
					cout <<"Actual output file size is " << intToString(std::filesystem::file_size(output_filename+".gz")/1000) << "KB" << endl;
					//out_file_res << intToString((((4+4+84*(ss.selected_superkmer_number/ss.actual_minimizer_number))*ss.actual_minimizer_number)/8)/1000) << "KB " << to_string(std::filesystem::file_size(output + ".gz")/1000) << "KB.\n";
					cout << "Minimizer number: " << ss.actual_minimizer_number << " Skmer/minimizer: " << ss.selected_superkmer_number/ss.actual_minimizer_number << endl;
					cout << "I have stored " << intToString(ss.selected_superkmer_number * (2*k-m1))<< " nucleotides in the output file containing the superkmers" << endl;

				}

			}
		}
        
    }
	    
}