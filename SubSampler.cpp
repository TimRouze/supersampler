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
#include "Decycling.h"
#include "utils.h"

#include "include/zstr.hpp"
#include "include/robin_hood.h"
#include "include/unordered_dense.h"



using namespace std;



void Subsampler::updateK(kmer& min, char nuc) {
	min <<= 2;
	min += nuc2int(nuc);
	min &= offsetUpdateAnchor;
}



void Subsampler::updateM(uint64_t& min, char nuc) {
	min <<= 2;
	min += nuc2int(nuc);
	min &= offsetUpdateMinimizer;
}



void Subsampler::updateRCK(kmer& min, char nuc) {
	min >>= 2;
	min += (nuc2intrc(nuc) << (2 * k - 2));
}



void Subsampler::updateRCM(uint64_t& min, char nuc) {
	min >>= 2;
	min += (nuc2intrc(nuc) << (2 * minimizer_size - 2));
}



//~ uint64_t Subsampler::unrevhash(uint64_t x){
	//~ uint64_t result = XXHash64::hash(&x, sizeof(x), 1312);
	//~ if(velo->mem(x)){
		//~ return result&offsetUpdateMinimizer;
	//~ }else{
		//~ return result|minimizer_number;
	//~ }
//~ }

uint64_t Subsampler::unrevhash(uint64_t x){
	uint64_t result = XXHash64::hash(&x, sizeof(x), 1312);
	return result;
	result&=mask;
	uint classe(velo->memDouble(x));
	if(classe==2){
		return result;
	}else if(classe==1){
		return result|second1;
	}
	return result|first1;
}






uint64_t Subsampler::regular_minimizer_pos(kmer seq, uint64_t& position, bool& is_rev, bool& multiple_min) {
	multiple_min=false;
	is_rev = false;
	uint64_t mini, mmer, canon_mmer;
	mmer = seq & offsetUpdateMinimizer;
	mini        = canonize(mmer, minimizer_size);
	if(mini != mmer){
		is_rev = true;
	}
	mmer = mini;
	uint64_t hash_mini = (unrevhash(mmer));
	position = 0;
	for (uint64_t i(1); i <= k - minimizer_size; i++) {
		bool local_rev;
		seq >>= 2;
		mmer = seq & offsetUpdateMinimizer;
		canon_mmer = canonize(mmer, minimizer_size);
		if(canon_mmer == mmer){
			local_rev = false;
		}else{
			local_rev = true;
		}
		mmer = canon_mmer;
		uint64_t hash = (unrevhash(mmer));
		if (hash_mini > hash) {
			position  = k - minimizer_size - i;
			mini      = mmer;
			is_rev = local_rev;
			hash_mini = hash;
			multiple_min = false;
		}else if(mmer == mini){
			multiple_min=true;
		}
	}
    return mini;
}



string extract_name(const string& str){
    string result;
    uint64_t begin(0);
    for(uint i(0);i<str.size();++i){
        if(str[i]=='/'){
            begin=i+1;
        }
    }
    for(uint i(begin);i<str.size();++i){
        if(str[i]=='.'){
            return result;
        }else{
            result.push_back(str[i]);
        }
    }
    return result;
}



string get_out_name(const string& str, const string& prefix){
    string result, path;
    uint64_t begin(0);
	path = "";
    for(uint i(0);i<str.size();++i){
        if(str[i]=='/'){
            begin=i+1;
			path = str.substr(0, i+1);
        }
    }
    for(uint i(begin);i<str.size();++i){
        if(str[i]=='.'){
            return prefix+result;
        }else{
            result.push_back(str[i]);
        }
    }
    return prefix + result;
}



void Subsampler::handle_superkmer(string& superkmer,map<uint32_t,pair<vector<bool>,string>>& sketch_max,kmer input_minimizer, bool inputrev, bool inputmultiple){
    selected_superkmer_number++;

    if(inputrev){
        superkmer=revComp(superkmer);
    }
    if(inputmultiple && type == 3){
        //MULTIPLE MINIMIZER SUPERKMER
        string kmerstr;
        for(uint i(0);i+k<=superkmer.size();++i){
            bool multiple;
            kmerstr=superkmer.substr(i, k);
            kmer seq(str2num(kmerstr));
            uint64_t position_min;
            kmer minimizer= regular_minimizer_pos(seq, position_min, inputrev,multiple);

            //REGULAR KMER FOUND INSIDE MULTIPLE MINIMIZER SUPERKMER
                uint64_t p =kmerstr.find(num2str(input_minimizer,minimizer_size));
                if(p==string::npos){
                    kmerstr=revComp(kmerstr);
                    p =kmerstr.find(num2str(input_minimizer,minimizer_size));
                }
                kmerstr = kmerstr.substr(0, p)+"\n"+kmerstr.substr(p+minimizer_size)+"\n";
                sketch_max[input_minimizer].second+=kmerstr;
                selected_kmer_number++;
				cursed_kmer_number++;
        }
    }else{
        selected_kmer_number+=superkmer.size()-k+1;
        if(superkmer.size()==2*k-minimizer_size){
            //Maximal SUPERKMER
            vector<bool> skmer;
            count_maximal_skmer++;
            skmer = str2boolv(superkmer.substr(0, k-minimizer_size) + superkmer.substr(0 + k, k-minimizer_size));
            sketch_max[input_minimizer].first.insert(sketch_max[input_minimizer].first.end(), skmer.begin(), skmer.end());
        }else if(type != 1){
            uint64_t p =superkmer.find(num2str(input_minimizer,minimizer_size));
            string skmerstr = superkmer.substr(0, p)+"\n"+superkmer.substr(p+minimizer_size)+"\n";
            sketch_max[input_minimizer].second+=skmerstr;
        }
    }
}



// RENAME
void Subsampler::parse_fasta_test(const string& input_file, const string& output_prefix) {
    total_kmer_number=actual_minimizer_number=total_superkmer_number=selected_kmer_number=selected_superkmer_number=count_maximal_skmer=0;
    uint64_t read_kmer(0);
	string tmp;
	array<mutex, 1024> mutexArray;
	zstr::ifstream* input_stream = openFile(input_file);
    if(input_stream==NULL){
        cout<<"Can't open file: "<<input_file<<endl;
        return;
    }
    if(not input_stream->good()){
        cout<<"Can't open file: "<<input_file<<endl;
        return;
    }
    //string clean_input_file=extract_name(input_file);
    //subsampled_file=output_prefix +clean_input_file+".gz";
	subsampled_file = get_out_name(input_file, output_prefix)+".gz";
	zstr::ofstream* out_file_skmer = (new zstr::ofstream(subsampled_file,21,9));
	map<uint32_t,pair<vector<bool>,string>> sketch_max;
	nb_mmer_selected = 0;
	{
		string ref, useless,superstr,skmerstr;
		uint32_t old_minimizer, minimizer;
		mutex mutexRead;
        vector<bool> skmer;
		while (not input_stream->eof()) {
			ref = useless = "";
			{
                Biogetline(input_stream,ref,'A',k);
				if (ref.size() < k) {
					ref = "";
				} else {
					read_kmer += ref.size() - k + 1;
				}
			}
			// FOREACH sequence
			if (not ref.empty()) {
				uint64_t skmer_size(0);
				bool is_rev, old_rev,multiple_min,old_multiple_min;
				old_minimizer = minimizer = minimizer_number;
				uint64_t last_position(0);
				kmer seq(str2num(ref.substr(0, k)));
				uint64_t position_min;
				uint64_t min_seq = (str2num(ref.substr(k-minimizer_size, minimizer_size))),
                min_rcseq(rcbc(min_seq, minimizer_size)),
				min_canon(min(min_seq, min_rcseq));
				minimizer         = regular_minimizer_pos(seq, position_min, old_rev,old_multiple_min);
				old_minimizer     = minimizer;
				uint64_t hash_min = unrevhash(minimizer);
				uint64_t i(0);
				uint64_t pos_end(0);
				// FOREACH KMER
				for (; i + k < ref.size(); ++i) {
					updateK(seq, ref[i + k]);
					updateM(min_seq, ref[i + k]);
					updateRCM(min_rcseq, ref[i + k]);
					min_canon      = (min(min_seq, min_rcseq));
					uint64_t new_h = unrevhash(min_canon);
					if (new_h < hash_min) {
						// THE NEW mmer is a MINIMIZER
						minimizer    = min_canon;
						multiple_min=false;
						hash_min     = new_h;
						position_min = i + k - minimizer_size + 1;
						if(min_canon == min_seq){
							is_rev = false;
						}else{
							is_rev = true;
						}
					//MAX_SUPERKMER_SIZE != 1 -> PREVENTS SKMERS SIZES GREATER THAN MAX SKMER SIZE WHEN K = M
					}else if (old_minimizer == min_canon && max_superkmer_size != 1){
						// DUPLICATE MINIMIZER
						old_multiple_min=true;
					}else {
						if (i >= position_min) {
							// the previous minimizer is outdated
							minimizer = regular_minimizer_pos(seq, position_min, is_rev, multiple_min);
							hash_min  = unrevhash(minimizer);
							position_min += (i + 1);
						}
					}

					if(old_minimizer != minimizer){
						//THE MINIMIZER CHANGED WE MUST HANDLE THE ASSOCIATED SUPERKMER
                    	if(unrevhash(old_minimizer) <= selection_threshold){
                            // pos_end: end of previous skmer (inclusive)
                            // last_position: start of current skmer
                            // i: start of last k-mer in current skmer
							if(last_position + minimizer_size - 2 > pos_end){
                                // at least one m-mer is missing
                                if (pos_end > 0) {
                                    nb_mmer_selected -= minimizer_size - 1; // S - m + 1
                                }
								nb_mmer_selected += i + k - last_position; // add curr skmer size
                                nb_mmer_selected -= k - minimizer_size; // remove start of curr skmer
							}
							else{
                                // add skmer size without overlap
								nb_mmer_selected += i + k - (pos_end + 1);
							}
						//if(unrevhash(old_minimizer)%subsampling_rate == 0){
                            //THE MINIMIZER IS SELECTED WE MUST OUTPUT THE ASSOCIATED SUPERKMER
                            handle_superkmer(superstr=ref.substr(last_position, i+k-last_position),sketch_max,old_minimizer,old_rev,old_multiple_min);
							pos_end = i + k - 1;
                        }
						total_kmer_number += (i - last_position + 1);
                        total_superkmer_number++;
						last_position = i + 1;
						old_minimizer = minimizer;
						old_rev = is_rev;
                        old_multiple_min=multiple_min;
						skmer_size = k;
					}else{
						skmer_size++;
					}
				}

				if (ref.size() - last_position > k - 1) {
					if(unrevhash(old_minimizer) <= selection_threshold){
						nb_mmer_selected -= minimizer_size - 1;
					//if(unrevhash(old_minimizer)%subsampling_rate == 0){
                        //THE MINIMIZER IS SELECTED WE MUST OUTPUT THE ASSOCIATED SUPERKMER
						handle_superkmer(superstr=ref.substr(last_position, i+k-last_position),sketch_max,old_minimizer,old_rev,old_multiple_min);
                        pos_end = i + k - 1;
					}
                    total_kmer_number += (i - last_position + 1);
                    total_superkmer_number++;
					skmer_size = k;
				}
			}
		}
	}
	nb_mmer_selected -= minimizer_size-1;
	tmp = to_string(k-1+max_superkmer_size) + " " + to_string(minimizer_size) + " " + to_string(selected_kmer_number) + " " + to_string(subsampling_rate) + "\n";
	out_file_skmer->write(tmp.c_str(), tmp.size());
	//HERE i is the bucket number
	for(auto const& [i,val]: sketch_max){
		string minstr=num2str(i,minimizer_size);
		out_file_skmer->write(minstr.c_str(), minimizer_size);
		string compressed(strCompressor(bool2strv(val.first)));
        uint32_t size_compressed(compressed.size());//TODO RISKY 16bit int
		out_file_skmer->write((const char*)&size_compressed, sizeof(size_compressed));
		out_file_skmer->write(compressed.c_str(), compressed.size());
        out_file_skmer->write(val.second.c_str(), val.second.size());
		out_file_skmer->write("\n\n", 2);

	}
    actual_minimizer_number=sketch_max.size();
	delete input_stream;
    delete out_file_skmer;
}



uint64_t Subsampler::compute_threshold(double sampling_rate){
	//UNCOMMENT THIS LINE TO HAVE THE SAME THRESHOLD AS SOURMASH
	//return (uint64_t)((pow(2, 64)-1)/sampling_rate);
    uint64_t mmerinkmer(k-minimizer_size+1);
    long double fraction_sampling=(long double)1/sampling_rate;
    long double root(pow((long double)1-fraction_sampling,(long double)1/mmerinkmer));
    long double result=((long double)1-root)*((uint64_t)1<<63);
    return (uint64_t) result*2;
}



void Subsampler::print_stat(){
    if(selected_kmer_number!=0){
        cout<<"I have seen "<<intToString(total_kmer_number)<<" kmers and I selected "<<intToString(selected_kmer_number)<<" kmers"<<endl;
        cout<<"This means a practical subsampling rate of "<<(double)total_kmer_number/selected_kmer_number<<endl;
        cout<<"I have seen "<<intToString(total_superkmer_number)<<" superkmers and I selected "<<intToString(selected_superkmer_number)<<" superkmers"<<endl;
        cout<<"This means a practical subsampling rate of "<<(double)total_superkmer_number/selected_superkmer_number<<endl;
        cout<<"This means a mean superkmer size of "<<(double)total_kmer_number/total_superkmer_number<<" kmer per superkmer in the input"<<endl;
        cout<<"This means a mean superkmer size of "<<(double)selected_kmer_number/selected_superkmer_number<<" kmer per superkmer in the output"<<endl;

        cout <<"Actual output file size is " << intToString(std::filesystem::file_size(subsampled_file)/1000) << "KB" << endl;
        cout <<"This mean " << ((double)std::filesystem::file_size(subsampled_file)*8/selected_kmer_number) << "bits per kmer" << endl;
        cout << "Minimizer number: " << actual_minimizer_number << " Skmer/minimizer: " << selected_superkmer_number/actual_minimizer_number << endl;
        cout << "Cursed kmers: " << cursed_kmer_number<<endl;
		cout << "Density is: " << (((double)selected_superkmer_number/nb_mmer_selected)*(k-minimizer_size+2)) << endl; // d * (w+1)
		cout << "Number of maximal skmer is: " << count_maximal_skmer << endl;
		cout << "Proportion of max skmers: " << ((double)count_maximal_skmer/selected_superkmer_number)*100 << "%" << endl;
        }else{
        cout<<"No kmer selected ***Crickets noise***"<<endl;
        }
}



int main(int argc, char** argv) {
	char ch;
	string input, inputfof, query, output("subsampled_");
	uint k(31);
	uint m1(11);
	uint c(8);
    double s(1000);
    bool verbose=true;
	uint type(3);

	while ((ch = getopt(argc, argv, "hdag:q:k:m:n:s:t:b:e:f:i:p:v:x:")) != -1) {
		switch (ch) {
			case 'i': input = optarg; break;
			case 'f': inputfof = optarg; break;
			case 'k': k = stoi(optarg); break;
			case 'm': m1 = stoi(optarg); break;
			case 't': c = stoi(optarg); break;
            case 's': s = stof(optarg); break;
			case 'p': output = optarg; break;
            case 'v': verbose = stoi(optarg); break;
            case 'x': type = stoi(optarg); break;
		}
	}
	if ((input == "" && inputfof == "")) {
		cout << "Core arguments:" << endl
		     << "	-i Input file" << endl
			 << "	-f Input file of file" << endl
			 << "	-p Output prefix (subsampled)" << endl
		     << "	-k Kmer size used  (31) " << endl
             << "	-s Subsampling used  (1000) " << endl
             << "	-t Threads used  (8) " << endl
             << "	-m Minimizer size used  (11) " << endl
             << "	-v Verbose level (1) " << endl
			 << "	-3/2/1 respectively Max skmers + any sized skmers + cursed skmers OR Max skmers and any sized skmers OR max skmers only. (default 3) " << endl
             ;
		return 0;
	}else{
        if(m1%2==0){
            cout<<"Minimizer size must be odd"<<endl;
            m1++;
        }
        if(k%2==0){
            cout<<"Kmer size must be odd"<<endl;
            k++;
        }
        cout<<" I use k="<<k<<" m="<<m1<<" s="<<s<<endl;
        cout<<"Maximal super kmer are of length "<<2*k-m1<<" or "<<k-m1+1<<" kmers" <<endl;
		if(input != ""){
            Subsampler ss(k,m1,s,c,type);
			ss.parse_fasta_test(input, output);
			if(verbose){
                #pragma omp critical (cout)
                {
                    ss.print_stat();
                }
            }
		}else{
			zstr::ifstream* fof = openFile(inputfof);
			string out_fof_name = get_out_name(inputfof, output);
			ofstream out_fof;
			out_fof.open(out_fof_name+".txt");
            if(fof==NULL){
                cout<<"Can't open file of file "<<inputfof<<endl;
            }
			#pragma omp parallel num_threads(c)
			{
				string curr_filename;
				while(not fof ->eof()){
					#pragma omp critical (fof)
					{
						getline(*fof, curr_filename);
					}
					if(curr_filename.size()>3){
                        #pragma omp critical (cout)
                        {
                            cout<<curr_filename<<endl;
							out_fof << get_out_name(curr_filename, output)+".gz\n";
                        }
                        Subsampler ss(k,m1,s,c,type);
						ss.parse_fasta_test(curr_filename, output);
                        if(verbose){
                            #pragma omp critical (cout)
                            {
                                ss.print_stat();
                            }
                        }
					}
				}
			}
            delete fof;
			out_fof.close();
		}
    }
}
