#include "Comparator.h"
#include "utils.h"



void Comparator::compare_buckets(const string& fileofile){
    uint64_t i(0), cpt(0), skmer_size;
    vector<istream*> input_files;
    vector<uint64_t> indices;
    vector<string> lines;
    string header, curr_filename;

    istream* fof = openFile(fileofile);
    while(not fof ->eof()){
        getline(*fof, curr_filename);
        if(curr_filename != ""){
            input_files.push_back(openFile(curr_filename));
            
        }
    }
    cout << input_files.size() << endl;
    nb_kmer_seen_infile.resize(input_files.size(), 0);
    minimizers.resize(input_files.size());
    vector<uint64_t> min_vector(input_files.size());
    for(int f = 0; f < input_files.size(); ++f){
        getline(*input_files[f], header);
        skmer_size = stoi(header.substr(0, header.find(" ")));
        header.erase(0, header.find(" ") + 1);
        m = stoi(header.substr(0, header.find(" ")));
        header.erase(0, header.find(" ") + 1);
        nb_kmer_tot = stoi(header);
        skmerWoM_size = skmer_size - m;
        k = ((skmer_size + m)/2) - m;
        cout <<  " There are " << intToString(nb_kmer_tot) << " " << k+m << "-mers in total in file " << f << endl;

    }
    
    cout << "kmers evaluated are of length: " << k << " minimizer size is " << m << endl; 
    indices.clear();
    increment_files(input_files, indices);
    while(run){
        indices = findMin(minimizers);
        if(indices.size() == 1){
            skip_bucket(input_files, indices);
            increment_files(input_files, indices);
        }
        else{
            update_colormap(input_files, indices);
            increment_files(input_files, indices);
        }
        indices.clear();
    }
    int file_n(0);
    for(auto elem: nb_kmer_seen_infile){
        cout << "unique kmers in file " << file_n << ": "<< intToString(elem) << endl;
        file_n++;
    }
}

void Comparator::skip_bucket(const vector<istream*>& files, vector<uint64_t> indices){
    robin_hood::unordered_map<uint64_t, uint64_t> skip_map;
    string skip;
    for(auto ind : indices){
        getline(*files[ind], skip);
        if (skip.size() < k) {
            skip = "";
        }
        if (not skip.empty()){
            uint64_t seq(str2num(skip.substr(0, k)));
            uint64_t i(0);
            while(i < skip.size()){
                while((i + k)%skmerWoM_size != 1) {
                    if(skip_map.count(seq) == 0){
                        nb_kmer_seen++;
                        skip_map[seq] = ((uint64_t)1 << ind);
                    }
                    updateK(seq, skip[i + k], k);
                    i++;
                }
                i = i+k;
            }
            nb_kmer_seen_infile[ind] += skip_map.size();
            score_A[ind*files.size()+ind] += skip_map.size();
            skip_map.clear();
        }
    }
}

void Comparator::update_colormap(const vector<istream*>& files, vector<uint64_t> indices){
    string ref;
    for(auto ind : indices){
        getline(*files[ind], ref);
        if (ref.size() < k) {
            ref = "";
        }
        if (not ref.empty()){
            uint64_t seq(str2num(ref.substr(0, k)));
            uint64_t i(0);
            while(i < ref.size()){
                while((i + k)%skmerWoM_size != 1) {
                    if(color_map.count(seq) == 0){
                        color_map[seq] = ((uint64_t)1 << ind);
                        nb_kmer_seen++;
                    }
                    else{
                        color_map[seq] |= ((uint64_t)1 << ind);
                    }
                    updateK(seq, ref[i + k], k);
                    i++;
                }
                i = i+k;
            }
            //cout << "indice: " << ind << " nb kmers ajoutés: " << color_map.size() << endl;
            //nb_kmer_seen_infile[ind] += color_map.size();
        }
    }
    compute_scores(files);
}

void Comparator::compute_scores(const vector<istream*>& files){
    string key = "";
    vector<uint64_t> score_v;
    uint64_t tmp;
    for (auto it : color_map){
        score_v.clear();
        tmp = it.second;
        for(int i(0);i<files.size();++i){
            score_v.push_back(tmp%2);
            if(tmp%2 == 1){
                nb_kmer_seen_infile[i]++;
            }
            tmp>>=1;
        }
        //score array.
        for(uint64_t i = 0; i < files.size(); ++i){
            for(uint64_t j = 0; j < files.size(); ++j){
                if(score_v[i] == score_v[j] && score_v[i] == 1){
                    score_A[i*files.size()+j] +=1;
                    
                }
            }
        }
    }
    color_map.clear();
}

string print_color(uint64_t c, int n){
	string res="";
	for(int i(0);i<n;++i){
		cout<<c%2;
		res += c%2;
		c>>=1;
	}
	cout<<" "<<flush;
	return res;
}

void Comparator::increment_files(const vector<istream*>& files, vector<uint64_t> indices){
    vector<string> lines;
    string minim, bucket_size;
    char sep = ' ';
    if(indices.size() != 0){
        for(uint64_t i(0); i < indices.size(); ++i){
            if(not files[indices[i]]->eof()){
                getline(*files[indices[i]], minim, sep);
                if(not minim.empty()){
                    minimizers[indices[i]] = stoi(minim);
                    getline(*files[indices[i]], bucket_size, sep);
                }
            }
            else{
                cout << minimizers << endl;
                minimizers[indices[i]] = -1;
                nb_files_eof++;
                cout << nb_files_eof << endl;

            }
        }
        if(nb_files_eof == files.size()){
            run = false;
        }
    }
    else{
        for(uint64_t i(0); i < files.size(); ++i){
            getline(*files[i], minim, sep);
            minimizers[i] = stoi(minim);
            getline(*files[i], bucket_size, sep);
        }
    }
}

vector<uint64_t> Comparator::findMin(const vector<uint64_t>& minims){
    vector<uint64_t> min_vector(minims.size(), -1);
    uint64_t min(-1);
    uint64_t i(0);
    for(; i < minims.size(); ++i){
        if(minims[i] < min){
            min = minims[i];
            min_vector.clear();
            min_vector.push_back(i);
            //min_vector.resize(minims.size(), -1);
            //min_vector[i] = min;
        }
        else if(minims[i] == min){
            min_vector.push_back(i);
        }
    }
    if(min == -1){
        run = false;
        min_vector.clear();
    }
    return min_vector;
}

int main(int argc, char** argv) {
	char ch;
	string inputfof, query;
	uint k(31);
	uint m1(10);
	uint c(1);
    uint64_t s(8); 

	while ((ch = getopt(argc, argv, "hdag:q:k:m:n:s:t:b:e:f:i:")) != -1) {
		switch (ch) {
			case 'f': inputfof = optarg; break;
			case 'k': k = stoi(optarg); break;
			case 'm': m1 = stoi(optarg); break;
			case 't': c = stoi(optarg); break;
            case 's': s = stoi(optarg); break;
		}
	}


	if ((inputfof == "")) {
		cout << "Core arguments:" << endl
		     << "	-i input files" << endl
		     << "	-k kmer size used  (31) " << endl
             << "	-s subsampling used  (8) " << endl
             << "	-m minimize size used  (10) " << endl;
		return 0;
	}else{
        Comparator comp = Comparator();
        comp.compare_buckets(inputfof);

        cout << "Containement index: " << endl;
        for(int i = 0; i < comp.nb_files_eof; ++i){
            cout << "file "<< i << "\t[";
            for(int j = 0; j < comp.nb_files_eof; ++j){
                cout << (double)comp.score_A[i*comp.nb_files_eof+j]/max(comp.nb_kmer_seen_infile[i], comp.nb_kmer_seen_infile[j]) << "  ";
            }
            cout << "]" << endl;
        }
        cout << "Jackard index:" << endl;
        for(int i = 0; i < comp.nb_files_eof; ++i){
            cout << "file "<< i << "\t[";
            for(int j = 0; j < comp.nb_files_eof; ++j){
                cout << (double)comp.score_A[i*comp.nb_files_eof+j]/(comp.nb_kmer_seen_infile[i] + comp.nb_kmer_seen_infile[j] - comp.score_A[i*comp.nb_files_eof+j])<< "  ";
            }
            cout << "]" << endl;
        }
        cout << "for " << intToString(comp.nb_kmer_seen) << " kmers evaluated" << endl;
    }
	    
}