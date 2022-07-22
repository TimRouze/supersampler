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
    robin_hood::unordered_map<uint64_t, uint64_t> color_map;
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
                    if(color_map.count(seq) == 0){
                        nb_kmer_seen++;
                        color_map[seq] = 1;
                    }
                    updateK(seq, skip
                    [i + k], k);
                    i++;
                }
            }
            nb_kmer_seen_infile[ind] += color_map.size();
            color_map.clear();
        }
    }
}

void Comparator::update_colormap(const vector<istream*>& files, vector<uint64_t> indices){
    robin_hood::unordered_map<uint64_t, uint64_t> color_map;
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
                        //nb_kmer_seen_infile[ind]++;
                    }
                    else{
                        color_map[seq] |= ((uint64_t)1 << ind);
                        //nb_kmer_seen_infile[ind]++;
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
    string key = "";
    vector<uint64_t> score_v;
    uint64_t tmp;
    for (auto it : color_map){
        score_v.clear();
        tmp = it.second;
        for(int i(0);i<files.size();++i){
            score_v.push_back(tmp%2);
            tmp>>=1;
        }
        //score Matrix creation, both loops iterate through files.
        for(uint64_t i = 0; i < files.size()-1; ++i){
            for(uint64_t j = i+1; j < files.size(); ++j){
                if(score_v[i] == score_v[j] && score_v[i] == 1){
                    key = to_string(i) + ";" + to_string(j);
                    if(score_map.count(key) == 0){
                        score_map[key] = 1;
                        nb_kmer_seen_infile[i]++;
                        nb_kmer_seen_infile[j]++;
                    }
                    else{
                        score_map[key] += 1;
                        nb_kmer_seen_infile[i]++;
                        nb_kmer_seen_infile[j]++;
                    }
                }
            }
        }
    }
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
                    // nb_kmer_seen_infile[indices[i]] += (stoi(bucket_size)-2*k);// + 1)-k-1;
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
            // nb_kmer_seen_infile[i] += (stoi(bucket_size)-2*k);//+1)-k-1;
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

/*void Comparator::compare_files(const string& fileofile) {
    uint64_t read_kmer(0), skmerWoM_size2(0), skmer_size2(0), minimizer_size2(0), tmp_size(0), nb_skmer1(0), nb_skmer2(0), nb_minim(0);
    string ref, header, skmers2;
    uint32_t minimizer, minimizer2, curr_minim;
    vector<string> skmers1, lines;
    string sep = " ";
    string curr_filename, first;
    vector<istream*> input_files;
    istream* fof = openFile(fileofile);

    
    while (not fof->eof()){
        getline(*fof, curr_filename);
        if(curr_filename != ""){
            input_files.push_back(openFile(curr_filename));
        }
    }
    istream* first_f = input_files.front();
    input_files.erase(input_files.begin());
    //File 1 opening and header reading
    getline(*first_f, header);
    skmer_size = stoi(header.substr(0, header.find(sep)));
    header.erase(0, header.find(sep) + 1);
    minimizer_size = stoi(header.substr(0, header.find(sep)));
    header.erase(0, header.find(sep) + 1);
    nb_skmer1 = stoi(header);
    skmerWoM_size = skmer_size - minimizer_size;
    uint64_t k = ((skmer_size + minimizer_size)/2) - minimizer_size;

    for(int i = 0;  i < input_files.size(); ++i){
        getline(*input_files[i], header);
        getline(*input_files[i], first);
        lines.push_back(first);
    }

    cout << "Skmer size is: " << skmer_size << " Minimizer size is: " << minimizer_size << " Size of skmers read is: "  << skmerWoM_size << endl;
    cout << "k used: " << k << endl;
    // REVOIR LA LOGIQUE DES BOUCLES ça VA PAS
    while (not first_f->eof()) {
        ref = "";
        getline(*first_f, ref);
        if(not ref.empty()){
            minimizer = stoi(ref.substr(0, ref.find(sep)));
            //A VIRER
            ref.erase(0, ref.find(sep) + 1);
            ref.erase(0, ref.find(sep) + 1);
            for (int i = 0;  i < input_files.size(); ++i){
                while(not input_files[i]->eof()){
                    nb_minim++;
                    if(stoi(lines[i].substr(0, lines[i].find(sep))) == minimizer){
                        lines[i].erase(0, lines[i].find(sep) + 1);
                        lines[i].erase(0, lines[i].find(sep) + 1);
                        for (int j = 0; j+k <= ref.length(); j+=k) {
                            nb_kmer_seen++;
                            if(lines[i].length() >= j+k){
                                
                                if(ref.substr(j, k) == lines[i].substr(j, k)){
                                    if(color_map.count(ref.substr(j, k)) == 0){
                                        color_map[ref.substr(j, k)] = 0;
                                        color_map[ref.substr(j, k)] += 1 << i;
                                    }
                                    else{
                                        color_map[ref.substr(j, k)] |= 1 << i;
                                    }
                                }
                            }
                        }
                        lines[i] = "";
                        getline(*input_files[i], lines[i]);
                        break;
                    }
                    else if(stoi(lines[i].substr(0, lines[i].find(sep))) > minimizer & stoi(lines[i].substr(0, lines[i].find(sep))) != 4294967295){
                        break;
                    }
                    else{
                        lines[i] = "";
                        getline(*input_files[i], lines[i]);
                    }
                }
            }
            
        }
    }
    string key = "";
    for (auto it : color_map) { 
        for(int i = 0; i < input_files.size(); ++i){
            if(it.second&1 << i+1){
                key == to_string(i) + "," + to_string(i+1);
                if(score_map.count(key) == 0){
                    score_map[key] = 1;
                }
                else{
                    score_map[key] += 1;
                }
            }
        }
    } 
    cout << "nb minimizers: " << nb_minim;
}

string Comparator::find(istream* file2, uint32_t minimizer){
    string skmers, ref;
    string sep = " ";
    while(not file2->eof()){
        if(curr_minim == minimizer){
            ref.erase(0, ref.find(sep) + 1);
            ref.erase(0, ref.find(sep) + 1);
            return ref;
        }
        else if(curr_minim > minimizer & curr_minim != 4294967295){
            return "";
        }
        ref = "";
        getline(*file2, ref);
        if(not ref.empty()){
            curr_minim = stoi(ref.substr(0, ref.find(sep)));
        }
    }
    return "EOF";
}*/

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
        double AUB(0);
        //cout<<" I use k="<<k<<" m="<<m1<<" s="<<s<<endl;
        //cout<<"Maximal super kmer are of length "<<2*k-m1<<" or "<<k-m1+1<<" kmers" <<endl;
        //cout << "Comparing files \'" + input1 << "\' and \'" + input2 + "\'" << endl;
        Comparator comp = Comparator();
        comp.compare_buckets(inputfof);

        cout << "score is: " << endl;
        for(auto it : comp.score_map){
            cout << intToString(comp.nb_kmer_seen_infile[stoi(it.first.substr(0, it.first.find(";")))]);
            cout << " + " << intToString(comp.nb_kmer_seen_infile[stoi(it.first.substr(it.first.find(";")+1, 1))]);
            cout << " - " << intToString(it.second) << endl;
            AUB = (double)(comp.nb_kmer_seen_infile[stoi(it.first.substr(0, it.first.find(";")))] + comp.nb_kmer_seen_infile[stoi(it.first.substr(it.first.find(";")+1, 1))]) - it.second;
            cout << "(AUB)" << AUB << endl;
            cout << it.first << " " << (double)it.second/AUB << endl;
            cout << "Score " << intToString(it.second) << endl;
            cout << endl;
        } 
        cout << "for " << intToString(comp.nb_kmer_seen) << " kmers evaluated" << endl;
    }
	    
}