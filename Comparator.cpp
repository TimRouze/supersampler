#include "Comparator.h"
#include "utils.h"
#include "include/zstr.hpp"



void Comparator::getfilesname(const string& fileofile, vector<string>& result){
    istream* fof = openFile(fileofile);
    if(fof==NULL){
        cout<<"Can't open "<<fileofile<<endl;
        return;
    }
    string current;
    while(not fof->eof()){
        current.clear();
        getline(*fof,current);
        if(current.size()>2){
            result.push_back(current);
        }
    }
}



void Comparator::get_header_info(const vector<istream*>& files){
    string header;
    for(uint64_t f = 0; f < files.size(); ++f){
        getline(*files[f], header);
        skmer_size = stoi(header.substr(0, header.find(" ")));
        header.erase(0, header.find(" ") + 1);
        m = stoi(header.substr(0, header.find(" ")));
        header.erase(0, header.find(" ") + 1);
        nb_kmer_tot = stoi(header.substr(0, header.find(" ")));
        header.erase(0, header.find(" ") + 1);
        sub_rate = stoi(header.substr(0, header.find(" ")));
        header.erase(0, header.find(" ") + 1);
        jaccard_only = stoi(header.substr(0, header.find(" ")));
        header.erase(0, header.find(" ") + 1);
        super_abundance = stoi(header.substr(0, header.find(" ")));
        // if(super_abundance){
        //     cout<<"Super abundance mode"<<endl;
        // }
        header.erase(0, header.find(" ") + 1);
        log_abundance = stoi(header);
        // if(log_abundance){
        //     cout<<"Log abundance mode"<<endl;
        // }
        k = ((skmer_size + m)/2); //-m
        skmer_size-=m;
    }
}



void Comparator::compare_sketches(uint size_query){
    query_size=size_query;
    vector<istream*> input_files;
    vector<uint64_t> indices;
    vector<string> lines;
    string header;
    for(uint i(0);i<files_names.size();++i){
        auto in=openFile(files_names[i]);
        if(in!=NULL){
            input_files.push_back(in);
        }
    }
    nb_files=input_files.size();
    nb_kmer_seen_infile.resize(input_files.size(), 0);
    nb_kmer_seen_infile_ab.resize(input_files.size(), 0);
    square_counts.resize(input_files.size(), 0);
    minimizers.resize(input_files.size());
    vector<uint64_t> min_vector(input_files.size());
    get_header_info(input_files);
    // cout << "kmers evaluated are of length: " << k << " minimizer size is " << m << endl;
    increment_files(input_files, indices);
    while(run){
        bool queryfound=findMin(minimizers,indices);
        if(indices.size() == 1 or  not queryfound){
            //NO COMPARISON TO DO
            skip_bucket(input_files, indices, num2str(minimizers[indices[0]],m));//, out_kmer);
            increment_files(input_files, indices);
        }else{
            count_intersection(input_files, indices,num2str(minimizers[indices[0]],m), out_kmer);
            increment_files(input_files, indices);
        }
    }
    // cout<<"Comparisons done"<<endl;
    for(uint64_t f = 0; f < input_files.size(); ++f){
        delete input_files[f];
    }
    delete out_kmer;
}



string Comparator::inject_minimizer(const string* str, const string& minimizerstr){
    string result;
    if(str->size() != 0){
        for(uint i(0);i<str->size();i+=skmer_size/2){
            result+=str->substr(i,skmer_size/2);
            i+=skmer_size/2;
            if(i < str->size()){
                result+=minimizerstr;
                result+=str->substr(i,skmer_size/2);
            }
        }
    }
    else{
        result = minimizerstr;
    }
    return result;
}



//TODO ADD NB KMER IN SKETCH HEADER + ONLY READ THE LINES WITHOUT ACTUALLY READING KMERS
//NO COMPARISON TO DO IN THIS BUCKET BUT WE STILL HAVE TO COUNT THE AMOUNT OF KMERS
void Comparator::skip_bucket(const vector<istream*>& files, const vector<uint64_t>& indices,const string& strminimizer){
    string skip1,skip2;
    string * skip = new string();
    kmer curr_kmer;
    string tmp, val;
    for(uint indice(0);indice<indices.size();++indice){
        uint64_t ind=indices[indice];
        uint32_t size_buffer=0;
        skip->clear();
        files[ind]->read((char*)&size_buffer,sizeof(size_buffer));
        skip->resize(size_buffer);
        files[ind]->read(skip->data(),size_buffer);
        *skip=strDecompressor(skip);
        *skip=inject_minimizer(skip,strminimizer);
        
        if (skip->size() < k-m) {
            *skip = "";
        }
        if (not skip->empty()){
            uint64_t kmers_in_skmer = ((k-m+1)*(skip->size()))/(2*k-m);
            
            nb_kmer_seen_infile[ind] += kmers_in_skmer;
            if(not jaccard_only){
                if(super_abundance){
                    for(uint64_t i(0); i<(skip->size())/(2*k-m); i++){
                        uint16_t val(0);
                        files[ind]->read((char*)&val, sizeof(val));
                        if(log_abundance){
                            val=(uint64_t)1 << val;
                        }
                        for(uint64_t j(0);j<k-m+1;++j){
                            nb_kmer_seen_infile_ab[ind] += val*val;
                        }
                    }
                }else{
                    for(uint64_t i(0); i<kmers_in_skmer; i++){
                        uint16_t val(0);
                        files[ind]->read((char*)&val, sizeof(val));
                        if(log_abundance){
                            val=(uint64_t)1 << val;
                        }
                        nb_kmer_seen_infile_ab[ind] += val*val;
                    }
                }
            }
        }
        //HERE WE READ THE NONMAXIMAL SUPERKMERS
        bool go=true;
        while (go){
            skip1=skip2="";
            getline(*files[ind], skip1);
            getline(*files[ind], skip2);
            if(skip1.empty() and skip2.empty()){
                go=false;
            }else{
                getline(*files[ind], tmp);
                stringstream abundance(tmp);
                nb_kmer_seen_infile[ind] += skip1.size()+skip2.size()-k+m+1;//TODO CHECK
                if(not jaccard_only){
                    if(super_abundance){
                        uint64_t ival=stoi(tmp);
                        if(log_abundance){
                            ival=(uint64_t)1 << ival;
                        }
                        nb_kmer_seen_infile_ab[ind] += ival*ival*(skip1.size()+skip2.size()-k+m+1);//TODO CHECK
                    }else{
                        while (getline(abundance, val, ' ')) {
                            uint64_t ival=stoi(val);
                            if(log_abundance){
                                ival=(uint64_t)1 << ival;
                            }
                            nb_kmer_seen_infile_ab[ind] += ival*ival;
                        }
                    }
                }
            }

        }
    }
    delete skip;
}



string kmer2str(kmer num,uint k){
	string str;
	int nuc;
	for(uint i(0);i<k;++i){
		nuc=num%4;
		switch (nuc){
			case 0:str.push_back('A');break;
			case 1:str.push_back('C');break;
			case 2:str.push_back('T');break;
			case 3:str.push_back('G');break;
		}
		num>>=2;
	}
	reverse( str.begin(), str.end());
	return str;
}



void Comparator::count_intersection(const vector<istream*>& files, const vector<uint64_t>& indices,const string& strminimizer, zstr::ofstream* out_kmer){
    string kmer_canon,skip,skip2;
    string *ref = new std::string();
    kmer curr_kmer;
    ankerl::unordered_dense::map<kmer, vector<uint32_t> > color_map;
    vector<kmer> interesting_hits;
    string tmp, val;
    vector<uint16_t> v_abundances;
    string fasta(">\n");
    string k_mer;
    for(auto ind : indices){
        v_abundances.clear();
        uint nbsuperkmer(0);
        //HERE WE READ THE MAXIMAL SUPERKMERS
        uint32_t size_buffer=0;
        ref->clear();
        files[ind]->read((char*)&size_buffer,sizeof(size_buffer));
        ref->resize(size_buffer);
        files[ind]->read(ref->data(),size_buffer);
        *ref=strDecompressor(ref); 
        *ref=inject_minimizer(ref,strminimizer);
        uint64_t cpt_abundance(0);
        if (ref->size() < k-m){
            *ref = "";
        }
        if (not ref->empty()){
            uint64_t i(0);
            uint64_t kmers_in_skmer(0);
            kmers_in_skmer = ((k-m+1)*(ref->size()))/(2*k-m);
            nb_kmer_seen_infile[ind] += kmers_in_skmer;
            if(not jaccard_only){
                if(super_abundance){
                    for(uint64_t i(0); i<(ref->size())/(2*k-m); i++){
                        uint16_t val(0);
                        files[ind]->read((char*)&val, sizeof(val));
                        if(log_abundance){
                            val=(uint64_t)1 << val;
                        }
                        for(uint64_t j(0);j<k-m+1;++j){
                            v_abundances.push_back(val);
                            nb_kmer_seen_infile_ab[ind] += val*val;
                        }
                    }
                }else{
                    for(uint64_t i(0); i<kmers_in_skmer; i++){
                        uint16_t val(0);
                        files[ind]->read((char*)&val, sizeof(val));
                        if(log_abundance){
                            val=(uint64_t)1 << val;
                        }
                        v_abundances.push_back(val);
                        nb_kmer_seen_infile_ab[ind] += val*val;
                    }
                }
            }
            while((i + k) <= ref->size()){
                curr_kmer = str2num(ref->substr(i, k-1));
                //Quand la longueur du skmer est atteinte, on saute de k-1 nucleotide pour aller au skmer suivant.
                for(uint j(0);j<k-m+1;++j){
                    updateK(curr_kmer, ref->at(i + k-1), k);
                    kmer canon=canonize(curr_kmer,k);
                    // Si c'est la première fois que l'on voit le kmer
                    if(color_map.count(canon) == 0){
                        // On aggrandi la color map et on ajoute a la case du fichier courant l'abondance du kmer pour ce fichier.
                        color_map[canon].resize(files.size()+1,0);
                        if(jaccard_only){
                            color_map[canon][ind]=1;
                        }else{
                            color_map[canon][ind]=v_abundances[cpt_abundance];
                        }
                    // Sinon
                    }else{
                        //Si c'est la première fois que l'on voit ce k-mer pour ce fichier
                        if(color_map[canon][ind]==0){
                            // On ajoute à la case du fichier l'abondance correspondante.
                            if(jaccard_only){
                                color_map[canon][ind]=1;
                            }else{
                                color_map[canon][ind]=v_abundances[cpt_abundance];
                            }
                            // Si c'est au moins la seconde fois qu'on voit ce kmer
                            if(not color_map[canon][files.size()]){
                                // On a un match donc on flag le kmer
                                interesting_hits.push_back(canon);
                                if(jaccard_only){
                                    color_map[canon][files.size()]=1;
                                }else{
                                    color_map[canon][files.size()]=v_abundances[cpt_abundance];
                                }
                            }
                        }
                    }
                    i++;
                    cpt_abundance++;
                }
                nbsuperkmer++;
                i +=k-1;
            }
        }
        // HERE WE READ THE NONMAXIMAL SUPERKMERS
        bool go=true;
        while (go){
            skip=skip2="";
            getline(*files[ind], skip);
            getline(*files[ind], skip2);
            if(skip.empty() and skip2.empty()){
                go=false;
            }else{
                nbsuperkmer++;
                skip+=strminimizer+skip2;
                uint64_t i(0);
                getline(*files[ind], tmp);
                
                uint64_t cpt_abundance(0), nb_kmer(0);
                stringstream abundance(tmp);
                v_abundances.clear();
                nb_kmer_seen_infile[ind] += skip.size()-k+1;
                if(not jaccard_only){
                    if(super_abundance){
                        uint16_t val(stoi(tmp));
                        
                        if(log_abundance){
                            val=(uint64_t)1 << val;
                        }
                        for(uint64_t j(0);j<skip.size()-k+1;++j){
                            v_abundances.push_back(val);
                            nb_kmer_seen_infile_ab[ind] += val*val;
                        }
                    }else{
                        while (getline(abundance, val, ' ')) {
                            uint32_t ival(stoi(val));
                            v_abundances.push_back(ival);
                            nb_kmer_seen_infile_ab[ind] += ival*ival;
                            //abund_tot += stoi(val);
                        }
                    }
                }
                curr_kmer = str2num(skip.substr(0, k-1));
                while(i+k <= skip.size()){
                    updateK(curr_kmer, skip[i + k-1 ], k);
                    kmer canon=canonize(curr_kmer,k);
                    // Si c'est la première fois que l'on voit le kmer
                    if(color_map.count(canon) == 0){
                        // On aggrandi la color map et on ajoute a la case du fichier courant l'abondance du kmer pour ce fichier.
                        color_map[canon].resize(files.size()+1,0);
                        if(jaccard_only){
                            color_map[canon][ind]=1;
                        }else{
                            color_map[canon][ind]=v_abundances[cpt_abundance];
                        }
                    }else{ 
                        // On ajoute à la case du fichier l'abondance correspondante.
                        if(not color_map[canon][ind]){
                            if(jaccard_only){
                                    color_map[canon][ind]=1;
                                }else{
                                    color_map[canon][ind]=v_abundances[cpt_abundance];
                                }
                            // Si c'est au moins la seconde fois qu'on voit ce kmer
                            if(not color_map[canon][files.size()]){
                                // On a un match donc on flag le kmer
                                interesting_hits.push_back(canon);
                                if(jaccard_only){
                                    color_map[canon][ind]=1;
                                }else{
                                    color_map[canon][ind]=v_abundances[cpt_abundance];
                                }
                            }
                        }
                    }
                    i++;
                    cpt_abundance++;
                }
            }
        }
    }
    delete ref;
    compute_scores(color_map,interesting_hits);
}


void Comparator::compute_scores(const ankerl::unordered_dense::map<kmer,vector<uint32_t> >& color_map,vector<kmer>& interesting_hits){
    vector<uint32_t> score_v;
    vector<uint32_t> ones;
    // Pour chaque k-mer interessant (un match dans au moins 2 sketchs)
    for(uint32_t z=(0);z<interesting_hits.size();++z){
        // On récupère l'abondance du k-mer pour chaque fichier
        score_v = color_map.at(interesting_hits[z]);
        ones.clear();
        // Simplement pour parcourir la map
        for(uint32_t x(0);x<nb_files;++x){
            if(score_v[x]!=0){
                ones.push_back(x);
            }
        }
        for(uint64_t i(0); i < ones.size(); ++i){
            for(uint64_t j(i+1); j < ones.size(); ++j){
                score_A[ones[i]*nb_files+ones[j]].first += (score_v[i]*score_v[j]);
                score_A[ones[i]*nb_files+ones[j]].second++;
                //PB HERE, SCORE_V[I] DOES NOT CONTAIN EVERY ABUNDANCE QUESTION IS WHY
                abund_tot += score_v[i];
            }
        }
    }
}



void Comparator::increment_files(const vector<istream*>& files, const vector<uint64_t>& indices){
    vector<string> lines;
    string minim, bucket_size;
    string buffer(m,'A');
    if(indices.size() != 0){
        //THE REAL USECASE OF THIS FUNCTION
        for(uint64_t i(0); i < indices.size(); ++i){
            if(not files[indices[i]]->eof()){
                files[indices[i]]->read(buffer.data(),m);
                if(not files[indices[i]]->eof()){
                    minimizers[indices[i]] = (uint64_t)str2num(buffer);
                }else{
                    minimizers[indices[i]] = -1;
                    nb_files_eof++;
                }
            }else{
                minimizers[indices[i]] = -1;
                nb_files_eof++;
            }
        }
        if(nb_files_eof == nb_files){
            //ALL FILES ARE CONSUMED STOP EVERYTHING
            run = false;
        }
    }else{
        //INITIALIZATION
        for(uint64_t i(0); i < files.size(); ++i){
            files[i]->read(buffer.data(),m);
            minimizers[i] = (uint64_t)str2num(buffer);
        }
    }
}



//FIND THE SET OF SMALLEST BUCKETS TO COMPARE
bool Comparator::findMin(const vector<uint64_t>& minims,vector<uint64_t>& min_vector){
    uint64_t min(-1);
    uint64_t i(0);
    uint64_t cpt(0);
    bool result(false);
    min_vector.clear();
    for(; i < minims.size(); ++i){
        if(minims[i] < min){
            //NEW MINIMUM FLUSH PREVIOUS RESULT
            min = minims[i];
            min_vector.clear();
            min_vector.push_back(i);
            if(i<query_size){
                result=true;            
            }
            else{
                result = false;
            }
        }else if(minims[i] == min){
            //NEW FILE THAT SHARE THE MINIMUM
            min_vector.push_back(i);
            if(i<query_size){
                result=true;
            }
        }
    }
    if(min == (uint64_t)-1){
        run = false;
        min_vector.clear();
    }
    return result;
}



void Comparator::print_containment(const string& outfile){
    zstr::ofstream out(outfile,21,1);
    for(uint32_t i = 0; i < nb_files; ++i){
        out << files_names[i];
        if(i!=(nb_files-1)){
            out<<',';
        }else{
            out<<'\n';
        }
    }
    out<<"\n";
    for(uint32_t i = 0; i < nb_files and i< query_size; ++i){
        for(uint32_t j = 0; j < nb_files; ++j){
            if(i==j){
                out<<"1";
            }else if(i<j){
                if(score_A.count(i*nb_files+j)==0){
                    out<<"0";
                }else{
                    double score=(double)score_A[i*nb_files+j].second/nb_kmer_seen_infile[i];
                    if(score<min_threshold){
                        out<<'0';
                    }else{
                        out << setprecision(precision) << score;
                    }
                }
            }else{
                if(score_A.count(j*nb_files+i)==0){
                    out<<"0";
                }else{
                    double score=(double)score_A[j*nb_files+i].second/nb_kmer_seen_infile[i];
                    if(score<min_threshold){
                        out<<'0';
                    }else{
                        out << setprecision(precision) << score;
                    }
                }
            }
            if(j!=(nb_files-1)){
                out<<',';
            }else{
                out<<'\n';
            }
        }
    }
}



void Comparator::print_jaccard(const string& outfile){
    zstr::ofstream out(outfile,21,1);
    cout << "Jackard index dump" << endl;
    for(uint32_t i = 0; i < nb_files; ++i){
        out << files_names[i];
        if(i!=(nb_files-1)){
            out<<',';
        }else{
            out<<'\n';
        }
    }
    for(uint32_t i = 0; i < nb_files and i< query_size; ++i){
        for(uint32_t j = 0; j < nb_files; ++j){
            if(i==j){
                out<<'1';
            }else if(i<j){
                if(score_A.count(i*nb_files+j)==0){
                    out<<"0";
                }else{
                    cout<<"Inter:"<<intToString(score_A[i*nb_files+j].second)<<" Union:"<<intToString((nb_kmer_seen_infile[i]+nb_kmer_seen_infile[j]-score_A[i*nb_files+j].second))<<" A:"<<intToString(nb_kmer_seen_infile[i])<<" B:"<<intToString(nb_kmer_seen_infile[j])<<endl;
                    // cout << intToString(abund_tot) << endl;
                    double score = (double)score_A[i*nb_files+j].second/(nb_kmer_seen_infile[i]+nb_kmer_seen_infile[j]-score_A[i*nb_files+j].second);
                    if(score<min_threshold){
                        out<<'0';
                    }else{
                        out << setprecision(precision) << score;
                    }
                }
            }else{
                if(score_A.count(j*nb_files+i)==0){
                    out<<"0";
                }else{
                    double score = (double)score_A[j*nb_files+i].second/(nb_kmer_seen_infile[i]+nb_kmer_seen_infile[j]-score_A[j*nb_files+i].second);
                    if(score<min_threshold){
                        out<<'0';
                    }else{
                        out << setprecision(precision) << score;
                    }
                }
            }
            if(j!=(nb_files-1)){
                out<<',';
            }else{
                out<<'\n';
            }
        }
    }
}


void Comparator::print_weighted_jaccard(const string& outfile){
    zstr::ofstream out(outfile,21,1);
    cout << "Cosine index dump" << endl;
    for(uint32_t i = 0; i < nb_files; ++i){
        out << files_names[i];
        if(i!=(nb_files-1)){
            out<<',';
        }else{
            out<<'\n';
        }
    }
    for(uint32_t i = 0; i < nb_files and i< query_size; ++i){
        for(uint32_t j = 0; j < nb_files; ++j){
            if(i==j){
                out<<'1';
            }else if(i<j){
                if(score_A.count(i*nb_files+j)==0){
                    out<<"0";
                }else{
                    // cout<<"Inter:"<<intToString(score_A[i*nb_files+j])<<" Union:"<<intToString(sqrt((double)nb_kmer_seen_infile_ab[i]*nb_kmer_seen_infile_ab[j]))<<" A:"<<intToString(nb_kmer_seen_infile[i])<<" B:"<<intToString(nb_kmer_seen_infile[j])<<endl;
                    // cout << intToString(abund_tot) << endl;
                    double score = (double)score_A[i*nb_files+j].first/sqrt((double)nb_kmer_seen_infile_ab[i]*nb_kmer_seen_infile_ab[j]);
                    if(score<min_threshold){
                        out<<'0';
                    }else{
                        out << setprecision(precision) << (1-(2*acos(score))/3.14159265358979323846);
                    }
                }
            }else{
                if(score_A.count(j*nb_files+i)==0){
                    out<<"0";
                }else{
                    double score = (double)score_A[j*nb_files+i].first/sqrt((double)nb_kmer_seen_infile_ab[j]*nb_kmer_seen_infile_ab[i]);
                    if(score<min_threshold){
                        out<<'0';
                    }else{
                        out << setprecision(precision) << (1-(2*acos(score))/3.14159265358979323846);
                    }
                }
            }
            if(j!=(nb_files-1)){
                out<<',';
            }else{
                out<<'\n';
            }
        }
    }
}



int main(int argc, char** argv) {
	char ch;
	string inputfof, query;
    uint64_t p(6);
    double min_threshold(0);
    string output_name("results");

	while ((ch = getopt(argc, argv, "hdag:q:k:m:n:s:t:b:e:f:i:p:o:")) != -1) {
		switch (ch) {
			case 'f': inputfof = optarg; break;
            case 'q': query = optarg; break;
            case 'p': p = stoi(optarg); break;
            case 'm': min_threshold = stod(optarg); break;
            case 'o': output_name = optarg; break;
		}
	}

	if ((inputfof == "")) {
		cout << "Core arguments:" << endl
            << "-f Index file of files (mandatory)" << endl
            << "-q Query file of files (\"\" for all versus all comparison of the index)"<<endl

            << "Ouput arguments:" << endl
            << "-m Minimum value to be output (0.0)"<<endl
            << "-p Required precision to be output in the CSV (6)"<<endl
            << "-o output prefix (results)"<<endl
            ;
		return 0;
	}else{
        if(query==""){
            cout<<"No query file, I will perform a all versus all comparison"<<endl;
            Comparator comp = Comparator(p,min_threshold);
            comp.getfilesname(inputfof,comp.files_names);

            cout<<"I found "<<comp.files_names.size()<<" documents"<<endl;
            auto start = std::chrono::system_clock::now();
            comp.compare_sketches(comp.files_names.size());
            auto middle = std::chrono::system_clock::now();
            std::chrono::duration<double> elapsed_seconds = middle - start;
            cout<<"Comparisons lasted "<<elapsed_seconds.count()<<" sec"<<endl;
            comp.print_containment(output_name+"_containment.csv.gz");//TODO NAME
            comp.print_jaccard(output_name+"_jaccard.csv.gz");
            if(not comp.jaccard_only){
                comp.print_weighted_jaccard(output_name+"_angular.csv.gz");
            }

            auto end = std::chrono::system_clock::now();
            std::chrono::duration<double> elapsed_seconds2 = end- middle;
            cout<<"Jaccard output lasted "<<elapsed_seconds2.count()<<" sec"<<endl;
        }else{
            Comparator comp = Comparator(p,min_threshold);
            comp.getfilesname(query,comp.files_names);
            uint query_size=comp.files_names.size();
            cout<<"I query "<<query_size<<" file(s) against the bank"<<endl;
            comp.getfilesname(inputfof,comp.files_names);
            comp.compare_sketches(query_size);
            comp.print_containment(output_name+"_containment.csv.gz");//TODO NAME
            comp.print_jaccard(output_name+"_jaccard.csv.gz");
        }
    }
}
