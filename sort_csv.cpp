#include "utils.h"
#include <map>


using namespace std;

uint32_t get_num(const string& str){
    uint64_t begin(0), end(0);
    for(uint i(0);i<str.size();++i){
        if(str[i]=='_'){
            begin=i+1;
        }
        if(str[i] == '.'){
            end = i;
            if(end - begin != 0 && begin != 0){

                return stoi(str.substr(begin, end-begin));
            }else{
                return 0;
            }
        }
    }
    return 0;
}

void sort_CSV(const string& filename, const string& outfilename, const string& fof_name){
    zstr::ifstream in(filename);
    ifstream fof(fof_name);
    ofstream out(outfilename);
    string line, name;
    vector<string> names_ordered;
    if(not in.good() && fof.good()){
        cout<<"cant open file"<<endl;
        return;
    }
    while(not fof.eof()){
        getline(fof, name);
        names_ordered.push_back(name);
    }
    getline(in,line);
    vector<string> files_names(split(line,','));

    uint N(files_names.size());
    double* matrix= new double[N*N];
    //cout<<N<<" files found"<<endl;
    map<uint32_t,uint32_t> sorted_names;
    map<uint32_t,string> names;
    map<uint32_t,uint32_t> old2new;
    uint initial_id(0), pos(0);
    for(uint i(0);i<files_names.size();++i){
        //cout<<"B    B"<<files_names[i]<<"E    E"<<endl;
        auto it = find(names_ordered.begin(), names_ordered.end(), files_names[i]);
        pos = it - names_ordered.begin();
        sorted_names[pos]=initial_id;
        names[pos] = files_names[i];
        initial_id++;
    }
    //cout<<"Map filled"<<endl;
    

    uint id(0);
    for (auto const& [key, val] : sorted_names){
        old2new[val]=id;
        //cout<<"Old: "<<val<<" new: "<<id<<endl;
        id++;
    }
    id = 0;
    for (auto const& [key, val] : names){
        out<<val;
        id++;
        if(id!=N){out<<',';}
    }

    out<<endl;
    //cout<<"Got sorted filenames"<<endl;

    uint line_id(0);
    while(not in.eof()){
        getline(in,line);
        if(line.size()<N){break;}
        vector<string> values(split(line,','));
        for(uint i(0);i<N;++i){
            matrix[old2new[i]*N+old2new[line_id]]=stod(values[i]);
        }
        line_id++;
    }
    //cout<<"Loaded matrice"<<endl;
    //~ cin.get();

    for(uint i(0);i<N;++i){
        for(uint j(0);j<N;++j){
            out<<matrix[i*N+j];
            if(matrix[i*N+j]!=matrix[j*N+i]){
                //cout<<i<<" "<<j<<"          ";
                //cout<<matrix[i*N+j]<<" "<<matrix[j*N+i]<<endl;
                cout<<"bug1 OR you are sorting a containment file"<<endl;
                //cin.get();
            }
            if(i==j){
                if(matrix[i*N+j]!=1){
                 cout<<matrix[i*N+j]<<endl;
                    cout<<"bug2"<<endl;
                    cin.get();
                }
            }
            if(j!=N-1){out<<',';}
        }
        out<<endl;
    }
    cout<<"The end"<<endl;
}



int main(int argc, char** argv){
    if(argc<4){
        cout<<"Need input, output filename and original fof"<<endl;
        return 0;
    }
    sort_CSV(argv[1],argv[2], argv[3]);
    return 0;
}
