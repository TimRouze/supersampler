#include "Comparator.h"
#include "utils.h"

void compare_files(const string& file1, const string& file2) {
	uint64_t total_nuc_number(0);
    uint64_t read_kmer(0);
    fstream input_stream;
    string ref, first, skmer_size, minimizer_size;
    uint32_t old_minimizer, minimizer;
    vector<uint32_t> keys;
    input_stream.open(file1);
    input_stream >> skmer_size;
    input_stream >> minimizer_size;
    cout << skmer_size << " " << minimizer_size << endl;
    uint64_t k = (stoi(skmer_size) + stoi(minimizer_size))/2;
    while (not input_stream.eof()) {
        ref = "";
        getline(input_stream, ref);
        if (ref.size() < k) {
            ref = "";
        } else {
            read_kmer += ref.size() - k + 1;
        }
        
    }
}
	

int main(int argc, char** argv) {
	char ch;
	string input, inputfof, query;
	uint k(31);
	uint m1(10);
	uint c(1);
    uint64_t s(8); 

	while ((ch = getopt(argc, argv, "hdag:q:k:m:n:s:t:b:e:f:i:")) != -1) {
		switch (ch) {
			case 'i': input = optarg; break;
			case 'k': k = stoi(optarg); break;
			case 'm': m1 = stoi(optarg); break;
			case 't': c = stoi(optarg); break;
            case 's': s = stoi(optarg); break;
		}
	}


	if ((input == "" )) {
		cout << "Core arguments:" << endl
		     << "	-i input file" << endl
		     << "	-k kmer size used  (31) " << endl
             << "	-s subsampling used  (8) " << endl
             << "	-m minimize size used  (10) " << endl;
		return 0;
	}else{
        cout<<" I use k="<<k<<" m="<<m1<<" s="<<s<<endl;
        cout<<"Maximal super kmer are of length "<<2*k-m1<<" or "<<k-m1+1<<" kmers" <<endl;
    }
	    
}