#include "utils.h"
#include <memory>
#include <map>
#include <algorithm>
#include <string>
#include <sys/stat.h>






kmer nuc2int(char c) {
    //~ cout<<c<<" "<<(c / 2) % 4<<endl;
	return (c / 2) % 4;
}



kmer nuc2intrc(char c) {
	return ((c / 2) % 4) ^ 2;
}



char int2nuc(unsigned char n){
	switch (n)
	{
	case 0:
		return 'A';
	case 1:
		return 'C';
	case 2:
		return 'T';
	case 3:
		return 'G';


	default:
		cout <<"wtf int2nuc"<<endl;
		cout<<(int)n<<endl;
		cin.get();
	}
	return 'A';
}


string strCompressor(const string& str){
	string result;
	if(str.empty()){
		return result;
	}
	char mod=str.size()%4;
	result+=mod;
	unsigned char c;
	for(uint64_t i(0);i<str.size();++i){
		c+=nuc2int(str[i]);
		if((i+1)%4==0){
			result+=c;
			c=0;
		}
        c<<=2;
	}
	if(mod!=0){
		result+=c;
	}
	return result;
}


string strDecompressor(const string* str){
	string result;
	if(str->empty()){
		return result;
	}
	char mod=str->at(0);
	unsigned char Packed_nuc;
	uint64_t last;
	if(mod==0){
		last=str->size();
	}else{
		last=str->size()-1;
	}
	char fchar[4];
	for(uint64_t i(1);i<last;++i){
		Packed_nuc=str->at(i);
		fchar[3]=int2nuc(Packed_nuc%4);
		Packed_nuc>>=2;
		fchar[2]=int2nuc(Packed_nuc%4);
		Packed_nuc>>=2;
		fchar[1]=int2nuc(Packed_nuc%4);
		Packed_nuc>>=2;
		fchar[0]=int2nuc(Packed_nuc%4);
		Packed_nuc>>=2;
		result+=fchar[0];
		result+=fchar[1];
		result+=fchar[2];
		result+=fchar[3];
	}
	if(mod!=0){
		Packed_nuc=str->at(last);
		for(uint64_t i(0);i<mod+1;++i){
			fchar[mod-i]=int2nuc(Packed_nuc%4);
			Packed_nuc>>=2;
		}
		for(uint64_t i(0);i<mod;++i){
			result+=fchar[i];
		}
	}
	return result;
}



string intToString(uint64_t n) {
	if (n < 1000) {
		return to_string(n);
	}
	string end(to_string(n % 1000));
	if (end.size() == 3) {
		return intToString(n / 1000) + "," + end;
	}
	if (end.size() == 2) {
		return intToString(n / 1000) + ",0" + end;
	}
	return intToString(n / 1000) + ",00" + end;
}



char revCompChar(char c) {
	switch (c) {
		case 'A': return 'T';
		case 'C': return 'G';
		case 'G': return 'C';
	}
	return 'A';
}



string revComp(const string& s) {
	string rc(s.size(), 0);
	for (int i((int)s.length() - 1); i >= 0; i--) {
		rc[s.size() - 1 - i] = revCompChar(s[i]);
	}
	return rc;
}



string getCanonical(const string& str) {
	return (min(str, revComp(str)));
}



kmer str2num(const string& str) {
	kmer res(0);
	for (uint64_t i(0); i < str.size(); i++) {
		res <<= 2;
		res += (str[i] / 2) % 4;
	}
	return res;
}


string num2str(kmer num,uint k){
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





uint64_t xorshift64star(uint64_t x) {
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    return x * 0x2545F4914F6CDD1DULL;
}



uint64_t xorshift64(uint64_t x){
	x ^= x << 13;
	x ^= x >> 7;
	x ^= x << 17;
	return x;
}


uint64_t murmur64(uint64_t h) {
  h ^= h >> 33;
  h *= 0xff51afd7ed558ccdL;
  h ^= h >> 33;
  h *= 0xc4ceb9fe1a85ec53L;
  h ^= h >> 33;
  return h;
}


static kmer hashtest( kmer u ){
  kmer v = u * 3935559000370003845 + 2691343689449507681;
  v ^= v >> 21;
  v ^= v << 37;
  v ^= v >>  4;
  v *= 4768777513237032717;
  v ^= v << 20;
  v ^= v >> 41;
  v ^= v <<  5;
  return v;
}



long hash64shift(long key)
{
  key = (~key) + (key << 21); // key = (key << 21) - key - 1;
  key = key ^ (key >> 24);
  key = (key + (key << 3)) + (key << 8); // key * 265
  key = key ^ (key >> 14);
  key = (key + (key << 2)) + (key << 4); // key * 21
  key = key ^ (key >> 28);
  key = key + (key << 31);
  return key;
}



uint64_t unrevhash(uint64_t x) {
	//return murmur64(x);
	//hash<uint64_t> my_hash;
    //return my_hash(x);
	uint64_t result2 = XXHash64::hash(&x, sizeof(x), 1312);
	return result2;
	x = ((x >> 32) ^ x) * 0xCFEE444D8B59A89B;
	x = ((x >> 32) ^ x) * 0xCFEE444D8B59A89B;
	x = ((x >> 32) ^ x);
	return x;
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



vector<bool> str2boolv(const string& str) {
	vector<bool> res;
	for (uint64_t i(0); i < str.size(); ++i) {
		if (str[i] == 'G' or str[i] == 'T') {
			res.push_back(true);
		} else {
			res.push_back(false);
		}
		if (str[i] == 'C' or str[i] == 'G') {
			res.push_back(true);
		} else {
			res.push_back(false);
		}
	}
	return res;
}



string bool2strv(const vector<bool>& v) {
	string res;
	for (uint64_t i(0); i < v.size(); i += 2) {
		if (v[i]) {
			if (v[i + 1]) {
				res += 'G';
			} else {
				res += 'T';
			}
		} else {
			if (v[i + 1]) {
				res += 'C';
			} else {
				res += 'A';
			}
		}
	}
	return res;
}






kmer hash64shift(kmer key) {
	key = (~key) + (key << 21); // key = (key << 21) - key - 1;
	key = key ^ (key >> 24);
	key = (key + (key << 3)) + (key << 8); // key * 265
	key = key ^ (key >> 14);
	key = (key + (key << 2)) + (key << 4); // key * 21
	key = key ^ (key >> 28);
	key = key + (key << 31);
	return key;
}



void cat_stream(istream& is, ostream& os) {
	const streamsize buff_size = 1 << 16;
	char* buff = new char[buff_size];
	while (true) {
		is.read(buff, buff_size);
		streamsize cnt = is.gcount();
		if (cnt == 0) break;
		os.write(buff, cnt);
	}
	delete[] buff;
}



void decompress_file(const string& file, const string& output_file) {
	unique_ptr<ofstream> ofs_p;
	ostream* os_p = &cout;
	if (not output_file.empty()) {
		ofs_p = unique_ptr<ofstream>(new strict_fstream::ofstream(output_file));
		os_p = ofs_p.get();
	}
	unique_ptr<istream> is_p(new zstr::ifstream(file));
	cat_stream(*is_p, *os_p);
}


struct stat STATbuffer;


zstr::ifstream* openFile(const string& input_file){
	zstr::ifstream* input_stream = new zstr::ifstream(input_file);
	if(not input_stream-> good()){
		cout << "Problem with file opening" << endl;
        return NULL;
	}
	return input_stream;
}



// It's quite complex to bitshift mmx register without an immediate (constant) count
// See: https://stackoverflow.com/questions/34478328/the-best-way-to-shift-a-m128i
__m128i mm_bitshift_left(__m128i x, unsigned count) {
	//~ assume(count < 128, "count=%u >= 128", count);
	__m128i carry = _mm_slli_si128(x, 8);
	if (count >= 64)                              // TODO: bench: Might be faster to skip this fast-path branch
		return _mm_slli_epi64(carry, count - 64); // the non-carry part is all zero, so return early
	// else
	carry = _mm_srli_epi64(carry, 64 - count);

	x = _mm_slli_epi64(x, count);
	return _mm_or_si128(x, carry);
}



__m128i mm_bitshift_right(__m128i x, unsigned count) {
	//~ assume(count < 128, "count=%u >= 128", count);
	__m128i carry = _mm_srli_si128(x, 8);
	if (count >= 64) return _mm_srli_epi64(carry, count - 64); // the non-carry part is all zero, so return early
	// else
	carry = _mm_slli_epi64(carry, 64 - count);

	x = _mm_srli_epi64(x, count);
	return _mm_or_si128(x, carry);
}



__uint128_t rcb(const __uint128_t& in, uint64_t n) {

	// assume(n <= 64, "n=%u > 64", n);

	union kmer_u {
		__uint128_t k;
		__m128i m128i;
		uint64_t u64[2];
		uint8_t u8[16];
	};

	kmer_u res = {.k = in};

	static_assert(sizeof(res) == sizeof(__uint128_t), "kmer sizeof mismatch");

	// Swap byte order

	kmer_u shuffidxs = {.u8 = {15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0}};

	res.m128i = _mm_shuffle_epi8(res.m128i, shuffidxs.m128i);

	// Swap nuc order in bytes

	const uint64_t c1 = 0x0f0f0f0f0f0f0f0f;

	const uint64_t c2 = 0x3333333333333333;

	for (uint64_t& x : res.u64) {

		x = ((x & c1) << 4) | ((x & (c1 << 4)) >> 4); // swap 2-nuc order in bytes

		x = ((x & c2) << 2) | ((x & (c2 << 2)) >> 2); // swap nuc order in 2-nuc

		x ^= 0xaaaaaaaaaaaaaaaa; // Complement;
	}

	// Realign to the right

	res.m128i = mm_bitshift_right(res.m128i, 128 - 2 * n);

	return res.k;
}



bool exists_test(const string& name) {
	struct stat buffer;
	return (stat(name.c_str(), &buffer) == 0);
}



uint64_t rcbc(uint64_t in, uint64_t n) {
	// assume(n <= 32, "n=%u > 32", n);
	// Complement, swap byte order
	uint64_t res = __builtin_bswap64(in ^ 0xaaaaaaaaaaaaaaaa);
	// Swap nuc order in bytes
	const uint64_t c1 = 0x0f0f0f0f0f0f0f0f;
	const uint64_t c2 = 0x3333333333333333;
	res = ((res & c1) << 4) | ((res & (c1 << 4)) >> 4); // swap 2-nuc order in bytes
	res = ((res & c2) << 2) | ((res & (c2 << 2)) >> 2); // swap nuc order in 2-nuc

	// Realign to the right
	res >>= 64 - 2 * n;
	return res;
}


uint64_t canonize(uint64_t x, uint64_t n) {
	return min(x, rcbc(x, n));
}


__uint128_t canonize(__uint128_t x, uint64_t n) {
	return min(x, rcb(x, n));
}



void print_bin(kmer n) {
	kmer mask = 1;
	mask <<= 63;
	for (uint64_t i(0); i < 64; ++i) {
		cout << (uint64_t)(n / mask);
		if (n / mask == 1) {
			n -= mask;
		}
		mask >>= 1;
	}
	cout << "\n";
}



kmer min_k(const kmer& k1, const kmer& k2) {
	if (k1 <= k2) {
		return k1;
	}
	return k2;
}



uint64_t asm_log2(const uint64_t x) {
	uint64_t y;
	asm("\tbsr %1, %0\n" : "=r"(y) : "r"(x));
	return y;
}



uint64_t mylog2(uint64_t val) {
	if (val == 0) return 0;
	if (val == 1) return 0;
	uint64_t ret = 0;
	while (val > 1) {
		val >>= 1;
		ret++;
	}
	return ret;
}



bool kmer_in_superkmer(const kmer canon, const vector<kmer>& V) {
	for (uint64_t i(0); i < V.size(); i++) {
		if (canon == V[i]) {
			return true;
		}
	}
	return false;
}



void dump_vector_bool(const vector<bool>& V, ostream* out) {
	int cmp = 0;
	uint8_t output = 0;
	vector<uint8_t> buf;
	for (uint64_t i(0); i < V.size(); ++i) {
		output = output | ((V[i] ? 1 : 0) << cmp);
		cmp++;
		if (cmp == 8) {
			buf.push_back(output);
			if (buf.size() >= 8000) {
				out->write((char*)buf.data(), buf.size());
				//~ *out<<flush;
				buf.clear();
			}
			cmp = 0;
			output = 0;
		}
	}
	if (V.size() % 8 != 0) {
		buf.push_back(output);
	}
	out->write((char*)buf.data(), buf.size());
}



void read_vector_bool(vector<bool>& V, zstr::ifstream* out, uint64_t n_bits) {
	uint64_t size_buffer(8000);
	uint64_t n_bytes(n_bits / 8 + (n_bits % 8 == 0 ? 0 : 1));
	uint64_t position(0);
	vector<uint8_t> buf(size_buffer, 0);
	while (position + size_buffer < n_bytes) {
		out->read((char*)buf.data(), size_buffer);
		for (uint64_t i(0); i < buf.size(); ++i) {
			V.push_back(buf[i] & 1);
			V.push_back(buf[i] & 2);
			V.push_back(buf[i] & 4);
			V.push_back(buf[i] & 8);
			V.push_back(buf[i] & 16);
			V.push_back(buf[i] & 32);
			V.push_back(buf[i] & 64);
			V.push_back(buf[i] & 128);
		}
		position += size_buffer;
	}
	buf.resize(n_bytes - position, 0);
	out->read((char*)buf.data(), n_bytes - position);
	for (uint64_t i(0); i < buf.size(); ++i) {
		V.push_back(buf[i] & 1);
		V.push_back(buf[i] & 2);
		V.push_back(buf[i] & 4);
		V.push_back(buf[i] & 8);
		V.push_back(buf[i] & 16);
		V.push_back(buf[i] & 32);
		V.push_back(buf[i] & 64);
		V.push_back(buf[i] & 128);
	}
}


vector<string> splitSTR(const string& s, char delim) {
	vector<string> res;
	uint pred(0);
	for (uint i(0); i < s.size(); ++i) {
		if (s[i] == delim) {
			res.push_back(s.substr(pred, i - pred));
			pred = i + 1;
		}
	}
    string laststr,end(s.substr(pred));

    cout<<"G"<<end<<"G"<<endl;
	res.push_back(end);
	return res;
}


vector<string> split(const string& s, char delim) {
	vector<string> res;
	uint pred(0);
	for (uint i(0); i < s.size(); ++i) {
		if (s[i] == delim) {
			res.push_back(s.substr(pred, i - pred));
			pred = i + 1;
		}
	}
    string laststr,end(s.substr(pred));
    for (uint i(0); i < end.size(); ++i) {
        if(isprint(end[i])){
            laststr+=end[i];
        }else{
            break;
        }
    }
    cout<<"g"<<laststr<<"g"<<endl;
	res.push_back(laststr);
	return res;
}



void split(const string& s, char delim, vector<string>& res) {
	res.clear();
	uint pred(0);
	for (uint i(0); i < s.size(); ++i) {
		if (s[i] == delim) {
			res.push_back(s.substr(pred, i - pred));
			pred = i + 1;
		}
	}
	res.push_back(s.substr(pred));
}



void Biogetline(zstr::ifstream* in,string& result,char type,uint K) {
  string discard;
  result.clear();
  switch(type){
    case 'Q':
      getline(*in,discard);
      getline(*in,result);
      getline(*in,discard);
      getline(*in,discard);
      break;
    case 'A':
      getline(*in,discard);
      char c=in->peek();
      while(c!='>' and c!=EOF){
        getline(*in,discard);
        transform(discard.begin(),discard.end(),discard.begin(),::toupper);
        result+=discard;
        c=in->peek();
      }
      break;
  }
  if(result.size()< K){
    result.clear();
  }
}



void clean_dna(string& str){
	vector<int> positions;
	for(uint i(0); i< str.size(); ++i){
		
		switch(str[i]){
			case 'a':break;
			case 'A':break;
			case 'c':break;
			case 'C':break;
			case 'g':break;
			case 'G':break;
			case 't':break;
			case 'T':break;
			default: 
			positions.push_back(i);

		}
	}
	if(positions.size() == str.size()){
		str = "";
	}else{
		sort(positions.begin(), positions.end(), greater{});
		for(uint i(0); i< positions.size(); ++i){
			str.erase(positions[i], 1);
		}
	}
	transform(str.begin(), str.end(), str.begin(), ::toupper);
}



string getLineFasta(zstr::ifstream* in) {
	string line, result;
	getline(*in, line);
	char c = static_cast<char>(in->peek());
	while (c != '>' and c != EOF) {
		getline(*in, line);
		result += line;
		c = static_cast<char>(in->peek());
	}

	clean_dna(result);
	return result;
}




void split2(const string& s, char delim, vector<string>& res) {
	res.clear();
	string word;
	uint siz(s.size());
	for (uint i(0); i < siz; ++i) {
		if (s[i] == delim) {
			res.push_back(word);
			word.clear();
		} else {
			word.push_back(s[i]);
		}
	}
	if (word.size() > 0) {
		res.push_back((word));
	}
}

string bool2str(vector<bool> V) {
	string result;
	for (uint64_t i(0); i < V.size(); ++i) {
		result += (V[i] ? '1' : '0');
	}
	return result;
}




//TODO UPDATEK AS A METHOD
void updateK(kmer& min, char nuc, uint64_t& k) {
	min <<= 2;
	min += nuc2int(nuc);
	min %= (kmer)1<<(2*k);
}





