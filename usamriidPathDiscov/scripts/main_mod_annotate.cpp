#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <cstdio>
#include <limits>
#include <iomanip>

using namespace std;

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//******************************************//
// constants								//
//******************************************//

// this is now a parameter
//const int CONTEXT_SIDE = 20;				// length of context on each side of variant

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//******************************************//
// feature type constants					//
//******************************************//

const int EXON = 0;
const int LSS = 1;	
const int RSS = 2;
const int INTRON = 3;
const int SNP = 5;
const int MIRNA = 6;
const int PROMOTER = 9;
const int TEXON = 10;
const int TLSS = 11;
const int TRSS = 12;
const int TINTRON = 13;
const int RNAEXON = 30;
const int VARIANT = 100;
const int THOU_LOW = 1001;
const int THOU_EXON = 1002;
const int THOU_TRIO = 1003;

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//******************************************//
// feature structs 							//
//******************************************//

struct base	
{
	// a base class from which all features will inherit
	
	// position, feature type, typeid
	//unsigned pos, typ, tid;
	unsigned pos, typ;
	
	void set0 (unsigned a1, unsigned a2) 
	{
  		pos = a1;
  		typ = a2;
	}	
};

struct variant : base	
{	
	unsigned tid; 						// type id
	string ref, var;					// ref nt, var nt
	
	void set (istringstream &iss) 
	{
		iss >> tid >> ws >> ref >> ws >> var;
	}	
};

struct base_exon : base	
{
	unsigned tid; 						// type id		 
	unsigned strt, fin, sen, num;		// feat_begin, feat_end, sense, exon num
};

struct rnaexon : base_exon	
{	 
	string txid;						// transcript id
	string dnachr;						// the corresponding dna chr
	unsigned dnastrt, dnafin;			// the corresponding dna start and end
	
	void set1 (string a1) 
	{
  		txid = a1;
	}	
	
	void set (istringstream &iss) 
	{
		iss >> tid >> ws >> strt >> ws >> fin >> ws >> sen >> ws >> num >> ws >> dnachr >> ws >> dnastrt >> ws >> dnafin;
	}	
};

struct exon : base_exon	
{
	unsigned ofst;						// exon offset w/in ccds
	
	void set (istringstream &iss) 
	{
		iss >> tid >> ws >> strt >> ws >> fin >> ws >> sen >> ws >> num >> ws >> ofst;
	}
};

struct texon : base_exon	
{	 
	unsigned ofst;						// offset
	
	void set (istringstream &iss) 
	{
		iss >> tid >> ws >> strt >> ws >> fin >> ws >> sen >> ws >> num >> ws >> ofst;
	}
};

struct intron : base_exon	
{
		
	void set (istringstream &iss) 
	{
		iss >> tid >> ws >> strt >> ws >> fin >> ws >> sen >> ws >> num;
	}	
};

struct tintron : base_exon	
{	

	void set (istringstream &iss) 
	{
		iss >> tid >> ws >> strt >> ws >> fin >> ws >> sen >> ws >> num;
	}	
};

struct splice : base_exon	
{	
	
	void set (istringstream &iss) 
	{
		iss >> tid >> ws >> strt >> ws >> fin >> ws >> sen >> ws >> num;
	}	
};

struct tsplice : base_exon	
{	
	
	void set (istringstream &iss) 
	{
		iss >> tid >> ws >> strt >> ws >> fin >> ws >> sen >> ws >> num;
	}	
};

struct snp : base	
{	
	unsigned tid; 					// type id
	unsigned strt, fin, sen;			// feat_begin, feat_end, sense
	string ntsnp;						// the SNP nucleotide change
	string bitfield;					// the SNP bitfield (not included in hg18)
	
	void set (istringstream &iss) 
	{
		iss >> tid >> ws >> strt >> ws >> fin >> ws >> sen >> ws >> ntsnp;
	}	

	// second method for dbSNP132 and on
	void set2 (istringstream &iss) 
	{
		iss >> tid >> ws >> strt >> ws >> fin >> ws >> sen >> ws >> ntsnp >> ws >> bitfield;
	}	
};

struct mirna : base	
{	
	unsigned tid; 					// type id
	unsigned strt, fin, sen;			// feat_begin, feat_end, sense
	
	void set (istringstream &iss) 
	{
		iss >> tid >> ws >> strt >> ws >> fin >> ws >> sen;
	}		
};

struct promoter : base	
{	
	unsigned tid; 					// type id
	unsigned strt, fin, sen;			// feat_begin, feat_end, sense

	void set (istringstream &iss) 
	{
		iss >> tid >> ws >> strt >> ws >> fin >> ws >> sen;
	}	
};

struct thousand : base	
{	
	string ref, var;					// nucleotides ref and var
	unsigned ac;						// AC allele count in genotypes, for each ALT allele, in the same order as listed
	unsigned an;						// AN total number of alleles in called genotypes	
	unsigned dp;						// DP combined depth across samples	

	void set (istringstream &iss, int what_type) 
	{
		if (what_type == THOU_TRIO)
		{
			iss >> ref >> ws >> var >> ws >> ac;
			an=0;
			dp=0;		
		}
		else
		{
			iss >> ref >> ws >> var >> ws >> ac >> ws >> an >> ws >> dp;
		}
	}
};

struct row
{	
	vector<variant> v_variant;			// vector of variants
	vector<exon> v_exon;				// vector of exons
	vector<intron> v_intron;			// vector of introns
	vector<splice> v_splice;			// vector of splice sites
	vector<texon> v_texon;				// vector of transcript exons
	vector<tintron> v_tintron;			// vector of transcript introns		
	vector<tsplice> v_tsplice;			// vector of transcript splice
	vector<snp> v_snp;					// vector of SNPs
	vector<mirna> v_mirna;				// vector of miRNAs
	vector<promoter> v_promoter;		// vector of promoters
	vector<thousand> v_thousand;		// vector of 1000 genome
	vector<rnaexon> v_rnaexon;			// vector of 1000 genome
	
	void clear_row ()
	{
		v_variant.clear();
		v_exon.clear();
		v_intron.clear();
		v_splice.clear();
		v_texon.clear();
		v_tintron.clear();
		v_tsplice.clear();
		v_snp.clear();
		v_mirna.clear();
		v_promoter.clear();
		v_thousand.clear();
		v_rnaexon.clear();
	}	
};

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//******************************************//
// other structs 							//
//******************************************//

// a structure containing 2 strings
struct twostr	 
{
	string s1;
	string s2;
};

// a structure containing 3 strings
struct threestr	 
{
	string s1;
	string s2;
	string s3;
};

// a structure containing 5 strings
struct fivestr	 
{
	string s1;
	string s2;
	string s3;
	string s4;
	string s5;	
};

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//******************************************//
// file reading functions					//
//******************************************//

// read humangencode.txt (codon	\t AA)
void read_gene_code(const char *f, map<string, string> &x) 
{
	ifstream is(f);
	string s1, s2;
	
	x.clear();
	while (1) 
	{
		is >> s1;
		if (is.eof()) break;
		is >> s2;

		x[s1] = s2; 
	}
}

// read a one column file of strings and store it in a vector of strings
void read_stringfile(char *s, vector<string> &m0) 
{	
	string line;
	istringstream iss;

	ifstream myfile (s);
	if (myfile.is_open())
  	{
    	while (!myfile.eof())
    	{
         	getline(myfile,line);
         	iss.clear();
         	iss.str(line);
         	if (!line.empty())    
         	{     	
      			m0.push_back(line);
         	}
    	}
    	myfile.close();
    }
  	else cout << "Unable to open file";    		
};

// read a 2 col file of strings, and store it in a vector of twostr structs
void read_2string(char *s, vector<twostr> &m1) 
{	
	string line;
	istringstream iss;
	twostr m;

	ifstream myfile (s);
	if (myfile.is_open())
  	{
    	while (!myfile.eof())
    	{
         	getline(myfile,line);
         	iss.clear();
         	iss.str(line);
         	if (!line.empty())    
         	{     	
				iss >> m.s1 >> ws >> m.s2;
				m1.push_back(m);
         	}
    	}
    	myfile.close();
    }
  	else cerr << "Unable to open file";    		
};

// read a 3 col file of strings, and store it in a vector of threestr structs
void read_3string(char *s, vector<threestr> &m1) 
{	
	string line;
	istringstream iss;
	threestr m;

	ifstream myfile (s);
	if (myfile.is_open())
  	{
    	while (!myfile.eof())
    	{
         	getline(myfile,line);
         	iss.clear();
         	iss.str(line);
         	if (!line.empty())    
         	{     	
				iss >> m.s1 >> ws >> m.s2 >> ws >> m.s3;
				m1.push_back(m);
         	}
    	}
    	myfile.close();
    }
  	else cerr << "Unable to open file";    		
};

// read a 5 col file of strings, and store it in a vector of fivestr structs
void read_5string(char *s, vector<fivestr> &m1) 
{	
	string line;
	istringstream iss;
	fivestr m;

	ifstream myfile (s);
	if (myfile.is_open())
  	{
    	while (!myfile.eof())
    	{
         	getline(myfile,line);
         	iss.clear();
         	iss.str(line);
         	if (!line.empty())    
         	{     	
				iss >> m.s1 >> ws >> m.s2 >> ws >> m.s3 >> ws >> m.s4 >> ws >> m.s5;
				m1.push_back(m);
         	}
    	}
    	myfile.close();
    }
  	else cerr << "Unable to open file";    		
};

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//******************************************//
// annotate functions						//
//******************************************//

void nt2aa(string &s, map<string, string> &g, string &r) 
{
	// nucleotide to amino acid

	string e;
	
	r = "";
	for(unsigned i = 0; i < s.length(); i += 3) {
		e = s.substr(i, 3);
		if (g.find(e) != g.end()) r = r + g[e];
		else r = r + 'X';
	}
}

string reverse_complement_seq(const string &str) 
{
	// return the rev comp of the string, but preserve dashes
	  	
  	string s;
  	
  	s.resize(str.size());
  	
  	int r;
  	
  	for (int i = str.size() - 1; i >= 0; --i)
  	{
  		r = str.size() - i - 1;
  		
  		switch(str[i]) 
  		{
  			case 'A':
  				s[r] = 'T';
  				break;
  			case 'T':
  				s[r] = 'A';
  				break;
  			case 'C':
  				s[r] = 'G';
  				break;
  			case 'G':
  				s[r] = 'C';
  				break;
  			case '-':
  				s[r] = '-';
  				break;  				  				
  			case 'a':
  				s[r] = 't';
  				break;
  			case 't':
  				s[r] = 'a';
  				break;
  			case 'c':
  				s[r] = 'g';
  				break;
  			case 'g':
  				s[r] = 'c';
  				break;  			
  		}
  	}
  	return s;
};

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void get_context(int start, int finish, string &full_seq, string &cntxt, int context_side)
{
	/*
	get the context around a variant 
	cntxt is passed by reference and modified

	full_seq is the chr seq on a single line
	CONTEXT_SIDE is the number of nucleotides on each side of the variant
	*/
		
	int a;									// a is left position
	int c;									// c is length of left flank
	int b;									// b is length of right flank	 						 
	
	cntxt="";								// reset value of context 
  	
 	// must convert variant positions to zero based counting
 	 		
 	a = max(0, start - 1 - context_side); 										// don't want to go back further than 0 	(a is left position, c is length of left flank)
 	if (a==0) 
 		c = start - 1; 
 	else 
 		c = context_side;
 		
 	b = min(context_side, ((int)full_seq.length() - 1) - (finish - 1)); 		// don't want to go past the end			(b is length of right flank)
 																				// must cast to int or get error: no matching function for call to ‘min(const int&, long unsigned int)’
 		
	cntxt = full_seq.substr(a, c) + "[" + full_seq.substr(start - 1, finish - start + 1) + "]" + full_seq.substr(finish, b); 			

}

bool get_rna_codon(rnaexon &f, variant &v, string &tx_full_seq, map<string, string> &g, string &tx_nam, string &cdn, string &aa, int strtcode, int endcode, int &variant_row) 
{
	
	// this doesnt output anything for indels
	// ASSUME UPPER CASE
	
	string refcdn, varcdn;		 										// a string to for the reference and variant codons
	string refaa, varaa;												// strings for the amino acids
	string refnt = v.ref;												// ref nt
	string varnt = v.var;												// var nt
	int ntpos, aapos;													// ntpos is (pos in the ccds seq = offset of variant in exon + offset of exon within ccds seq)
	ostringstream ossaa, osscdn;
	bool nsyn;
	
	// get the coding portion only
	string tx_seq = tx_full_seq.substr(strtcode - 1, endcode - strtcode + 1);

	// calculate ntpos position in coding seq							// ntpos is 0-based counting

	ntpos = v.pos - strtcode;

 	aapos = ntpos/3;													// assume zero based counting here

	nsyn = false;

	if (ntpos < 0)
		cerr << "error at v " << variant_row << ": " << tx_nam << " ntpos negative" << endl;
	else if (refnt != "-" && tx_seq.substr(ntpos, 1) != refnt) 			// check to make sure the ref is what it's stated to be			
	  	cerr << "error at v " << variant_row << ": " << tx_nam << " ref doesnt match seq at " << ntpos << ", stated ref: " << refnt << ", actual ref: " << tx_seq.substr(ntpos, 1) << endl;
	else if (tx_seq.size() % 3 > 0) 									// check to make sure coding seq divisible by 3
	  	cerr << "error at v " << variant_row << ": " << tx_nam << " coding seq not divisible by 3" << endl;
	else 
	{															
		refcdn = tx_seq.substr(ntpos - (ntpos % 3),3); 			    	// get the ref codon
		nt2aa(refcdn, g, refaa);										// get the ref AA	
		
		if (refnt == "-") 
			if (varnt == "-") 
			{
				varcdn = "?";
				varaa = "?";
			} 
			else 
			{
				varcdn = "+";
				varaa = "+";
			}
		else if (varnt == "-") 
		{
			varcdn = "-";
			varaa = "-";
		} 
		else 
		{
			varcdn = refcdn;
			varcdn.replace(ntpos % 3, 1, varnt);		 				// get the var codon - replace the appropriate position with the variant nucleotide
			nt2aa(varcdn, g, varaa);									// get the var AA		
		}
		
		osscdn << refcdn << (1+ntpos % 3) << varcdn;
		ossaa << refaa << (aapos+1) << varaa;							// need to make 1 based counting

		nsyn = (refaa != varaa);
  	}
	
	cdn = osscdn.str();
	aa = ossaa.str();
	return nsyn;
}

bool get_codon(exon &f, variant &v, string &ccds_seq, map<string, string> &g, string &ccds_nam, string &cdn, string &aa, int &variant_row) 
{
	
	// this doesnt output anything for indels
	// ASSUME UPPER CASE
	
	string refcdn, varcdn;		 										// a string to for the reference and variant codons
	string refaa, varaa;												// strings for the amino acids
	string refnt = v.ref;												// ref nt
	string varnt = v.var;												// var nt
	int ntpos, aapos;													// ntpos is (pos in the ccds seq = offset of variant in exon + offset of exon within ccds seq)
	ostringstream ossaa, osscdn;
	bool nsyn;

	// calculate ntpos position in ccds seq								// ntpos is 0-based counting
	// + : (exon ofst win ccds seq) + (pos - exon_start)
	// - : (exon ofst win ccds seq) + (exon_end - pos)
	if (f.sen == 1)
		ntpos = f.ofst + (f.pos - f.strt);
	else if (f.sen == 0)	
		ntpos = f.ofst + (f.fin - f.pos);

 	aapos = ntpos/3;													// assume zero based counting here

	if (f.sen == 0)														// account for sense - if negative strand, reverse complement
	{
  		refnt = reverse_complement_seq(refnt);
  		varnt = reverse_complement_seq(varnt);
  	}

	nsyn = false;
	if (refnt != "-" && ccds_seq.substr(ntpos, 1) != refnt) 			// check to make sure the ref is what it's stated to be
	{			
	  	cerr << "error at v " << variant_row << ": " << ccds_nam << " ref doesnt match seq at " << ntpos << ", stated ref: " << refnt << ", actual ref: " << ccds_seq.substr(ntpos, 1) << endl;
  		osscdn << "";
		ossaa << "";
	} 
	else if (ccds_seq.size() % 3 > 0) 									// check to make sure ccds seq divisible by 3
	{				 				
	  	cerr << "error at v " << variant_row << ": " << ccds_nam << " seq not divisible by 3" << endl;
  		osscdn << "";
		ossaa << "";
	} 
	else 
	{															
		refcdn = ccds_seq.substr(ntpos - (ntpos % 3),3); 			    // get the ref codon
		nt2aa(refcdn, g, refaa);										// get the ref AA	
		
		if (refnt == "-") 
			if (varnt == "-") 
			{
				varcdn = "?";
				varaa = "?";
			} 
			else 
			{
				varcdn = "+";
				varaa = "+";
			}
		else if (varnt == "-") 
		{
			varcdn = "-";
			varaa = "-";
		} 
		else 
		{
			varcdn = refcdn;
			varcdn.replace(ntpos % 3, 1, varnt);		 				// get the var codon - replace the appropriate position with the variant nucleotide
			nt2aa(varcdn, g, varaa);									// get the var AA		
		}
		
		osscdn << refcdn << (1+ntpos % 3) << varcdn;
		ossaa << refaa << (aapos+1) << varaa;							// need to make 1 based counting

		nsyn = (refaa != varaa);
  	}
	
	cdn = osscdn.str();
	aa = ossaa.str();
	return nsyn;
}

template <typename T>
void get_names(vector<T> &my_vector, vector<threestr> &m1, stringstream &ccds_exon_stream) 
{
	// get gene name, ccds names, exon numbers, and senses; and concatenate them onto comma-delimited strings

	string strccds, strgene, strsen;
	stringstream strexonum;

	// current & previous elt - to avoid concating duplicates onto strgene
	string c_gene, p_gene;		
	string c_ccds;	
	string c_sen;
	
	T ex; // a base_exon
		
	for(int i = 0; i < my_vector.size(); i++)	
	{	
		ex = my_vector[i];			// do this so dont have to repeatedly lookup
						
		if (i==0)
		{
			strccds = m1[ex.tid-1].s1;
			strgene = m1[ex.tid-1].s2;
			if (m1[ex.tid-1].s3 == "+") {strsen="1";} else {strsen="0";}
			strexonum << ex.num;
					
			// set prev
			p_gene=strgene;
		}
		else
		{
			// get current
			c_ccds = m1[ex.tid-1].s1;
			c_gene = m1[ex.tid-1].s2;
			if (m1[ex.tid-1].s3 == "+") {c_sen="1";} else {c_sen="0";}		
			
			// ccds should be unique - not nec to check for duplicates
			if (c_gene != p_gene) {strgene = strgene + "," + c_gene;}
			strccds = strccds + "," + c_ccds;
			strsen = strsen + "," + c_sen;
			strexonum << "," << ex.num;		

			p_gene = c_gene;						
		}
	}
	
	// set the output
	ccds_exon_stream << strgene << "\t" << strccds << "\t" << strsen << "\t" << strexonum.str();	
}

void get_t_names(row &r, vector<threestr> &m1, stringstream &t_exon_stream) 
{
	
	// get transcript accession name, transcript exon numbers, and transcript senses; and concatenate them onto comma-delimited strings
	
	string stracc, strgene, strsen; 
	stringstream strexonum, str_t_pos;
	
	string c_sen;
	unsigned t_pos;				// the position in the transcript in 1-based counting
	
	// current & previous elt - to avoid concating duplicates onto strgene
	string c_gene, p_gene;	
	
	texon tex;
		
	for(int i = 0; i < r.v_texon.size(); i++)	
	{
		tex = r.v_texon[i];		// do this so dont have to repeatedly lookup
		
		// calculate position within transcript
		if (tex.sen == 1)
			t_pos = tex.ofst + (tex.pos - tex.strt) + 1;
		else	
			t_pos = tex.ofst + (tex.fin - tex.pos) + 1;		
					
		if (i==0)
		{
			stracc = m1[tex.tid-1].s1;
			strgene = m1[tex.tid-1].s2;			
			if (m1[tex.tid-1].s3 == "+") {strsen="1";} else {strsen="0";}
			strexonum << tex.num;	
			str_t_pos << t_pos;	
			
			// set prev
			p_gene=strgene;			
		}
		else
		{
			c_gene = m1[tex.tid-1].s2;
			if (c_gene != p_gene) {strgene = strgene + "," + c_gene;}
			p_gene = c_gene;							
			
			stracc = stracc + "," + m1[tex.tid-1].s1;
			strexonum << "," << tex.num;		
			if (m1[tex.tid-1].s3 == "+") {c_sen="1";} else {c_sen="0";}			
			strsen = strsen + "," + c_sen;				
			str_t_pos << "," << t_pos;	
		}
	}
	
	t_exon_stream << strgene << "\t" << stracc << "\t" << strsen << "\t" << strexonum.str() << "\t" << str_t_pos.str();	
}

void get_p_names(row &r, vector<threestr> &m1, stringstream &t_promoter_stream) 
{	
	
	// get promoter names onto comma-delimited strings

	string stracc, strgene, strsen;
	stringstream str_p_pos;

	promoter pro;
	unsigned p_pos;		// position - distance from transcription start site
	
	// current & previous elt - to avoid concating duplicates onto strgene
	string c_gene, p_gene, c_sen;	
		
	for(int i = 0; i < r.v_promoter.size(); i++)	
	{
		pro = r.v_promoter[i];	// do this so dont have to repeatedly lookup
		
		// calculate distance from transcription start site
		if (pro.sen == 1)
			p_pos = (pro.pos - pro.strt) + 1;	
		else
			p_pos = (pro.fin - pro.pos) + 1;

						
		if (i==0)
		{
			stracc = m1[pro.tid-1].s1;
			strgene = m1[pro.tid-1].s2;
			if (m1[pro.tid-1].s3 == "+") {strsen="1";} else {strsen="0";}
			str_p_pos << p_pos;
			// set prev
			p_gene=strgene;					
		}
		else
		{
			c_gene = m1[pro.tid-1].s2;
			if (c_gene != p_gene) {strgene = strgene + "," + c_gene;}
			p_gene = c_gene;
			
			if (m1[pro.tid-1].s3 == "+") {c_sen="1";} else {c_sen="0";}			
			strsen = strsen + "," + c_sen;			
						
			stracc = stracc + "," + m1[pro.tid-1].s1;
			str_p_pos << "," << p_pos;			
		}
	}
	
	t_promoter_stream << strgene << "\t" << stracc << "\t" << strsen << "\t" << str_p_pos.str();	
}

void get_mirna_names(row &r, vector<string> &m4, stringstream &mirna_stream) 
{	
	
	// get mirna names onto comma-delimited strings

	string strgene, strsen;
	string c_gene, c_sen;
	
	mirna mir;
	
	for(int i = 0; i < r.v_mirna.size(); i++)	
	{
		mir = r.v_mirna[i];	// do this so dont have to repeatedly lookup
						
		if (i==0)
		{
			strgene = m4[mir.tid-1];
			if (mir.sen == 1) {strsen="1";} else {strsen="0";}			
				
		}
		else
		{			
			c_gene = m4[mir.tid-1];
			if (mir.sen == 1) {c_sen="1";} else {c_sen="0";}
			strgene = strgene + "," + c_gene;
			strsen = strsen + "," + c_sen;	
		}
	}
	
	mirna_stream << strgene << "\t" << strsen;
}

bool is_match(string &nt1, string &nt2, string &nt3, string &nt4, unsigned sense)
{
	bool match=false;
	
	//cout << "match_test " << nt1 << "," << nt2 << " " << nt3 << "," << nt4 << endl;
	
	if (sense==1)
	{
		//cout << "match " << nt1 << "," << nt2 << " " << nt3 << "," << nt4 << endl;
		
		if ( (nt3 == nt1 && nt4 == nt2) || (nt3 == nt2 && nt4 == nt1) )
			match=true;
	}
	else										
	{
		// neg sense => must reverse complement

		//cout << "match_rc " << nt1 << "," << nt2 << " " << reverse_complement_seq(nt3) << "," << reverse_complement_seq(nt4) << endl;

		if ( (reverse_complement_seq(nt3) == nt1 && reverse_complement_seq(nt4) == nt2) || (reverse_complement_seq(nt3) == nt2 && reverse_complement_seq(nt4) == nt1) )
			match=true;		
	}

	return match;

}

short is_thou(thousand &f, variant &v, stringstream &strac, stringstream &stran, stringstream &strdp, int &n_thousand) 
{	
	
	// get thou genome data ids if match etc

	short val = 0;																	// return value - starts off false	

	if (is_match(v.ref, v.var, f.ref, f.var, 1))									// test for match
	{	
		val = 1;
		if (n_thousand==0)
		{
			strac << f.ac;
			
			// ugly patch - account for differences in format for trio data
			if (f.typ == THOU_TRIO)
			{
				stran << "-";
				strdp << "-";
			}
			else
			{
				stran << f.an;
				strdp << f.dp;
			}			
		}
		else
		{
			strac << "," << f.ac;
			
			if (f.typ == THOU_TRIO)
			{
				stran << ",-";
				strdp << ",-";
			}
			else
			{
				stran << "," << f.an;
				strdp << "," << f.dp;
			}									
		}
	}
	return val;
}

short is_snp(snp &f, variant &v, vector<string> &m2, string &strsnp, string &strbit, short &typ_snp, bool is_medical) 
{
	// typ_snp: 1 - variant is a point mutation that exactly matches the SNP, 2 - variant is on exactly SNP position, 0 - otherwise
	// return: 1 if typ_snp==1 || typ_snp==2 
	// SO FAR just checking single point mutations SNPs
	
	// FIX THIS FUNCTION LATER
	
	string nt3, nt4;																// strings for the nucleotides of the SNP
		
	/*
		Because the reference has not yet been corrected, must account for the following phenomenon:
	 	UCSC: "Question: I am confused about the start coordinates for items in the refGene table. It looks like you need to add "1" to the starting point in order to get 
	 	the same start coordinate as is shown by the Genome Browser. Why is this the case?"
	 	Response:
	 	Our internal database representations of coordinates always have a zero-based start and a one-based end. We add 1 to the start before displaying coordinates in the Genome Browser. 
	 	Therefore, they appear as one-based start, one-based end in the graphical display. The refGene.txt file is a database file, and consequently is based on the internal representation.
	 	We use this particular internal representation because it simplifies coordinate arithmetic, i.e. it eliminates the need to add or subtract 1 at every step. 
	 	Unfortunately, it does create some confusion when the internal representation is exposed or when we forget to add 1 before displaying a start coordinate. 
	 	However, it saves us from much trickier bugs. In summary, if you use a database dump file but would prefer to see the one-based start coordinates, 
		you will always need to add 1 to each start coordinate."
	*/
	
	short val = 0;																	// return value - starts off false
	
	if (f.fin - f.strt == 1 && f.pos == f.fin)										// if at correct position (only consider single point mutations, not indels or MNPs)
	{
		val = 1;																	// it intersects
		
		if (typ_snp==0)															// concat SNPs
		{
			strsnp=m2[f.tid-1];
			if (is_medical) {strbit=f.bitfield;}
		} 
		else 
		{
			strsnp=strsnp+","+m2[f.tid-1];
			if (is_medical) {strbit=strbit+","+f.bitfield;}
		}

		
		if (typ_snp!=1)																// if perfect, dont re-test
		{
			typ_snp=2;																// we know it's at least on the right position			
			
			if (f.ntsnp.length() == 3)												// now test for nucleotide match (string looks like "A/G")
			{
				nt3 = f.ntsnp.substr(0,1);											// get first nt
				nt4 = f.ntsnp.substr(2,1);											// get second nt
				if (is_match(v.ref, v.var, nt3, nt4, f.sen))
					typ_snp=1;
			}
		}
	}
		
	return val;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

// this function does the printing
void print_row(row &r, map<string, string> &g, string &chromosome, vector<threestr> &m0, vector<threestr> &m1, vector<string> &m2, vector<string> &m3, vector<string> &m4, map<string, bool> &flag_map, vector<string> &flag_order, int &variant_row)
{
	
	stringstream ccds_exon_stream;		// stringstream for output of -exon flag
	stringstream ccds_intron_stream;	// stringstream for output of -intron flag
	stringstream ccds_splice_stream;	// stringstream for output of -splice flag
	
	stringstream t_exon_stream;			// stringstream for output of -texon flag
	stringstream t_intron_stream;		// stringstream for output of -tintron flag
	stringstream t_splice_stream;		// stringstream for output of -tsplice flag
	
	stringstream t_promoter_stream;		// stringstream for output of -tpromoter flag
	
	stringstream mirna_stream;			// stringstream for output of -mirna flag	
		
	string straa = "";					// concat aa_change(s) 
	string strcdn = "";					// concat codons
		
	string cdn;						// string for codon to send into get_codon function
	string aa;						// string for aa
	
	//string snp_hit;					// an intersecting SNP
	string strsnp = "";				// concat SNPs 
	string strbit = "";				// concat bitfield 
			
	stringstream strac, stran, strdp;	// concat thousand genome allele count, allele number, depths
	
	stringstream row_flags; 			// the flags portion of the row
	stringstream row_out; 				// the non-flags portion of the row
	
	variant v; 						// variant holder so dont have to do repeated lookups

	// flags
	bool flag_nsyn;					// 1 if nsyn, 0 if not
	
	// number of times variant intersects w each feature
	int n_exon = 0;
	int n_intron = 0;
	int n_splice = 0;
	int n_texon = 0;
	int n_tintron = 0;
	int n_tsplice = 0;
	int n_promoter = 0;
	// these are more complicated - have to check for intersection
	int n_thousand;
	int n_snp;
	int n_mirna;
	int n_tot;							// total number of all features
	short typ_snp; 						// SNP feature type: 0 no intersect, 1 is_true_snp, 2 is_on_snp	
		
	// --------- deal with coding part ---------- //	
	// ***these just depend on position, not the specific variant***
	
	// concat coding gene names, ccds names, senses, and exon numbers onto strings
	if (flag_map["exon"])
	{
		n_exon = r.v_exon.size();
		get_names(r.v_exon, m1, ccds_exon_stream);
	}
	
	if (flag_map["intron"])
	{
		n_intron = r.v_intron.size();
		get_names(r.v_intron, m1, ccds_intron_stream);
	}
	
	if (flag_map["splice"])
	{
		n_splice = r.v_splice.size();
		get_names(r.v_splice, m1, ccds_splice_stream);		
	}	
	// concat transcript accession numbers, senses, and exon numbers onto strings
	if (flag_map["texon"])
	{
		n_texon = r.v_texon.size();
		get_t_names(r, m0, t_exon_stream);
	}
	
	if (flag_map["tintron"])
	{
		n_tintron = r.v_tintron.size();
		get_names(r.v_tintron, m0, t_intron_stream);
	}
	
	if (flag_map["tsplice"])
	{
		n_tsplice = r.v_tsplice.size();
		get_names(r.v_tsplice, m0, t_splice_stream);
	}	
	
	// concat promoter transcript accession numbers onto string
	if (flag_map["tpromoter"])
	{
		n_promoter = r.v_promoter.size();
		get_p_names(r, m0, t_promoter_stream);
	}

	if (flag_map["mirna"])
	{
		n_mirna = r.v_mirna.size();
		get_mirna_names(r, m4, mirna_stream);
	}			
	// i : loop over all variants at a given position - each variant gets a row
	for(int i = 0; i < r.v_variant.size(); i++)	
	{	
		// ***these depend on the specific variant***		
		v = r.v_variant[i];

		n_tot = n_exon + n_intron + n_splice + n_texon + n_tintron + n_tsplice + n_promoter + n_mirna;	// reset this value for each variant. add SNPs & thous below b/c we dont know the number yet
		n_thousand = 0;																			// reset
		n_snp = 0;																				// reset
		typ_snp = 0;																			// reset
		flag_nsyn = false;																		// reset
		row_out.clear();																		// reset
		row_out.str("");
		row_flags.clear();																		// reset
		row_flags.str("");
		
//		if (!flag_map["position_only"])
//		{
		if (flag_map["snp"])
		{
			// k : loop over all SNPs for a given variant
			for(int k = 0; k < r.v_snp.size(); k++)	
			{	
				// sets strsnp and typ_snp
				n_snp = n_snp + is_snp(r.v_snp[k], v, m2, strsnp, strbit, typ_snp, flag_map["medical"]);
			}	
			//n_tot = n_tot + n_snp;
			if (n_snp > 0) {n_tot = n_tot + 1;}		// for now, treat SNPs as a 0 or 1 thing
		}
		
		if (flag_map["thousand"])
		{
			// l : loop over all thous for a given variant
			for(int l = 0; l < r.v_thousand.size(); l++)	
			{	
				n_thousand = n_thousand + is_thou(r.v_thousand[l], v, strac, stran, strdp, n_thousand);
			}	
			n_tot = n_tot + n_thousand;
		}
			
		if (flag_map["translate_aa"])
		{
			// j : loop over all ccds's (coding exons) for a given variant
			for(int j = 0; j < r.v_exon.size(); j++)	
			{			
				// do a bunch of things simultaneously - get the codon (cdn), the AA (aa), and set the nonsyn flag
				flag_nsyn = (get_codon(r.v_exon[j], v, m3[r.v_exon[j].tid-1], g, m1[r.v_exon[j].tid-1].s1, cdn, aa, variant_row) || flag_nsyn);
				if (j==0)
				{	
					strcdn = cdn;
					straa = aa;
				}
				else
				{
					strcdn = strcdn + "," + cdn;
					straa = straa + "," + aa;
				}
			}
		}
//		}				
			
		// set up to print output
		// accumulate flags in row flags 
		// and output in row out	
		for(int ii = 0; ii < flag_order.size(); ii++)
		{
	
			if (flag_order[ii] == "exon")
			{
				row_out << ccds_exon_stream.str() << "\t" << strcdn << "\t" << straa;
				if (ii != flag_order.size() - 1) {row_out << "\t";}
				row_flags << n_exon << "\t" << flag_nsyn << "\t";
			}
			else if (flag_order[ii] == "intron")
			{			
				row_out << ccds_intron_stream.str();
				if (ii != flag_order.size() - 1) {row_out << "\t";}
				row_flags << n_intron << "\t";
			}
			else if (flag_order[ii] == "splice")
			{
				row_out << ccds_splice_stream.str();
				if (ii != flag_order.size() - 1) {row_out << "\t";}
				row_flags << n_splice << "\t";							
			}						
			else if (flag_order[ii] == "texon")
			{
				row_out << t_exon_stream.str();
				if (ii != flag_order.size() - 1) {row_out << "\t";}
				row_flags << n_texon << "\t";							
			}						
			else if (flag_order[ii] == "tintron")
			{
				row_out << t_intron_stream.str();
				if (ii != flag_order.size() - 1) {row_out << "\t";}
				row_flags << n_tintron << "\t";							
			}						
			else if (flag_order[ii] == "tsplice")
			{			
				row_out << t_splice_stream.str();
				if (ii != flag_order.size() - 1) {row_out << "\t";}
				row_flags << n_tsplice << "\t";
			}						
			else if (flag_order[ii] == "tpromoter")
			{			
				row_out << t_promoter_stream.str();
				if (ii != flag_order.size() - 1) {row_out << "\t";}
				row_flags << n_promoter << "\t";
			}						
			else if (flag_order[ii] == "snp")
			{			
				row_out << strsnp;
				if (flag_map["medical"]) {row_out << "\t" << strbit;}
				if (ii != flag_order.size() - 1) {row_out << "\t";}
				row_flags << typ_snp << "\t";				
			}
			else if (flag_order[ii] == "mirna")
			{			
				row_out << mirna_stream.str();
				if (ii != flag_order.size() - 1) {row_out << "\t";}
				row_flags << n_mirna << "\t";
			}									
			else if (flag_order[ii] == "thousand")
			{	
				row_out << strac.str() << "\t" << stran.str() << "\t" << strdp.str();
				if (ii != flag_order.size() - 1) {row_out << "\t";}
				row_flags << n_thousand << "\t";						
			}						
		}

		// print output
		cout << row_flags.str() << n_tot << "\t";
		cout << "chr" << chromosome << ":" << v.pos << "\t";
		cout << v.ref << "/" << v.var << "\t";
		cout << row_out.str() << endl;

	}
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

// this function does the printing for rna rows
void print_rna_row(row &r, map<string, string> &g, vector<twostr> &m0, vector<threestr> &m1, vector<fivestr> &m2, int &variant_row)
{
	// m0 is seq_transcript (txid \t seq)
	// m1 is seg_transcript_coding (txid \t coding start \t coding end)
	// m2 is name_transcript (txid \t chr \t acc id \t gene \t sen)
			
	string cdn = "";				// string for codon
	string aa = "";					// string for aa
	
	rnaexon this_rnaexon;			// the current rna exon
	unsigned rnaid;					// id
	
	unsigned strtcode, endcode;		// coding region begin end
		
	// flags
	bool flag_nsyn;					// 1 if nsyn, 0 if not
	bool flag_coding;				// 1 if coding, 0 if not
	
	variant v; 						// variant holder so dont have to do repeated lookups
	
	for(int i = 0; i < r.v_variant.size(); i++)	
	{		
		v = r.v_variant[i];
		
		flag_nsyn = false;			// reset
		flag_coding = false;		// reset			
		
		// make sure it's in a known transcript
		if (r.v_rnaexon.size()==0)
			cerr << "error - transcript not found at " << v.pos << "\n";
		else
		{
			// transcript variants must only intersect with a single uniq transcript
			this_rnaexon=r.v_rnaexon[0];
			
			// cast string to unsigned
			strtcode = atoi(m1[this_rnaexon.tid-1].s2.c_str());
			endcode = atoi(m1[this_rnaexon.tid-1].s3.c_str());
			
			//m1
			
			// check if coding and, if so, get codon
			if (v.pos >= strtcode && v.pos <= endcode)
			{
				flag_coding = true;
				// do a bunch of things simultaneously - get the codon (cdn), the AA (aa), and set the nonsyn flag
				flag_nsyn = (get_rna_codon(this_rnaexon, v, m0[this_rnaexon.tid-1].s2, g, this_rnaexon.txid, cdn, aa, strtcode, endcode, variant_row) || flag_nsyn);				
			}	 
			
			// print output
			cout << flag_coding << "\t" << flag_nsyn << "\t";
			cout << m2[this_rnaexon.tid-1].s3 << ":" << v.pos << "\t";
			cout << v.ref << "/" << v.var << "\t";
			cout << m2[this_rnaexon.tid-1].s2 << "\t" << m2[this_rnaexon.tid-1].s4 << "\t" << m2[this_rnaexon.tid-1].s5 << "\t" << this_rnaexon.num << "\t";
			cout << cdn << "\t" << aa;
			cout << endl;
		}
	}	
	
	
}

// this function prints context only
void print_context(row &r, string &full_seq, int context_side)
{
	
	string cntxt;	// a string for the context
	
	// i : loop over all variants at a given position - each variant gets a row
	for(int i = 0; i < r.v_variant.size(); i++)	
	{	
		get_context(r.v_variant[i].pos, r.v_variant[i].pos, full_seq, cntxt, context_side);
		cout << cntxt << endl;
	}	
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void main_annotate(int argc, char **argv) 
{	
	/*
	arguments:
	
	$1 humangencode.txt (file: codon \t AA)
	$2 chromosome name
	$3 ccds/name_transcript (file: refseq accession id \t gene \t strand)
	$4 ccds/name (file: CCDSid \t gene \t strand)
	$5 snp/name (file: SNP name) 
	$6 CCDS sequences (file: CCDS seq)
	$7 chr sequeces (file: chr seq)
	$8 context boolean: 0 for standard output, 1 for contexts only
	*/
	
	// to compile: g++ -O2 main_mod_annotate.cpp -o modannotate
	
	// to run, e.g.:
	// c=1; (output of segsone 100, e.g., 4a) | modannotate -h -t /ifs/home/c2b2/rr_lab/shares/scripts/humangencode.txt /ifs/scratch/c2b2/rr_lab/shares/ref/hg18/results/chr/${c}/ccds/seq -c ${c} -transcript /ifs/scratch/c2b2/rr_lab/shares/ref/hg18/results/chr/${c}/ccds/name_transcript -ccds /ifs/scratch/c2b2/rr_lab/shares/ref/hg18/results/chr/${c}/ccds/name -exon -intron -splice -texon -tintron -tsplice -tpromoter -thousand -snp /ifs/scratch/c2b2/rr_lab/shares/ref/hg18/results/chr/${c}/snp/name
	
	// to run for rna:
	// (output of segsone 100 -i) | modannotate -h -rna /ifs/scratch/c2b2/rr_lab/shares/ref/hg19/results/refseq_transcriptome/seq_transcript /ifs/scratch/c2b2/rr_lab/shares/ref/hg19/results/refseq_transcriptome/seg_transcript_coding
	
	row r;												// a row structure
	string line;										// a line piped in
	istringstream iss;									// input stream
	unsigned vid, vprev, pos, empty, typ, tid;			// variant_id, prev variant_id, position, an empty holder, feature type
	bool fl = true;										// first line flag (first line starts as true)
	string chromosome;									// the chr 
	string transcript_id;								// the transcript id in the case of RNA
	int context_side;			 						// length of the context on one side
	int variant_row = 1;								// a variant row counter, which starts as one (i.e., everytime a new variant comes up, this is incremented)
	
	// these are a set of boolean flags set by user input to modularize the output of this program
	map<string, bool> flag_map;	
	flag_map["header"] = false;							// turns on header
	flag_map["context"] = false;	
	flag_map["exon"] = false;
	flag_map["intron"] = false;
	flag_map["splice"] = false;	
	flag_map["texon"] = false;
	flag_map["tintron"] = false;	
	flag_map["tsplice"] = false;	
	flag_map["tpromoter"] = false;	
	flag_map["snp"] = false;		
	flag_map["mirna"] = false;	
	flag_map["thousand"] = false;		
	flag_map["translate_aa"] = false;	
	flag_map["position_only"] = false;					// this flag turns on the option of just giving nucleotides (not ref and var nucleotides)
	flag_map["medical"] = false;						// this flag should be ON for hg19 SNPs that have an extra column - the bitfield specifying clinical assoc
	flag_map["rna"] = false;							// this flag goes hi when using rna data
	
	// this vector tracks the order of the flags so when we print the output we know what order they came in
	vector<string> flag_order;
	
	// features
	variant my_variant;				
	exon my_exon;					
	intron my_intron;					
	splice my_splice;
	texon my_texon;	
	tintron my_tintron;	
	tsplice my_tsplice;				
	snp my_snp;						
	promoter my_promoter;			
	thousand my_thousand;
	rnaexon my_rnaexon;
	mirna my_mirna;
	
	// args
	map<string, string> g;						// codon to amino acid mapping
	vector<threestr> m0;						// a vector of 3 strings (refseq accession id, genename, strand)
	vector<threestr> m1;						// a vector of 3 strings (ccdsname, genename, strand)
	vector<string> m2;							// a vector of strings of snpname
	vector<string> m3;							// a vector of strings of ccds_seq
	vector<string> m4;							// a vector of strings of mirna name
	string full_seq;							// full chr seq for the context		

	// for rna:
	vector<twostr> refseq_ref;					// a vector of 2 strings for the transcript transcript ref (transcript id, seq)
	vector<threestr> coding_region;				// a vector of 3 strings for the coding region (transcript id, start, end)
	vector<fivestr> rna_name;					// a vector of 5 strings for the transcript names (transcript id, chr, assession id, gene, sense)
	
	stringstream row_header_flags; 				// the flags portion of the row header
	stringstream row_header_out; 				// the non-flags portion of the row	header
	
	// process flags
	// ADD INTERNAL FLAGS
	for(int i = 1; i < argc;)
	{	
		if (strcmp(argv[i],"-c") == 0)
		{
			// chr flag (e.g., ${c})
			chromosome = argv[i+1];				// chr name
			i=i+2;			
		}
		else if (strcmp(argv[i],"-h") == 0)
		{
			// header flag
 			flag_map["header"] = true;
			i++;
		}
		else if (strcmp(argv[i],"-t") == 0)
		{
			// translate flag - so read gene code (e.g., /ifs/home/c2b2/rr_lab/shares/scripts/humangencode.txt) and ccds seq (e.g., /ifs/scratch/c2b2/rr_lab/shares/ref/hg18/results/chr/${c}/ccds/seq)
			read_gene_code(argv[i+1], g); 		// read humangencode.txt
 			read_stringfile(argv[i+2], m3);		// ccds/seq
 			flag_map["translate_aa"] = true;
			i=i+3;
		}
		else if (strcmp(argv[i],"-p") == 0)
		{
			// position only flag
			flag_map["position_only"] = true;
			i++;
		}
		else if (strcmp(argv[i],"-context") == 0)
		{
			
			// context flag (e.g., /ifs/scratch/c2b2/rr_lab/shares/ref/hg18/results/chr/${c}/seq)
			
		   	ifstream myfile (argv[i+1]);					// chr/seq
		  	if (myfile.is_open())
		  	{
		     	getline (myfile,full_seq);
		    	myfile.close();
		  	}
			else cerr << "Unable to open chr seq file";
			
			context_side = atoi(argv[i+2]); 				// length of the context on one side
			
			flag_map["context"] = true; 	// flag			
			
			// add context length here
			
			i=i+3;
		}
		else if (strcmp(argv[i],"-rna") == 0)
		{
			
			// rna flag
			
			flag_map["rna"] = true; 	// flag	
					
			read_gene_code(argv[i+1], g); 					// read humangencode.txt
			read_2string(argv[i+2], refseq_ref);			// read the refseq ref seqs		
			read_3string(argv[i+3], coding_region);			// read the coding region
			read_5string(argv[i+4], rna_name);				// read the rna names file		
			
			i=i+5;
		}		
		else if (strcmp(argv[i],"-ccds") == 0)
		{
			// get ccds names (e.g., /ifs/scratch/c2b2/rr_lab/shares/ref/hg18/results/chr/${c}/ccds/name)	
			read_3string(argv[i+1], m1);					// ccds/name
			i=i+2;
		}
		else if (strcmp(argv[i],"-exon") == 0)
		{
			flag_map["exon"] = true;
			flag_order.push_back("exon");
			row_header_flags << "#exon" << "\t" << "nsyn" << "\t";
			row_header_out << "exon_gene" << "\t" << "exon_CCDS" << "\t" << "exon_sense" << "\t" << "exon_#" << "\t" << "codon" << "\t" << "AA" << "\t";
			i++;
		}
		else if (strcmp(argv[i],"-intron") == 0)
		{
			flag_map["intron"] = true;
			flag_order.push_back("intron");
			row_header_flags << "#intron" << "\t";
			row_header_out << "intr_gene" << "\t" << "intr_CCDS" << "\t" << "intr_sense" << "\t" << "intr_#" << "\t";			
			i++;
		}
		else if (strcmp(argv[i],"-splice") == 0)
		{
			flag_map["splice"] = true;
			flag_order.push_back("splice");
			row_header_flags << "#spl" << "\t";			
			row_header_out << "spl_gene" << "\t" << "spl_CCDS" << "\t" << "spl_sense" << "\t" << "spl_exon_#" << "\t";
			i++;
		}
		else if (strcmp(argv[i],"-transcript") == 0)
		{
			// get ccds names (e.g., /ifs/scratch/c2b2/rr_lab/shares/ref/hg18/results/chr/${c}/ccds/name_transcript)	
			read_3string(argv[i+1], m0);
			i=i+2;
		}		
		else if (strcmp(argv[i],"-texon") == 0)
		{
			flag_map["texon"] = true;
			flag_order.push_back("texon");
			row_header_flags << "#tr_ex" << "\t";
			row_header_out << "texon_gene" << "\t" << "texon_acc" << "\t" << "texon_sense" << "\t" << "texon_#" << "\t" << "texon_position" << "\t";					
			i++;
		}
		else if (strcmp(argv[i],"-tintron") == 0)
		{
			flag_map["tintron"] = true;
			flag_order.push_back("tintron");
			row_header_flags << "#tr_in" << "\t";
			row_header_out << "tintron_gene" << "\t" << "tintron_acc" << "\t" << "tintron_sense" << "\t" << "tintron_#" << "\t";						
			i++;
		}
		else if (strcmp(argv[i],"-tsplice") == 0)
		{
			flag_map["tsplice"] = true;
			flag_order.push_back("tsplice");
			row_header_flags << "#tr_spl" << "\t";
			row_header_out << "tsplice_gene" << "\t" << "tsplice_acc" << "\t" << "tsplice_sense" << "\t" << "tsplice_exon_#" << "\t";
			i++;
		}		
		else if (strcmp(argv[i],"-tpromoter") == 0)
		{
			flag_map["tpromoter"] = true;
			flag_order.push_back("tpromoter");			
			row_header_flags << "#tr_pro" << "\t";
			row_header_out << "pro_gen" << "\t" << "pro_acc" << "\t" << "pro_sen" << "\t" << "pro_pos" << "\t";		
			i++;
		}
		else if (strcmp(argv[i],"-snp") == 0)
		{
			// cerr << "reading snp names" << endl;
			read_stringfile(argv[i+1], m2);				// snp/name
			flag_map["snp"] = true;
			flag_order.push_back("snp");
			row_header_flags << "typeSNP" << "\t";
			row_header_out << "SNP" << "\t";												
			i=i+2;
		}
		else if (strcmp(argv[i],"-mirna") == 0)
		{
			read_stringfile(argv[i+1], m4);				// mirna/name
			flag_map["mirna"] = true;
			flag_order.push_back("mirna");
			row_header_flags << "#mirna" << "\t";
			row_header_out << "mirna_gen" << "\t" << "mirna_sen" << "\t";												
			i=i+2;
		}
		else if (strcmp(argv[i],"-medsnp") == 0)
		{
			// cerr << "reading snp names" << endl;
			read_stringfile(argv[i+1], m2);				// snp/name
			flag_map["snp"] = true;
			flag_map["medical"] = true;
			flag_order.push_back("snp");
			row_header_flags << "typeSNP" << "\t";
			row_header_out << "SNP" << "\t" << "bitfield" << "\t";												
			i=i+2;
		}
		else if (strcmp(argv[i],"-thousand") == 0)
		{
			flag_map["thousand"] = true;	
			flag_order.push_back("thousand");
			row_header_flags << "#1000" << "\t";
			row_header_out << "1000_AC" << "\t" << "1000_AN" << "\t" << "1000_DP" << "\t";				
			i++;
		}
		else
			i++;		
	}
	
	// if not context, print header
	if (flag_map["context"] && flag_map["header"])
	{
		cout << "context" << endl;
	}
	if (flag_map["rna"] && flag_map["header"])
	{
		cout << "code" << "\t" << "nsyn" << "\t" << "tx:pos" << "\t" << "ref/var" << "\t" << "chr" << "\t" << "gene" << "\t" << "exon_sense" << "\t" << "exon_#" << "\t" << "codon" << "\t" << "AA" << "\t" << endl;
	}	
	else if (flag_map["header"])
	{
		// remove trailing tab (last character) of row_header_out
		cout << row_header_flags.str() << "#tot" << "\t" << "chr:pos" << "\t" << "ref/var" << "\t" << row_header_out.str().erase(row_header_out.str().size()-1,1) << endl;	
		// FIX THIS LATER (UGLY!)
		
	}
	
	while (true)
	{
				
		// read line, get current vid
		std::getline(std::cin, line);
		
		// read file piped in until hit an empty line				
		if (line.empty())
		{
			// print last 
			if (flag_map["context"])
				print_context(r, full_seq, context_side);
			else if (flag_map["rna"])
				print_rna_row(r, g, refseq_ref, coding_region, rna_name, variant_row);
			else
				print_row(r, g, chromosome, m0, m1, m2, m3, m4, flag_map, flag_order, variant_row);
						
			break;
		}
		else 
		{
			iss.clear();
			iss.str(line);	

			if (flag_map["rna"])
				iss >> vid >> ws >> transcript_id >> pos >> ws >> empty >> ws >> typ;		
			else
				iss >> vid >> ws >> pos >> ws >> empty >> ws >> typ;
						
			// print the previous row as soon as vid changes, then flush row
			if (!fl && vid!=vprev) 																		
			{
				if (flag_map["context"])
					print_context(r, full_seq, context_side);
				else if (flag_map["rna"])
					print_rna_row(r, g, refseq_ref, coding_region, rna_name, variant_row);					
				else
					print_row(r, g, chromosome, m0, m1, m2, m3, m4, flag_map, flag_order, variant_row);
					
				r.clear_row();
				
				variant_row++; // increment variant counter
			}	
			vprev = vid;
			fl = false;					

  			switch(typ) 
  			{
  				case VARIANT:
  					my_variant.set0(pos, typ);
  					my_variant.set(iss);
  					r.v_variant.push_back(my_variant);
  					break;    					
  				case EXON:
  					my_exon.set0(pos, typ);
  					my_exon.set(iss);
  					r.v_exon.push_back(my_exon);
  					break;
  				case RNAEXON:
  					my_rnaexon.set0(pos, typ);
  					my_rnaexon.set1(transcript_id);
  					my_rnaexon.set(iss);
  					r.v_rnaexon.push_back(my_rnaexon);
  					break;   						
  				case TEXON:
  					my_texon.set0(pos, typ);
  					my_texon.set(iss);
  					r.v_texon.push_back(my_texon);
  					break;
  				case INTRON:
  					my_intron.set0(pos, typ);
  					my_intron.set(iss);
  					r.v_intron.push_back(my_intron);
  					break;  					
  				case TINTRON:  				
  					my_tintron.set0(pos, typ);
  					my_tintron.set(iss);
  					r.v_tintron.push_back(my_tintron);
  					break;
  				case LSS:
  				case RSS:
  					my_splice.set0(pos, typ);
  					my_splice.set(iss);
  					r.v_splice.push_back(my_splice);
  					break;  					
  				case TLSS:
  				case TRSS:  				  				
  					my_tsplice.set0(pos, typ);
  					my_tsplice.set(iss);
  					r.v_tsplice.push_back(my_tsplice);
  					break;  					
  				case SNP:
  					my_snp.set0(pos, typ);
  					if (flag_map["medical"])
  					{
  						my_snp.set2(iss);
  					}
  					else
  					{
  						my_snp.set(iss);
  					}
  					r.v_snp.push_back(my_snp);
  					break;
  				case MIRNA:
  					my_mirna.set0(pos,typ);
  					my_mirna.set(iss);
  					r.v_mirna.push_back(my_mirna);
  					break;
  				case PROMOTER:
  					my_promoter.set0(pos, typ);
  					my_promoter.set(iss);
  					r.v_promoter.push_back(my_promoter);
  					break;
  				case THOU_LOW:
  				case THOU_EXON:
  				case THOU_TRIO:
  					my_thousand.set0(pos, typ);
  					my_thousand.set(iss, typ);
  					r.v_thousand.push_back(my_thousand);
  					break;	
  			}
		}		
	}	
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    
int main(int argc, char **argv) 
{    
	main_annotate(argc, argv);
  	return 0;
}
