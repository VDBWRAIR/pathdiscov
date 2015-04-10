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
// #include <regex>
// #include <boost/regex.hpp>

// turn debug on:
// #define DEBUG true

using namespace std;

/*
Name	Description	Size*	Range*
char		Character or small integer.	1byte	signed: -128 to 127 unsigned: 	0 to 255
short int 	(short)	Short Integer.		2bytes	signed: -32768 to 32767 unsigned: 	0 to 65535
int			Integer.					4bytes	signed: -2147483648 to 2147483647 unsigned: 0 to 4294967295
long int (long)	Long integer.			4bytes	signed: -2147483648 to 2147483647 unsigned: 0 to 4294967295
bool	Boolean value. It can take one of two values: true or false.	1byte	true or false
float	Floating point number.	4bytes	+/- 3.4e +/- 38 (~7 digits)
double	Double precision floating point number.	8bytes	+/- 1.7e +/- 308 (~15 digits)
long double	Long double precision floating point number.	8bytes	+/- 1.7e +/- 308 (~15 digits)
wchar_t	Wide character.	2 or 4 bytes	1 wide character
*/

/*
// preprocessor directives
#define #error #import	#undef #elif #if #include #using #else #ifdef #line #endif #ifndef #pragma
*/

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//******************************************//
// file reading functions					//
//******************************************//

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
  	else cout << "Unable to open file" << endl;   		
}

/* - - - - - - - - - - - - - - - - - - - - */

// read a fastq file and store IDs as a vector of strings (every 4 lines, get ID)
void read_fastq_id(char *s, vector<string> &fastq_id) 
{	
	string line;
	string lineid;	
	istringstream iss;

	ifstream myfile (s);
	if (myfile.is_open())
  	{
  		
  		int nr = 1; // a variable for row number
  		
    	while (!myfile.eof())
    	{    		
         	getline(myfile,line);
         	iss.clear();
         	iss.str(line);
         	if (!line.empty())    
         	{   
				if (nr%4==1)
				{
					lineid = line.substr (1,string::npos);	// skip first char and go to end of line
					fastq_id.push_back(lineid);
				}
				nr++;
         	}
    	}
    	myfile.close();
    }
  	else cout << "Unable to open file" << endl;   		
}

/* - - - - - - - - - - - - - - - - - - - - */

// read a fasta file and store IDs as a vector of strings (ASSUME IDs begin with ">")
void read_fasta_id(char *s, vector<string> &fasta_id) 
{	
	string line;
	string lineid;	
	istringstream iss;

	ifstream myfile (s);
	if (myfile.is_open())
  	{
  		
  		int nr = 1; // a variable for row number
  		
    	while (!myfile.eof())
    	{    		
         	getline(myfile,line);
         	iss.clear();
         	iss.str(line);
         	if (!line.empty())    
         	{   
				if (line.compare(0,1,">") == 0)
				{
					lineid = line.substr (1,string::npos);	// skip first char and go to end of line
					fasta_id.push_back(lineid);
				}
				nr++;
         	}
    	}
    	myfile.close();
    }
  	else cout << "Unable to open file" << endl;   		
}

/* - - - - - - - - - - - - - - - - - - - - */

// read a fastq file and store IDs as a map (every 4 lines, get ID)
void read2map_fastq_id(char *s, std::map<string, bool> &fastq_id) 
{	
	string line;
	string lineid;	
	istringstream iss;

	ifstream myfile (s);
	if (myfile.is_open())
  	{
  		
  		int nr = 1; // a variable for row number
  		
    	while (!myfile.eof())
    	{    		
         	getline(myfile,line);
         	iss.clear();
         	iss.str(line);
         	if (!line.empty())    
         	{   
				if (nr%4==1)
				{
					lineid = line.substr (1,string::npos);	// skip first char and go to end of line
					fastq_id[lineid] = true;
				}
				nr++;
         	}
    	}
    	myfile.close();
    }
  	else cout << "Unable to open file" << endl;   		
}

/* - - - - - - - - - - - - - - - - - - - - */

// read a fasta file and store IDs as a vector of strings (ASSUME IDs begin with ">")
void read2map_fasta_id(char *s, std::map<string, bool> &fasta_id) 
{	
	string line;
	string lineid;	
	istringstream iss;

	ifstream myfile (s);
	if (myfile.is_open())
  	{
  		
  		int nr = 1; // a variable for row number
  		
    	while (!myfile.eof())
    	{    		
         	getline(myfile,line);
         	iss.clear();
         	iss.str(line);
         	if (!line.empty())    
         	{   
				if (line.compare(0,1,">") == 0)
				{
					lineid = line.substr (1,string::npos);	// skip first char and go to end of line
					fasta_id[lineid] = true;
				}
				nr++;
         	}
    	}
    	myfile.close();
    }
  	else cout << "Unable to open file" << endl;   		
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//******************************************//
// functions								//
//******************************************//

// print IDs
void print_id(std::map<string, bool> &file_id) 
{	
    for(std::map<string, bool>::iterator iter = file_id.begin(); iter != file_id.end(); ++iter)
    {
    	// http://stackoverflow.com/posts/1443798/edit
		// map is associative container. Hence, iterator is a pair of key,val. IF you need only keys, you can ignore the value part from the pair
	    
	    string key = iter->first;
	    //ignore value
	    //bool v = iter->second;
	    
		cout << key << endl;
    }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

// does 3 things: (1) change fastq IDs to numbers; (2) truncate every 3rd line to "+" if nec; (3) replace dots in seq w Ns
// output 2 files: (1) file with new IDs (2) mapping of former numeric IDs to new numberical IDs 
void change_fastqid_2_numeric(char *s) 
{	
	string line;
	istringstream iss;

	ifstream myfile (s);
	if (myfile.is_open())
  	{
  		
  		int nr = 1; // a variable for row number
  		
    	while (!myfile.eof())
    	{    		
         	getline(myfile,line);
         	iss.clear();
         	iss.str(line);
         	if (!line.empty())    
         	{   
				if (nr%4==1)
				{
					cout << "@";
					cout << (nr+3)/4 << endl;
				}
				else
				{
					cout << line << endl; 
				}
				nr++;
         	}
    	}
    	myfile.close();
    }
  	else cout << "Unable to open file" << endl; 		
}

/* - - - - - - - - - - - - - - - - - - - - */

// print IDs in file1 but not in file2
void get_file1_not_file2_id(std::map<string, bool> &file1_id, std::map<string, bool> &file2_id) 
{	
    for(std::map<string, bool>::iterator iter = file1_id.begin(); iter != file1_id.end(); ++iter)
    {
    	// http://stackoverflow.com/posts/1443798/edit
		// map is associative container. Hence, iterator is a pair of key,val. IF you need only keys, you can ignore the value part from the pair
	    
	    string key = iter->first;
	    //ignore value
	    //bool v = iter->second;
	    
	    if (!file2_id[key])
	    {
			cout << key << endl;
	    }	    
    }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//******************************************//
// get options								//
//******************************************//

void get_opt(int argc, char **argv) 
{
	// process flags
	;
}

void check_args(int argc, char **argv) 
{

	if (argc == 1 || strcmp(argv[1],"--help") == 0 || strcmp(argv[1],"-help") == 0 || strcmp(argv[1],"-h") == 0)
	{
		cout << endl;
		cout << "NAME" << endl;
		cout << "\t" << "fastaq_tools_diff.exe - print the ID difference between two fastq/fasta files - to be precise, print what's in first file and not second file" << endl;
		cout << endl;
		cout << "SYNOPSIS" << endl;
		cout << "\t" << "fastaq_tools_diff.exe [--fastq fastq_file] [--fasta fasta_file]" << endl;
		cout << endl;
		cout << "EXAMPLES" << endl; 
		cout << "\t" << "# argument order matters!" << endl;
		cout << "\t" << "fastaq_tools_diff.exe --fastq 1.fastq --fasta 1.fasta" << endl;
		cout << "\t" << "fastaq_tools_diff.exe --fastq 1.fastq --fastq 2.fastq" << endl;
		cout << "\t" << "fastaq_tools_diff.exe --fasta 1.fasta --fasta 2.fasta" << endl;
		cout << endl;
		cout << "COMPILATION" << endl;
		cout << "\t" << "g++ -O2 fastaq_tools.cpp -o fastaq_tools_diff.exe" << endl;
		cout << endl;
		exit; 
	}

}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void main_program(int argc, char **argv) 
{
	// Declare variables

//	string file1;						// string for file 1 name
//	string file2;						// string for file 2 name

	// containers for all possible file types:
	std::map<string, bool> file1_id;		// hold file 1 fastq IDs
	std::map<string, bool> file2_id;		// hold file 2 fasta IDs
	
	bool other_options = false; 			// this is a bool flag to go hi in the case youre doing something OTHER THAN a file difference
	
	// store IDs as vectors	
	/*
	
	// containers for all possible file types:
	vector<string> file1_id;		// hold file 1 fastq IDs
	vector<string> file2_id;		// hold file 2 fasta IDs	
	
	bool file1_isfastq;					// boolean for fastq
	bool file2_isfastq;					// boolean for fastq	
	
	// process flags
	for(int i = 1; i < argc;)
	{	
		if (strcmp(argv[i],"--fastq") == 0)
		{
			// if first arg set file1 else set file2
			if (i == 1)
			{
				file1_isfastq = 1;
				read_fastq_id(argv[i+1], file1_id);				
			}
			else
			{
				file2_isfastq = 1;
				read_fastq_id(argv[i+1], file2_id);
			}
			
			i=i+2;			
		}
		else if (strcmp(argv[i],"--fasta") == 0)
		{
			if (i == 1)
			{
				file1_isfastq = 0;
				read_fasta_id(argv[i+1], file1_id);				
			}
			else
			{
				file2_isfastq = 0;
				read_fasta_id(argv[i+1], file2_id);
			}
			
			i=i+2;
		}
	}
	
	cout << "*** file1 ***" << endl;	
	// i : loop over IDs
	for(int i = 0; i < file1_id.size(); i++)	
	{	
		// cout << i << endl;
		cout << file1_id[i] << endl;		
	}

	cout << "*** file2 ***" << endl;	
	// i : loop over IDs
	for(int i = 0; i < file2_id.size(); i++)	
	{	
		// cout << i << endl;
		cout << file2_id[i] << endl;		
	}

	get_id_diff(file1_id, file2_id);	
	*/

	check_args(argc, argv);	
	
	// store IDs as maps
	// process flags
	for(int i = 1; i < argc;)
	{	
		
		// cout << i << endl;
				
		// if fastq for difference function
		if (strcmp(argv[i],"--fastq") == 0)
		{
			// if first arg set file1 else set file2
			if (i == 1)
			{
				read2map_fastq_id(argv[i+1], file1_id);				
			}
			else
			{
				read2map_fastq_id(argv[i+1], file2_id);
			}
			
			i=i+2;			
		}
		// if fasta	for difference function	
		else if (strcmp(argv[i],"--fasta") == 0)
		{
			// if first arg set file1 else set file2
			if (i == 1)
			{
				read2map_fasta_id(argv[i+1], file1_id);				
			}
			else
			{
				read2map_fasta_id(argv[i+1], file2_id);
			}
			
			i=i+2;
		}
		// convert fastq IDs to numbers
		else if (strcmp(argv[i],"--fastq_id2num") == 0)
		{
			other_options = true;
			change_fastqid_2_numeric(argv[i+1]);			
			i++;
		}
		// convert fastq IDs to numbers
		else 
		{
			i++;
		}

	}

	if (!other_options)
	{
		#ifdef DEBUG
			cout << "*** file1 ***" << endl;	
			print_id(file1_id);
			cout << "*** file2 ***" << endl;	
			print_id(file2_id);	
		#endif

		get_file1_not_file2_id(file1_id,file2_id);	
	}
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    
int main(int argc, char **argv) 
{
	// g++ -O2 fastaq_tools.cpp -o fastaq_tools_diff.exe
	
	main_program(argc, argv);
  	return 0;
}
