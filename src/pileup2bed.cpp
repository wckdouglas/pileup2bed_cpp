#include <iostream>
#include <vector>
#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <cstdlib>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <unordered_map>
#include "stringManipulation.h"
#include "pileupFix.h"
#include "stringVector.h"

typedef unordered_map<string , int> count_dict;

//print usage
int usage(char *argv[])
{
	cerr << "usage: "<< argv[0] << " <filename>|<stdin> ";
	cerr << "<quality threshold> <coverage threshold>" << endl;
	cerr << endl;
	cerr << "Takes in mpileup result format:"<<endl;
	cerr << "samtools mpileup -f <ref.fa> <bamFile> | ";
	cerr << argv[0] << " - ";
	cerr << "<quality threshold> <coverage threshold>" << endl;
	cerr << endl;
	return 0;
}

string parseBases(string bases, string ref, int &deletion, int &insertion, string &insertion_bases, string &deletion_bases)
{
    int i_max = bases.length(), i = 0, j, insert_count, del_count;
    string _bases = "";
	string reverse_ref;
	string insert_digit, del_digit;
	char c;
    while (i < i_max)
	{
        c = bases.at(i);
        if (c == '.')
		{
            _bases +=  ref;
		}
        else if( c == ',')
		{
			reverse_ref = complementBases(ref);
			transform(reverse_ref.begin(), reverse_ref.end(), reverse_ref.begin(), ::tolower);
            _bases += reverse_ref;
		}
        else if (c == '+')
		{
            j = i + 1;
            insert_count = 0;
			insert_digit = "";
            while ( isdigit(bases.at(j)) )
			{
                insert_digit +=  bases.at(j);
                j ++;
			}
			stringstream ss;
			insert_count = atoi(insert_digit.c_str());
			cout << insert_count << endl;
            i = j + insert_count - 1;
            insertion_bases += bases.substr( j , i+1 );
            insertion += insert_count;
		}
        else if (c == '-')
		{
            j = i + 1;
            del_count = 0;
			del_digit = "";
            while ( isdigit(bases.at(j)) )
			{
                del_digit +=  bases.at(j);
                j ++;
			}
			stringstream ss;
			del_count = atoi(del_digit.c_str());
            i = j + del_count - 1;
        	deletion_bases += bases.substr( j , i+1 );
            deletion += del_count;
		}
		else if (c == '^')
		{
		 i ++;
		}
		else if (c!='$')
		{
			_bases += c;
		}
        i ++;
	}
    return _bases;
}

string qualityBases(string bases, string quals, int qual_threshold)
{
	stringstream ss;
	string high_qual_bases = "";
	int no_of_base = bases.length(), i = 0;
	int base_qual;
	string hq_base;
	for (i = 0; i < no_of_base; i++)
	{
		base_qual = quals[i] - 33;
		if (base_qual >= qual_threshold)
		{
			ss << bases.at(i);
			ss >> hq_base;
			high_qual_bases += hq_base;
		}
	}
	return high_qual_bases;
}

void countBase( string high_qual_bases, string positive_strand_bases, string negative_strand_bases)
{
	
}

// extract from each line different columns
// and give them to further processing
void processLine( stringList columns, int qualThreshold, int coverageThreshold, string line)
{
	if (columns[2] != "N" && columns[2] != "." && columns[2] != "_")
	{
	    string transcriptID, pos, ref, bases, quals;
		string extracted_bases, high_qual_bases;
		string positive_strand_bases = "", negative_strand_bases= "";
		int insertion = 0, deletion = 0;
    	string insertion_bases = "", deletion_bases = "";
	    int cov;
	    if (columns.size() == 6)
	    {
	        cov = atoi(columns[3].c_str());
	        if (cov > coverageThreshold)
	        {
	            transcriptID = columns[0];
	            pos = columns[1];
	            ref = columns[2];
	            bases = columns[4];
	    		quals = columns[5];
	            assert( quals.length() == cov ) ;
				extracted_bases = parseBases(bases, ref, deletion, insertion , insertion_bases, deletion_bases);
				cout << extracted_bases << '\n' << quals << endl;
				cout << extracted_bases.length() << '\n' << cov	 << endl;
				assert(extracted_bases.length() == cov);
				high_qual_bases = qualityBases(extracted_bases, quals, qualThreshold);
            }
        }
    }
}


// if lines are read from file,
// this function takes in and open the file and
// parse it line by line
void readFile(const char* filename, int qualThreshold, int coverageThreshold)
{
    ifstream myfile(filename);
    for (string line; getline(myfile, line);)
    {
        stringList columns = split(line,'\t');
        processLine(columns, qualThreshold, coverageThreshold, line);
    }
}

// if lines are read from stdin,
// this function takes in and open the file and
// parse it line by line
void readStream(int qualThreshold, int coverageThreshold)
{
    for (string line; getline(cin, line);)
    {
        stringList columns = split(line,'\t');
        processLine(columns, qualThreshold, coverageThreshold, line);
    }
}

void printHeader()
{
	cout << "chrom\tstart\tend\tref_base\tcoverage\tstrand\tA\tC\tT\tG\tinsertion\tdeletion\n";
}

// main function
int main(int argc, char *argv[])
{
    ios::sync_with_stdio(false);
    // warnings
    if (argc != 4)
    {
        usage(argv);
	return 0;
    }

    // create modified RNA index
    int qualThreshold, coverageThreshold;
    qualThreshold = atoi(argv[2]);
    coverageThreshold = atoi(argv[3]);
    cerr << "Using quality threshold:  "<< qualThreshold << '\n';
    cerr << "Using coverage threshold: "<< coverageThreshold << '\n';
    // read lines
    printHeader();
    if (strcmp(argv[1],"-") == 0)
    {
	cerr << "Reading from stdin" << endl;
        readStream(qualThreshold, coverageThreshold);
    }
    else
    {
        const char* filename = argv[1];
		cerr << "Reading from: " << filename << endl;
        readFile(filename, qualThreshold, coverageThreshold);
    }
    return 0;
}
