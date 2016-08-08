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

string parseBases(string bases, string ref, int &deletion, int &insertion)
{
    int i_max = bases.size(), i = 0;
    string _bases = "";
	char c;
    while (i < i_max)
	{
        c = bases[i];
        if (c == '.')
		{
            _bases = _bases +  ref;
		}
        else if( c == ',')
		{
            _bases = _bases + reverseStrandcomplement(ref);
		}
/*
        elif c == '^':
            i += 1
        elif c == '+':
            j = i + 1
            indel_count = 0
            while bases[j].isdigit():
                indel_count += int(bases[j])
                j += 1
            i = j + indel_count - 1
            insertion += indel_count
        elif c == '-':
            j = i + 1
            indel_count = 0
            while bases[j].isdigit():
                indel_count += int(bases[j])
                j += 1
            i = j + indel_count - 1
            deletion += indel_count
        elif c != '$':
            _bases += c
*/
        i ++;
	}
	cout << _bases << endl;
//    return _bases, insertion, deletion
}

// extract from each line different columns
// and give them to further processing
void processLine( stringList columns, int qualThreshold, int coverageThreshold, string line)
{
	if (columns[2] != "N" && columns[2] != "." && columns[2] != "_")
	{
	    string transcriptID, pos, ref, bases, quals;
		int insertion = 0, deletion = 0;
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
				parseBases(bases, ref, deletion, insertion );
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
