#include <string>
#include "stringVector.h"

using namespace std;

int countDigits(int number)
{
	if (number < 10)
	{
		return 1;
	}
	int count = 0;
	while (number > 0)
	{
		number /= 10;
		count++;
	}
	return count;
}

string complementBases(string base)
{
	string basepair;
	if (base.compare("A")==0)
	{
		basepair = "T";
	}
	else if (base.compare("T")==0)
	{
		basepair = "A";
	}
	else if (base.compare("C") == 0)
	{
		basepair = "G";
	}
	else if (base.compare("G")==0)
	{
		basepair = "C";
	}
	else if (base.compare("N")==0)
	{
		basepair = "N";
	}
	return basepair;
}
