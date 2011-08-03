/*
 * The contents of this file were provided by Quang Le.
*/

#include <limits>
#include "CVCFT.h"

CVCFT :: CVCFT (std::istream &i_ins, const char* key)
{
	keyid = key;
	ins = &i_ins;

	while (std::getline((*ins), line))
		if (strstr(line.c_str(), "#CHROM") == line.c_str()) break; //find the #CHROM


	bf = (char*)calloc(line.size() + 10, sizeof(char));
	strcpy(bf, line.c_str());

	vector<char*> vc;
	split(bf, '\t',  vc);
	for (int i = 9; i != (int) vc.size(); ++i) names.push_back(vc[i]);
	ns = (int) names.size();

	//fprintf(stderr, "ns = %d\n", ns);
}


CVCFT :: ~CVCFT()
{
	free(bf);
}



int CVCFT :: split(char*bf, char det, vector<char*>  &vc)
{
	vc.clear();
	vc.push_back(bf);
	int len = strlen(bf);
	for (int i = 1; i < len; ++i)
		if (bf[i-1] == det) {
			bf[i-1] = '\0';
			vc.push_back(&bf[i]);
		}
	return (int) vc.size();
}



void CVCFT :: GetProb(vector<double> & prob)
{
	prob.clear();
	vector<char*> ts, ts1;

	for (int i = 0; i != ns; ++i)
	{
		//      printf("%s\n", data[i+9]);
		if( ikeyid == std::numeric_limits< int >::max() ) {
			// Key not in format.  Don't set anything.
		}
		else if (split(data[i+9], ':', ts) <= ikeyid || split(ts[ikeyid], ',', ts1) < 3)
		{
			prob.push_back(0);
			prob.push_back(0);
			prob.push_back(0);
			//	  fprintf(stderr, "no data for sample: %s at site %s\n", names[i].c_str(), data[2]);
		}
		else {
			for (int k = 0; k != 3; ++k) prob.push_back( atof (ts1[k]));
		}
	}
}


int CVCFT :: ReadNext(void) {

	vector<char*> ts;
	string line;


	if (getline((*ins), line))
	{
		free(bf);
		bf = (char*)calloc(line.size() + 10, sizeof(char));
		strcpy(bf, line.c_str());

		ikeyid = std::numeric_limits< int >::max() ;

		split(bf, '\t', data);
		//      printf("key = %s, data = %s\n", keyid.c_str(), data[8]);
		split(data[8], ':', ts);



		for (int i = 0; i != (int) ts.size(); ++i)
			if (strcmp(ts[i], keyid.c_str()) == 0) ikeyid = i;
		//if (ikeyid == -1)  return 0;
		if (strlen(data[3]) > 1 || strlen(data[4]) > 1)	return 0;  //not biallele
		else       return 1;
	}
	else return 0;
}
