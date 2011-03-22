/*
 * The contents of this file were provided by Quang Le.
*/

#ifndef __CVCFT__
#define  __CVCFT__

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <iostream>

using namespace std;



#define M20 20000000


class CVCFT
{
private:

	string keyid;
	int   ikeyid;
	int split(char* bf, char det,  vector<char*> &vc);

	vector<string> names;
	int ns;
	vector<char*> data;

	string line;

	istream *ins;

	char *bf;

public:

	CVCFT ( std::istream &i_ins, const char* key);
	~CVCFT();

	int ReadNext(void);

	vector<string>* GetNames(void) {
		return &names;
	}

	const char* GetRSID(void) {
		return data[2];
	}
	char GetRef(void) {
		return data[3][0];
	}
	char GetAlt(void) {
		return data[4][0];
	}
	void GetChrPos(string &chr, int &pos) {
		chr = data[0];
		pos = atoi(data[1]) ;
	}
	void GetProb(vector<double>& prob);



};
#endif
