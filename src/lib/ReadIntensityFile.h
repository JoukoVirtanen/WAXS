#ifndef _ReadIntensityFile_included_
#define _ReadIntensityFile_included_

# include <iostream>
# include <fstream>
# include <string>
# include <cstring>
# include <cmath>
# include <sstream>
# include <time.h>
# include <vector>
# include <iomanip>
# include <complex>
# include <numeric>
# include <ctype.h>

# include "IOUtils.h"
# include "LinkedList.h"
# include "StringUtils.h"
# include "TypeDef.h"

void ReadIntensityFile(string IntensityFile, vector<Real> &s, vector<Real> &expi, vector<Real> &error)
{
	bool valid;
	char charexperiment[1000], test[1000];
	string line, str;
        vector<string> word;
	int wordcount, wordcount2;
        s.clear();
        expi.clear();
        error.clear();
	strcpy(charexperiment, IntensityFile.c_str());
	valid=exist(charexperiment);
	if (!valid)
	{
		cout <<"Unable to open file "<<charexperiment<<endl;
		exit(0);
	}


	fstream ExperimentFile;
	ExperimentFile.open(charexperiment,ios::in);
	getline(ExperimentFile, line);
	while (true)
	{
		if (ExperimentFile.eof()) break;
		wordcount=count(line.begin(),line.end(),' ');
                wordcount2=count(line.begin(), line.end(), '\t');
                if (wordcount2>wordcount) wordcount=wordcount2;
		strcpy(test, GetWord2(line, 1).c_str());
		if (isdigit(test[0]))
		{
			if (wordcount>1)
			{
				s.push_back(StrToFloat(GetWord2(line, 1)));
				expi.push_back(StrToFloat(GetWord2(line, 2)));
				error.push_back(StrToFloat(GetWord2(line, 3)));
				//cout <<"one "<<"s= "<<s[s.size()-1]<<" "<<expi[expi.size()-1]<<" "<<error[error.size()-1]<<endl;
			}
			else
			{
				s.push_back(StrToFloat(GetWord2(line, 1)));
				expi.push_back(StrToFloat(GetWord2(line, 2)));
				error.push_back((expi[expi.size()-1]+expi[0]*0.01)*0.02);
				//cout <<"two "<<"s= "<<s[s.size()-1]<<" "<<expi[expi.size()-1]<<" "<<error[error.size()-1]<<endl;
				//cout <<"s.size()= "<<s.size()<<endl;
			}
		}
		getline(ExperimentFile, line);
	}
	for (int n=0;n<int(s.size());n++)
	{
		//cout <<"s["<<n<<"]= "<<s[n]<<" expi= "<<expi[n]<<" error= "<<error[n]<<endl;
	}
	cout <<"error[0]="<<error[0]<<" error[1]="<<error[1]<<" error[2]="<<error[2]<<endl;
}
#endif

