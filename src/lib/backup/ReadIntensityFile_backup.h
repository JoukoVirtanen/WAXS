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
	char charexperiment[100], test[1000];
	string line, str;
	int position, wordcount;

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
		strcpy(test, GetWord2(line, 1).c_str());
		if (isdigit(test[0]))
		{
			if (wordcount>0)
			{
				s.push_back(StrToFloat(GetWord2(line, 1)));
				expi.push_back(StrToFloat(GetWord2(line, 2)));
				error.push_back(StrToFloat(GetWord2(line, 3)));
				//cout <<"one "<<"s["<<t<<"]= "<<s[t]<<" "<<expi[t]<<" "<<error[t]<<endl;
			}
			else
			{
				s.push_back(StrToFloat(GetWord2(line, 1)));
				expi.push_back(StrToFloat(GetWord2(line, 2)));
				error.push_back((expi[expi.size()-1]+expi[0]*0.01)*0.02);
				//cout <<"two "<<"s["<<t<<"]= "<<s[t]<<" "<<expi[t]<<" "<<error[t]<<endl;
			}
		}
		getline(ExperimentFile, line);
	}

	cout <<"error[0]="<<error[0]<<" error[1]="<<error[1]<<" error[2]="<<error[2]<<endl;
}
#endif

