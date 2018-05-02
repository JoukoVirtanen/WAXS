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

using namespace std;

#include "/home2/jouko/project/HeaderFiles/StringUtils.h"
#include "/home2/jouko/project/HeaderFiles/LinkedList.h"
#include "/home2/jouko/project/HeaderFiles/TypeDef.h"

string GetWord(string line, int WordNum)
{
        int n, position, WordCount;
        string str;

	line=" "+line+" ";

        WordCount=count(line.begin(),line.end(),' ');

        if (WordNum<WordCount)
        {
                for (n=0;n<WordNum;n++)
                {

                        while (true)
                        {
                                str=line[0];

                                if (str==" ") line=line.substr(1,100);
                                else break;

                                WordCount=count(line.begin(),line.end(),' ');
                                if (WordCount==0) break;
                        }

                        if (WordCount>0)
                        {
                                position=line.find(" ");
                                str=line.substr(0,position);
                                line=line.substr(position, 100);
                        }
                        else str=line.substr(0,100);
                }
        }

        return str;
}

string GetWord2(string line, int WordNum)
{
        int n, position, WordCount;
        string str;

	line="\t"+line+"\t";

        WordCount=count(line.begin(),line.end(),'\t');

        if (WordNum<WordCount)
        {
                for (n=0;n<WordNum;n++)
                {

                        while (true)
                        {
                                str=line[0];

                                if (str=="\t") line=line.substr(1,100);
                                else break;

                                WordCount=count(line.begin(),line.end(),'\t');
                                if (WordCount==0) break;
                        }

                        if (WordCount>0)
                        {
                                position=line.find("\t");
                                str=line.substr(0,position);
                                line=line.substr(position, 100);
                        }
                        else str=line.substr(0,100);
                }
        }

        return str;
}

bool exist(char charfile[100])
{
	bool n;

	fstream file;
	file.open(charfile,ios::in);

	if (!file) n=false;
	else n=true;

	return n;
}

void FileNotValid(char CharFile[], string FileType)
{
	bool valid;
	string File;

	cout <<"Error "<<FileType<<" "<<CharFile<<" cannot be opened."<<endl;
	exit(1);
}

void ReadCalculatedData(string experiment, Real i[], Real s[], int points)
{
	bool valid;
	char charexperiment[1000];
        string line, str;
        int t;
        int position, wordcount;

	strcpy(charexperiment, experiment.c_str());
	valid=exist(charexperiment);
	if (!valid) FileNotValid(charexperiment, "calculated data file");
        //cout <<"Checked for file"<<endl;
	if (experiment=="")
	{
	        experiment="/home2/jouko/WAXS/LysozymeDiluted3.txt";
	        strcpy(charexperiment, experiment.c_str());
	}
        //cout <<"After strcpy"<<endl;
	fstream ExperimentFile;
	ExperimentFile.open(charexperiment,ios::in);
        //cout <<"Opened file"<<endl;
        t=0;
	while (t<points)
	{
                //cout <<"t= "<<t<<" points= "<<points<<endl;
		getline(ExperimentFile, line);

		//line="  "+line+"   ";
                //cout <<"GetWord2(line,1)[0]= "<<GetWord2(line,1)[0]<<endl;
                //cout <<"GetWord2(line,1)= "<<GetWord2(line,1)<<endl;
                //cout <<"line= "<<line<<endl;
                if (isdigit(GetWord2(line,1)[0]))
                {
		        s[t]=StrToFloat(GetWord2(line,1));
		        i[t]=StrToFloat(GetWord2(line,2));
		        t++;
                }
                //cout <<line<<endl;
		//cout <<GetWord(line,1)<<"       g              "<<GetWord(line,2)<<endl;
		//cout <<s[t]<<"  "<<i[t]<<endl;
                if (ExperimentFile.eof()) break;
	}

}

void ReadExperimentalData(string experiment, Real expi[], Real error[], int points)
{
	bool valid;
	char charexperiment[100];
        string line, str;
        int t;
        int position, wordcount;

	strcpy(charexperiment, experiment.c_str());
	valid=exist(charexperiment);
	if (!valid) FileNotValid(charexperiment, "experimental data file");

	if (experiment=="")
	{
	        experiment="/home2/jouko/WAXS/LysozymeDiluted3.txt";
	        strcpy(charexperiment, experiment.c_str());
	}

	fstream ExperimentFile;
	ExperimentFile.open(charexperiment,ios::in);

	for (t=0;t<points;t++)
	{
		getline(ExperimentFile, line);
		wordcount=count(line.begin(),line.end(),' ');

		if (wordcount>0)
		{
			position=line.find(" ");
			str=line.substr(0,position);
			std::istringstream ee (str);
			ee >> expi[t];

			str=line.substr(position,12);
			std::istringstream ff (str);
			ff >> error[t];
		}
		else
		{
			std::istringstream ee (line);
			ee >> expi[t];
			error[t]=(expi[t]+expi[0]*0.01)*0.02;
		}
	}
	for (t=0;t<points;t++)
	{
	//	cout <<"expi["<<t<<"]= "<<expi[t]<<endl;
	}
        //cout <<"error[0]="<<error[0]<<" error[1]="<<error[1]<<" error[2]="<<error[2]<<endl;
}

Real RFactor(Real i[], Real expi[], int beginfit, int endfit, int points, Real error[])
{
	int n;
	Real rfactor, scale, SumDiff, sum;
	SumDiff=0;
	sum=0;
	scale=expi[0]/i[0];
	//cout <<"scale= "<<scale<<endl;
	for (n=0;n<points;n++) i[n]*=scale;
	for (n=beginfit;n<endfit;n++)
	{
		SumDiff+=abs(expi[n]-i[n]);
		sum+=(expi[n]+i[n]);
		//cout <<"SumDiff= "<<SumDiff<<" sum= "<<sum<<" expi["<<n<<"]= "<<expi[n]<<" i["<<n<<"]= "<<i[n]<<endl;
	}
	rfactor=SumDiff/sum;
	return rfactor;
}

Real MatchI0(Real i[], Real expi[], int beginfit, int endfit, int points, Real error[])
{
        int t;
        Real chisqr;
        Real isqr;
        Real crossterm;
        Real scale;
        Real sumsqr;
        Real * inverror2;

        inverror2=new Real[endfit+1];

        for (t=0;t<endfit+1;t++)
        {
                inverror2[t]=0;
        }

        //cout <<"error[0]="<<error[0]<<" error[1]="<<error[1]<<endl;

        for (t=0;t<endfit+1;t++)
        {
                //cout <<"In normalize.  expi["<<t<<"]="<<expi[t]<<endl;
                if (error[t]!=0)
                {
                        inverror2[t]=1.0/(error[t]*error[t]);
                }
                else
                {
                        if (expi[t]!=0)
                        {
                                inverror2[t]=1.0/(expi[t]*expi[t]*0.02*0.02);
                        }
                        else
                        {
                                inverror2[t]=0.01;
                                //cout <<"Error.  Experimental intensity is 0"<<endl;
                        }
                }
        }

        isqr=0;
        crossterm=0;
        sumsqr=0;

	scale=expi[0]/i[0];
	//cout <<"scale= "<<scale<<endl;
        for (t=0;t<points;t++)
        {
                i[t]=i[t]*scale;
        }

        for (t=beginfit;t<endfit+1;t++)
        {
                sumsqr=sumsqr+(i[t]-expi[t])*(i[t]-expi[t])*inverror2[t];
        }
	
        chisqr=sumsqr/float(endfit-beginfit+1);

        delete [] inverror2;

        return chisqr;
}

Real NoScaling(Real i[], Real expi[], int beginfit, int endfit, int points, Real error[])
{
        int t;
        Real chisqr;
        Real isqr;
        Real crossterm;
        Real scale;
        Real sumsqr;
        Real * inverror2;

        inverror2=new Real[endfit+1];

        for (t=0;t<endfit+1;t++)
        {
                inverror2[t]=0;
        }

        //cout <<"error[0]="<<error[0]<<" error[1]="<<error[1]<<endl;

        for (t=0;t<endfit+1;t++)
        {
                //cout <<"In normalize.  expi["<<t<<"]="<<expi[t]<<endl;
                if (error[t]!=0)
                {
                        inverror2[t]=1.0/(error[t]*error[t]);
                }
                else
                {
                        if (expi[t]!=0)
                        {
                                inverror2[t]=1.0/(expi[t]*expi[t]*0.02*0.02);
                        }
                        else
                        {
                                inverror2[t]=0.01;
                                //cout <<"Error.  Experimental intensity is 0"<<endl;
                        }
                }
        }

        isqr=0;
        crossterm=0;
        sumsqr=0;

        for (t=beginfit;t<endfit+1;t++)
        {
                sumsqr=sumsqr+(i[t]-expi[t])*(i[t]-expi[t])*inverror2[t];
        }
	
        chisqr=sumsqr/float(endfit-beginfit+1);

        delete [] inverror2;

        return chisqr;
}

Real log(Real i[], Real expi[], int beginfit, int endfit, int points, Real error[])
{
        int t;
        Real chisqr;
        Real isqr;
        Real crossterm;
        Real scale;
        Real sumsqr;

        isqr=0;
        crossterm=0;
        sumsqr=0;

        for (t=beginfit;t<endfit+1;t++)
        {
                sumsqr=sumsqr+(log(i[t]/expi[t]))*(log(i[t]/expi[t]));
        }
	
        chisqr=sumsqr/Real(endfit-beginfit+1);

        return chisqr;
}

Real CalcLogNormalization(Real expi[], Real i[], int beginfit, int endfit)
{
        Real sum=0;

        for (int t=beginfit;t<endfit;t++)
        {
                //cout <<"sum= "<<sum<<" expi["<<t<<"]= "<<expi[t]<<" i= "<<i[t]<<endl;
                sum+=log(expi[t]/i[t]);
        }
        sum/=(endfit+1-beginfit);
        return exp(sum);
}

void ScaleIntensity(Real scale, Real i[], int points)
{
        for (int t=0;t<points;t++)
        {
                //cout <<"i["<<t<<"]= "<<i[t]<<endl;
                i[t]*=scale;
                //cout <<"i["<<t<<"]= "<<i[t]<<endl;
        }
}

Real log_normalize(Real i[], Real expi[], int beginfit, int endfit, int points, Real error[])
{
        int t, numpoints=0;
        Real chisqr;
        Real isqr;
        Real crossterm;
        Real scale;
        Real sumsqr;
        //cout <<"In log_normalize"<<endl;
        isqr=0;
        crossterm=0;
        sumsqr=0;
        scale=CalcLogNormalization(expi, i, beginfit, endfit);
        //cout <<"scale= "<<scale<<endl;
        ScaleIntensity(scale, i, points);
        for (t=beginfit;t<endfit;t++)
        {
                //cout <<"expi["<<t<<"]= "<<expi[t]<<" i= "<<i[t]<<" sqr= "<<(log(i[t]/expi[t]))*(log(i[t]/expi[t]))<<endl;
                sumsqr=sumsqr+(log(i[t]/expi[t]))*(log(i[t]/expi[t]));
                numpoints++;
        }
	//cout <<"numpoints= "<<numpoints<<endl;
        //cout <<"endfit-beginfit= "<<endfit-beginfit<<endl;
        chisqr=sumsqr/Real(endfit-beginfit);

        return chisqr;
}

Real normalize(Real i[], Real expi[], int beginfit, int endfit, int points, Real error[])
{
        int t;
        Real chisqr;
        Real isqr;
        Real crossterm;
        Real scale;
        Real sumsqr;
        Real * inverror2;

        inverror2=new Real[endfit+1];

        for (t=0;t<endfit+1;t++)
        {
                inverror2[t]=0;
        }

        //cout <<"error[0]="<<error[0]<<" error[1]="<<error[1]<<endl;

        for (t=0;t<endfit+1;t++)
        {
                //cout <<"In normalize.  expi["<<t<<"]="<<expi[t]<<endl;
                if (error[t]!=0)
                {
                        inverror2[t]=1.0/(error[t]*error[t]);
                }
                else
                {
                        if (expi[t]!=0)
                        {
                                inverror2[t]=1.0/(expi[t]*expi[t]*0.02*0.02);
                        }
                        else
                        {
                                inverror2[t]=0.01;
                                //cout <<"Error.  Experimental intensity is 0"<<endl;
                        }
                }
        }

        isqr=0;
        crossterm=0;
        sumsqr=0;

        for (t=beginfit;t<endfit+1;t++)
        {
                //cout <<"In normalize.  i["<<t<<"]="<<i[t]<<endl;
                isqr+=i[t]*i[t]*inverror2[t];
                crossterm+=i[t]*expi[t]*inverror2[t];
		//cout <<"i["<<t<<"]= "<<i[t]<<" expi["<<t<<"]= "<<expi[t]<<endl;
        }
	//cout <<"crossterm= "<<crossterm<<" isqr= "<<isqr<<endl;
        if (isqr!=0)
        {
                scale=crossterm/isqr;
        }
        else
        {
                scale=1.0;
                cout <<"Error scalling\n";
        }
	//cout <<"scale= "<<scale<<endl;
        for (t=0;t<points;t++)
        {
                i[t]=i[t]*scale;
        }

        for (t=beginfit;t<endfit+1;t++)
        {
                sumsqr=sumsqr+(i[t]-expi[t])*(i[t]-expi[t])*inverror2[t];
        }
	
        chisqr=sumsqr/float(endfit-beginfit+1);

        delete [] inverror2;

        return chisqr;
}

int main (int argc, char *argv[])
{
	string CalculatedDataFile, ExperimentalDataFile, option;
	int n;
	int beginfit, endfit, points;
	Real chisqr;
	Real *i, *error, *expi, *s;
	//cout <<"Starting"<<endl;
        //cout <<"starting"<<endl;
        if (argc<6)
        {
                cout <<"Usage: "<<argv[0]<<" CalculatedDataFile ExperimentalDataFile beginfit endfit points option"<<endl;
                exit(EXIT_FAILURE);
        }
	CalculatedDataFile=argv[1];
	ExperimentalDataFile=argv[2];
	beginfit=StrToInt(argv[3]);
	endfit=StrToInt(argv[4]);
	points=StrToInt(argv[5]);
	option=(argv[6]);
	//cout <<"Got arguments"<<endl;
        i=new Real[points];
	error=new Real[points];
	expi=new Real[points];
	s=new Real[points];
        //cout <<"Allocated memory"<<endl;
	for (n=0;n<points;n++)
	{
		i[n]=0;
		error[n]=0;
		expi[n]=0;
		s[n]=0;
	}
        //cout <<"Zeroed arrays"<<endl;
	//cout <<"About to read calculated data"<<endl;
	ReadCalculatedData(CalculatedDataFile, i, s, points);
	//cout <<"ReadCalculated Data"<<endl;
	ReadCalculatedData(ExperimentalDataFile, expi, s, points);
	//cout <<"Read experimental data"<<endl;
        //ReadExperimentalData(ExperimentalDataFile, expi, error, points);
	//cout <<"ReadExperimental Data"<<endl;
	if (option=="normalize") chisqr=normalize(i, expi, beginfit, endfit, points, error);
	else if (option=="NoScaling") chisqr=NoScaling(i, expi, beginfit, endfit, points, error);
	else if (option=="RFactor") chisqr=RFactor(i, expi, beginfit, endfit, points, error);
	else if (option=="MatchI0") chisqr=MatchI0(i, expi, beginfit, endfit, points, error);
	else if (option=="log") chisqr=log(i, expi, beginfit, endfit, points, error);
	else if (option=="log_normalize") chisqr=log_normalize(i, expi, beginfit, endfit, points, error);
        else
        {
                cout <<"ERROR: Unrecognized option "<<option<<endl;
                exit(EXIT_FAILURE);
        }
        //chisqr=NoScaling(i, expi, beginfit, endfit, points, error);
	cout <<"chisqr= "<<chisqr<<endl;
	
        delete [] i;
        delete [] error;
        delete [] expi;
        delete [] s;
        
        return 0;
}

	



