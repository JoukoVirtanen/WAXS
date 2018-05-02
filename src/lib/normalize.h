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

#include "StringUtils.h"
#include "LinkedList.h"
#include "IOUtils.h"
#include "TypeDef.h"

Real RFactor(vector<Real> &i, vector<Real> &expi, int beginfit, int endfit, int points, vector<Real> &error)
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

Real MatchI0(vector<Real> &i, vector<Real> &expi, int beginfit, int endfit, int points, vector<Real> &error)
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

Real NoScaling(vector<Real> &i, vector<Real> &expi, int beginfit, int endfit, int points, vector<Real> &error)
{
        int t;
        Real chisqr;
        Real isqr;
        Real crossterm;
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

Real log(vector<Real> &i, vector<Real> &expi, int beginfit, int endfit, int points, vector<Real> &error)
{
        int t;
        Real chisqr;
        Real isqr;
        Real crossterm;
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

Real CalcLogNormalization(vector<Real> &expi, vector<Real> &i, int beginfit, int endfit)
{
        int nexp=expi.size();
        int npoint=i.size();
        Real sum=0;

        if (nexp<endfit) endfit=nexp;
        if (npoint<endfit) endfit=npoint;

        for (int t=beginfit;t<endfit;t++)
        {
                cout <<"sum= "<<sum<<" expi["<<t<<"]= "<<expi[t]<<" i= "<<i[t]<<endl;
                if (expi[t]>0 && i[t]>0) sum+=log(expi[t]/i[t]);
        }
        sum/=(endfit+1-beginfit);
        return exp(sum);
}

void ScaleIntensity(Real scale, vector<Real> &i)
{
        int points=i.size();
        for (int t=0;t<points;t++)
        {
                //cout <<"i["<<t<<"]= "<<i[t]<<endl;
                i[t]*=scale;
                //cout <<"i["<<t<<"]= "<<i[t]<<endl;
        }
}

Real log_normalize(vector<Real> &i, vector<Real> &expi, int beginfit, int endfit, vector<Real> &error)
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
        cout <<"scale= "<<scale<<endl;
        ScaleIntensity(scale, i);
        for (t=beginfit;t<endfit;t++)
        {
                if (expi[t]>0)
                {
                        cout <<"expi["<<t<<"]= "<<expi[t]<<" i= "<<i[t]<<" sqr= "<<(log(i[t]/expi[t]))*(log(i[t]/expi[t]))<<endl;
                        sumsqr=sumsqr+(log(i[t]/expi[t]))*(log(i[t]/expi[t]));
                        numpoints++;
                }
        }
        //cout <<"numpoints= "<<numpoints<<endl;
        //cout <<"endfit-beginfit= "<<endfit-beginfit<<endl;
        chisqr=sumsqr/Real(endfit-beginfit);

        return chisqr;
}

Real normalize(vector<Real> &i, vector<Real> &expi, int beginfit, int endfit, int points, vector<Real> &error)
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

Real normalize(vector<Real> &i, vector<Real> &expi, int beginfit, int endfit)
{
        int t;
        Real chisqr;
        Real isqr;
        Real crossterm;
        Real scale;
        Real sumsqr;

        //cout <<"error[0]="<<error[0]<<" error[1]="<<error[1]<<endl;


        isqr=0;
        crossterm=0;
        sumsqr=0;

        for (t=beginfit;t<=endfit;t++)
        {
                //cout <<"In normalize.  i["<<t<<"]="<<i[t]<<endl;
                isqr+=i[t]*i[t];
                crossterm+=i[t]*expi[t];
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
        for (t=0;t<int(i.size());t++)
        {
                //i[t]=i[t]*scale;
        }

        for (t=beginfit;t<endfit+1;t++)
        {
                sumsqr=sumsqr+(i[t]-expi[t])*(i[t]-expi[t]);
        }

        chisqr=sumsqr/Real(endfit-beginfit+1);


        return scale;
}

Real normalize(vector<Real> &i, vector<Real> &expi, int beginfit, int endfit, int points, vector<Real> &error, string option)
{
        Real chisqr;
        if (option=="normalize") chisqr=normalize(i, expi, beginfit, endfit, points, error);
        else if (option=="NoScaling") chisqr=NoScaling(i, expi, beginfit, endfit, points, error);
        else if (option=="RFactor") chisqr=RFactor(i, expi, beginfit, endfit, points, error);
        else if (option=="MatchI0") chisqr=MatchI0(i, expi, beginfit, endfit, points, error);
        else if (option=="log") chisqr=log(i, expi, beginfit, endfit, points, error);
        else if (option=="log_normalize") chisqr=log_normalize(i, expi, beginfit, endfit, error);
        else
        {
                cout <<"ERROR: Unrecognized option "<<option<<endl;
                exit(EXIT_FAILURE);
        }
        //chisqr=NoScaling(i, expi, beginfit, endfit, points, error);
        cout <<"chisqr= "<<chisqr<<endl;

        return chisqr;
}





