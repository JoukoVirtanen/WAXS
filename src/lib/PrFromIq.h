#ifndef _PrFromScattering_included_
#define _PrFromScattering_included_

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

# include "IOUtils.h"
# include "MathUtils.h"
# include "minimize.h"

struct ArgStruct
{
        int BeginFit, EndFit;
        Real Dmax, PenaltyCoefficient;
        Real **imn;
        vector<Real> i, expi, error;
};

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

        for (t=0;t<endfit+1;t++) inverror2[t]=0;

        //cout <<"error[0]="<<error[0]<<" error[1]="<<error[1]<<endl;

        for (t=0;t<endfit+1;t++)
        {
                //cout <<"In normalize.  expi["<<t<<"]="<<expi[t]<<endl;
                if (error[t]!=0) inverror2[t]=1.0/(error[t]*error[t]);
                else
                {
                        if (expi[t]!=0) inverror2[t]=1.0/(expi[t]*expi[t]*0.02*0.02);
                        else
                        {
                                inverror2[t]=0.01;
                                cout <<"Error.  Experimental intensity is 0"<<endl;
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
        }

        if (isqr!=0) scale=crossterm/isqr;
        else
        {
                scale=1.0;
                cout <<"Error scalling"<<endl;
        }

        for (t=0;t<points;t++) i[t]=i[t]*scale;

        for (t=beginfit;t<endfit+1;t++) sumsqr+=(i[t]-expi[t])*(i[t]-expi[t])*inverror2[t];
	
        chisqr=sumsqr/float(endfit-beginfit+1);

        delete [] inverror2;

        return chisqr;
}
/************************************************************************
Real fitpr3(Real s[], Real expi[], int beginfit, int endfit)
{

        int m, n, p, t, u;
        int points;
	int start, pick, NumIncrease;
	char CharPrFile[100], CharPrFile2[100], CharIntensityFile[100];
	string IntensityFile, PrFile2;
        Real a, b, c;
        Real bestchisqr;
        Real chisqr, chisqr1, chisqr2, old_chisqr;
        Real TotalPenalty;
	Real dChisqr, move, MoveSize, Temp, Temp2, NumAccepted, StnDev, delta;
	Real NumRejected, TotalMoves;
	Real TotalChisqr, TotalChisqrSqr, Sum;
        Real RandNum;
        Real dist;
        Real dist2;
        Real inc;
	Real TotalPr;
        Real PR[45];
        Real PR2[5000], FitPr2[5000];
        Real i[10000];
	Real InitTempArray[100];
        Real derivative2;
        Real parameter[100], parameter2[100], parameter3[100], fitparameter[100];
        Real step[100];

        bestchisqr=1000000000;
	Temp=20.0;
	delta=1.0;
	NumAccepted=0;
	NumIncrease=0;

	ExtrapolateIntensity(s, expi, endfit);

        for (t=0;t<45;t++)
        {
                PR[t]=0;
        }

        for (t=0;t<5000;t++)
        {
                PR2[t]=0;
		FitPr2[t]=0;
        }

        for (t=0;t<1000;t++)
        {
                i[t]=0;
        }

        for (t=0;t<100;t++)
        {
                fitparameter[t]=0;
                parameter[t]=0;
                parameter2[t]=0;
                parameter3[t]=0;
                step[t]=0;
		InitTempArray[t]=0;
        }

        dist=0;
        inc=1.3;
	MoveSize=0.001;

        for (t=0;t<40;t++)
        {
                PR[t]=dist*dist*(2.0-3.0*dist+dist*dist*dist);
                dist=dist+1/40.0;
        }

        step[0]=0.1;

        points=endfit+1;

        while (step[0]*step[0]>0.0000001)                   //This optimizes inc.  It optimizes dmax.
        {
                dist=0;

                a=((PR[0]-PR[1])*(inc-2*inc)-(PR[1]-PR[2])*(0-inc))/((0-inc*inc)*(inc-2*inc)-(inc*inc-2*inc*2*inc)*(0-inc));
                b=(PR[0]-PR[1]-a*(0-inc*inc))/(0-inc);
                c=0;

                dist2=0;
                for (m=0;m<50;m++)
                {
                        PR2[m]=a*dist2*dist2+b*dist2;
                        dist2=dist2+inc/101.0;
                }

                dist=inc;

		start=50;
		for (t=1;t<44;t++)
		{
			cout <<"t= "<<t<<endl;
			if (PR[t]==0 && PR[t-1]!=0)
			{
				//cout <<"PR["<<t<<"]= "<<PR[t]<<" PR["<<t-1<<"]= "<<PR[t-1]<<endl;
			        a=((PR[t+1]-PR[t])*(inc)-(PR[t]-PR[t-1])*(inc))/(((dist+inc)*(dist+inc)-dist*dist)*(inc)-(dist*dist-(dist-inc)*(dist-inc))*(inc));
			        b=(PR[t+1]-PR[t]-a*((dist+inc)*(dist+inc)-dist*dist))/(inc);
			        c=PR[t]-b*dist-a*dist*dist;
			        dist2=dist-inc*0.5+inc*0.5/101.0;
				for (m=start;m<start+50;m++)
				{
				        PR2[m]=a*dist2*dist2+b*dist2+c;
				        dist2=dist2+inc/101.0;
					//cout  <<"PR2["<<m<<"]= "<<PR2[m]<<" dist2= "<<dist2<<endl;
				}

				for (m=start+50;m<start+101;m++)
				{
				        PR2[m]=0;
				}
			}

			if (PR[t]==0 && PR[t-1]==0)
			{
				for (m=start;m<start+101;m++)
				{
				        PR2[m]=0;
				}
			}

			if (PR[t]!=0)
			{
				cout <<"PR["<<t<<"]= "<<PR[t]<<endl;
			        a=((PR[t+1]-PR[t])*(inc)-(PR[t]-PR[t-1])*(inc))/(((dist+inc)*(dist+inc)-dist*dist)*(inc)-(dist*dist-(dist-inc)*(dist-inc))*(inc));
			        b=(PR[t+1]-PR[t]-a*((dist+inc)*(dist+inc)-dist*dist))/(inc);
			        c=PR[t]-b*dist-a*dist*dist;
			        dist2=dist-inc*0.5+inc*0.5/101.0;
				for (m=start;m<start+101;m++)
				{
				        PR2[m]=a*dist2*dist2+b*dist2+c;
				        dist2=dist2+inc/101.0;
				}
			}
		        start=start+101;
		        dist=dist+inc;
		}

                for (m=0;m<endfit+1;m++)
                {
                        dist=inc/101.0;
                        i[m]=0;
                        for (n=1;n<5000;n++)
                        {
                                i[m]+=PR2[n]*sin(s[m]*dist)/(s[m]*dist);
                                dist=dist+inc/101.0;
                        }
                }

                chisqr2=chisqr1;
                chisqr1=normalize(i, expi, 0, endfit, points, error);

                if (chisqr2>chisqr1)
                {
                        step[0]=1.2*step[0];
                }

                if (chisqr2<chisqr1)
                {
                        step[0]=-0.5*step[0];
                }

                inc=inc+step[0];		
        }

        cout <<inc<<endl;

        for (t=0;t<endfit+1;t++)
        {
                ofstream intensity("/home2/jouko/WAXS/Intensity_fitpr2.txt", ios::app);
                intensity << s[t] <<" "<<i[t]<<endl;
        }
        ofstream intensity("/home2/jouko/WAXS/Intensity_fitpr2.txt", ios::app);
        intensity << endl;

	dist=0.0;
        for (t=0;t<4500;t++)
        {
                ofstream pofr2("/home2/jouko/WAXS/pr2lysozyme_fitpr2.txt", ios::app);
                pofr2 << dist <<" "<<PR2[t]<<endl;
		dist+=inc/101.0;
        }
        ofstream pofr2("/home2/jouko/WAXS/pr2lysozyme_fitpr2.txt", ios::app);
        pofr2 << endl;

	inc=1.0;
        for (t=0;t<40;t++)
        {
                parameter[t]=PR[t];
        }

        points=endfit+1;
	old_chisqr=chisqr1;

        for (u=0;u<25000;u++)
        {
		pick=int(Real(rand())*45.0/Real(RAND_MAX));
		move=MoveSize*(Real(rand())/Real(RAND_MAX)-0.5);
		if (pick==0) pick=1;
		if (pick==44) pick=43;

		PR[pick]+=move;

		cout <<"chisqr= "<<chisqr<<endl;
		p=0;
		TotalPr=0.0;
		for (n=0;n<45;n++)
		{
			TotalPr+=PR[n];
		}

		for (n=0;n<45;n++)
		{
			PR[n]=PR[n]/TotalPr;
		}

		TotalPenalty=0;
		for (n=1;n<44;n++)
		{
			derivative2=PR[n+1]-2.0*PR[n]+PR[n-1];
			TotalPenalty+=derivative2*derivative2;
		}

		TotalPenalty+=PR[44]*PR[44];

		cout <<"PR[42]= "<<PR[42]<<" PR[43]= "<<PR[43]<<" PR[44]= "<<PR[44]<<" derivative2= "<<PR[44]-2.0*PR[43]+PR[42]<<endl;

		a=((PR[0]-PR[1])*(inc-2*inc)-(PR[1]-PR[2])*(0-inc))/((0-inc*inc)*(inc-2*inc)-(inc*inc-2*inc*2*inc)*(0-inc));
		b=(PR[0]-PR[1]-a*(0-inc*inc))/(0-inc);
		c=0;

                dist2=0;
                for (m=0;m<50;m++)
                {
                        PR2[m]=a*dist2*dist2+b*dist2;
                        dist2=dist2+inc/101.0;
                }

                dist=inc;

		start=50;
		for (t=1;t<44;t++)
		{
			cout <<"t= "<<t<<endl;
			if (PR[t]==0 && PR[t-1]!=0)
			{
				//cout <<"PR["<<t<<"]= "<<PR[t]<<" PR["<<t-1<<"]= "<<PR[t-1]<<endl;
			        a=((PR[t+1]-PR[t])*(inc)-(PR[t]-PR[t-1])*(inc))/(((dist+inc)*(dist+inc)-dist*dist)*(inc)-(dist*dist-(dist-inc)*(dist-inc))*(inc));
			        b=(PR[t+1]-PR[t]-a*((dist+inc)*(dist+inc)-dist*dist))/(inc);
			        c=PR[t]-b*dist-a*dist*dist;
			        dist2=dist-inc*0.5+inc*0.5/101.0;
				for (m=start;m<start+50;m++)
				{
				        PR2[m]=a*dist2*dist2+b*dist2+c;
				        dist2=dist2+inc/101.0;
					//cout  <<"PR2["<<m<<"]= "<<PR2[m]<<" dist2= "<<dist2<<endl;
				}

				for (m=start+50;m<start+101;m++)
				{
				        PR2[m]=0;
				}
			}

			if (PR[t]==0 && PR[t-1]==0)
			{
				for (m=start;m<start+101;m++)
				{
				        PR2[m]=0;
				}
			}

			if (PR[t]!=0)
			{
				cout <<"PR["<<t<<"]= "<<PR[t]<<endl;
			        a=((PR[t+1]-PR[t])*(inc)-(PR[t]-PR[t-1])*(inc))/(((dist+inc)*(dist+inc)-dist*dist)*(inc)-(dist*dist-(dist-inc)*(dist-inc))*(inc));
			        b=(PR[t+1]-PR[t]-a*((dist+inc)*(dist+inc)-dist*dist))/(inc);
			        c=PR[t]-b*dist-a*dist*dist;
			        dist2=dist-inc*0.5+inc*0.5/101.0;
				for (m=start;m<start+101;m++)
				{
				        PR2[m]=a*dist2*dist2+b*dist2+c;
				        dist2=dist2+inc/101.0;
				}
			}
		        start=start+101;
		        dist=dist+inc;
		}

		dist=0;

		for (m=0;m<endfit+1;m++)
		{
			dist=inc/101.0;
			i[m]=0;
			for (n=1;n<4500;n++)
			{
				i[m]+=PR2[n]*sin(s[m]*dist)/(s[m]*dist);
				dist+=inc/101.0;
			}
		}
		chisqr=normalize(i, expi, 0, endfit, points, error);
		cout <<"chisqr= "<<chisqr<<" TotalPenalty= "<<TotalPenalty*PenaltyCoefficient<<endl;
		chisqr+=TotalPenalty*PenaltyCoefficient;

		RandNum=Real(rand())/Real(RAND_MAX);
		dChisqr=chisqr-old_chisqr;
		cout <<"RandNum= "<<RandNum<<" pick= "<<pick<<" move= "<<move<<" dChisqr= "<<dChisqr<<" exp= "<<exp(-dChisqr/Temp)<<" Temp= "<<Temp<<endl;

		if ( RandNum<exp(-dChisqr/Temp) )
		{
			old_chisqr=chisqr;
			cout <<"Move accepted"<<endl;
			NumAccepted+=1.0;
		}
		else
		{
			PR[pick]-=move;
			cout <<"Move rejected"<<endl;
			NumRejected+=1.0;
		}

		if (NumIncrease<=100 && dChisqr>0)
		{
			if (NumIncrease<100)
			{
				InitTempArray[NumIncrease]=dChisqr;
			}
			NumIncrease++;
			if (NumIncrease==100)
			{
				Temp2=0.1;
				while (true)
				{
					Sum=0;
					for (n=0;n<100;n++)
					{
						Sum+=exp(-InitTempArray[n]/Temp2);
					}
					cout <<"Sum= "<<Sum<<endl;
					if (Sum<1.0) Temp2=Temp2*1.01;
					else
					{
						Temp=Temp2;
						cout <<"Temp= "<<Temp<<" Sum= "<<Sum<<endl;
						break;
					}
				}
			}
		}
		

		TotalChisqr+=chisqr;
		TotalChisqrSqr+=chisqr*chisqr;

		if (chisqr<bestchisqr)
		{
			cout <<"chisqr= "<<chisqr<<endl;
			bestchisqr=chisqr;
			for (n=0;n<44;n++)
			{
				fitparameter[n]=PR[n];
			}

			for (n=0;n<4500;n++)
			{
				FitPr2[n]=PR2[n];
			}
		}

		if (NumAccepted>1000)
		{
			TotalMoves=NumAccepted+NumRejected;
			cout <<"TotalChisqrSqr= "<<TotalChisqrSqr<<" TotalMoves= "<<TotalMoves<<" TotalChisqr= "<<TotalChisqr<<endl;
			StnDev=sqrt( TotalChisqrSqr/TotalMoves-TotalChisqr*TotalChisqr/(TotalMoves*TotalMoves) );
			Temp=Temp/(1+Temp*delta/(3.0*StnDev));
			cout <<"Fraction accepted= "<<NumAccepted/(TotalMoves);
			TotalChisqrSqr=0;
			TotalChisqr=0;
			NumAccepted=0;
			NumRejected=0;
		}
        }

	IntensityFile=PrFile+"Intensity.txt";
	PrFile2=PrFile+"_2.txt";
	strcpy(CharPrFile, PrFile.c_str());
	strcpy(CharPrFile2, PrFile2.c_str());
	strcpy(CharIntensityFile, IntensityFile.c_str());
	
	dist=0;
        for (t=0;t<45;t++)
        {
                ofstream pofr(CharPrFile, ios::app);
                pofr << dist <<"\t"<<fitparameter[t]<<endl;
		dist+=inc;
        }
        ofstream pofr(CharPrFile, ios::app);
        pofr << endl;

	dist=0;
	for (t=0;t<4500;t++)
	{
		ofstream pofr(CharPrFile2, ios::app);
		pofr << dist <<"\t"<<FitPr2[t]<<endl;
		dist+=inc/101;
	}

	for (t=0;t<endfit;t++)
	{
		ofstream intensity2(CharIntensityFile, ios::app);
		intensity2 <<s[t]<<"\t"<<i[t]<<"\t"<<expi[t]<<endl;
	}

	ofstream intensity2(CharIntensityFile, ios::app);
	intensity2 << endl;

        return inc;
}
********************************************************************/
/*******************************************************************
void FitPr4(int BeginFit, int EndFit, int NumVariables, int points, Real s[], Real expi[])
{
	int n, m, t;
	int NumPrPoints, NumIterations, pick;
	int MaxIterations;
	Real BestChisqr, chisqr, chisqr_old;
	Real ConvergenceCriterion;
	Real dChisqr, delta;
	Real derivative2;
	Real dist, Dmax, inc;
	Real move, MaxMove;
	Real RandNum, StnDev, Temp;
	Real TotalChisqr, TotalChisqrSqr;
	Real TotalPenalty, TotalPr;
	Real *coefficient, *BestCoefficient;
	Real **imn;
	Real *PR, *BestPR, *i;

	NumPrPoints=5000;
	BestChisqr=10000000;
	MaxIterations=500000;
	ConvergenceCriterion=0.1;
	delta=1.0;
	MaxMove=0.1;
	Temp=200.0;
	Dmax=45.0;
	TotalPenalty=1000000.0;
	
	inc=Dmax/Real(NumPrPoints);

	coefficient=new Real[NumVariables];
	BestCoefficient=new Real[NumVariables];
	PR=new Real[NumPrPoints];
	BestPR=new Real[NumPrPoints];
	i=new Real[points];

	imn=new Real *[NumVariables];
	for (n=0;n<NumVariables;n++)
	{
		*(imn+n)=new Real[points];
	}

	cout <<"Allocated memory"<<endl;
	for (n=0;n<NumPrPoints;n++)
	{
		PR[n]=0;
		BestPR[n]=0;
	}

	for (n=0;n<NumVariables;n++)
	{
		coefficient[n]=0;
		BestCoefficient[n]=0;
	}

	coefficient[0]=1.0;
	cout <<"Initialized arrays"<<endl;

	for (m=0;m<NumVariables;m++)
	{
		for (n=0;n<points;n++)
		{
			dist=inc;
			for (t=1;t<NumPrPoints;t++)
			{
				imn[m][n]+=sin(pi*Real((m+1)*t)/Real(NumPrPoints-1))*sin(s[n]*dist)/(s[n]*dist);
				dist+=inc;
			}
		}
	}
	
	for (m=0;m<NumVariables;m++)
	{
		for (n=0;n<EndFit;n++)
		{
			i[n]+=coefficient[m]*imn[m][n];
		}
	}	

	cout <<"BeginFit= "<<BeginFit<<" EndFit= "<<EndFit<<endl;
	chisqr_old=normalize(i, expi, BeginFit, EndFit, points, error);

	for (n=0;n<NumVariables;n++);
	{
		cout <<"coefficeint["<<n<<"]= "<<coefficient[n]<<endl;
	}
	cout <<endl;

	NumIterations=0;
	while (NumIterations<MaxIterations)
	{
		RandNum=rand();
		pick=int(Real(RandNum*NumVariables)/Real(RAND_MAX));
		move=(Real(rand())/Real(RAND_MAX)-0.5)*MaxMove;
		NumIterations++;

		coefficient[pick]+=move;

		for (n=0;n<NumPrPoints;n++)
		{
			PR[n]=0;
		}

		for (m=0;m<NumVariables;m++)
		{
			for (n=beginfit;n<endfit;n++)
			{
				i[n]+=coefficient[m]*imn[m][n];
			}
		}

		chisqr=normalize(i, expi, beginfit, endfit, points, error);

		for (m=0;m<NumVariables;m++)
		{
			for (n=0;n<NumPrPoints;n++)
			{
				PR[n]+=coefficient[m]*sin(pi*Real((m+1)*n)/Real(NumPrPoints-1));
			}
		}
		TotalPr=0;
		for (n=0;n<NumPrPoints;n++)
		{
			TotalPr+=PR[n];
		}

		for (n=0;n<NumPrPoints;n++)
		{
			PR[n]=PR[n]/TotalPr;
		}
		
		TotalPenalty=0;
		for (n=1;n<NumPrPoints-1;n++)
		{
			derivative2=(PR[n-1]-2.0*PR[n]+PR[n+1])/(inc*inc);
			TotalPenalty+=derivative2*derivative2;
		}
		Real calc=0;
		Real product, product2;
		for (n=0;n<NumVariables;n++)
		{
			product=pi*Real(n+1)/Dmax;
			product2=product*product;
			calc+=0.5*coefficient[n]*coefficient[n]*product2*product2*Dmax;
		}

		cout <<"chisqr= "<<chisqr<<" TotalPenalty= "<<TotalPenalty*PenaltyCoefficient<<" calc= "<<calc*PenaltyCoefficient<<endl;
		chisqr+=TotalPenalty*PenaltyCoefficient;

		TotalChisqr+=chisqr;
		TotalChisqrSqr+=chisqr*chisqr;
		

		if (chisqr<BestChisqr)
		{
			cout <<"New best chisqr= "<<chisqr<<endl;
			BestChisqr=chisqr;
			for (n=0;n<NumVariables;n++)
			{
				BestCoefficient[n]=coefficient[n];
			}
		}
		
		dChisqr=chisqr-chisqr_old;
		RandNum=Real(rand())/Real(RAND_MAX);
		if ( RandNum < exp(-dChisqr/Temp) )
		{
			chisqr_old=chisqr;
		}
		else
		{
			coefficient[pick]-=move;
		}

		if (NumIterations%1000==0)
		{
			StnDev=sqrt(TotalChisqrSqr/1000-TotalChisqr*TotalChisqr/(1000*1000));
			Temp=Temp/(1.0+Temp*delta/(3.0*StnDev));
			cout <<"StnDev= "<<StnDev<<" Criterion= "<<ConvergenceCriterion*Temp<<endl;
			if (StnDev<ConvergenceCriterion*Temp) break;
			TotalChisqr=0;
			TotalChisqrSqr=0;
			cout <<"Temp= "<<Temp<<endl;
		}
	}

	for (n=0;n<points;n++)
	{
		i[n]=0;
	}
	
	for (m=0;m<NumVariables;m++)
	{
		for (n=0;n<points;n++)
		{
			i[n]+=BestCoefficient[m]*imn[m][n];
		}
	}

	normalize(i, expi, beginfit, endfit, points, error);

	for (m=0;m<NumVariables;m++)
	{
		for (n=0;n<NumPrPoints;n++)
		{
			PR[n]+=BestCoefficient[m]*sin(pi*Real((m+1)*n)/Real(NumPrPoints-1));
		}
	}

	dist=0;
	for (n=0;n<NumPrPoints;n++)
	{
		ofstream pofr("/home2/jouko/WAXS/pofr_ubiquitin.txt", ios::app);
		pofr << dist<<"\t"<<PR[n]<<endl;
		dist+=inc;
	}

	for (n=0;n<points;n++)
	{
		ofstream intensity("/home2/jouko/WAXS/pofr_intensity.txt", ios::app);
		intensity <<s[n]<<"\t"<<i[n]<<"\t"<<expi[n]<<endl;
	}
}
**************************************************************************/

void PrintIntensity(string IntensityFile, Real s[], Real i[], int points)
{
	int n;
	char CharIntensityFile[100];
	
	AddIndexToFile(IntensityFile);

	strcpy(CharIntensityFile, IntensityFile.c_str());

        for (n=0;n<points;n++)
        {
                ofstream pofr(CharIntensityFile, ios::app);
                pofr << s[n]<<"\t"<<i[n]<<endl;
        }
}

Real FitPr5(int NumVariables, int BeginFit, int EndFit, int NumPrPoints, Real Dmax, vector<Real> &s, vector<Real> &i, vector<Real> &expi, vector<Real> &PR, Real PenaltyCoefficient)
{
	//cout <<"In FitPr5"<<endl;
	int n, m, t;
	Real chisqr;
	Real dist, inc;
	Real **imn;
        vector<long double> variable;
        vector< vector<long double> > Coefficient;

        Safe2DAlloc(Coefficient, NumVariables+1, NumVariables, "Coefficient");
        SafeAlloc(variable, NumVariables, "variables");

	cout <<"In FitPr5"<<endl;	
	inc=Dmax/Real(NumPrPoints-1);
	//cout <<"BeginFit= "<<BeginFit<<" EndFit= "<<EndFit<<" s= "<<s[EndFit]<<endl;

	imn=new Real *[NumVariables];
	for (n=0;n<NumVariables;n++)
	{
		*(imn+n)=new Real[EndFit];
	}
	cout <<"Allocated memory"<<endl;
	for (n=0;n<NumVariables;n++)
	{
		for (m=0;m<EndFit;m++)
		{
			imn[n][m]=0;
		}
	}
	//cout <<"Initialized Coefficient"<<endl;
	//cout <<"Allocated memory"<<endl;
	for (n=0;n<NumPrPoints;n++)
	{
		PR[n]=0;
	}
	cout <<"Initialized PR"<<endl;

	//cout <<"Initialized arrays"<<endl;

	for (m=0;m<NumVariables;m++)
	{
		for (n=0;n<EndFit;n++)
		{
			dist=inc;
			for (t=1;t<NumPrPoints;t++)
			{
				if (m==3 && n==33)
				{
					//cout <<"imn["<<m<<"]["<<n<<"]= "<<imn[m][n]<<" s= "<<s[n]<<" dist= "<<dist<<endl;
					//cout <<"s*dist= "<<s[n]*dist<<endl;
					//cout <<"sin(s*dist)= "<<s[n]*dist<<endl;
					//cout <<"NumPrPoints= "<<NumPrPoints<<endl;
					//cout <<"sin(pi...= "<<sin(pi*Real((m+1)*t)/Real(NumPrPoints-1))<<endl;
				}
				imn[m][n]+=sin(pi*Real((m+1)*t)/Real(NumPrPoints-1))*sin(s[n]*dist)/(s[n]*dist);
				dist+=inc;
			}
		}
	}
	cout <<"Calculated imn"<<endl;
	for (m=0;m<NumVariables;m++)
	{
		for (n=0;n<EndFit;n++)
		{
			imn[m][n]*=inc;
		}
	}

	for (m=0;m<NumVariables;m++)
	{
		for (n=0;n<NumVariables;n++)
		{
			for (t=BeginFit;t<EndFit;t++)
			{
				Coefficient[m][n]+=imn[m][t]*imn[n][t];
				//cout <<"imn["<<m<<"]["<<t<<"]= "<<imn[m][t]<<" imn["<<n<<"]["<<t<<"]= "<<imn[n][t]<<endl;
			}	
		}
	}
	cout <<"Calculated coefficient1"<<endl;
	//PrintMatrix(Coefficient);
	//for (n=0;n<NumVariables;n++)
	//{
	//	cout <<"Coefficient["<<n<<"]["<<n<<"]= "<<Coefficient[n][n]<<endl;
	//}
	//cout <<endl;

	for (n=0;n<NumVariables;n++)
	{
		for (t=BeginFit;t<EndFit;t++)
		{
			Coefficient[NumVariables][n]+=expi[t]*imn[n][t];
		}
	}
	cout <<"Calculated coefficient2"<<endl;
	//PrintMatrix(Coefficient);
	
	for (n=0;n<NumVariables;n++)
	{
		//Coefficient[n][n]+=Real((n+1)*(n+1)*(n+1)*(n+1))*PenaltyCoefficient;
		Coefficient[n][n]+=Real((n+1)*(n+1))*PenaltyCoefficient*Dmax*Dmax;
	}

	cout <<"Calculated coefficient3"<<endl;
        //for (n=0;n<NumVariables;n++)
        //{
        //        cout <<"Coefficient["<<n<<"]["<<n<<"]= "<<Coefficient[n][n]<<endl;
        //}
        //cout <<endl;

	//PrintMatrix(Coefficient);
	Equation(Coefficient, variable, NumVariables);
	cout <<"Solvented Equation"<<endl;
	for (n=0;n<EndFit;n++)
	{
		i[n]=0;
	}
	
	for (m=0;m<NumVariables;m++)
	{
		//cout <<"variable["<<m<<"]= "<<variable[m]<<endl;
		for (n=0;n<EndFit;n++)
		{
			i[n]+=variable[m]*imn[m][n];
		}
	}
	//PrintIntensity("/home2/jouko/WAXS/PrFromScatteringTestIntensity.txt", s, i, EndFit);	
	chisqr=0;
        //cout <<"BeginFit= "<<BeginFit<<" EndFit= "<<EndFit<<endl;
	for (n=BeginFit;n<EndFit;n++)
	{
		chisqr+=(expi[n]-i[n])*(expi[n]-i[n]);
                //cout <<"expi["<<n<<"]= "<<expi[n]<<" i["<<n<<"]= "<<i[n]<<endl;
	}
	chisqr/=(EndFit-BeginFit+1);
	//chisqr=normalize(i, expi, beginfit, endfit, points, error);
	Real calc=0;
	Real TotalVariable2=0;
	for (n=0;n<NumVariables;n++)
	{
		TotalVariable2=variable[n]*variable[n];
	}
	for (n=0;n<NumVariables;n++)
	{
                //cout <<"variable["<<n<<"]= "<<variable[n]<<endl;
		calc+=variable[n]*variable[n]*Real((n+1)*(n+1))*PenaltyCoefficient*Dmax*Dmax/(TotalVariable2);
	}

	//cout <<"chisqr= "<<chisqr<<" penalty= "<<calc<<" total= "<<chisqr+calc<<endl;
	chisqr+=calc;

	for (m=0;m<NumVariables;m++)
	{
		for (n=0;n<NumPrPoints;n++)
		{
			PR[n]+=variable[m]*sin(pi*Real((m+1)*n)/Real(NumPrPoints-1));
			if (m==1)
			{
				//cout <<"PR["<<n<<"]= "<<PR[n]<<" variable[1]= "<<variable[m]<<" NumPrPoints= "<<NumPrPoints<<endl;
			}
		}
	}
	//for (n=0;n<NumPrPoints;n++) cout <<"PR["<<n<<"]= "<<PR[n]<<endl;
	for (n=0;n<NumVariables;n++)
	{
                delete [] imn[n];
	}
        delete [] imn;
	return chisqr;

}

void CalculatePrFromVariable(vector<Real> &variable, vector<Real> &PR)
{
        int NumPrPoints=PR.size();
	int NumVariables=variable.size();
        for (int m=0;m<NumVariables;m++)
	{
		for (int n=0;n<NumPrPoints;n++)
		{
			PR[n]+=variable[m]*sin(pi*Real((m+1)*n)/Real(NumPrPoints-1));
		}
	}
}

Real CalculateNegativeArea(vector<Real> &pr)
{
        int NumPrPoints=pr.size();
        Real NegativeArea=0;
        for (int i=0;i<NumPrPoints;i++)
        {
                if (pr[i]<0) NegativeArea+=pr[i];
        }
        return -NegativeArea;
}

Real CalculateAbsValueArea(vector<Real> &pr)
{
        int NumPrPoints=pr.size();
        Real AbsArea=0;
        for (int i=0;i<NumPrPoints;i++)
        {
                AbsArea+=abs(pr[i]);
        }
        return AbsArea;
}

Real CalculateNegativePenalty(vector<Real> &variable)
{
        int NumPrPoints;
        Real NegativeArea, AbsArea;
        vector<Real> pr;
        NumPrPoints=5000;
        SafeAlloc(pr, 5000, "pr in CalculateNegativePenalty");
        CalculatePrFromVariable(variable, pr);
        NegativeArea=CalculateNegativeArea(pr);
        AbsArea=CalculateAbsValueArea(pr);
        return NegativeArea/AbsArea;
}

Real EvaluateVariables(vector<Real> &variable, Real **&imn, vector<Real> &i, vector<Real> &expi, vector<Real> &error, int EndFit, int BeginFit, Real Dmax, Real PenaltyCoefficient)
{
        bool UseErrors=false;
	int NumVariables=variable.size();
	Real chisqr;
	for (int n=0;n<EndFit;n++)
	{
		i[n]=0;
	}
		
	for (int m=0;m<NumVariables;m++)
	{
		//cout <<"variable["<<m<<"]= "<<variable[m]<<endl;
		for (int n=0;n<EndFit;n++)
		{
			i[n]+=variable[m]*imn[m][n];
		}
	}
	//PrintIntensity("/home2/jouko/WAXS/PrFromScatteringTestIntensity.txt", s, i, EndFit);	
	chisqr=0;
        //cout <<"BeginFit= "<<BeginFit<<" EndFit= "<<EndFit<<endl;
	for (int n=BeginFit;n<EndFit;n++)
	{
                if (UseErrors) chisqr+=(expi[n]-i[n])*(expi[n]-i[n])/(error[n]*error[n]);
                else chisqr+=(expi[n]-i[n])*(expi[n]-i[n]);
                //cout <<"expi["<<n<<"]= "<<expi[n]<<" i["<<n<<"]= "<<i[n]<<endl;
	}
	chisqr/=(EndFit-BeginFit+1);
	//chisqr=normalize(i, expi, beginfit, endfit, points, error);
	Real calc=0;
	Real TotalVariable2=0;
	for (int n=0;n<NumVariables;n++)
	{
		TotalVariable2+=variable[n]*variable[n];
	}
	for (int n=0;n<NumVariables;n++)
	{
                //cout <<"variable["<<n<<"]= "<<variable[n]<<endl;
		calc+=variable[n]*variable[n]*Real((n+1)*(n+1))*PenaltyCoefficient*Dmax*Dmax/(TotalVariable2);
	}
        Real NegativePenalty;
        NegativePenalty=CalculateNegativePenalty(variable)*100000*PenaltyCoefficient;
	//cout <<"chisqr= "<<chisqr<<" penalty= "<<calc<<" total= "<<chisqr+calc<<endl;
	chisqr+=calc+NegativePenalty;
	return chisqr;
}

Real EvaluateVariables(vector<Real> variable, ArgStruct args)
{
        bool UseErrors=false;
        int EndFit, BeginFit;
	int NumVariables=variable.size();
        Real Dmax, PenaltyCoefficient;
        Real **imn;
        vector<Real> i, expi, error;
	Real chisqr;

        EndFit=args.EndFit;
        BeginFit=args.BeginFit;
        Dmax=args.Dmax;
        PenaltyCoefficient=args.PenaltyCoefficient;
        i=args.i;
        expi=args.expi;
        error=args.error;
        imn=args.imn;

	for (int n=0;n<EndFit;n++)
	{
		i[n]=0;
	}
		
	for (int m=0;m<NumVariables;m++)
	{
		//cout <<"variable["<<m<<"]= "<<variable[m]<<endl;
		for (int n=0;n<EndFit;n++)
		{
			i[n]+=variable[m]*imn[m][n];
		}
	}
	//PrintIntensity("/home2/jouko/WAXS/PrFromScatteringTestIntensity.txt", s, i, EndFit);	
	chisqr=0;
        //cout <<"BeginFit= "<<BeginFit<<" EndFit= "<<EndFit<<endl;
	for (int n=BeginFit;n<EndFit;n++)
	{
		if (UseErrors) chisqr+=(expi[n]-i[n])*(expi[n]-i[n])/(error[n]*error[n]);
		else chisqr+=(expi[n]-i[n])*(expi[n]-i[n]);
                //cout <<"error["<<n<<"]= "<<error[n]<<endl;
                //cout <<"expi["<<n<<"]= "<<expi[n]<<" i["<<n<<"]= "<<i[n]<<endl;
	}
	chisqr/=(EndFit-BeginFit+1);
	//chisqr=normalize(i, expi, beginfit, endfit, points, error);
	Real calc=0;
	Real TotalVariable2=0;
	for (int n=0;n<NumVariables;n++)
	{
		TotalVariable2+=variable[n]*variable[n];
	}
	for (int n=0;n<NumVariables;n++)
	{
                //cout <<"variable["<<n<<"]= "<<variable[n]<<endl;
		calc+=variable[n]*variable[n]*Real((n+1)*(n+1))*PenaltyCoefficient*Dmax*Dmax/(TotalVariable2);
	}
        Real NegativePenalty;
        NegativePenalty=CalculateNegativePenalty(variable)*100000*PenaltyCoefficient;

	//cout <<"chisqr= "<<chisqr<<" penalty= "<<calc<<" total= "<<chisqr+calc<<endl;
	chisqr+=calc+NegativePenalty;
	return chisqr;
}

Real MCOptVariable(vector<Real> &variable, Real **&imn, vector<Real> &i, vector<Real> &expi, vector<Real> &error, int EndFit, int BeginFit, Real Dmax, Real PenaltyCoefficient)
{
	int pick, NumVariables=variable.size();
	int NumIterations=0, MaxIterations=1000;
	Real Temp=1000, Rand;
	Real chisqr, chisqr_old, BestChisqr, dChisqr, move, MaxMove;
	vector<Real> best_variable;
	chisqr_old=EvaluateVariables(variable, imn, i, expi, error, EndFit, BeginFit, Dmax, PenaltyCoefficient);
	best_variable=variable;
        BestChisqr=chisqr_old;
	while (NumIterations<MaxIterations)
	{
		pick=RandInt(NumVariables-1);
                cout <<"pick= "<<pick<<endl;
                cout <<"variable["<<pick<<"]= "<<variable[pick]<<endl;
		MaxMove=variable[pick]*0.1;
		move=MaxMove*(Real(rand())/Real(RAND_MAX)-0.5);
                cout <<"move= "<<move<<" MaxMove= "<<MaxMove<<" RAND_MAX= "<<RAND_MAX<<endl;
		variable[pick]+=move;
		chisqr=EvaluateVariables(variable, imn, i, expi, error, EndFit, BeginFit, Dmax, PenaltyCoefficient);
		
		dChisqr=chisqr-chisqr_old;
		Rand=Real(rand())/Real(RAND_MAX);

		if (dChisqr<0) chisqr_old=chisqr;
		else
		{
			if (Rand < exp(-dChisqr/Temp)) chisqr_old=chisqr;
			else variable[pick]-=move;
		}

		if (chisqr<BestChisqr)
		{
			BestChisqr=chisqr;
			for (int n=0;n<NumVariables;n++)
			{
				best_variable[n]=variable[n];
			}
		}
	
		NumIterations++;
	}

	for (int n=0;n<NumVariables;n++)
	{
		variable[n]=best_variable[n];
	}

}

Real SteepestDescent(vector<Real> &variable, Real **&imn, vector<Real> &i, vector<Real> &expi, vector<Real> &error, int EndFit, int BeginFit, Real Dmax, Real PenaltyCoefficient)
{
        Real (*func_to_minimize)(vector<Real>, ArgStruct)=EvaluateVariables;
        Real ConvergenceCriteria=10e-8;
        ArgStruct args;
        
        args.imn=imn;
        args.i=i;
        args.expi=expi;
        args.error=error;
        args.EndFit=EndFit;
        args.BeginFit=BeginFit;
        args.Dmax=Dmax;
        args.PenaltyCoefficient=PenaltyCoefficient;

        OptimizeCoefficients(variable, args, func_to_minimize, ConvergenceCriteria);
}

Real FitPr6(int NumVariables, int BeginFit, int EndFit, int NumPrPoints, Real Dmax, vector<Real> &s, vector<Real> &i, vector<Real> &expi, vector<Real> &error, vector<Real> &PR, Real PenaltyCoefficient)
{
	//cout <<"In FitPr6"<<endl;
        bool UseErrors=false;
	int n, m, t;
	Real chisqr;
	Real dist, inc;
	Real **imn;
        vector<Real> variable;
	vector< vector<Real> > Coefficient;
        SafeAlloc(variable, NumVariables, "variables");
	Safe2DAlloc(Coefficient, NumVariables+1, NumVariables, "Coefficient");

	cout <<"In FitPr6"<<endl;	
	inc=Dmax/Real(NumPrPoints-1);
	//cout <<"BeginFit= "<<BeginFit<<" EndFit= "<<EndFit<<" s= "<<s[EndFit]<<endl;

	imn=new Real *[NumVariables];
	for (n=0;n<NumVariables;n++)
	{
		*(imn+n)=new Real[EndFit];
	}
	cout <<"Allocated memory"<<endl;
	for (n=0;n<NumVariables;n++)
	{
		for (m=0;m<EndFit;m++)
		{
			imn[n][m]=0;
		}
	}
	//cout <<"Initialized Coefficient"<<endl;
	//cout <<"Allocated memory"<<endl;
	for (n=0;n<NumPrPoints;n++)
	{
		PR[n]=0;
	}
	cout <<"Initialized PR"<<endl;

	//cout <<"Initialized arrays"<<endl;

	for (m=0;m<NumVariables;m++)
	{
		for (n=0;n<EndFit;n++)
		{
			dist=inc;
			for (t=1;t<NumPrPoints;t++)
			{
				if (m==3 && n==33)
				{
					//cout <<"imn["<<m<<"]["<<n<<"]= "<<imn[m][n]<<" s= "<<s[n]<<" dist= "<<dist<<endl;
					//cout <<"s*dist= "<<s[n]*dist<<endl;
					//cout <<"sin(s*dist)= "<<s[n]*dist<<endl;
					//cout <<"NumPrPoints= "<<NumPrPoints<<endl;
					//cout <<"sin(pi...= "<<sin(pi*Real((m+1)*t)/Real(NumPrPoints-1))<<endl;
				}
				imn[m][n]+=sin(pi*Real((m+1)*t)/Real(NumPrPoints-1))*sin(s[n]*dist)/(s[n]*dist);
				dist+=inc;
			}
		}
	}
	cout <<"Calculated imn"<<endl;
	for (m=0;m<NumVariables;m++)
	{
		for (n=0;n<EndFit;n++)
		{
			imn[m][n]*=inc;
		}
	}

	for (m=0;m<NumVariables;m++)
	{
		for (n=0;n<NumVariables;n++)
		{
			for (t=BeginFit;t<EndFit;t++)
			{
				Coefficient[m][n]+=imn[m][t]*imn[n][t];
				//cout <<"imn["<<m<<"]["<<t<<"]= "<<imn[m][t]<<" imn["<<n<<"]["<<t<<"]= "<<imn[n][t]<<endl;
			}	
		}
	}
	cout <<"Calculated coefficient1"<<endl;
	//PrintMatrix(Coefficient);
	//for (n=0;n<NumVariables;n++)
	//{
	//	cout <<"Coefficient["<<n<<"]["<<n<<"]= "<<Coefficient[n][n]<<endl;
	//}
	//cout <<endl;

	for (n=0;n<NumVariables;n++)
	{
		for (t=BeginFit;t<EndFit;t++)
		{
			Coefficient[NumVariables][n]+=expi[t]*imn[n][t];
		}
	}
	cout <<"Calculated coefficient2"<<endl;
	//PrintMatrix(Coefficient);
	
	for (n=0;n<NumVariables;n++)
	{
		//Coefficient[n][n]+=Real((n+1)*(n+1)*(n+1)*(n+1))*PenaltyCoefficient;
		//Coefficient[n][n]+=Real((n+1)*(n+1))*PenaltyCoefficient*Dmax*Dmax;
	}

	cout <<"Calculated coefficient3"<<endl;
        //for (n=0;n<NumVariables;n++)
        //{
        //        cout <<"Coefficient["<<n<<"]["<<n<<"]= "<<Coefficient[n][n]<<endl;
        //}
        //cout <<endl;

	//PrintMatrix(Coefficient);
	Equation(Coefficient, variable, NumVariables);

	MCOptVariable(variable, imn, i, expi, error, EndFit, BeginFit, Dmax, PenaltyCoefficient);
	cout <<"About to enter SteepestDescent"<<endl;
        SteepestDescent(variable, imn, i, expi, error, EndFit, BeginFit, Dmax, PenaltyCoefficient);
	cout <<"Solvented Equation"<<endl;
	for (n=0;n<EndFit;n++)
	{
		i[n]=0;
	}
	
	for (m=0;m<NumVariables;m++)
	{
		//cout <<"variable["<<m<<"]= "<<variable[m]<<endl;
		for (n=0;n<EndFit;n++)
		{
			i[n]+=variable[m]*imn[m][n];
		}
	}
	//PrintIntensity("/home2/jouko/WAXS/PrFromScatteringTestIntensity.txt", s, i, EndFit);	
	chisqr=0;
        //cout <<"BeginFit= "<<BeginFit<<" EndFit= "<<EndFit<<endl;
	for (n=BeginFit;n<EndFit;n++)
	{
		if (UseErrors) chisqr+=(expi[n]-i[n])*(expi[n]-i[n])/(error[n]*error[n]);
		else chisqr+=(expi[n]-i[n])*(expi[n]-i[n]);
                //cout <<"expi["<<n<<"]= "<<expi[n]<<" i["<<n<<"]= "<<i[n]<<endl;
	}
	chisqr/=(EndFit-BeginFit+1);
	//chisqr=normalize(i, expi, beginfit, endfit, points, error);
	Real calc=0;
	Real TotalVariable2=0;
	for (n=0;n<NumVariables;n++)
	{
		TotalVariable2=variable[n]*variable[n];
	}
	for (n=0;n<NumVariables;n++)
	{
                //cout <<"variable["<<n<<"]= "<<variable[n]<<endl;
		calc+=variable[n]*variable[n]*Real((n+1)*(n+1))*PenaltyCoefficient*Dmax*Dmax/(TotalVariable2);
	}
        Real NegativePenalty;
        NegativePenalty=CalculateNegativePenalty(variable)*100000*PenaltyCoefficient;

	//cout <<"chisqr= "<<chisqr<<" penalty= "<<calc<<" total= "<<chisqr+calc<<endl;
	chisqr+=calc+NegativePenalty;

	for (m=0;m<NumVariables;m++)
	{
		for (n=0;n<NumPrPoints;n++)
		{
			PR[n]+=variable[m]*sin(pi*Real((m+1)*n)/Real(NumPrPoints-1));
			if (m==1)
			{
				//cout <<"PR["<<n<<"]= "<<PR[n]<<" variable[1]= "<<variable[m]<<" NumPrPoints= "<<NumPrPoints<<endl;
			}
		}
	}
	//for (n=0;n<NumPrPoints;n++) cout <<"PR["<<n<<"]= "<<PR[n]<<endl;
	for (n=0;n<NumVariables;n++)
	{
                delete [] imn[n];
	}
        delete [] imn;
	return chisqr;

}

void PrintPr(string PrFile, Real Dmax, int NumPrPoints, vector<Real> &PR)
{
	char CharPrFile[1000];
	int n; 
	Real dist, inc;
	AddIndexToFile(PrFile);

	strcpy(CharPrFile, PrFile.c_str());

	inc=Dmax/Real(NumPrPoints-1);

        dist=0;
        for (n=0;n<NumPrPoints;n++)
        {
                ofstream pofr(CharPrFile, ios::app);
                pofr << dist<<"\t"<<PR[n]<<endl;
                dist+=inc;
        }
}

void PrFromScattering(int NumVariables, int BeginFit, int EndFit, int NumPrPoints, vector<Real> &s, vector<Real> &expi, vector<Real> &error, vector<Real> &PR, Real &Dmax, Real PenaltyCoefficient, string IntensityFile)
{
	//Finds P(r) which optimizes the fit to experimental data and is smooth.  Also can optimize Dmax.
	int n;
	int NumIterations, MaxSteps=100000, steps=0;;
	char CharIntensityFile[1000], CharPrFile[1000];
	Real BestChisqr, BestDmax;
	Real chisqr, chisqr_old, dChisqr;
	Real dist, inc;
	Real MaxMove, move, Temp, step;
	Real Rand;
	Real derivative, criterion=0.0001;
	vector<Real> i;

	NumIterations=0;
	BestChisqr=10000000;
        BestDmax=Dmax;
        SafeAlloc(i, EndFit, "i in PrFromScattering");

	step=0.5;
	MaxMove=Dmax*0.1;
	Temp=10000.0;
	cout <<"About to enter FitPr5"<<endl;
	//cout <<" Dmax= "<<Dmax<<" s[100]= "<<s[100]<<" i[100]= "<<i[100]<<" expi[100]= "<<expi[100]<<" PR[100]= "<<PR[100]<<endl;
	cout <<"NumVariables= "<<NumVariables<<endl;	
	BestChisqr=FitPr6(NumVariables, BeginFit, EndFit, NumPrPoints, Dmax, s, i, expi, error, PR, PenaltyCoefficient);
	cout <<"NumVariables= "<<NumVariables<<endl;
	cout <<"Finished with FirPr5"<<endl;
	cout <<"Entering MC"<<endl;
	cout <<"BestChisqr= "<<BestChisqr<<endl;
	chisqr_old=BestChisqr;
	while (NumIterations<7)
	{
		move=MaxMove*(Real(rand()/Real(RAND_MAX)-0.5));
		Dmax+=move;
		
		chisqr=FitPr6(NumVariables, BeginFit, EndFit, NumPrPoints, Dmax, s, i, expi, error, PR, PenaltyCoefficient);
		dChisqr=chisqr-chisqr_old;
		cout <<"chisqr= "<<chisqr<<" Dmax= "<<Dmax<<" BestDmax= "<<BestDmax<<" BestChisqr= "<<BestChisqr<<endl;
		cout <<"PR[10]= "<<PR[10]<<endl;
		Rand=Real(rand())/Real(RAND_MAX);

		if (dChisqr<0) chisqr_old=chisqr;
		else
		{
			if (Rand < exp(-dChisqr/Temp)) chisqr_old=chisqr;
			else Dmax-=move;
		}

		if (chisqr<BestChisqr)
		{
			BestChisqr=chisqr;
			BestDmax=Dmax;
		}
	
		NumIterations++;
	}

	Dmax=BestDmax;
	chisqr=BestChisqr;
	cout <<"Entering steepest descent"<<endl;
	/*
	while (step*step>0.00001 && steps<MaxSteps)
	{
		Dmax+=step;

		chisqr_old=chisqr;
		chisqr=FitPr5(NumVariables, BeginFit, EndFit, NumPrPoints, Dmax, s, i, expi, PR, PenaltyCoefficient);
		cout <<"PR[100]= "<<PR[100]<<endl;
		cout <<"chisqr= "<<chisqr<<" chisqr_old= "<<chisqr_old<<" Dmax= "<<Dmax<<endl;

		if (chisqr>chisqr_old)
		{
			Dmax=BestDmax;
			step=-step*0.5;
			chisqr=chisqr_old;
		}
		else 
		{
			step=step*1.2;
			BestDmax=Dmax;
		}
		derivative=(chisqr-chisqr_old)/step;
		if (derivative*derivative<criterion) break;
                steps++;
	}
	*/
	strcpy(CharIntensityFile, IntensityFile.c_str());

        for (n=0;n<EndFit;n++)
        {
                ofstream intensity(CharIntensityFile, ios::app);
                intensity <<s[n]<<"\t"<<i[n]<<"\t"<<expi[n]<<endl;
        }
}
#endif
