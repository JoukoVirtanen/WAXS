#ifndef _minimize_included_
#define _minimize_included_

//template<class Args>
//Real OptimizeCoefficients(vector<Real> &coefficient, Args args, Real (*func_to_minimize)(vector<Real>, Args), Real dScoreCriteria=10e-15);

template <class Args>
bool AtMaxima(vector<Real> coefficient, Real &TestStep, int index, Args args, Real (*func_to_minimize)(vector<Real>, Args))
{
	Real score1, score2, score3;

	coefficient[index]-=TestStep;
	score1=(*func_to_minimize)(coefficient, args);
	coefficient[index]+=TestStep;
	score2=(*func_to_minimize)(coefficient, args);
	coefficient[index]+=TestStep;
	score3=(*func_to_minimize)(coefficient, args);
	coefficient[index]-=TestStep;

	if (score2>score1 && score2>score3) return true;
	return false;
}

template <class Args>
bool AtMinima(vector<Real> coefficient, Real &TestStep, int index, Args args, Real (*func_to_minimize)(vector<Real>, Args))
{
	Real score1, score2, score3;

	coefficient[index]-=TestStep;
	score1=(*func_to_minimize)(coefficient, args);
	coefficient[index]+=TestStep;
	score2=(*func_to_minimize)(coefficient, args);
	coefficient[index]+=TestStep;
	score3=(*func_to_minimize)(coefficient, args);
	coefficient[index]-=TestStep;

	if (score2<score1 && score2<score3) return true;
	return false;
}
template <class Args>
Real CalcDerivative(vector<Real> coefficient, Real &TestStep, int index, Args args, Real (*func_to_minimize)(vector<Real>, Args))
{
	Real derivative, dScore;
	Real score, score_old;
        Real MinTestStep=10e-30;

	if (TestStep==0) TestStep=0.001;
	score_old=(*func_to_minimize)(coefficient, args);
	//cout <<"coefficient["<<index<<"]= "<<coefficient[index]<<"\t";
	coefficient[index]+=TestStep;
	score=(*func_to_minimize)(coefficient, args);
	//cout <<"In CalcDerivative"<<endl;
	//cout <<"score_old= "<<score_old<<"\t"<<"score= "<<score<<endl;
	dScore=score-score_old;
	if (dScore>0)
	{
		coefficient[index]-=TestStep;
		//cout <<"coefficient["<<index<<"]= "<<coefficient[index]<<endl;
		TestStep=-TestStep*0.05;
                if (abs(TestStep)<MinTestStep) return 0;
		derivative=CalcDerivative(coefficient, TestStep, index, args, func_to_minimize);
	}
	else
	{
		coefficient[index]-=TestStep;
		//cout <<"coefficient["<<index<<"]= "<<coefficient[index]<<endl;
		derivative=dScore/TestStep;
	}
	return derivative;
}

template <class Args>
Real CalcGradient(vector<Real> coefficient, Args args, Real (*func_to_minimize)(vector<Real>, Args), vector<Real> &step)
{
	int n, Size, min;
	Real MinTestStep, TestStep, TotalDerivative, TotalTestStep;
	vector<Real> derivative;
        
	Size=coefficient.size();
	CreateVector(derivative, Size);
	TotalDerivative=0;
	TotalTestStep=0;
	cout <<"In CalcGradient"<<endl;
	//for (n=0;n<Size;n++) cout <<"coefficient["<<n<<"]= "<<coefficient[n]<<endl;
	for (n=0;n<Size;n++)
	{
                cout <<"n= "<<n<<endl;
		//TestStep=0.1;
		TestStep=coefficient[n]*0.001;
		if (TestStep==0) TestStep=0.0001;
		derivative[n]=CalcDerivative(coefficient, TestStep, n, args, func_to_minimize);
		if (n==0) MinTestStep=TestStep;
		if (TestStep<MinTestStep) 
		{
			MinTestStep=TestStep;
			min=n;
		}
		//if (n==7) cout <<"TestStep= "<<TestStep<<" coefficient[7]= "<<coefficient[7]<<endl;
		TotalDerivative+=derivative[n]*derivative[n];
		TotalTestStep+=TestStep*TestStep;
	}

	for (n=0;n<Size;n++)
	{
		if (TotalDerivative!=0) step[n]=derivative[n]*MinTestStep/sqrt(TotalDerivative);
		else step[n]=coefficient[n]*10e-5;
		//cout <<"derivative["<<n<<"]= "<<derivative[n]<<endl;
		//cout <<"step["<<n<<"]= "<<step[n]<<endl;
	}
	//cout <<"derivative[0]= "<<derivative[0]<<endl;
	//cout <<"TotalTestStep= "<<TotalTestStep<<" TotalDerivative= "<<TotalDerivative<<endl;
	
	//for (n=0;n<Size;n++) cout <<"coefficient["<<n<<"]= "<<coefficient[n]<<endl;
	//cout <<"TotalDerivative= "<<TotalDerivative<<endl;
	return TotalDerivative;
}

template<class Args>
int CountMaxima(vector<Real> coefficient, Args args, Real (*func_to_minimize)(vector<Real>, Args))
{
	int n, NumMaxima, Size;
	Real TestStep;
	bool Maxima;
        cout <<"In CountMaxima"<<endl;
	Size=coefficient.size();
	NumMaxima=0;
	for (n=0;n<Size-1;n++)
	{
                cout <<"n= "<<n<<endl;
		TestStep=coefficient[n]*0.001;
		Maxima=AtMaxima(coefficient, TestStep, n, args, func_to_minimize);
		if (Maxima) NumMaxima++;
	}
	return NumMaxima;
}

template<class Args>
int CountMinima(vector<Real> coefficient, Args args, Real (*func_to_minimize)(vector<Real>, Args))
{
	int n, NumMinima, Size;
	Real TestStep;
	bool Minima;

	Size=coefficient.size();
	NumMinima=0;
	for (n=0;n<Size;n++)
	{
		TestStep=coefficient[n]*1e-30;
		Minima=AtMinima(coefficient, TestStep, n, args, func_to_minimize);
		if (Minima) NumMinima++;
	}
	return NumMinima;
}

template<class Args>
Real MinimizeAlongVector(vector<Real> &coefficient, vector<Real> step, Args args, Real (*func_to_minimize)(vector<Real>, Args))
{
        bool verbose=false;
	int m, n, NumMinima, Size;
	Real criteria, criteria2, dScoreCriteria;
	Real ConvergedScore;
	Real score, score_old, dScore;
	Real TotalDerivative;
	Real TestStep, Total, TotalStepSqr;
	vector<Real> BestCoefficient;
	if (verbose) cout <<"In MinimizeAlongVector"<<endl;
	dScoreCriteria=1e-15;
	Size=coefficient.size();
	CreateVector(BestCoefficient, Size);
	for (n=0;n<Size;n++) BestCoefficient[n]=coefficient[n];
	score_old=(*func_to_minimize)(coefficient, args);
	if (verbose) cout <<"score_old= "<<score_old<<endl;
	ConvergedScore=score_old;	
	while (true)
	{
		for (n=0;n<Size;n++) coefficient[n]+=step[n];

		//for (n=0;n<Size;n++) cout <<"step["<<n<<"]= "<<step[n]<<"\t"<<"coefficient= "<<"\t"<<coefficient[n]<<endl;
		//cout <<endl;
		//for (n=0;n<Size;n++) cout <<"step["<<n<<"]= "<<step[n]<<"\t"<<"coefficient= "<<"\t"<<coefficient[n]<<endl;
		//cout <<endl;
		//cout <<"Incremented coefficient"<<endl;
		//cout <<"Left SetY"<<endl;
		score=(*func_to_minimize)(coefficient, args);
		dScore=score-score_old;

		//for (n=0;n<Size;n++) cout <<"score= "<<score<<"\t"<<"coefficient["<<n<<"]= "<<coefficient[n]<<"\t"<<" step["<<n<<"]= "<<step[n]<<endl;
		//cout <<"score= "<<score<<"\t"<<"coefficient[0]= "<<coefficient[0]<<"\t"<<" step[0]= "<<step[0]<<endl;
		//cout <<"score= "<<score<<endl;
		if (score>score_old)
		{
			for (n=0;n<Size;n++)
			{
				coefficient[n]-=step[n];
				step[n]=-step[n]*0.5;
				coefficient[n]=BestCoefficient[n];
			}
		}
		else
		{
			score_old=score;
			for (n=0;n<Size;n++)
			{
				step[n]=step[n]*1.2;
				BestCoefficient[n]=coefficient[n];
			}
		}
		//cout <<"Updated step"<<endl;
		TotalStepSqr=0;
		for (n=0;n<Size;n++)
		{
			TotalStepSqr+=step[n]*step[n];
		}
		//if (verbose) cout <<"TotalStepSqr= "<<TotalStepSqr<<"\t"<<"dScore2= "<<dScore*dScore<<"\t"<<"eval= "<<dScore*dScore/TotalStepSqr<<"\t"<<"crit= "<<dScoreCriteria<<endl;
		if (TotalStepSqr==0) break;
		//dScoreSqr=dScore*dScore;
		if (dScore*dScore/TotalStepSqr<dScoreCriteria)
		{
			//cout <<"About to enter CalcGradient"<<endl;
			for (n=0;n<Size;n++) coefficient[n]=BestCoefficient[n];
			break;
		}
	}
}

template<class Args>
Real IterateOverCoefficients(vector<Real> &coefficient, Args args, Real (*func_to_minimize)(vector<Real>, Args))
{
        bool verbose;
	int Size=coefficient.size();
	Real TestStep;
	vector<Real> step;
	if (verbose) cout <<"In IterateOverCoefficients"<<endl;
	CreateVector(step, Size);
	for (int i=0;i<Size;i++)
	{
		TestStep=coefficient[i]*0.0001;
		for (int j=0;j<Size;j++) step[j]=0;
		CalcDerivative(coefficient, TestStep, i, args, func_to_minimize);
		step[i]=TestStep;
		MinimizeAlongVector(coefficient, step, args, func_to_minimize);
	}
}

template<class Args>
void MinimizeOneAtATime(vector<Real> &coefficient, Args args, Real (*func_to_minimize)(vector<Real>, Args))
{
        bool verbose=false;
	int Size, NumMinima, n;
	Real score, score_old;
	Real TotalGradient, dScoreCriterion, GradientCriterion;
	vector<Real> derivative;

	n=0;
	dScoreCriterion=1e-5;
	GradientCriterion=1e-30;	
	Size=coefficient.size();
	CreateVector(derivative, Size);
	score=(*func_to_minimize)(coefficient, args);
	if (verbose) cout <<"In MinimizeOneAtATime"<<endl;
	while (true)
	{
		score_old=score;
		IterateOverCoefficients(coefficient, args, func_to_minimize);
		TotalGradient=CalcGradient(coefficient, args, func_to_minimize, derivative);
		score=(*func_to_minimize)(coefficient, args);
		NumMinima=CountMinima(coefficient, args, func_to_minimize);
		if (n%1000==0 && verbose)
		{
			//cout <<"score_old= "<<score_old<<" score= "<<score<<" score_old-score= "<<score_old-score<<endl;
			//cout <<"TotalGradient= "<<TotalGradient<<endl;
			//cout <<"NumMinima= "<<NumMinima<<endl;
		}
		if (score_old-score<dScoreCriterion) break;
		if (TotalGradient<GradientCriterion) break;
		if (NumMinima==Size) break;
		n++;
	}
}

template<class Args>
Real OptimizeCoefficients(vector<Real> &coefficient, Args &args, Real (*func_to_minimize)(vector<Real>, Args), Real dScoreCriteria)
{
        bool verbose=true;
	int m, n, minSteps=10, NumMinima, numSteps=0, Size;
	Real criteria, criteria2;
	Real ConvergedScore;
	Real score, score_old, dScore;
	Real TotalDerivative;
	Real TestStep, Total, TotalStepSqr;
	vector<Real> step, BestCoefficient;

	criteria=0.00001;
	criteria2=1e-5;
	TestStep=0.001;
	Size=coefficient.size();
	CreateVector(step, Size);
	CreateVector(BestCoefficient, Size);
	for (n=0;n<Size;n++) BestCoefficient[n]=coefficient[n];
	score_old=(*func_to_minimize)(coefficient, args);
	if (verbose) cout <<"score_old= "<<score_old<<endl;
	CalcGradient(coefficient, args, func_to_minimize, step);
	if (verbose) cout <<"Done with CalcGradient"<<endl;
	ConvergedScore=score_old;	
	while (true)
	{
		for (n=0;n<Size;n++) coefficient[n]+=step[n];

		//for (n=0;n<Size;n++) cout <<"step["<<n<<"]= "<<step[n]<<"\t"<<"coefficient= "<<"\t"<<coefficient[n]<<endl;
		//cout <<endl;
		score=(*func_to_minimize)(coefficient, args);
		dScore=score-score_old;
                numSteps++;
		//for (n=0;n<Size;n++) cout <<"score= "<<score<<"\t"<<"coefficient["<<n<<"]= "<<coefficient[n]<<"\t"<<" step["<<n<<"]= "<<step[n]<<endl;
		cout <<"score= "<<score<<"\t"<<score_old<<"\t"<<"coefficient[0]= "<<coefficient[0]<<"\t"<<" step[0]= "<<step[0]<<endl;
		//cout <<"score= "<<score<<endl;
		if (score>score_old)
		{
			for (n=0;n<Size;n++)
			{
				coefficient[n]-=step[n];
				step[n]=-step[n]*0.5;
				coefficient[n]=BestCoefficient[n];
			}
		}
		else
		{
			score_old=score;
			for (n=0;n<Size;n++)
			{
				step[n]=step[n]*1.2;
				BestCoefficient[n]=coefficient[n];
			}
		}
		//cout <<"Updated step"<<endl;
		TotalStepSqr=0;
		for (n=0;n<Size;n++)
		{
			TotalStepSqr+=step[n]*step[n];
		}
		//cout <<"TotalStepSqr= "<<TotalStepSqr<<"\t"<<"dScore2= "<<dScore*dScore<<"\t"<<"eval= "<<dScore*dScore/TotalStepSqr<<"\t"<<"crit= "<<dScoreCriteria<<endl;
		if (TotalStepSqr==0) break;
		//dScoreSqr=dScore*dScore;
		if (dScore*dScore/TotalStepSqr<dScoreCriteria && numSteps>minSteps)
		{
			//cout <<"About to enter CalcGradient"<<endl;
                        numSteps=0;
			for (n=0;n<Size;n++) coefficient[n]=BestCoefficient[n];
			TotalDerivative=CalcGradient(coefficient, args, func_to_minimize, step);
			NumMinima=CountMinima(coefficient, args, func_to_minimize);
			if (verbose)
                        {
                        cout <<"NumMinima= "<<NumMinima<<endl;
			cout <<"TotalDerivative= "<<TotalDerivative<<endl;
			cout <<"score_old= "<<score_old<<" ConvergedScore= "<<ConvergedScore<<endl;
			cout <<"score_old-ConvergedScore= "<<score_old-ConvergedScore<<endl;
			cout <<endl;
                        }
                        if (TotalDerivative<criteria2) break;
			if (NumMinima==Size) break;
			if (abs(score_old-ConvergedScore)<ConvergedScore*1e-8) break;
			else ConvergedScore=score_old;
		}
	}
		
	for (n=0;n<Size;n++)
	{
		if (verbose) cout <<"coefficient["<<n<<"]= "<<coefficient[n]<<endl;
		coefficient[n]=BestCoefficient[n];
	}
}
/*
Real NegativeFunc(vector<Real> &coefficient, Args args, Real (*func_to_maximize)(vector<Real>, Args))
{
	return -func_to_maximize(coefficient, args);
}

Real MaximizeFunc(vector<Real> &coefficient, Args args, Real (*func_to_maximize)(vector<Real>, Args))
{
	Real (*func_to_minimize)(vector<Real>, Args)=Negav
}
*/
#endif
