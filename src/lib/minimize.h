#ifndef _minimize_included_
#define _minimize_included_

# include "MathUtils.h"

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

vector<Real> solveHessian(vector<Real> &gradient, Matrix &hessian)
{
        int ncoeff=gradient.size();
        vector<Real> step;
        vector< vector<Real> > coefficient;

        SafeAlloc(step, ncoeff, "step");
        Safe2DAlloc(coefficient, ncoeff+1, ncoeff, "coefficient");

        for (int i=0;i<ncoeff;i++)
        {
                //cout <<"coefficient["<<i<<"]["<<ncoeff<<"]= "<<coefficient[i][ncoeff]<<endl;
                coefficient[ncoeff][i]=+gradient[i];
                for (int j=0;j<ncoeff;j++)
                {
                        coefficient[i][j]=hessian[i][j];
                }
        }
        Equation(coefficient, step);
        return step;
}

        template <class Args>
Real CalcDerivative(vector<Real> coefficient, Real &TestStep, int index, Args args, Real (*func_to_minimize)(vector<Real>, Args))
{
        Real derivative, dScore;
        Real score, score_old;
        Real MinTestStep=10e-30;

        if (TestStep==0) TestStep=10e-6;
        score_old=(*func_to_minimize)(coefficient, args);
        //cout <<"coefficient["<<index<<"]= "<<coefficient[index]<<"\t";
        coefficient[index]+=TestStep;
        cout <<"coefficient["<<index<<"]= "<<coefficient[index]<<" TestStep= "<<TestStep<<endl;
        score=(*func_to_minimize)(coefficient, args);
        //cout <<"In CalcDerivative"<<endl;
        //cout <<"score_old= "<<score_old<<"\t"<<"score= "<<score<<endl;
        dScore=score-score_old;
        if (dScore>0)
        {
                coefficient[index]-=TestStep;
                cout <<"coefficient["<<index<<"]= "<<coefficient[index]<<" TestStep= "<<TestStep<<endl;
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
        cout <<"derivative= "<<derivative<<endl;
        return derivative;
}

        template <class Args>
Real CalcGradient(vector<Real> coefficient, Args args, Real (*func_to_minimize)(vector<Real>, Args), vector<Real> &derivative)
{
        int Size, min;
        Real MinTestStep, TestStep;

        Size=coefficient.size();
        //cout <<"In CalcGradient"<<endl;
        //for (n=0;n<Size;n++) cout <<"coefficient["<<n<<"]= "<<coefficient[n]<<endl;
        for (int n=0;n<Size;n++)
        {
                //cout <<"n= "<<n<<endl;
                TestStep=coefficient[n]*10e-6;
                if (TestStep==0) TestStep=10e-6;
                derivative[n]=CalcDerivative(coefficient, TestStep, n, args, func_to_minimize);
                //cout <<"derivative["<<n<<"]= "<<derivative[n]<<endl;
                if (n==0) MinTestStep=abs(TestStep);
                if (abs(TestStep)<MinTestStep) 
                {
                        MinTestStep=abs(TestStep);
                        min=n;
                }
        }
        return 0.1;
}

        template <class Args>
Real calcHessianElement(vector<Real> coefficient, Args args, Real (*func_to_minimize)(vector<Real>, Args), int index1, int index2)
{
        Real derivative1, derivative2;
        Real TestStep1, TestStep2;
        TestStep1=coefficient[index1]*10e-6;
        TestStep2=TestStep1*0.1;
        derivative1=CalcDerivative(coefficient, TestStep1, index1, args, func_to_minimize);
        coefficient[index2]+=TestStep2;
        derivative2=CalcDerivative(coefficient, TestStep1, index1, args, func_to_minimize);
        return (derivative2-derivative1)/TestStep2;
}

        template <class Args>
void calcHessian(vector<Real> coefficient, Args args, Real (*func_to_minimize)(vector<Real>, Args), Matrix &hessian)
{
        int Size=coefficient.size();

        for (int i=0;i<Size;i++)
        {
                for (int j=0;j<Size;j++)
                {
                        hessian[i][j]=calcHessianElement(coefficient, args, func_to_minimize, i, j);
                        hessian[j][i]=hessian[i][j];
                }
        }
}
        
        template <class Args>
void calcStepForConjugateGradient(vector<Real> coefficient, Args args, Real (*func_to_minimize)(vector<Real>, Args), vector<Real> &step)
{
        int Size=step.size();
        vector<Real> derivative;
        //Matrix hessian;
        vector< vector<Real> > hessian;
        SafeAlloc(derivative, Size, "derivative");
        Safe2DAlloc(hessian, Size, Size, "hessian");
        CalcGradient(coefficient, args, func_to_minimize, derivative);
        calcHessian(coefficient, args, func_to_minimize, hessian);
        //PrintMatrix(hessian);
        for (int i=0;i<Size;i++)
        {
                for (int j=0;j<Size;j++)
                {
                        cout <<setiosflags(ios::left)<<setw(10)<<setprecision(5)<<hessian[i][j]<<"\t";
                }
                cout <<endl;
        }
        
        step=solveHessian(derivative, hessian);
}

Real calcVectorMagnitude(vector<Real> &v)
{
        int Size=v.size();
        Real sum=0;

        for (int i=0;i<Size;i++)
        {
                sum+=v[i]*v[i];
        }
        return sqrt(sum);
}
        template <class Args>
Real calcGradientMagnitude(vector<Real> coefficient, Args args, Real (*func_to_minimize)(vector<Real>, Args), vector<Real> &derivative)
{
        CalcGradient(coefficient, args, func_to_minimize, derivative);
        return calcVectorMagnitude(derivative);
}

void scaleVector(vector<Real> &v, Real length)
{
        int Size=v.size();
        Real magnitude;
        magnitude=calcVectorMagnitude(v);
        for (int n=0;n<Size;n++)
        {
                v[n]*=length/magnitude;
        }
}

        template <class Args>
void calcStepForSteepestDescent(vector<Real> coefficient, Args args, Real (*func_to_minimize)(vector<Real>, Args), vector<Real> &step)
{
        Real length;
        length=CalcGradient(coefficient, args, func_to_minimize, step);
        scaleVector(step, length);
}

        template<class Args>
int CountMaxima(vector<Real> coefficient, Args args, Real (*func_to_minimize)(vector<Real>, Args))
{
        int NumMaxima, Size;
        Real TestStep;
        bool Maxima;
        //cout <<"In CountMaxima"<<endl;
        Size=coefficient.size();
        NumMaxima=0;
        for (int n=0;n<Size-1;n++)
        {
                //cout <<"n= "<<n<<endl;
                TestStep=coefficient[n]*0.001;
                Maxima=AtMaxima(coefficient, TestStep, n, args, func_to_minimize);
                if (Maxima) NumMaxima++;
        }
        return NumMaxima;
}

        template<class Args>
int CountMinima(vector<Real> coefficient, Args args, Real (*func_to_minimize)(vector<Real>, Args))
{
        int NumMinima, Size;
        Real TestStep;
        bool Minima;

        Size=coefficient.size();
        NumMinima=0;
        for (int n=0;n<Size;n++)
        {
                TestStep=coefficient[n]*1e-3;
                Minima=AtMinima(coefficient, TestStep, n, args, func_to_minimize);
                if (Minima) NumMinima++;
        }
        return NumMinima;
}

        template<class Args>
Real MinimizeAlongVector(vector<Real> &coefficient, vector<Real> step, Args args, Real (*func_to_minimize)(vector<Real>, Args))
{
        bool verbose=false;
        int Size, numSteps=0, minSteps=10, maxSteps=1000000;
        int skip=1;
        Real dScoreCriteria;
        Real ConvergedScore;
        Real score, score_old, dScore;
        Real TotalStepSqr;
        vector<Real> BestCoefficient;
        if (verbose) cout <<"In MinimizeAlongVector"<<endl;
        dScoreCriteria=1e-15;
        Size=coefficient.size();
        CreateVector(BestCoefficient, Size);
        BestCoefficient=coefficient;
        score_old=(*func_to_minimize)(coefficient, args);
        if (verbose) cout <<"score_old= "<<score_old<<endl;
        ConvergedScore=score_old;	
        while (true)
        {
                for (int n=0;n<Size;n++) coefficient[n]+=step[n];

                //for (int n=0;n<Size;n++) cout <<"step["<<n<<"]= "<<step[n]<<"\t"<<"coefficient= "<<"\t"<<coefficient[n]<<endl;
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
                        for (int n=0;n<Size;n++)
                        {
                                coefficient[n]-=step[n];
                                //cout <<"step["<<n<<"]= "<<step[n]<<" about to decrease step"<<endl;
                                step[n]=-step[n]*0.5;
                                //cout <<"step["<<n<<"]= "<<step[n]<<" just decreased step"<<endl;
                                coefficient[n]=BestCoefficient[n];
                        }
                }
                else
                {
                        score_old=score;
                        for (int n=0;n<Size;n++)
                        {
                                //cout <<"step["<<n<<"]= "<<step[n]<<" about to increase step"<<endl;
                                step[n]=step[n]*1.2;
                                //cout <<"step["<<n<<"]= "<<step[n]<<" just to increased step"<<endl;
                                BestCoefficient[n]=coefficient[n];
                        }
                }
                //cout <<"Updated step"<<endl;
                TotalStepSqr=0;
                for (int n=0;n<Size;n++)
                {
                        TotalStepSqr+=step[n]*step[n];
                }
                //if (verbose) cout <<"TotalStepSqr= "<<TotalStepSqr<<"\t"<<"dScore2= "<<dScore*dScore<<"\t"<<"eval= "<<dScore*dScore/TotalStepSqr<<"\t"<<"crit= "<<dScoreCriteria<<endl;
                if (TotalStepSqr==0) break;
                //dScoreSqr=dScore*dScore;
                if (numSteps%skip==0 && verbose) cout <<"dScore= "<<dScore<<" dScore2/TotalStepSqr= "<<dScore*dScore/TotalStepSqr<<" dScoreCriteria= "<<dScoreCriteria<<" numSteps= "<<numSteps<<endl;
                //if (dScore*dScore/TotalStepSqr<dScoreCriteria && numSteps>minSteps && dScore!=0)
                if (dScore*dScore/TotalStepSqr<dScoreCriteria && numSteps>minSteps)
                {
                        //cout <<"About to enter CalcGradient"<<endl;
                        break;
                }
                if (numSteps>maxSteps) break;
                if (dScore==0) dScoreCriteria*=2.0;
                numSteps++;
        }
        coefficient=BestCoefficient;
        return score_old;
}

        template<class Args>
void IterateOverCoefficients(vector<Real> &coefficient, Args args, Real (*func_to_minimize)(vector<Real>, Args))
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
                TotalGradient=calcGradientMagnitude(coefficient, args, func_to_minimize, derivative);
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

Real calcFletcherReeves(Real f2, Real f1)
{
        return f2*f2/(f1*f1);
}

void updateStep(vector<Real> &derivative, vector<Real> &derivative_old, vector<Real> &step, vector<Real> &step_old)
{
        int nstep=step.size();
        Real gamma;

        for (int i=0;i<nstep;i++)
        {
                if (derivative_old[i]!=0)
                {
                        gamma=calcFletcherReeves(derivative[i], derivative_old[i]);
                        step[i]=derivative[i]+gamma*step_old[i];
                }
                else step[i]=derivative[i];
        }
}

        template<class Args>
Real nonlinearConjugateGradient(vector<Real> &coefficient, Args &args, Real (*func_to_minimize)(vector<Real>, Args), Real dScoreCriteria)
{
        //Algorithm from en.wikipedia.org/wiki/Energy_minimization
        int Size=coefficient.size();
        Real dScore, score, score_old, length;
        vector <Real> step, step_old, derivative, derivative_old;
        score_old=(*func_to_minimize)(coefficient, args);
        SafeAlloc(step, Size, "step");
        SafeAlloc(derivative, Size, "derivative");
        SafeAlloc(step_old, Size, "step_old");
        CalcGradient(coefficient, args, func_to_minimize, derivative);
        step=derivative;
        while (true)
        {
                score=MinimizeAlongVector(coefficient, step, args, func_to_minimize);
                //cout <<"score= "<<score<<endl;
                dScore=score-score_old;
                if (abs(dScore)<abs(score)*dScoreCriteria && score!=0)
                {
                        ZeroVector(step_old);
                        CalcGradient(coefficient, args, func_to_minimize, derivative);
                        step=derivative;
                        score=MinimizeAlongVector(coefficient, step, args, func_to_minimize);
                }
                dScore=score-score_old;
                score_old=score;
                if (abs(dScore)<abs(score)*dScoreCriteria && score!=0) break;
                else if (abs(dScore)<dScoreCriteria) break;
                derivative_old=derivative;
                length=CalcGradient(coefficient, args, func_to_minimize, derivative);
                updateStep(derivative, derivative_old, step, step_old);
                step_old=step;
                scaleVector(step, length);
        }
        return score;
}
        
        template<class Args>
Real conjugateGradient(vector<Real> &coefficient, Args &args, Real (*func_to_minimize)(vector<Real>, Args), Real dScoreCriteria)
{
        bool verbose=false;
        int minSteps=10, NumMinima, numSteps=0, Size;
        Real criteria, criteria2;
        Real ConvergedScore;
        Real score, score_old, dScore;
        Real TotalDerivative;
        Real TestStep, Total, TotalStepSqr;
        vector<Real> step;

        criteria=0.00001;
        criteria2=1e-5;
        TestStep=0.001;
        Size=coefficient.size();
        SafeAlloc(step, Size, "step");
        score_old=(*func_to_minimize)(coefficient, args);
        if (verbose) cout <<"score_old= "<<score_old<<endl;
        calcStepForConjugateGradient(coefficient, args, func_to_minimize, step);
        if (verbose) cout <<"Done with CalcGradient"<<endl;
        while (true)
        {

                score=MinimizeAlongVector(coefficient, step, args, func_to_minimize);
                if (verbose) cout <<"score= "<<score<<endl;
                dScore=score-score_old;
                score_old=score;
                if (abs(dScore)<abs(score)*dScoreCriteria && score!=0) break;
                else if (abs(dScore)<dScoreCriteria) break;
                calcStepForConjugateGradient(coefficient, args, func_to_minimize, step);
        }

        for (int n=0;n<Size;n++)
        {
                if (verbose) cout <<"coefficient["<<n<<"]= "<<coefficient[n]<<endl;
        }
        return score;
}

        template<class Args>
Real OptimizeCoefficients(vector<Real> &coefficient, Args &args, Real (*func_to_minimize)(vector<Real>, Args), Real dScoreCriteria)
{
        bool verbose=true;
        int minSteps=10, NumMinima, numSteps=0, Size;
        Real criteria, criteria2;
        Real ConvergedScore;
        Real score, score_old, dScore;
        Real TotalDerivative;
        Real TestStep, Total, TotalStepSqr;
        vector<Real> step;

        criteria=0.00001;
        criteria2=1e-5;
        TestStep=0.001;
        Size=coefficient.size();
        SafeAlloc(step, Size, "step");
        score_old=(*func_to_minimize)(coefficient, args);
        if (verbose) cout <<"score_old= "<<score_old<<endl;
        calcStepForSteepestDescent(coefficient, args, func_to_minimize, step);
        if (verbose) cout <<"Done with CalcGradient"<<endl;
        while (true)
        {

                score=MinimizeAlongVector(coefficient, step, args, func_to_minimize);
                cout <<"score= "<<score<<endl;
                dScore=score-score_old;
                score_old=score;
                if (abs(dScore)<abs(score)*dScoreCriteria && score!=0) break;
                else if (abs(dScore)<dScoreCriteria) break;
                calcStepForSteepestDescent(coefficient, args, func_to_minimize, step);
        }

        for (int n=0;n<Size;n++)
        {
                if (verbose) cout <<"coefficient["<<n<<"]= "<<coefficient[n]<<endl;
        }
        return score;
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
