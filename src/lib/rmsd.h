#ifndef _rmsd_included_
#define _rmsd_included_

# include <numeric>

# include "AtomUtils.h"
# include "Constants.h"
# include "Structures.h"
# include "TypeDef.h"

Real RMSDNoAlign(vector<AtomStruct> &Atoms1, vector<AtomStruct> &Atoms2)
{
        Real dist, rmsd=0, NumAve=0;
        int i, natom;
        natom=Atoms1.size();
        //if (natom!=Atoms2.size()) 
        //{
        //	cout <<"Error the two molecules do not have the same number of atoms.";
        //	exit(0);
        //}

        for (i=0;i<natom;i++)
        {
                if (Atoms1[i].AtomName=="CA" || Atoms1[i].AtomName=="O1P")
                {
                        dist=AtomSquareDistance(Atoms1[i], Atoms2[i]);
                        rmsd+=dist;
                        NumAve+=1.0;
                }
        }
        rmsd/=NumAve;
        rmsd=sqrt(rmsd);
        return rmsd;
}

Real RMSDNoAlignAll(vector<AtomStruct> &Atoms1, vector<AtomStruct> &Atoms2)
{
        Real dist, rmsd=0, NumAve=0;
        int i, natom;
        natom=Atoms1.size();
        for (i=0;i<natom;i++)
        {
                //cout <<endl;
                //PrintAtomInfo(Atoms1[i]);
                //PrintAtomInfo(Atoms2[i]);
                //cout <<endl;
                dist=AtomSquareDistance(Atoms1[i], Atoms2[i]);
                rmsd+=dist;
                NumAve+=1.0;
        }
        rmsd/=NumAve;
        rmsd=sqrt(rmsd);
        return rmsd;
}

Real EvaluateOrientation(vector<AtomStruct> Atoms1, vector<AtomStruct> Atoms2, Real theta, Real phi, Real psi)
{
        RotateAtoms(Atoms2, theta, phi, psi);
        return RMSDNoAlign(Atoms1, Atoms2);
}

void Rotate(vector<AtomStruct> &Atoms, Real step, int RotationType)
{
        if (RotationType==X) RotateAtoms(Atoms, step, 0, 0);
        else if (RotationType==Y) RotateAtoms(Atoms, 0, step, 0);
        else if (RotationType==Z) RotateAtoms(Atoms, 0, 0, step);
        else 
        {
                cout <<"Error unallowed rotation type "<<RotationType<<endl;
                exit(0);
        }
}

Real CalcDerivative(vector<AtomStruct> Atoms1, vector<AtomStruct> Atoms2, Real step, int RotationType, Real theta, Real phi, Real psi)
{
        Real derivative, dRMSD, rmsd, rmsd_old;
        rmsd_old=EvaluateOrientation(Atoms1, Atoms2, theta, phi, psi);
        if (RotationType==X) theta+=step;
        else if (RotationType==Y) phi+=step;
        else if (RotationType==Z) psi+=step;
        rmsd=EvaluateOrientation(Atoms1, Atoms2, theta, phi, psi);
        if (RotationType==X) theta-=step;
        else if (RotationType==Y) phi-=step;
        else if (RotationType==Z) psi-=step;
        dRMSD=rmsd-rmsd_old;
        //cout <<"RotationType= "<<RotationType<<" rmsd= "<<rmsd<<" rmsd_old= "<<rmsd_old<<" dRMSD= "<<dRMSD<<endl;
        //cout <<"step= "<<step<<endl;
        //cout <<endl;
        if (dRMSD>0)
        {
                step*=-0.5;
                derivative=CalcDerivative(Atoms1, Atoms2, step, RotationType, theta, phi, psi);		
        }
        else derivative=dRMSD/step;
        return derivative;
}

Real CalcGradient(vector<AtomStruct> Atoms1, vector<AtomStruct> Atoms2, Real derivative[3], Real theta, Real phi, Real psi)
{
        int i;
        Real gradient=0, step=0.0001;
        for (i=0;i<3;i++) derivative[i]=CalcDerivative(Atoms1, Atoms2, step, i, theta, phi, psi);
        for (i=0;i<3;i++) gradient+=derivative[i]*derivative[i];
        return gradient;
}

Real GridSearch(vector<AtomStruct> Atoms1, vector<AtomStruct> Atoms2, Real &theta, Real &phi, Real &psi)
{
        Real BestTheta, BestPhi, BestPsi;
        Real BestRmsd, rmsd;
        Real ThetaInc, PhiInc, PsiInc;
        ThetaInc=2.0*pi/5.0;
        PhiInc=2.0*pi/5.0;
        PsiInc=2.0*pi/5.0;
        BestTheta=theta; BestPhi=phi; BestPsi=psi;
        BestRmsd=EvaluateOrientation(Atoms1, Atoms2, theta, phi, psi);
        for (theta=0;theta<2.0*pi;theta+=ThetaInc)
        {
                for (phi=0;phi<2.0*pi;phi+=PhiInc)
                {
                        for (psi=0;psi<2.0*pi;psi+=PsiInc)
                        {
                                rmsd=EvaluateOrientation(Atoms1, Atoms2, theta, phi, psi);
                                //cout <<"rmsd= "<<rmsd<<" psi= "<<psi<<endl;
                                if (rmsd<BestRmsd)
                                {
                                        BestRmsd=rmsd; BestTheta=theta; BestPhi=phi; BestPsi=psi;
                                        //cout <<"rmsd= "<<rmsd<<" theta= "<<theta<<" phi= "<<phi<<" psi= "<<psi<<endl;
                                }
                        }
                }
        }
        theta=BestTheta; phi=BestPhi, psi=BestPsi;
        return BestRmsd;
}

void GetStep(Real derivative[], Real step[], Real TotalStep)
{
        int i;
        Real TotalSqr=0;
        for (i=0;i<3;i++) TotalSqr+=derivative[i]*derivative[i];
        for (i=0;i<3;i++) step[i]=derivative[i]*TotalStep/sqrt(TotalSqr);
}

Real MinimizeAlongVector(vector<AtomStruct> Atoms1, vector<AtomStruct> Atoms2, Real derivative[], Real &theta, Real &phi, Real &psi)
{
        int i;
        Real BestTheta, BestPhi, BestPsi;
        Real dRMSD, rmsd, rmsd_old, TotalStep=0.001;
        Real dRMSDCriteria, GradientCriteria, TotalStepSqr, SqrGradient;
        Real step[3];
        GradientCriteria=1.0e-15; dRMSDCriteria=1.0e-15;
        BestTheta=theta; BestPhi=phi; BestPsi=psi;
        rmsd_old=EvaluateOrientation(Atoms1, Atoms2, theta, phi, psi);
        GetStep(derivative, step, TotalStep);
        while (true)
        {
                theta+=step[X];
                phi+=step[Y];
                psi+=step[Z];
                rmsd=EvaluateOrientation(Atoms1, Atoms2, theta, phi, psi);
                dRMSD=rmsd-rmsd_old;
                if (rmsd>rmsd_old)
                {
                        theta-=step[X];
                        phi-=step[Y];
                        psi-=step[Z];
                        for (i=0;i<3;i++) step[i]=-step[i]*0.5;
                        theta=BestTheta; phi=BestPhi; psi=BestPsi;
                }
                else
                {
                        rmsd_old=rmsd;
                        for (i=0;i<3;i++) step[i]*=1.2;
                        BestTheta=theta; BestPhi=phi; BestPsi=psi;
                }
                TotalStepSqr=0;
                for (i=0;i<3;i++) TotalStepSqr+=step[i]*step[i];
                SqrGradient=dRMSD*dRMSD/TotalStepSqr;
                if (SqrGradient<GradientCriteria || abs(dRMSD)<dRMSDCriteria) break;
        }
        theta=BestTheta; phi=BestPhi; psi=BestPsi;
        rmsd=EvaluateOrientation(Atoms1, Atoms2, BestTheta, BestPhi, BestPsi);
        return rmsd;
}

Real RMSD(vector<AtomStruct> &Atoms1, vector<AtomStruct> &Atoms2)
{
        int i;
        Real dRMSD, ConvergenceCriterion;
        Real gradient, GradientCriterion=10e-15;
        Real rmsd, rmsd_old;
        Real theta, phi, psi;
        Real XOrigin, YOrigin, ZOrigin;
        Real derivative[3];
        vector<AtomStruct> CopyAtoms1, CopyAtoms2;
        CopyVector(Atoms1, CopyAtoms1);
        CopyVector(Atoms2, CopyAtoms2);
        RemoveWaters(CopyAtoms1);
        RemoveWaters(CopyAtoms2);
        theta=phi=psi=0;
        ConvergenceCriterion=10e-10;
        for (i=0;i<3;i++) derivative[i]=0;
        CalcCenter(CopyAtoms2, XOrigin, YOrigin, ZOrigin);
        CalcCenter(CopyAtoms1, XOrigin, YOrigin, ZOrigin);
        center(CopyAtoms1);
        center(CopyAtoms2);
        rmsd_old=GridSearch(CopyAtoms1, CopyAtoms2, theta, phi, psi);
        rmsd=rmsd_old;
        //rmsd_old=RMSDNoAlign(Atoms1, Atoms2);
        while (true)
        {
                gradient=CalcGradient(CopyAtoms1, CopyAtoms2, derivative, theta, phi, psi);
                if (gradient<GradientCriterion) break;
                rmsd=MinimizeAlongVector(CopyAtoms1, CopyAtoms2, derivative, theta, phi, psi);
                dRMSD=rmsd-rmsd_old;
                //cout <<"rmsd= "<<rmsd<<" rmsd_old= "<<rmsd_old<<endl;
                if (-dRMSD<ConvergenceCriterion) break;
                rmsd_old=rmsd;
        }
        center(Atoms2);
        RotateAtoms(Atoms2, theta, phi, psi);
        MoveAtoms(Atoms2, XOrigin, YOrigin, ZOrigin);
        cout <<"rmsd= "<<rmsd<<endl;
        return rmsd;
}

Real RMSD_CA(vector<AtomStruct> Atoms1, vector<AtomStruct> Atoms2)
{
        RemoveNonCA(Atoms1);
        RemoveNonCA(Atoms2);
        cout <<"Atoms1.size()= "<<Atoms1.size()<<endl;
        cout <<"Atoms2.size()= "<<Atoms2.size()<<endl;
        return RMSD(Atoms1, Atoms2);
}
#endif
