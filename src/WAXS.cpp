/*******************************************************************************

  WAXS.cpp
  -------

  Copyright (C) 2007 The University of Chicago

Authors:
Jouko Virtanen

Description:
Calculates wide angle x-ray scattering given the cartesian coordinates of a
macromolecule.


Changes from version 91:  FFT for cube method implemented.  Working on
representing protein surface with spherical harmonics.

Changes from version 97: Fixed problems with solventtypes and cube
subroutines.  Working on calculating p(r) from scattering and atomic 
coordinates.

Changes from version 134: Added structures.
Changes from version 141: Made sure Debye is working.
Changes from version 149: Make sure water sphere is working.
Changes from version 153: Made FFT more memory efficient.
Changes from version 156: Added SWANS
 *******************************************************************************/


# include <iostream>
# include <algorithm>
# include <stdio.h>
# include <stdlib.h>
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

# include "lib/StringUtils.h"
# include "lib/LinkedList.h"
# include "WAXS.h"
# include "lib/AssignAtomIDs.h"
# include "lib/AtomCode.h"
# include "lib/AtomIDs.h"
# include "lib/AtomUtils.h"
# include "lib/Constants.h"
# include "lib/ConvertStructures.h"
# include "lib/GetAtomType.h"
# include "lib/GetNextFrame.h"
# include "lib/HydrationShell_fast4.h"
# include "lib/IOUtils.h"
# include "lib/MathUtils.h"
# include "lib/normalize.h"
# include "lib/PointsOnASphere.h"
# include "lib/ReadPdb.h"
# include "lib/ReadDcd.h"
# include "lib/ReadIntensityFile.h"
# include "lib/ResidueCode.h"
# include "lib/Structures.h"
# include "lib/TypeDef.h"
# include "lib/VectorManip.h"
# include "lib/WritePdb.h"
# include "ReadParameterFile.h"

bool verbose=false;
const string HomeDir="/home/jouko/project/WAXS/";
const string Version="WAXS 1.0";

const int HYDRATION_OPT=0;
const int EXCLUDED_VOLUME_SCALE_OPT=1;

vector < vector< vector<Real> > > pr;
Real solventatomr[NumAtomTypes];

bool CommandLine;

inline Real CalcGaussianSphereExcludedElectrons(Real density, Real atomr)
{
        //Calculates how many electrons a Gaussina excluded volume dummy atom 
        //has
        cout <<"GaussianSphereExcludedElectrons= "<<density*atomr*atomr*atomr*pi32;
        return density*atomr*atomr*atomr*pi32;
}

inline Real CalcHardSphereExcludedElectrons(Real density, Real atomr)
{
        //Calculates how many electrons a hard sphere excluded volume dummy atom
        //has.
        cout <<"HardSphereExcludedElectrons= "<<density*atomr*atomr*atomr*pi*4.0/3.0;
        return density*atomr*atomr*atomr*pi*4.0/3.0;
}

inline Real calcElectrons(Real contrast, Real atomr, ParamStruct &params)
{
        //The bulk solvent density needs to be subtracted from every
        //point in the protein.  This is done with excluded volume
        //dummy atoms.  This function specifies the number of electrons
        //an excluded volume dummy atom has given the bulk solvent density,
        //the radius of the atom, and the type of the atom.
        if (params.ExcludedVolumeSphereType == "GaussianSphere") 
        {
                return CalcGaussianSphereExcludedElectrons(contrast, atomr);
        }
        else if (params.ExcludedVolumeSphereType == "HardSphere") 
        {
                return CalcHardSphereExcludedElectrons(contrast, atomr);
        }
        else 
        {
                cout <<"Fatal Error.  Unknown ExcludedVolumeSphereType "
                        <<params.ExcludedVolumeSphereType<<endl;
                exit(EXIT_FAILURE);
        }
        return UNK_REAL;
}

Real calcExcludedElectrons(Real density, Real atomr, ParamStruct &params)
{
        return calcElectrons(density, atomr, params);
}

void CalculateExcludedElectrons(Real ExcludedElectrons[], Real density, Real atomr[], ParamStruct &params)
{
        //Calculates the number of electrons all excluded volume dummy atoms 
        //have.
        for (int n=0;n<NumAtomTypes;n++)
        {
                ExcludedElectrons[n]=calcElectrons(density, atomr[n], params);
                cout <<"density= "<<density<<endl;
                cout <<"atomr["<<n<<"]= "<<atomr[n]<<endl;
                cout <<"ExcludedElectrons["<<n<<"]= "<<ExcludedElectrons[n]<<endl;
        }
}

void CalculateSolventCorrectedElectrons(Real SolventCorrectedElectrons[], Real NumElectrons[], Real ExcludedElectrons[], Real atomr[], ParamStruct &params)
{
        //Subtracts the excluded volume electrons from the real electrons of 
        //an atom when the supperimpose method is used.
        for (int n=0;n<NumAtomTypes;n++)
        {
                SolventCorrectedElectrons[n]=NumElectrons[n]-ExcludedElectrons[n];
                cout <<"SolventCorrectedElectrons["<<n<<"]= "<<SolventCorrectedElectrons[n]<<endl;
                cout <<"NumElectrons["<<n<<"]= "<<NumElectrons[n]<<endl;
                cout <<"ExcludedElectrons["<<n<<"]= "<<ExcludedElectrons[n]<<endl;
        }
        SolventCorrectedElectrons[HydrationShell]=calcElectrons(params.hsdensity, atomr[HydrationShell], params);
}

inline void CalculateSolventCorrectedElectrons(Real SolventCorrectedElectrons[], Real contrast, Real atomr[], ParamStruct &params)
{
        //The bulk solvent density is subtracted from every point in the 
        //protein.  This can be accomplished with dummy excluded volume
        //atoms, which have negative electron density supperimposed on real 
        //atoms.  This function calculates the total number of electrons
        //that each atom type has when the bulk density is subtracted out.
        Real NumElectrons[NumAtomTypes], ExcludedElectrons[NumAtomTypes];
        SetNumElectrons(NumElectrons);
        for (int i=0;i<NumAtomTypes;i++) ExcludedElectrons[i]=0;
        CalculateExcludedElectrons(ExcludedElectrons, contrast, atomr, params);
        CalculateSolventCorrectedElectrons(SolventCorrectedElectrons, NumElectrons, ExcludedElectrons, atomr, params);
}

int ReadPdbCube(vector<CubeStruct> &cubes, ParamStruct &params)
{
        //Reads in a pdb and saves the info in cube format.  The BFactors are
        //the densities.  This is used when the cube method is used and 
        //AssignPdb is set to FromPdb.
        string line;
        fstream protein;
        CubeStruct cube;

        InitializeCube(cube);
        OpenFile(params.InputCubeDensityPdbFile, protein, "InputCubeDensityPdbFile");

        while (true)
        {
                //cout <<line<<endl;
                if (line.substr(0,4)=="CRYS")
                {
                        params.XBoxLength=StrToFloat(GetWord(line,2));
                        params.YBoxLength=StrToFloat(GetWord(line,3));
                        params.ZBoxLength=StrToFloat(GetWord(line,4));
                }
                if (protein.eof()) break;
                if (line.substr(0,3)=="END") break;

                if (line.substr(0,4)=="HETA" || line.substr(0,4)=="ATOM")
                {

                        cube.x=(StrToFloat(line.substr(31,7)));
                        cube.y=(StrToFloat(line.substr(39,7)));
                        cube.z=(StrToFloat(line.substr(47,7)));
                        cube.density=(StrToFloat(line.substr(61,10)));
                        SafePushBack(cubes, cube, "cubes");
                }

                getline(protein, line);
        }

        return cubes.size();
}

int ReadPdbCube(Real contrast, vector<CubeStruct> &cubes, ParamStruct &params)
{
        //Reads a pdb, saves the info in cube format and subtracts out the 
        //bulk solvent density.
        int CubeNum;

        CubeNum=ReadPdbCube(cubes, params);
        SubtractBulkDensity(cubes, contrast);

        return CubeNum;
}

void GetTinkerCharges(string keyfile, vector<AtomStruct> &Atoms)
{
        //Make sure this works
        //This reads in charges from a Tinker parameter file and uses them
        //to calculate scattering factors for the atoms.  This is probably
        //broken.
        string line;
        string str;
        vector<string> vstr;
        int n, natom=Atoms.size();
        Real *temp;

        SafeArrayAlloc(temp, natom, "natom");

        fstream charges;
        OpenFile(keyfile, charges, "Tinker keyfile");

        getline(charges, line);

        while (true)
        {
                str=line.substr(0,6);

                if (str=="charge")
                {
                        Tokenize2(line, " ", vstr);
                        temp[n]=StrToFloat(vstr[2]);
                        n++;
                }
        }

        for (n=0;n<natom;n++)
        {
                //Atoms[n].charge=(atomm[Atoms[n].atomid]-temp[Atoms[n].atomtype])/atomm[Atoms[n].atomid];
        }
}

//bool StrToBool2(string str)
//{
//	return ( str=="yes" || str == "y" || str == "true" || str == "t" );
//}

void cart2sphere(Real r[], Real phi[], Real stheta[], Real ctheta[], vector<AtomStruct> &Atoms)
{
        //This converts the cartesian coordinates of the atoms to an array
        //of polar coordinates.  This is used for multipole expansion.
        int natom;
        Real invr;
        Real xn, yn, zn;
        natom=Atoms.size();
        for (int n=0;n<natom;n++)
        {
                xn=Atoms[n].x;
                yn=Atoms[n].y;
                zn=Atoms[n].z;
                cout <<"xn= "<<xn<<endl;
                cout <<"yn= "<<yn<<endl;
                cout <<"zn= "<<zn<<endl;
                r[n]=sqrt(xn*xn+yn*yn+zn*zn);
                invr=1.0/r[n];
                ctheta[n]=(zn*invr);
                stheta[n]=(sqrt(xn*xn+yn*yn)*invr);
                phi[n]=( asin(yn/sqrt(xn*xn+yn*yn)) );

                if (xn<0) phi[n]=pi-phi[n];

                if (xn>0 && yn<0) phi[n]+=2.0*pi;
        }
}

Real bessel2(long double r, int l)
{
        //Calculates the bessel function using a recursive formula.  This is 
        //used for multipole expansion.  This can fail especially for small r.
        int i;
        vector<long double> j;
        long double invr;
        long double longi;

        j.resize(l+2);

        invr=1.0/r;
        j[0]=sin(r)*invr;
        j[1]=(j[0]-cos(r))*invr;

        longi=3.0;
        for (i=2;i<l+1;i++)
        {
                j[i]=j[i-1]*longi*invr-j[i-2];
                longi=longi+2.0;
        }

        if (abs(j[l])>1)
        {
                j[l]=0;
                cout <<"Bessel overflow\n";
        }

        return j[l];
}

Real bessel(long double r, int l, long double sterm[])
{
        //Calculates the bessel function using Taylor expansion.  It can give
        //bad answers for large r.  So sometimes it switches to bessel2.
        int m;
        long double total;
        long double term;
        long double n;
        long double s;
        long double r2;

        term=1.0;
        total=0.0;
        r2=r*r;
        m=1;

        for (n=1.0;n<=l;n+=1.0)
        {
                term=term*r/(2.0*n+1.0);
        }

        s=1.0;

        while (abs(term)>0.0000001)
        {
                total+=term;
                term=-term*sterm[m]*r2;
                s+=1.0;
                m++;
        }

        if (abs(total)>1) total=bessel2(r, l);

        return total;

}

void PrintAtomr(Real atomr[])
{
        //Prints the excluded volume radii.
        for (int n=0;n<NumAtomTypes;n++)
        {
                cout <<"aromr["<<n<<"]= "<<atomr[n]<<endl;
        }
}

void FindProteinCenter(Real atomm[], Real contrast, Real atomr[], vector<AtomStruct> &Atoms, ParamStruct &params, Real &avex, Real &avey, Real &avez)
{
        //Finds the center of mass of the protein corrected for excluded volume.
        int AtomID, natom=Atoms.size();
        Real xsum, ysum, zsum;
        Real masssum;
        Real atommass;
        avex=0, avey=0, avez=0;

        xsum=0;
        ysum=0;
        zsum=0;
        masssum=0;
        PrintAtomr(atomr);
        for (int n=0;n<natom;n++)
        {
                if (IsProtein(Atoms[n]) || isNA(Atoms[n]) || Atoms[n].atomid==HydrationShell || Atoms[n].atomid==ExcludedVolume)
                {
                        atommass=atomm[Atoms[n].atomid]-calcExcludedElectrons(contrast, atomr[Atoms[n].atomid], params);
                        atommass*=Atoms[n].weight;
                        //cout <<"atommass= "<<atommass<<endl;
                        xsum+=Atoms[n].x*atommass;
                        ysum+=Atoms[n].y*atommass;
                        zsum+=Atoms[n].z*atommass;
                        cout <<"xsum= "<<xsum<<" x= "<<Atoms[n].x<<" Electrons= "<<atommass<<endl;
                        masssum+=atommass;
                }
        }
        cout <<"xsum= "<<xsum<<endl;
        cout <<"avex= "<<avex<<" avey= "<<avey<<" avez= "<<avez<<endl;
        if (masssum!=0)
        {
                avex=xsum/masssum;
                avey=ysum/masssum;
                avez=zsum/masssum;
                cout <<"avex= "<<avex<<" avey= "<<avey<<" avez= "<<avez<<endl;
        }
        else
        {
                cout <<"ERROR: No protein atoms.  Unable to find center of protein."<<endl;
                avex=0;
                avey=0;
                avez=0;
        }

}

void center(Real atomm[], Real density, Real atomr[], vector<AtomStruct> &Atoms, ParamStruct &params)
{
        Real avex, avey, avez;
        //Moves the solvent corrected center of mass of the protein to the
        //origin.
        FindProteinCenter(atomm, density, atomr, Atoms, params, avex, avey, avez);
        MoveAtoms(Atoms, -avex, -avey, -avez);
}

void DistanceMatrix(vector<AtomStruct> Atoms)
{
        //Write
}

void findradii(int natom, Real numatom[], Real atomr[], vector<AtomStruct> Atoms)
{
        //Finds the appropriate excluded volume of the atoms from the protein
        //coordinates.  It finds all distances between neighboring atoms, 
        //and bases the radii on those distances.  The neighbors are determined
        //by the radii so this is an iterative proceedure.
        int nearest;
        int atomidn;
        Real mindist;
        Real dist;
        Real dist2;
        Real xn, yn, zn;
        Real dx, dy, dz;
        Real atomr2[NumAtomTypes];
        Real *avevolume;
        Real neighbordist[NumAtomTypes];
        Real neighbortype[NumAtomTypes][NumAtomTypes];

        SafeArrayAlloc(avevolume, natom, "avevolume");
        mindist=0;

        for (int n=0;n<NumAtomTypes;n++)
        {
                atomr2[n]=0;
                neighbordist[n]=0;
                for (int m=0;m<NumAtomTypes;m++)
                {
                        neighbortype[n][m]=0;
                }
        }

        for (int n=0;n<natom;n++) avevolume[n]=0;
        DistanceMatrix(Atoms);

        for (int n=0;n<natom;n++)
        {
                xn=Atoms[n].x;
                yn=Atoms[n].y;
                zn=Atoms[n].z;
                atomidn=Atoms[n].atomid;
                mindist=100000000;

                for (int m=0;m<natom;m++)
                {
                        dx=xn-Atoms[m].x;
                        dy=yn-Atoms[m].y;
                        dz=zn-Atoms[m].z;

                        if (m!=n)
                        {
                                dist=dx*dx+dy*dy+dz*dz;
                                if (dist<mindist) mindist=dist;
                        }
                }
                mindist=sqrt(mindist);
                avevolume[atomidn]=avevolume[atomidn]+mindist;
        }

        for (int n=0;n<NumAtomTypes;n++)
        {
                if (numatom[n]>0) atomr[n]=avevolume[n]*0.5/numatom[n];
        }

        for (int n=0;n<natom;n++)
        {
                xn=Atoms[n].x;
                yn=Atoms[n].y;
                zn=Atoms[n].z;
                atomidn=Atoms[n].atomid;
                mindist=1000;
                for (int m=0;m<natom;m++)
                {
                        if (m!=n)
                        {
                                dx=xn-Atoms[m].x;
                                dy=yn-Atoms[m].y;
                                dz=zn-Atoms[m].z;
                                dist=dx*dx+dy*dy+dz*dz;
                                if (dist<100)
                                {
                                        dist=sqrt(dist)-atomr[Atoms[m].atomid];
                                        if (dist<mindist)
                                        {
                                                dist2=dist+atomr[Atoms[m].atomid];
                                                mindist=dist;
                                                nearest=Atoms[m].atomid;
                                        }
                                }
                        }
                }
                neighbordist[atomidn]+=dist2;
                neighbortype[atomidn][nearest]+=1;
        }

        for (int n=0;n<6;n++)
        {
                atomr2[n]=atomr[n];
        }

        for (int k=0;k<1000;k++)
        {
                for (int n=0;n<6;n++)
                {
                        atomr[n]=neighbordist[n];
                        for (int m=0;m<6;m++)
                        {
                                cout <<"atomr["<<n<<"]="<<atomr[n]<<" neighbortype["<<n<<"]["<<m<<"]="<<neighbortype[n][m]<<" atomr2["<<m<<"]="<<atomr2[m]<<endl;
                                atomr[n]=atomr[n]-neighbortype[n][m]*atomr2[m];
                        }

                        if (numatom[n]>0)
                        {
                                atomr[n]=atomr[n]/numatom[n];
                                atomr[n]=(atomr[n]+atomr2[n])*0.5;
                                atomr2[n]=atomr[n];
                        }
                }
        }
}

void FindRadii2(int natom, Real atomr[], Real numatom[], vector<AtomStruct> &Atoms)
{
        //Optimizes the excluded volume radii by minimizing overlaps between
        //them.  This might not be working.
        int fit[50];
        Real FirstOverlap;
        Real Overlap1, Overlap2, Overlap3;
        Real total;
        Real volume;
        Real TotalVolume;
        Real fitatomr[100];
        Real solventatomr[100];
        Real PreviousAtomr[100];
        Real step[100];
        Real OverlapVolume[100];

        for (int t=0;t<NumAtomTypes;t++)
        {
                fitatomr[t]=0;
                solventatomr[t]=0;
        }

        for (int t=0;t<50;t++)
        {
                OverlapVolume[t]=0;
                fit[t]=0;
                PreviousAtomr[t]=0;
                step[t]=0;
        }

        for (int t=0;t<50;t++)
        {
                step[t]=0.0001;
        }

        for (int t=0;t<6;t++)
        {
                fit[t]=1;
        }

        for (int n=0;n<6;n++)
        {
                cout <<"atomr["<<n<<"]="<<atomr[n]<<endl;
        }

        //findradii(natom, numatom, atomr);

        atomr[HYDROGEN]=0.902255;
        atomr[CARBON]=1.43693;
        atomr[NITROGEN]=1.24527;
        atomr[OXYGEN]=1.22099;
        atomr[SULFUR]=2.19596;

        volume=findvolume(natom, atomr, numatom, Atoms);
        cout <<"volume="<<volume<<endl;

        for (int n=0;n<6;n++)
        {
                cout <<"atomr["<<n<<"]="<<atomr[n]<<endl;
        }

        for (int u=0;u<10;u++)
        {
                FirstOverlap=Overlap(natom, atomr, Atoms);

                cout <<"FirstOverlap="<<FirstOverlap<<endl;

                for (int t=0;t<50;t++)
                {
                        step[t]=0.00001;
                }

                for (int t=0;t<6;t++)
                {
                        if (fit[t]==1)
                        {
                                for (int n=0;n<6;n++)
                                {
                                        PreviousAtomr[n]=atomr[n];
                                }

                                atomr[t]=atomr[t]+step[t];

                                TotalVolume=0;

                                for (int n=0;n<6;n++)
                                {
                                        TotalVolume+=numatom[n]*4.0*pi*atomr[n]*atomr[n]*atomr[n]/3.0;
                                }

                                for (int n=0;n<6;n++)
                                {
                                        if (TotalVolume!=0) atomr[n]=atomr[n]*exp(log(volume/TotalVolume)/3.0);
                                }

                                for (int n=0;n<6;n++)
                                {
                                        cout <<"atomr["<<n<<"]="<<atomr[n]<<endl;
                                }

                                OverlapVolume[t]=Overlap(natom, atomr, Atoms);

                                for (int n=0;n<6;n++)
                                {
                                        atomr[n]=PreviousAtomr[n];
                                }

                                cout <<OverlapVolume[t]<<"\t"<<t<<endl;
                        }
                }

                total=0;

                for (int t=0;t<6;t++)
                {
                        if (fit[t]==1) total += (FirstOverlap-OverlapVolume[t])*(FirstOverlap-OverlapVolume[t])/step[t];
                }

                for (int t=0;t<6;t++)
                {
                        if (fit[t]==1) step[t]=(FirstOverlap-OverlapVolume[t])*0.00001/(sqrt(total*step[t]));
                }

                Overlap1=FirstOverlap;
                Overlap2=0;
                Overlap3=0;

                while (true)
                {
                        for (int n=0;n<6;n++)
                        {
                                PreviousAtomr[n]=atomr[n];
                        }

                        for (int t=0;t<6;t++)
                        {
                                if (fit[t]==1) atomr[t] += step[t];
                        }

                        TotalVolume=0;

                        for (int n=0;n<6;n++)
                        {
                                TotalVolume+=numatom[n]*4.0*pi*atomr[n]*atomr[n]*atomr[n]/3.0;
                        }

                        for (int n=0;n<6;n++)
                        {
                                atomr[n]=atomr[n]*exp(log(volume/TotalVolume)/3.0);
                        }

                        Overlap2=Overlap1;
                        Overlap1=Overlap(natom, atomr, Atoms);

                        cout <<Overlap1<<"\t"<<atomr[0]<<"\t"<<atomr[1]<<"\t"<<atomr[2]<<"\t"<<atomr[3]<<"\t"<<atomr[4]<<endl;

                        if (Overlap2<Overlap1)
                        {
                                for (int t=0;t<6;t++)
                                {
                                        atomr[t]=PreviousAtomr[t];
                                        step[t]=-step[t]*0.5;
                                }
                        }


                        if (Overlap1<Overlap2)
                        {
                                for (int t=0;t<6;t++)
                                {
                                        step[t]=step[t]*1.2;
                                }
                        }

                        if ((Overlap2-Overlap1)*(Overlap2-Overlap1)<0.00001) break;
                }
        }

}

void FitQuadratic(Real x1, Real y1, Real x2, Real y2, Real x3, Real y3, Real &a, Real &b, Real &c)
{
        //Given three x y pairs finds the coefficients of the equation 
        //ax^2+bx+c=y which passes through all of the points.
        a=( (y1-y2)*(x2-x3)-(y2-y3)*(x1-x2) )/( (x1*x1-x2*x2)*(x2-x3)-(x2*x2-x3*x3)*(x1-x2) );
        b=( y1-y2-a*(x1*x2-x2*x2))/(x1-x2);
        c=y1-a*x1*x1-b*x1;
}

void FindRadii3(int natom, Real atomr[], Real numatom[], vector<AtomStruct> &Atoms)
{
        //Optimizes the radii of the excluded volume dummy atoms  by 
        //minimizing the overalap between them.  This is probably not any better
        //than FindRadii2, and is overly complicated.  
        //Should just use minimize.h.
        bool parabolic;
        int fit[50];
        int u;
        Real Overlap1, Overlap2, Overlap3;
        Real volume, TotalVolume;
        Real a, b, c;
        Real BestOverlap;
        Real temp1, temp2, temp3;
        Real temppar1, temppar2, temppar3;
        Real fitatomr[100];
        Real fitatomr2[NumAtomTypes];
        Real solventatomr[100];
        Real parameter[100], parameter2[100], parameter3[100];
        Real BestAtomr[100];
        Real fitparameter[100];
        Real step[100];


        volume=findvolume(natom, atomr, numatom, Atoms);
        cout <<"volume="<<volume<<endl;


        for (int t=0;t<NumAtomTypes;t++)
        {
                fitatomr[t]=0;
                fitatomr2[t]=0;
                solventatomr[t]=0;
        }

        for (int t=0;t<50;t++)
        {
                parameter[t]=0;
                parameter2[t]=0;
                parameter3[t]=0;
                fitparameter[t]=0;
                BestAtomr[t]=0;
                step[t]=0;
        }

        for (int t=0;t<50;t++)
        {
                step[t]=0.0001;
                fit[t]=0;
        }

        for (int n=0;n<5;n++)
        {
                fit[n]=1;
        }

        for (int n=0;n<6;n++)
        {
                parameter[n]=atomr[n];
                BestAtomr[n]=atomr[n];
                cout <<atomr[n]<<endl;
        }

        u=0;
        Overlap1=0;
        BestOverlap=1000000;
        while (u==0)
        {
                for (int t=0;t<8;t++)
                {
                        step[t]=0.001;
                        if (Overlap1>BestOverlap) u=1;
                        Overlap1=0;
                        Overlap2=0;
                        Overlap3=0;

                        //cout <<"Optimizing overlap"<<endl;

                        if (fit[t]!=1)
                        {
                                step[t]=0;
                        }

                        parabolic=false;

                        for (int n=0;n<6;n++)
                        {
                                atomr[n]=BestAtomr[n];
                        }


                        while (step[t]*step[t]>0.00000001 && atomr[t]>0 && atomr[t]<2.4)
                        {

                                if (!parabolic)
                                {
                                        Overlap3=Overlap2;
                                        Overlap2=Overlap1;
                                        Overlap1=Overlap(natom, atomr, Atoms);
                                }

                                if (parabolic)
                                {
                                        temp1=Overlap(natom, atomr, Atoms);

                                        if (temp1==Overlap1 || temp1==Overlap2 || temp1==Overlap3) break;

                                        if (temp1<Overlap1)
                                        {
                                                Overlap3=Overlap2;
                                                Overlap2=Overlap1;
                                                Overlap1=temp1;
                                                parameter3[t]=parameter2[t];
                                                parameter2[t]=fitparameter[t];
                                                fitparameter[t]=atomr[t];
                                        }

                                        if (temp1>Overlap1 && temp1<Overlap2)
                                        {
                                                Overlap3=Overlap2;
                                                Overlap2=temp1;
                                                parameter3[t]=parameter2[t];
                                                parameter2[t]=atomr[t];
                                        }

                                        if (temp1>Overlap2 && temp1<Overlap3)
                                        {
                                                Overlap3=temp1;
                                                parameter3[t]=atomr[t];
                                        }

                                        atomr[t]=fitparameter[t];

                                        //cout <<"parameter["<<t<<"]="<<parameter[t]<<" fitparameter["<<t<<"]="<<fitparameter[t]<<endl;
                                }

                                if (Overlap1<BestOverlap)
                                {
                                        BestOverlap=Overlap1;
                                        for (int n=0;n<6;n++)
                                        {
                                                BestAtomr[n]=atomr[n];
                                        }
                                }

                                cout <<setw(9)<<left<<"Overlap"<<setw(9)<<Overlap1<<setw(9)<<left<<t<<setw(12)<<atomr[t]<<endl;
                                //cout <<"step["<<t<<"]="<<step[t]<<" chisqr1="<<chisqr1<<" chisqr2="<<chisqr2<<" chisqr3="<<chisqr3<<endl;

                                if (Overlap3>Overlap2 && Overlap1>Overlap2 && !parabolic)
                                {
                                        fitparameter[t]=parameter2[t];
                                        temp1=Overlap1;
                                        temp2=Overlap2;
                                        temp3=Overlap3;
                                        temppar1=atomr[t];
                                        temppar2=parameter2[t];
                                        temppar3=parameter3[t];
                                        Overlap1=temp2;

                                        if (temp1<temp3)
                                        {
                                                Overlap2=temp1;
                                                Overlap3=temp3;
                                                parameter2[t]=temppar1;
                                                parameter3[t]=temppar3;
                                        }
                                        else
                                        {
                                                Overlap2=temp3;
                                                Overlap3=temp1;
                                                parameter2[t]=temppar3;
                                                parameter3[t]=temppar1;
                                        }

                                        parabolic=true;
                                        atomr[t]=fitparameter[t];
                                        cout <<endl;
                                }

                                if (parabolic)
                                {
                                        FitQuadratic(atomr[t], Overlap1, parameter2[t], Overlap2, parameter3[t], Overlap3, a, b, c);
                                        atomr[t]=-b/(2*a);
                                        //cout <<"Assign parameter["<<t<<"]="<<parameter[t]<<" fitparameter["<<t<<"]="<<fitparameter[t]<<endl;
                                        step[t]=fitparameter[t]-atomr[t];
                                        if (temp1>Overlap3)
                                        {
                                                step[t]=0;
                                                atomr[t]=fitparameter[t];
                                        }
                                }


                                if (!parabolic)
                                {
                                        parameter3[t]=parameter2[t];
                                        parameter2[t]=atomr[t];

                                        if (Overlap2>Overlap1) step[t]=step[t]*1.2;
                                        if (Overlap2<Overlap1) step[t]=-step[t]*0.5;
                                        atomr[t]=atomr[t]+step[t];
                                }

                                if (atomr[t]<0) atomr[t]=0;

                                TotalVolume=0;

                                for (int n=0;n<6;n++)
                                {
                                        TotalVolume+=numatom[n]*4.0*pi*atomr[n]*atomr[n]*atomr[n]/3.0;
                                }

                                for (int n=0;n<6;n++)
                                {
                                        if (TotalVolume!=0)
                                        {
                                                atomr[n]=atomr[n]*exp(log(volume/TotalVolume)/3.0);
                                        }
                                }
                        }
                }
        }

}

Real Overlap(int TotalParticles, Real atomr[], vector<AtomStruct> &Atoms)
{
        //Calculates the total overlap volume of a collection of spheres.
        Real dx, dy, dz;
        Real xm, ym, zm;
        Real atomr1, atomr2;
        Real costheta1, costheta2;
        Real d, r1;
        Real height1, height2;
        Real cone, section;
        Real DVolume, OverlapVolume;

        OverlapVolume=0;

        for (int m=0;m<TotalParticles-1;m++)
        {
                xm=Atoms[m].x;
                ym=Atoms[m].y;
                zm=Atoms[m].z;
                atomr1=atomr[Atoms[m].atomid];

                for (int n=m+1;n<TotalParticles;n++)
                {
                        dx=Atoms[m].x-xm;
                        dy=Atoms[m].y-ym;
                        dz=Atoms[m].z-zm;
                        atomr2=atomr[Atoms[n].atomid];
                        d=sqrt(dx*dx+dy*dy+dz*dz);

                        if (d < (atomr1+atomr2) && ((d>atomr1 && d> atomr2) || (d<atomr1 && d< atomr2)) )
                        {
                                costheta1=(atomr1*atomr1+d*d-atomr2*atomr2)/(2.0*atomr1*d);
                                height1=atomr1*costheta1;

                                costheta2=(atomr2*atomr2+d*d-atomr1*atomr1)/(2.0*atomr2*d);
                                height2=atomr2*costheta2;

                                if (atomr1>height1)
                                {
                                        r1=sqrt(atomr1*atomr1-height1*height1);
                                        cone=r1*r1*(height1+height2);
                                        section=2.0*atomr1*atomr1*atomr1*(1-costheta1)+2.0*atomr2*atomr2*atomr2*(1-costheta2);
                                        DVolume=section-cone;
                                        OverlapVolume+=DVolume;
                                }
                                else
                                {
                                        //cout <<"Error calculating overlap volume. d="<<d<<" atomr1="<<atomr1<<" atomr2="<<atomr2<<" n="<<n<<" m="<<m<<endl;
                                        ofstream waxslog("OverlapError.txt", ios::app);
                                        waxslog <<"Error calculating overlap volume. d="<<d<<" atomr1="<<atomr1<<" atomr2="<<atomr2<<" n="<<n<<" m="<<m<<endl;
                                }
                        }

                        if ( d<atomr1 && d>atomr2 && (d+atomr1) < atomr2)
                        {
                                costheta1=(atomr1*atomr1+d*d-atomr2*atomr2)/(2.0*atomr1*d);
                                height1=atomr1*costheta1;
                                height2=d-height1;
                                costheta2=(atomr2*atomr2+d*d-atomr1*atomr1)/(2.0*atomr2*d);
                                if (atomr1>height1)
                                {
                                        r1=sqrt(atomr1*atomr1-height1*height1);
                                        DVolume=r1*r1*height2+2.0*atomr2*atomr2*atomr2*(1-costheta2)+2.0*atomr1*atomr1*atomr1*(1-costheta1)-r1*r1*height1;
                                }
                                else
                                {
                                        cout <<"Error calculating overlap volume. d="<<d<<" atomr1="<<atomr1<<" atomr2="<<atomr2<<" n="<<n<<" m="<<m<<endl;
                                        ofstream waxslog("OverlapError.txt", ios::app);
                                        waxslog <<"Error calculating overlap volume. d="<<d<<" atomr1="<<atomr1<<" atomr2="<<atomr2<<" n="<<n<<" m="<<m<<endl;
                                }
                        }

                        if ( d<atomr2 && d>atomr1 && (d+atomr1) > atomr2)
                        {
                                costheta2=(atomr2*atomr2+d*d-atomr1*atomr1)/(2.0*atomr2*d);
                                height2=atomr2*costheta1;
                                height1=d-height1;
                                costheta1=(atomr1*atomr1+d*d-atomr2*atomr2)/(2.0*atomr1*d);
                                if (atomr2>height2)
                                {
                                        r1=sqrt(atomr2*atomr2-height2*height2);
                                        DVolume=r1*r1*height1+2.0*atomr1*atomr1*atomr1*(1-costheta1)+2.0*atomr2*atomr2*atomr2*(1-costheta2)-r1*r1*height2;
                                }
                                else
                                {
                                        cout <<"Error calculating overlap volume. d="<<d<<" atomr1="<<atomr1<<" atomr2="<<atomr2<<" n="<<n<<" m="<<m<<endl;
                                        ofstream waxslog("OverlapError.txt", ios::app);
                                        waxslog <<"Error calculating overlap volume. d="<<d<<" atomr1="<<atomr1<<" atomr2="<<atomr2<<" n="<<n<<" m="<<m<<endl;
                                }
                        }

                        if ( (d+atomr1) < atomr2) OverlapVolume+=4.0*atomr1*atomr1*atomr1;
                        if ( (d+atomr2) < atomr1) OverlapVolume+=4.0*atomr2*atomr2*atomr2;
                }
        }

        OverlapVolume=OverlapVolume*pi/3.0;

        //cout <<"OverlapVolume="<<OverlapVolume<<endl;

        return OverlapVolume;
}

void MinMax(Real &XMin, Real &YMin, Real &ZMin, Real &XMax, Real &YMax, Real &ZMax, vector<AtomStruct> Atoms, int natom)
{
        //Finds the minimum and maximum coordinates of the first natom atoms
        //of a protein.  Similar functions are in MinMax.h.
        XMin=Atoms[0].x;
        YMin=Atoms[0].y;
        ZMin=Atoms[0].z;
        XMax=Atoms[0].x;
        YMax=Atoms[0].y;
        ZMax=Atoms[0].z;
        for (int n=0;n<natom;n++)
        {
                if (Atoms[n].x<XMin) XMin=Atoms[n].x;
                if (Atoms[n].y<YMin) YMin=Atoms[n].y;
                if (Atoms[n].z<ZMin) ZMin=Atoms[n].z;
                if (Atoms[n].x>XMax) XMax=Atoms[n].x;
                if (Atoms[n].y>YMax) YMax=Atoms[n].y;
                if (Atoms[n].z>ZMax) ZMax=Atoms[n].z;
        }
}

Real findvolume(int natom, Real atomr[], Real numatom[], vector<AtomStruct> Atoms)
{
        //Finds the total volume of a protein and scales the radii of the
        //excluded volume dummy atoms so that the total excluded volumes of the
        //excluded volume dummy atoms matches the volume of the protein.
        //This should be split into at least two functions.
        int nearest;
        Real inc, inc3;
        Real xgrid, ygrid, zgrid;
        Real value;
        Real volume, mindist;
        Real xmax, ymax, zmax;
        Real xmin, ymin, zmin;
        Real totalvolume;
        Real OverlapVolume;

        volume=0;
        inc=0.5;
        inc3=inc*inc*inc;

        MinMax(xmin, ymin, zmin, xmax, ymax, zmax, Atoms, natom);

        cout <<"In FindVolume"<<endl;

        for (xgrid=xmin-4;xgrid<xmax+4+inc;xgrid+=inc)
        {
                for (ygrid=ymin-4;ygrid<ymax+4+inc;ygrid+=inc)
                {
                        for (zgrid=zmin-4;zgrid<zmax+4+inc;zgrid+=inc)
                        {
                                value=location(xgrid, ygrid, zgrid, nearest, mindist, Atoms, 0);
                                if (value==InProtein)
                                {
                                        volume+=inc3;
                                }
                        }
                }
        }

        totalvolume=0;

        for (int n=0;n<6;n++)
        {
                //totalvolume+=numatom[n]*5.568328*atomr[n]*atomr[n]*atomr[n];
                totalvolume+=numatom[n]*4.0*pi*atomr[n]*atomr[n]*atomr[n]/3.0;
        }

        for (int n=0;n<6;n++)
        {
                if (totalvolume!=0)
                {
                        atomr[n]=atomr[n]*exp(log(volume/totalvolume)/3.0);
                }
        }
        cout <<"volume="<<volume<<endl;
        cout <<"About to enter OverlapVolume"<<endl;


        cout <<"volume="<<volume<<endl;
        OverlapVolume=Overlap(natom, atomr, Atoms);
        cout <<"OverlapVolume="<<OverlapVolume<<endl;

        return volume;

}


void FourGaussianParameters(Real a[][NumAtomTypes], Real b[][NumAtomTypes], Real c[])
{
        //The scattering factors of the atoms.  This is better than the 
        //scattering factors which are expressed as the sum of five Gaussians
        //and an off set, because this goes to zero as q goes to infinity.
        //Data obtained form:
        //Dirac-Fock calculations of X-ray scattering factors. Acta Cryst. (1994). A50, 481-497	

        for (int m=0;m<NumAtomTypes;m++)
        {
                c[m]=0;
                for (int n=0;n<5;n++)
                {
                        a[n][m]=0;
                        b[n][m]=0;
                }
        }

        a[0][HYDROGEN]=0.493002;                  
        b[0][HYDROGEN]=10.5109;                  
        a[1][HYDROGEN]=0.322912;
        b[1][HYDROGEN]=26.1257;
        a[2][HYDROGEN]=0.140191;
        b[2][HYDROGEN]=3.14236;
        a[3][HYDROGEN]=0.040810;
        b[3][HYDROGEN]=57.7997;

        a[0][CARBON]=2.6158;
        b[0][CARBON]=11.3632;
        a[1][CARBON]=0.2279;
        b[1][CARBON]=3.0830;
        a[2][CARBON]=1.5983;
        b[2][CARBON]=0.3787;
        a[3][CARBON]=1.5602;
        b[3][CARBON]=49.7088;

        a[0][NITROGEN]=0.4743;
        b[0][NITROGEN]=0.1041;
        a[1][NITROGEN]=2.9082;
        b[1][NITROGEN]=9.1890;
        a[2][NITROGEN]=2.2778;
        b[2][NITROGEN]=27.0869;
        a[3][NITROGEN]=1.3332;
        b[3][NITROGEN]=0.4612;

        a[0][OXYGEN]=3.4716;
        b[0][OXYGEN]=11.9964;
        a[1][OXYGEN]=1.8289;
        b[1][OXYGEN]=4.7941;
        a[2][OXYGEN]=1.7198;
        b[2][OXYGEN]=0.2372;
        a[3][OXYGEN]=0.9790;
        b[3][OXYGEN]=31.9917;

        a[0][SULFUR]=7.1301;
        b[0][SULFUR]=1.4247;
        a[1][SULFUR]=5.0712;
        b[1][SULFUR]=21.7545;
        a[2][SULFUR]=2.0611;
        b[2][SULFUR]=0.0994;
        a[3][SULFUR]=1.7369;
        b[3][SULFUR]=54.2128;

        a[0][IRON]=8.7267;
        b[0][IRON]=3.8443;
        a[1][IRON]=8.0070;
        b[1][IRON]=0.2319;
        a[2][IRON]=6.4440;
        b[2][IRON]=9.4874;
        a[3][IRON]=2.8058;
        b[3][IRON]=65.8373;

        a[0][PHOSPHORUS]=7.1583;
        b[0][PHOSPHORUS]=1.7276;
        a[1][PHOSPHORUS]=3.4832;
        b[1][PHOSPHORUS]=23.7455;
        a[2][PHOSPHORUS]=2.0781;
        b[2][PHOSPHORUS]=0.1140;
        a[3][PHOSPHORUS]=2.2770;
        b[3][PHOSPHORUS]=57.1173;
        /****************************************
          a[0][O_Minus]=3.106934;
          b[0][O_Minus]=19.86808;
          a[1][O_Minus]=3.235142;
          b[1][O_Minus]=6.960252;
          a[2][O_Minus]=1.148886;
          b[2][O_Minus]=0.170043;
          a[3][O_Minus]=0.783981;
          b[3][O_Minus]=65.693509;
          a[4][O_Minus]=0.676953;
          b[4][O_Minus]=0.630757;

          c[O_Minus]=0.046136;
         *****************************************/
        a[0][NA_Plus]=4.4278;
        b[0][NA_Plus]=5.2355;
        a[1][NA_Plus]=2.4274;
        b[1][NA_Plus]=2.2658;
        a[2][NA_Plus]=1.7182;
        b[2][NA_Plus]=0.1272;
        a[3][NA_Plus]=1.4258;
        b[3][NA_Plus]=12.6031;

        a[0][MG]=4.3643;
        b[0][MG]=2.1618;
        a[1][MG]=3.9083;
        b[1][MG]=5.9459;
        a[2][MG]=1.6872;
        b[2][MG]=0.0879;
        a[3][MG]=0.0382;
        b[3][MG]=0.0000;

        a[0][CL]=1.8908;
        b[0][CL]=0.0682;
        a[1][CL]=7.2032;
        b[1][CL]=1.1532;
        a[2][CL]=6.1851;
        b[2][CL]=18.8805;
        a[3][CL]=2.7190;
        b[3][CL]=51.3334;

        for (int m=0;m<NumAtomTypes;m++)
        {
                for (int n=0;n<5;n++)
                {
                        b[n][m]/=16.0*pi*pi;
                }
        }
}

void FiveGaussianParameters(Real a[][NumAtomTypes], Real b[][NumAtomTypes], Real c[])
{
        //Parameters for scattering factors for atoms.  This is not very good, 
        //because the off set means that the scattering factors never go to 
        //zero, which means obtaining p(r) by Fourier transform will not work.
        
        //Data obtained from:
        //Waasmaier, D., and Kirfel, A. 1995. New Analytical Scattering-Factor Functions for Free Atoms and Ions. Acta Cryst.A51: 416-431

        for (int m=0;m<NumAtomTypes;m++)
        {
                c[m]=0;
                for (int n=0;n<5;n++)
                {
                        a[n][m]=0;
                        b[n][m]=0;
                }
        }

        a[0][HYDROGEN]=0.493002;                     
        b[0][HYDROGEN]=10.5109;                      
        a[1][HYDROGEN]=0.322912;
        b[1][HYDROGEN]=26.1257;
        a[2][HYDROGEN]=0.140191;
        b[2][HYDROGEN]=3.14236;
        a[3][HYDROGEN]=0.040810;
        b[3][HYDROGEN]=57.7997;
        a[4][HYDROGEN]=0;
        b[4][HYDROGEN]=0;

        c[HYDROGEN]=0.003038;

        a[0][CARBON]=2.657506;
        b[0][CARBON]=14.780758;
        a[1][CARBON]=1.078079;
        b[1][CARBON]=0.7767775;
        a[2][CARBON]=1.490909;
        b[2][CARBON]=42.086843;
        a[3][CARBON]=-4.24107;
        b[3][CARBON]=-0.000294;
        a[4][CARBON]=0.713791;
        b[4][CARBON]=0.239535;

        c[CARBON]=4.297983;

        a[0][NITROGEN]=11.893780;
        b[0][NITROGEN]=0.000158;
        a[1][NITROGEN]=3.277479;
        b[1][NITROGEN]=10.232723;
        a[2][NITROGEN]=1.858092;
        b[2][NITROGEN]=30.34469;
        a[3][NITROGEN]=0.858927;
        b[3][NITROGEN]=0.656065;
        a[4][NITROGEN]=0.912985;
        b[4][NITROGEN]=0.217287;

        c[NITROGEN]=-11.804902;

        a[0][OXYGEN]=2.960427;
        b[0][OXYGEN]=14.182259;
        a[1][OXYGEN]=2.508818;
        b[1][OXYGEN]=5.936858;
        a[2][OXYGEN]=0.637853;
        b[2][OXYGEN]=0.112726;
        a[3][OXYGEN]=0.722838;
        b[3][OXYGEN]=34.958481;
        a[4][OXYGEN]=1.142756;
        b[4][OXYGEN]=0.39024;

        c[OXYGEN]=0.027014;

        a[0][SULFUR]=6.372157;
        b[0][SULFUR]=1.514347;
        a[1][SULFUR]=5.154568;
        b[1][SULFUR]=22.092528;
        a[2][SULFUR]=1.473732;
        b[2][SULFUR]=0.061373;
        a[3][SULFUR]=1.635073;
        b[3][SULFUR]=55.445176;
        a[4][SULFUR]=1.209372;
        b[4][SULFUR]=0.646925;

        c[SULFUR]=0.154722;

        a[0][IRON]=12.311098;
        b[0][IRON]=5.009415;
        a[1][IRON]=1.876623;
        b[1][IRON]=0.014461;
        a[2][IRON]=3.066177;
        b[2][IRON]=18.743041;
        a[3][IRON]=2.070451;
        b[3][IRON]=82.767874;
        a[4][IRON]=6.975185;
        b[4][IRON]=0.346506;

        c[IRON]=-0.304931;

        a[0][PHOSPHORUS]=1.950541;
        b[0][PHOSPHORUS]=0.908139;
        a[1][PHOSPHORUS]=4.14930;
        b[1][PHOSPHORUS]=27.044953;
        a[2][PHOSPHORUS]=1.494560;
        b[2][PHOSPHORUS]=0.071280;
        a[3][PHOSPHORUS]=1.522042;
        b[3][PHOSPHORUS]=67.520190;
        a[4][PHOSPHORUS]=5.729711;
        b[4][PHOSPHORUS]=1.981173;

        c[PHOSPHORUS]=0.155233;
        /*
           a[0][O_Minus]=3.106934;
           b[0][O_Minus]=19.86808;
           a[1][O_Minus]=3.235142;
           b[1][O_Minus]=6.960252;
           a[2][O_Minus]=1.148886;
           b[2][O_Minus]=0.170043;
           a[3][O_Minus]=0.783981;
           b[3][O_Minus]=65.693509;
           a[4][O_Minus]=0.676953;
           b[4][O_Minus]=0.630757;

           c[O_Minus]=0.046136;
           */
        a[0][NA_Plus]=3.148690;
        b[0][NA_Plus]=2.594987;
        a[1][NA_Plus]=4.073989;
        b[1][NA_Plus]=6.046925;
        a[2][NA_Plus]=0.767888;
        b[2][NA_Plus]=0.070139;
        a[3][NA_Plus]=0.995612;
        b[3][NA_Plus]=14.122657;
        a[4][NA_Plus]=0.968249;
        b[4][NA_Plus]=0.217037;

        c[NA_Plus]=0.0445300;

        
        a[0][MG]=3.062918;
        b[0][MG]=2.015803;
        a[1][MG]=4.135106;
        b[1][MG]=4.417941;
        a[2][MG]=0.853742;
        b[2][MG]=0.065307;
        a[3][MG]=1.036792;
        b[3][MG]=9.669710;
        a[4][MG]=0.852520;
        b[4][MG]=0.187818;

        c[MG]=0.058851;

        a[0][CL]=1.061802;
        b[0][CL]=0.144727;
        a[1][CL]=7.139886;
        b[1][CL]=1.171795;
        a[2][CL]=6.524271;
        b[2][CL]=19.467656;
        a[3][CL]=2.355626;
        b[3][CL]=60.320301;
        a[4][CL]=35.829404;
        b[4][CL]=0.000436;

        c[CL]=-34.916604;

        for (int m=0;m<NumAtomTypes;m++)
        {
                for (int n=0;n<5;n++)
                {
                        b[n][m]/=16.0*pi*pi;
                }
        }
}

void SANSParameters(Real a[][NumAtomTypes], Real b[][NumAtomTypes], Real c[])
{
        //Parmaters for the neutron scattering, scatting factors.
        for (int m=0;m<NumAtomTypes;m++)
        {
                c[m]=0;
                for (int n=0;n<5;n++)
                {
                        a[n][m]=0;
                        b[n][m]=0;
                }
        }

        c[HYDROGEN]=-0.3742;        //Data obtained from: en.wikipedia.org/wiki/Small-angle_neutron_scattering
        c[CARBON]=0.6651;
        c[NITROGEN]=0.940;
        c[OXYGEN]=0.5804;
        c[SULFUR]=0.2847;
        //c[IRON]=-0.304931;           Get value for this
        c[PHOSPHORUS]=0.517;

}
Real ScatteringFactor(Real a[][NumAtomTypes], Real b[][NumAtomTypes], Real c[], int atom, Real s)
{
        //Calculates the vacuum scattering factor of one atom.
        Real f=0;
        for (int j=0;j<5;j++)
        {
                f+=a[j][atom]/exp(s*s*b[j][atom]);
        }

        f+=c[atom];

        return f;
}

Real ScatteringFactor(int atom, Real s, string ScatteringType)
{
        //Gets the scattering factor parameters and calculates the 
        //scattering factor for one value of s for one atom.
        Real a[5][NumAtomTypes], b[5][NumAtomTypes], c[NumAtomTypes];
        //FiveGaussianParameters(a, b, c);
        if (ScatteringType=="X-ray") FourGaussianParameters(a, b, c);
        else if (ScatteringType=="Neutron") SANSParameters(a, b, c);
        else
        {
                cout <<"ERROR: Unknown ScatteringType "<<ScatteringType
                        <<".  The allowed options are X-ray and Neutron"<<endl;
                exit(EXIT_FAILURE);
        }
        return ScatteringFactor(a, b, c , atom, s);
}

void AtomScatteringFactor(IntensityStruct &i, string ScatteringType)
{
        //Gets the scattering factor parameters and calculates the 
        //scattering factors for all atom types and points.
        //This and the above functions cold be split up more.
        int points=i.calc.size();
        Real a[5][NumAtomTypes], b[5][NumAtomTypes], c[NumAtomTypes];
        //FiveGaussianParameters(a, b, c);
        if (ScatteringType=="X-ray") FourGaussianParameters(a, b, c);
        else if (ScatteringType=="Neutron") SANSParameters(a, b, c);
        else
        {
                cout <<"ERROR: Unknown ScatteringType "<<ScatteringType
                        <<".  The allowed options are X-ray and Neutron"<<endl;
                exit(EXIT_FAILURE);
        }
        for (int m=0;m<NumAtomTypes;m++)
        {
                for (int n=0;n<points;n++) i.f[n][m]=ScatteringFactor(a, b, c, m, i.s[n]);
        }
        //cout <<"HydrationShell atom scattering"<<endl;
        //for (int n=0;n<points;n++) cout <<i.s[n]<<"\t"<<i.f[n][HydrationShell]<<endl;
}

Real SphereScattering(Real radius, Real density, Real s)
{
        //Calculates the scattering factor of a sphere.
        Real f, satomr;
        satomr=s*radius;
        f=density*4.0*pi*(sin(satomr)-satomr*cos(satomr))/(s*s*s);
        return f;
}

Real GaussianScattering(Real radius, Real density, Real s)
{
        //Calculates the scattering factor of a sphere with a Gaussina density
        //profile.
        return density*radius*radius*radius*pi32/exp(s*s*radius*radius*0.25);
}

Real TrapazoidalScattering(Real radius, Real density, Real s)
{
        //Calculates the scattering factor of a sphere with a trapazoidal
        //density profile.
        Real f=0;
        Real atomr1=radius*0.407;
        Real atomr2=radius*1.593;
        Real cos1=cos(s*atomr1);
        Real cos2=cos(s*atomr2);
        Real sin1=sin(s*atomr1);
        Real sin2=sin(s*atomr2);
        f+=(sin1-s*atomr1*cos1)/(s*s*s);
        f+=(atomr2/(atomr2-atomr1))*(sin2-s*atomr2*cos2)/(s*s*s);
        f-=(atomr2/(atomr2-atomr1))*(sin1-s*atomr1*cos1)/(s*s*s);
        f-=(-atomr2*atomr2*cos2/s+2.0*atomr2*sin2/(s*s)+2.0*cos2/(s*s*s))/(s*(atomr2-atomr1));
        f+=(-atomr1*atomr1*cos1/s+2.0*atomr1*sin1/(s*s)+2.0*cos1/(s*s*s))/(s*(atomr2-atomr1));
        f*=4.0*pi*density;
        return f;
}

Real SphericalScattering(Real radius, Real density, Real s, string SphereType)
{
        //Calculates the scattering factor of a sphecally symettric dummy atom.
        //The atom has to be either a hard sphere, Gaussian sphere, or,
        //trapazoidal sphere.
        Real f;
        if (SphereType=="HardSphere") f=SphereScattering(radius, density, s);
        else if (SphereType=="GaussianSphere") f=GaussianScattering(radius, density, s);
        else if (SphereType=="TrapazoidalScattering") f=TrapazoidalScattering(radius, density, s);
        else
        {
                cout <<"Error.  Unrecognized SphereType "<<SphereType<<endl;
                exit(EXIT_FAILURE);
        }
        return f;
}

Real SolventCorrectedScattering(int n, Real radius, Real s, Real f, Real density, Real hsdensity, ParamStruct params)
{
        //Returns the scattering factor of a dummy excluded volume atom 
        //plus the scattering factor of the real atom at s.
        Real g;
        if (n!=WaterSphere && n!=HydrationShell)
        {
                g=SphericalScattering(radius, density, s, params.ExcludedVolumeSphereType);
                return f-g;
        }
        if (n==WaterSphere)
        {
                f=SphereScattering(radius, -density, s);
                if (params.ScatteringType=="X-ray") 
                {
                        f*=(ScatteringFactor(HYDROGEN, s, params.ScatteringType)*2.0+ScatteringFactor(OXYGEN, s, params.ScatteringType))/10.0;
                }
        }
        if (n==HydrationShell) f=SphericalScattering(radius, hsdensity, s, params.ExcludedVolumeSphereType);
        return f;
}

Real SolventCorrectedScattering(int n, Real radius, Real s, Real density, Real hsdensity, ParamStruct params)
{
        //Calculates the real scattering factor of an atom at s and subtracts
        //the scattering factor of a dummy excluded volume atom.
        Real f, g;
        f=ScatteringFactor(n, s, params.ScatteringType);
        g=SolventCorrectedScattering(n, radius, s, f, density, hsdensity, params);
        return g;
}

void MinMax2(int TotalParticles, Real &Min, Real &Max, vector<AtomStruct> Atoms)
{
        //Finds the minimum and maximum dimensions of the first TotalParticles
        //atoms of a protein.
        Min=Atoms[0].x;
        Max=Atoms[0].x;

        for (int n=0;n<TotalParticles;n++)
        {
                if (Atoms[n].x<Min) Min=Atoms[n].x;
                if (Atoms[n].y<Min) Min=Atoms[n].y;
                if (Atoms[n].z<Min) Min=Atoms[n].z;
                if (Atoms[n].x>Max) Max=Atoms[n].x;
                if (Atoms[n].y>Max) Max=Atoms[n].y;
                if (Atoms[n].z>Max) Max=Atoms[n].z;
        }
}

void VectorAverage(int natom, Real atomr[], Real contrast, Real hsdensity, int points, IntensityStruct &i, ParamStruct &params, vector<AtomStruct> Atoms)
{
        //Calculates the scattering from many different vectors and averages
        //to obtain the scattering pattern.
        bool TrigLookUp;
        char CharIntensityOutputFile[1000];
        string dir, base, IntensityOutputFile;
        int TotalParticles=Atoms.size();
        int IntMaxLookUp, IntDotProduct;
        Real angle;
        Real BinSize, InvBinSize;
        Real dx, dy, dz;
        Real Min, Max;
        Real MinX, MinY, MinZ;
        Real MaxX, MaxY, MaxZ;
        Real MinLookUp, MaxLookUp;
        Real mux, muy, muz;
        Real dotproduct;
        Real NumVectors;
        Real TotalIReal, TotalIImag;
        Real cosphi, nslices, sinphi, theta, thetainc;
        Real ireal[NumAtomTypes], iimag[NumAtomTypes];
        Real *CosLookUp, *SinLookUp;
        Real **f;

        dir=RemoveExtension(params.IntensityOutputFile) + "/";
        base=GetBase(params.IntensityOutputFile);
        IntensityOutputFile=dir + base;
        IntensityOutputFile=RemoveExtension(IntensityOutputFile)+"_phases.txt";
        AddIndexToFile(IntensityOutputFile);
        strcpy(CharIntensityOutputFile, IntensityOutputFile.c_str());
        ofstream file(CharIntensityOutputFile, ios::app);

        cout <<"In VectorAverage"<<endl;
        TrigLookUp=false;
        cout <<"TotalParticles= "<<TotalParticles<<endl;
        Safe2DArrayAlloc(f, points, NumAtomTypes, "f");

        for (int t=0;t<points;t++)
        {
                i.calc[t]=0;
                for (int n=0;n<NumAtomTypes;n++) f[t][n]=0;
        }

        for (int n=0;n<NumAtomTypes;n++)
        {
                for (int t=0;t<points;t++) f[t][n]=SolventCorrectedScattering(n, atomr[n], i.s[t], i.f[t][n], contrast, 0.0, params);
        }
        //Sets up the look up table for sin and cos.  This should be its own
        //function.
        if (TrigLookUp)
        {
                cout <<"TotalParticles= "<<TotalParticles<<endl;
                MinMax(MinX, MinY, MinZ, MaxX, MaxY, MaxZ, Atoms);
                Min=MinX+MinY+MinZ;
                Max=MaxX+MaxY+MaxZ;
                cout <<"Min= "<<Min<<" Max= "<<Max<<endl;
                if (abs(Min)>abs(Max)) Max=abs(Min);
                else Min=abs(Max);
                MinLookUp=-i.s[points-1]*Max;
                cout <<"MinLookUp= "<<MinLookUp<<endl;
                MaxLookUp=2.0*i.s[points-1]*Max;
                BinSize=0.001;
                InvBinSize=1.0/BinSize;
                IntMaxLookUp=int(floor(MaxLookUp*InvBinSize+0.5));
                cout <<"IntMaxLookUp= "<<IntMaxLookUp<<endl;
                SafeArrayAlloc(CosLookUp, IntMaxLookUp, "CosLookUp");
                SafeArrayAlloc(SinLookUp, IntMaxLookUp, "SinLookUp");

                angle=MinLookUp;
                for (int n=0;n<IntMaxLookUp;n++)
                {
                        CosLookUp[n]=cos(angle);
                        SinLookUp[n]=sin(angle);
                        if (n<100) cout <<"angle= "<<angle<<" CosLookUp["<<n<<"]= "<<CosLookUp[n]<<" SinLookUp["<<n<<"]= "<<SinLookUp[n]<<endl;
                        angle+=BinSize;
                }
        }


        nslices=Real(params.VectorsPerInclination);
        //nslices=5.0;
        cosphi=1.0/nslices-1.0;
        thetainc=Real(params.VectorsPerInclination);
        //thetainc=5.0;
        theta=pi/thetainc;

        NumVectors=0.0;

        while (true)
        {
                cout <<"cosphi= "<<cosphi<<" theta= "<<theta<<endl;
                sinphi=sqrt(1.0-cosphi*cosphi);
                //Finds the direction of the scattering vector.
                dx=sinphi*cos(theta);
                dy=sinphi*sin(theta);
                dz=cosphi;

                NumVectors+=1.0;
                for (int t=0;t<points;t++)
                {
                        //Everything in this loop could be a separate function.
                        for (int n=0;n<NumAtomTypes;n++)
                        {
                                ireal[n]=0;
                                iimag[n]=0;
                        }
                        TotalIReal=0;
                        TotalIImag=0;
                        //The components of the scattering vector.
                        mux=i.s[t]*dx;
                        muy=i.s[t]*dy;
                        muz=i.s[t]*dz;
                        //These could be there own functions.
                        if (TrigLookUp)
                        {
                                for (int n=0;n<TotalParticles;n++)
                                {
                                        //cout <<"n= "<<n<<endl;
                                        dotproduct=mux*Atoms[n].x+muy*Atoms[n].y+muz*Atoms[n].z;
                                        //if (n<100 && t==0) cout <<"n= "<<n<<" dotproduct= "<<dotproduct<<endl;
                                        dotproduct-=MinLookUp;
                                        IntDotProduct=int(floor(dotproduct*InvBinSize+0.5));
                                        //cout <<"dotproduct= "<<dotproduct<<" CosLookUp["<<IntDotProduct<<"]= "<<CosLookUp[IntDotProduct]<<" SinLookUp= "<<SinLookUp[IntDotProduct]<<endl;
                                        ireal[Atoms[n].atomid]+=CosLookUp[IntDotProduct]*Atoms[n].weight;
                                        iimag[Atoms[n].atomid]+=SinLookUp[IntDotProduct]*Atoms[n].weight;
                                        //if (n<100 && t==0)
                                        //{
                                        //	cout <<"CosLookUp["<<IntDotProduct<<"]= "<<CosLookUp[IntDotProduct]<<" SinLookUp= "<<SinLookUp[IntDotProduct]<<endl;
                                        //	cout <<"weight= "<<weight[n]<<endl;
                                        //}
                                }

                        }
                        else
                        {
                                for (int n=0;n<TotalParticles;n++)
                                {
                                        dotproduct=mux*Atoms[n].x+muy*Atoms[n].y+muz*Atoms[n].z;
                                        ireal[Atoms[n].atomid]+=cos(dotproduct)*Atoms[n].weight;
                                        iimag[Atoms[n].atomid]+=sin(dotproduct)*Atoms[n].weight;
                                }
                        }
                        //Multiply the real and imaginary parts by the scatteing factors.
                        for (int n=0;n<NumAtomTypes;n++) 
                        {
                                TotalIReal+=ireal[n]*f[t][n];
                                TotalIImag+=iimag[n]*f[t][n];
                                //if (t==0)
                                //{
                                //	cout <<"TotalIReal= "<<TotalIReal<<" TotalIImag= "<<TotalIImag<<
                                //	" ireal["<<n<<"]= "<<ireal[n]<<" iimag["<<n<<"]= "<<iimag[n]<<
                                //	" f["<<t<<"]["<<n<<"]= "<<f[t][n]<<endl;
                                //}
                        }
                        file <<t<<"\t"<<i.s[t]<<"\t"<<mux<<"\t"<<muy<<"\t"<<muz<<"\t"<<TotalIReal<<"\t"<<TotalIImag<<endl;
                        i.calc[t]=i.calc[t]+TotalIReal*TotalIReal+TotalIImag*TotalIImag;
                }
                //Change the direction.  Break out of loop if scattering has been calculated for all vectors.
                theta += 2.0*pi/thetainc;
                if (theta>2.0*pi)
                {
                        theta=pi/thetainc;
                        cosphi += 2.0/nslices;
                        if (cosphi>0.0-0.9/nslices) break;
                }
        }

        for (int n=0;n<points;n++) i.calc[n]=i.calc[n]*0.5/NumVectors;

}

void cube3(int natom, Real atomr[], Real numatom[], Real contrast, Real hsdensity, int points, IntensityStruct &i, ParamStruct &params, vector<AtomStruct> Atoms, bool UniformHydrationShell)
{
        //An early version of the cube method.  Saves time by joining cubes
        //with the same density together.  Probably broken.
        int TotalParticles=Atoms.size();
        int nearest, distance;
        int previous, value;
        int CubeNum;
        int MAX;
        vector<Real> CubeX, CubeY, CubeZ;
        vector<Real> CubeWeight;
        Real xgrid, ygrid, zgrid;
        Real xmax, xmin, ymax, ymin, zmax, zmin;
        Real dx, dy, dz;
        Real mux, muy, muz;
        Real ireal, iimag;
        Real dotproduct, product;
        Real sinmuz;
        vector<Real> cubelz;
        vector< vector<Real> > AtomTypesDensityMatrix, ElementsDensityMatrix;
        Real f[1000][NumAtomTypes];
        Real dist, mindist;
        Real inc;
        Real volume;
        Real thetainc;
        Real theta;
        Real stop;

        Real phibin[100];
        Real thetabin[100];
        Real phi2;
        Real cosphi;
        Real nslices;
        Real sinphi;

        MAX=TotalParticles*1000;

        CubeNum=0;
        MinMax(xmin, ymin, zmin, xmax, ymax, zmax, Atoms, natom);
        ReadGofRFile(params.AtomTypesGofRFile, contrast, AtomTypesDensityMatrix, params.RecBin);
        ReadGofRFile(params.ElementsGofRFile, contrast, ElementsDensityMatrix, params.RecBin);


        inc=0.50;

        cout <<"In cube3"<<endl;

        for (int t=0;t<points;t++) i.calc[t]=0;
        //The cube density should not be calculated at the same time that
        //the cubes are joined.  This should also be its own function.
        for (xgrid=xmin-4;xgrid<xmax+4+inc;xgrid=xgrid+inc)
        {
                for (ygrid=ymin-4;ygrid<ymax+4+inc;ygrid=ygrid+inc)
                {
                        previous=InSolution;
                        for (zgrid=zmin-4;zgrid<zmax+4+inc;zgrid=zgrid+inc)
                        {
                                value=location(xgrid, ygrid, zgrid, nearest, mindist, Atoms, params.maxhs);

                                if (value==InProtein)
                                {
                                        //cout <<"value="<<value<<" previous="<<previous<<endl;
                                        if (previous==InProtein)
                                        {
                                                cubelz[CubeNum]+=inc;
                                                CubeZ[CubeNum]+=inc*0.5;
                                                //cout <<"cubelz["<<CubeNum<<"]="<<cubelz[CubeNum]<<endl;
                                        }

                                        if (previous!=InProtein)
                                        {
                                                previous=InProtein;
                                                SafePushBack(CubeX, xgrid, "CubeX");
                                                SafePushBack(CubeY, ygrid, "CubeY");
                                                SafePushBack(CubeZ, zgrid, "CubeZ");
                                                SafePushBack(cubelz, inc, "cubelz");
                                        }
                                }

                                if ( value!=InProtein && previous==InProtein)
                                {
                                        previous=value;
                                        SafePushBack(CubeWeight, -contrast, "CubeWeight");
                                        CubeNum++;
                                }

                                if (value==InHydrationShell)
                                {
                                        SafePushBack(CubeX, xgrid, "CubeX");
                                        SafePushBack(CubeY, ygrid, "CubeY");
                                        SafePushBack(CubeZ, zgrid, "CubeZ");
                                        SafePushBack(cubelz, inc, "cubelz");

                                        if (UniformHydrationShell) SafePushBack(CubeWeight, hsdensity, "CubeWeight");
                                        else
                                        {
                                                dx=xgrid-Atoms[nearest].x;
                                                dy=ygrid-Atoms[nearest].y;
                                                dz=zgrid-Atoms[nearest].z;
                                                dist=sqrt(dx*dx+dy*dy+dz*dz);
                                                distance=int(floor(dist*10.0+0.5));
                                                SafePushBack(CubeWeight, AtomTypesDensityMatrix[Atoms[nearest].AtomType][distance], "CubeWeight");
                                        }
                                        CubeNum++;
                                }
                        }
                }
        }

        volume=0;
        int CubeLZ=cubelz.size();
        for (int t=0;t<CubeLZ;t++)
        {
                volume+=inc*inc*cubelz[t];
        }

        for (int t=0;t<1000;t++)
        {
                for (int n=0;n<NumAtomTypes;n++)
                {
                        f[t][n]=0;
                }
        }

        for (int n=0;n<NumAtomTypes;n++)
        {
                for (int t=0;t<points;t++)
                {
                        f[t][n]=SolventCorrectedScattering(n, 0, i.s[t], i.f[t][n], 0, hsdensity, params);
                }
        }

        cout <<"natom="<<natom<<endl;
        cout <<"f[0][0]="<<f[0][0]<<" f[0][1]="<<f[0][1]<<endl;

        cout <<"contrast="<<contrast<<endl;

        cout <<"CubeNum="<<CubeNum<<" cubelz.size="<<cubelz.size()<<"CubeWeight.size="<<CubeWeight.size()<<" CubeX.size="<<CubeX.size()<<endl;

        stop=0;
        nslices=10.0;
        cosphi=1.0/nslices-1.0;
        thetainc=10.0;
        theta=pi/thetainc;

        //cout <<"About to average"<<endl;
        //Calculates scattering for many different vectors and averages.
        //See void VectorAverage for more info.
        while (true)
        {
                //cout <<"cosphi="<<cosphi<<endl;
                sinphi=sqrt(1-cosphi*cosphi);
                dx=sinphi*cos(theta);
                dy=sinphi*sin(theta);
                dz=cosphi;

                //PrintPdbLine("c://Directions3.pdb", m, "O", "HOH", dx, dy, dz, 1.0, 1.0);

                nearest=int(floor(theta*10+0.5));
                thetabin[nearest]=thetabin[nearest]+1.0;
                phi2=atan(dz/dx);

                if (phi2<0)
                {
                        phi2=phi2+pi;
                }
                nearest=int(floor(phi2*10+0.5));
                phibin[nearest]=phibin[nearest]+1.0;
                theta += 2.0*pi/thetainc;
                if (theta>2*pi)
                {
                        theta=pi/thetainc;
                        cosphi += 2.0/nslices;
                        if (cosphi>1.0-1.0/nslices)
                        {
                                break;
                        }
                }

                for (int t=0;t<points;t++)
                {
                        ireal=0;
                        iimag=0;
                        mux=i.s[t]*dx;
                        muy=i.s[t]*dy;
                        muz=i.s[t]*dz;
                        product=8.0*sin(mux*inc*0.5)*sin(muy*inc*0.5)/(mux*muy*muz);
                        for (int n=0;n<natom;n++)
                        {
                                dotproduct=mux*Atoms[n].x+muy*Atoms[n].y+muz*Atoms[n].z;
                                ireal+=f[t][Atoms[n].atomid]*cos(dotproduct);
                                iimag+=f[t][Atoms[n].atomid]*sin(dotproduct);
                        }

                        for (int n=0;n<CubeNum;n++)
                        {
                                dotproduct=mux*CubeX[n]+muy*CubeY[n]+muz*CubeZ[n];
                                sinmuz=sin(muz*cubelz[n]*0.5)*CubeWeight[n]*product;
                                ireal+=sinmuz*cos(dotproduct);
                                iimag+=sinmuz*sin(dotproduct);
                        }

                        i.calc[t]=i.calc[t]+ireal*ireal+iimag*iimag;
                }
        }
}

void cube4b(int natom, Real atomr[], Real contrast, Real hsdensity, IntensityStruct &i, vector<AtomStruct> Atoms, ParamStruct &params)
{
        //Another version of the cube method.  The difference between
        //this and the previous cube3 seems to be that cavities are 
        //eliminated.
        char CharCubeDensityPdbFile[1000];
        int TotalParticles=Atoms.size();
        int nearest, distance;
        int previous, value;
        int CubeNum;
        int MAX;
        int *WithinCutOff;
        Real *CubeX, *CubeY, *CubeZ;
        Real *CubeWeight;
        vector< vector<Real> > AtomTypesDensityMatrix, ElementsDensityMatrix;
        Real mindist;
        Real xgrid, ygrid, zgrid;
        Real xmax, xmin, ymax, ymin, zmax, zmin;
        Real dx, dy, dz;
        Real mux, muy, muz;
        Real ireal, iimag, IrealCube, IimagCube;
        Real dotproduct, product;
        Real sinmuz;
        Real *cubelz;
        Real f[1000][NumAtomTypes];
        Real dist, r;
        Real inc;
        Real volume;
        Real thetainc;
        Real theta;
        Real stop;

        Real cosphi;
        Real nslices;
        Real sinphi;
        string PdbOutputFile; //Initialize
        PrintPdb(PdbOutputFile, Atoms);

        ReadGofRFile(params.AtomTypesGofRFile, contrast, AtomTypesDensityMatrix, params.RecBin);
        ReadGofRFile(params.ElementsGofRFile, contrast, ElementsDensityMatrix, params.RecBin);

        //center(atomm, natom, 0, atomr);

        MoveAtoms(Atoms, -params.XOrigin, -params.YOrigin, -params.ZOrigin);

        cout <<"natom= "<<natom<<" TotalParticles= "<<TotalParticles<<endl;

        PrintPdb(PdbOutputFile, Atoms);

        cout <<"In cube2"<<endl;

        MAX=TotalParticles*1000;

        cout <<"MAX= "<<MAX<<endl;

        CubeNum=0;

        SafeArrayAlloc(CubeX, MAX, "CubeX");
        cout <<"Allocated cubex"<<endl;
        SafeArrayAlloc(CubeY, MAX, "CubeX");
        cout <<"Allocated cubey"<<endl;
        SafeArrayAlloc(CubeZ, MAX, "CubeX");
        cout <<"Allocated cubez"<<endl;
        SafeArrayAlloc(CubeWeight, MAX, "CubeX");
        cout <<"Allocated cubeweight"<<endl;
        SafeArrayAlloc(cubelz, MAX, "CubeX");

        cout <<"Allocated memory"<<endl;

        for (int n=0;n<MAX;n++)
        {
                CubeX[n]=0;
                CubeY[n]=0;
                CubeZ[n]=0;
                CubeWeight[n]=0;
        }

        cout <<"Initialized CubeX CubeY CubeZ CubeWeight"<<endl;

        strcpy(CharCubeDensityPdbFile, params.OutputCubeDensityPdbFile.c_str());

        inc=0.50;
        for (int t=0;t<MAX;t++) cubelz[t]=0;

        cout <<"After initializing cubelz"<<endl;

        for (int t=0;t<params.points;t++) i.calc[t]=0;

        MinMax(xmin, ymin, zmin, xmax, ymax, zmax, Atoms);

        for (xgrid=xmin-10;xgrid<xmax+10+inc;xgrid=xgrid+inc)
        {
                cout <<"xgrid= "<<xgrid<<" CubeNum= "<<CubeNum<<" MAX= "<<MAX<<endl;
                for (ygrid=ymin-10;ygrid<ymax+10+inc;ygrid=ygrid+inc)
                {
                        previous=InSolution;
                        for (zgrid=zmin-10;zgrid<zmax+10+inc;zgrid=zgrid+inc)
                        {
                                value=location(xgrid, ygrid, zgrid, nearest, mindist, Atoms, params.maxhs);

                                if (value==InProtein)
                                {
                                        //cout <<"value="<<value<<" previous="<<previous<<endl;
                                        if (previous==InProtein)
                                        {
                                                cubelz[CubeNum]+=inc;
                                                CubeZ[CubeNum]+=inc*0.5;
                                                //cout <<"cubelz["<<CubeNum<<"]="<<cubelz[CubeNum]<<endl;
                                        }

                                        if (previous!=InProtein)
                                        {
                                                previous=InProtein;
                                                CubeX[CubeNum]=xgrid;
                                                CubeY[CubeNum]=ygrid;
                                                CubeZ[CubeNum]=zgrid;
                                                cubelz[CubeNum]=inc;
                                        }
                                }

                                if ( value!=InProtein && previous==InProtein)
                                {
                                        previous=value;
                                        CubeWeight[CubeNum]=-contrast;
                                        CubeNum++;
                                }

                                r=sqrt(xgrid*xgrid+ygrid*ygrid+zgrid*zgrid);
                                if (value==InHydrationShell && r<42.5)
                                {
                                        CubeX[CubeNum]=xgrid;
                                        CubeY[CubeNum]=ygrid;
                                        CubeZ[CubeNum]=zgrid;
                                        cubelz[CubeNum]=inc;

                                        if (params.UniformHydrationShell) CubeWeight[CubeNum]=hsdensity;
                                        else
                                        {
                                                //cout <<"nearest= "<<nearest<<endl;
                                                dx=xgrid-Atoms[nearest].x;
                                                dy=ygrid-Atoms[nearest].y;
                                                dz=zgrid-Atoms[nearest].z;

                                                dist=sqrt(dx*dx+dy*dy+dz*dz);
                                                distance=int(floor(dist/params.RecBin+0.5));
                                                //cout <<"distance= "<<distance<<endl;
                                                //cout <<"CubeNum="<<CubeNum<<" nearest="<<nearest<<" distance="<<distance<<" type="<<solventtype[nearest]<<" MAX="<<MAX<<endl;
                                                //cout <<"solventtype["<<nearest<<"]= "<<solventtype[nearest]<<endl;
                                                GetDensity(CubeWeight[CubeNum], Atoms[nearest], distance, AtomTypesDensityMatrix, ElementsDensityMatrix, 0, params);
                                        }
                                        CubeNum++;
                                }
                        }
                }
        }

        Real xn, yn, zn;

        SafeArrayAlloc(WithinCutOff, CubeNum, "WithinCutOff");

        for (int n=0;n<CubeNum;n++)
        {
                WithinCutOff[n]=0;
        }
        Real CutOff=5.0;
        Real CutOff2=CutOff*CutOff;
        for (int n=0;n<CubeNum-1;n++)
        {
                if (CubeWeight[n]!=-contrast)
                {
                        xn=CubeX[n];
                        yn=CubeY[n];
                        zn=CubeZ[n];
                        for (int m=n+1;m<CubeNum;m++)
                        {
                                if (CubeWeight[m]!=-contrast)
                                {
                                        dx=xn-CubeX[m];
                                        dy=yn-CubeY[m];
                                        dz=zn-CubeZ[m];
                                        r=dx*dx+dy*dy+dz*dz;
                                        if (r<CutOff2)
                                        {
                                                WithinCutOff[n]++;
                                                WithinCutOff[m]++;
                                        }
                                }
                        }
                }
        }

        for (int n=0;n<CubeNum;n++)
        {
                if (WithinCutOff[n]<500.0 && CubeWeight[n]!=-contrast)
                {
                        cout <<"WithinCutOff["<<n<<"]= "<<WithinCutOff[n]<<endl;
                }
                if (WithinCutOff[n]<500)
                {
                        CubeWeight[n]=-contrast;	
                }
        }

        Real HydrationShellMass=0;

        for (int n=0;n<CubeNum;n++)
        {
                HydrationShellMass+=CubeWeight[n]*cubelz[n];
        }

        HydrationShellMass=HydrationShellMass*inc*inc;
        cout <<"HydrationShellMass= "<<HydrationShellMass<<endl;

        for (int n=0;n<CubeNum;n++)
        {
                if (CubeWeight[n]!=-contrast)
                {
                        PrintPdbLine(params.OutputCubeDensityPdbFile, true, 1, "O", "HOH", 1, CubeX[n], CubeY[n], CubeZ[n], 1.0, CubeWeight[n]);
                }
        }

        ofstream pdb2(CharCubeDensityPdbFile, ios::app);
        pdb2  << "END"<<endl;

        cout <<"CubeNum="<<CubeNum<<endl;

        volume=0;
        for (int t=0;t<3000;t++)
        {
                volume+=inc*inc*cubelz[t];
        }

        for (int t=0;t<1000;t++)
        {
                for (int n=0;n<NumAtomTypes;n++)
                {
                        f[t][n]=0;
                }
        }

        for (int n=0;n<NumAtomTypes;n++)
        {
                for (int t=0;t<params.points;t++)
                {
                        f[t][n]=SolventCorrectedScattering(n, atomr[n], i.s[t], i.f[t][n], 0, 0, params);
                }
        }

        cout <<"f[0][0]="<<f[0][0]<<endl;

        stop=0;
        nslices=20.0;
        cosphi=1.0/nslices-1.0;
        thetainc=20.0;
        theta=pi/thetainc;

        //cout <<"About to average"<<endl;

        while (true)
        {
                cout <<"cosphi="<<cosphi<<" theta="<<theta<<endl;
                sinphi=sqrt(1-cosphi*cosphi);
                dx=sinphi*cos(theta);
                dy=sinphi*sin(theta);
                dz=cosphi;


                /*********************************************************************************************************
                  PrintPdbLine("/home/jouko/WAXS/Directions4.pdb", m, "O", "HOH", 1, dx, dy, dz, 1.0, 1.0);
                 ***************************************************************************************************/

                for (int t=0;t<params.points;t++)
                {
                        ireal=0;
                        iimag=0;
                        IrealCube=0;
                        IimagCube=0;
                        mux=i.s[t]*dx;
                        muy=i.s[t]*dy;
                        muz=i.s[t]*dz;
                        product=8.0*sin(mux*inc*0.5)*sin(muy*inc*0.5)/(mux*muy*muz);

                        for (int n=0;n<natom;n++)
                        {
                                dotproduct=mux*Atoms[n].x+muy*Atoms[n].y+muz*Atoms[n].z;
                                ireal+=f[t][Atoms[n].atomid]*cos(dotproduct);
                                iimag+=f[t][Atoms[n].atomid]*sin(dotproduct);
                        }

                        if (t==100)
                        {
                                cout <<"ireal="<<ireal<<" iimag="<<iimag<<" i[100]="<<i.calc[100]<<endl;
                        }

                        for (int n=0;n<CubeNum;n++)
                        {
                                dotproduct=mux*CubeX[n]+muy*CubeY[n]+muz*CubeZ[n];
                                sinmuz=sin(muz*cubelz[n]*0.5);
                                IrealCube+=cos(dotproduct)*CubeWeight[n]*sinmuz;
                                IimagCube+=sin(dotproduct)*CubeWeight[n]*sinmuz;
                        }

                        ireal+=IrealCube*product;
                        iimag+=IimagCube*product;

                        if (t==100)
                        {
                                cout <<"ireal="<<ireal<<" iimag="<<iimag<<" i[100]="<<i.calc[100]<<endl;
                                cout <<endl<<endl;
                        }

                        i.calc[t]=i.calc[t]+ireal*ireal+iimag*iimag;
                }

                theta += 2.0*pi/thetainc;
                if (theta>2.0*pi)
                {
                        theta=pi/thetainc;
                        cosphi += 2.0/nslices;
                        if (cosphi>0.0-0.9/nslices)
                        {
                                break;
                        }
                }
        }

        if (params.points>100) cout <<"i[100]="<<i.calc[100]<<endl;

}

Real FindProteinElectrons(vector<AtomStruct> &Atoms, vector<CubeStruct> &cubes, Real CubeSize)
{
        //Finds the total number of elextrons in the system, when the cube
        //method is used.
        int natom=Atoms.size();
        int ncube=cubes.size();
        Real TotalElectrons=0, AtomElectrons=0, ExcludedVolumeElectrons=0;
        Real HydrationShellElectrons=0, CubeVolume;
        Real NumElectrons[NumAtomTypes];
        CubeVolume=CubeSize*CubeSize*CubeSize;
        cout <<"CubeVolume= "<<CubeVolume<<endl;
        cout <<"ncube= "<<ncube<<endl;
        SetNumElectrons(NumElectrons);
        for (int j=0;j<natom;j++) AtomElectrons+=NumElectrons[Atoms[j].atomid];
        for (int j=0;j<ncube;j++)
        {
                if (cubes[j].LocatedIn==InProtein) ExcludedVolumeElectrons+=cubes[j].density*CubeVolume;
                if (cubes[j].LocatedIn==InHydrationShell) 
                {
                        //cout <<"cubes["<<j<<"].density= "<<cubes[j].density<<endl;
                        HydrationShellElectrons+=cubes[j].density*CubeVolume;
                }
        }
        TotalElectrons=AtomElectrons+ExcludedVolumeElectrons+HydrationShellElectrons;
        cout <<"TotalElectrons= "<<TotalElectrons<<endl
                <<"AtomElectrons= "<<AtomElectrons<<endl
                <<"ExcludedVolumeElectrons= "<<ExcludedVolumeElectrons<<endl
                <<"HydrationShellElectrons= "<<HydrationShellElectrons<<endl;
        return TotalElectrons;
}

void CalcCubeIntensity(vector<CubeStruct> &cubes, IntensityStruct &i, vector<AtomStruct> &Atoms, ParamStruct &params)
{
        //Calculates scattering from a set of cubes using Fourier transform.
        //Calculates scattering for many different vectors and averages.
        int CubeNum=cubes.size();
        int points=i.calc.size();
        int natom;
        Real cosphi, sinphi, theta, thetainc;
        Real dotproduct;
        Real dx, dy, dz;
        Real ireal, iimag;
        Real IrealCube, IimagCube;
        Real mux, muy, muz;
        Real nslices, NumVectors;
        Real product;
        Real **f;

        cout <<"In CalcCubeIntensity"<<endl;
        Safe2DArrayAlloc(f, points, NumAtomTypes, "f");
        cout <<"Allocated memory for f"<<endl;
        for (int m=0;m<points;m++)
        {
                for (int n=0;n<NumAtomTypes;n++)
                {
                        f[m][n]=0;
                }
        }

        for (int m=0;m<NumAtomTypes;m++)
        {
                for (int n=0;n<points;n++)
                {
                        f[n][m]=SolventCorrectedScattering(m, 0, i.s[n], i.f[n][m], 0, 0, params);
                }
        }
        cout <<"About to enter FindProteinElectrons"<<endl;
        FindProteinElectrons(Atoms, cubes, params.CubeSize);
        nslices=params.VectorsPerInclination;
        cosphi=1.0/nslices-1.0;
        thetainc=params.VectorsPerInclination;
        theta=pi/thetainc;
        NumVectors=0.0;
        cout <<"About to enter FindNumProteinAtoms"<<endl;
        natom=FindNumProteinAtoms(Atoms);
        cout <<"natom= "<<natom<<endl;
        while (true)
        {
                cout <<"cosphi="<<cosphi<<" theta="<<theta<<endl;
                sinphi=sqrt(1-cosphi*cosphi);
                dx=sinphi*cos(theta);
                dy=sinphi*sin(theta);
                dz=cosphi;
                NumVectors+=1.0;

                /*********************************************************************************************************
                  PrintPdbLine("/home/jouko/WAXS/Directions4.pdb", m, "O", "HOH", 1, dx, dy, dz, 1.0, 1.0);
                 ***************************************************************************************************/

                for (int t=0;t<points;t++)
                {
                        ireal=0;
                        iimag=0;
                        IrealCube=0;
                        IimagCube=0;
                        mux=i.s[t]*dx;
                        muy=i.s[t]*dy;
                        muz=i.s[t]*dz;
                        product=8.0*sin(mux*params.CubeSize*0.5)*sin(muy*params.CubeSize*0.5)*sin(muz*params.CubeSize*0.5)/(mux*muy*muz);

                        for (int n=0;n<natom;n++)
                        {
                                dotproduct=mux*Atoms[n].x+muy*Atoms[n].y+muz*Atoms[n].z;
                                ireal+=f[t][Atoms[n].atomid]*cos(dotproduct);
                                iimag+=f[t][Atoms[n].atomid]*sin(dotproduct);
                                if (t==10)
                                {
                                        //PrintAtomInfo(Atoms[t]);
                                        //cout <<"dotproduct= "<<dotproduct<<endl;
                                }
                        }

                        if (t==0)
                        {
                                cout <<"dx= "<<dx<<" dy= "<<dy<<" dz= "<<dz<<endl;
                                cout <<"mux= "<<mux<<" muy= "<<muy<<" muz= "<<muz<<endl;
                                cout <<"f["<<t<<"][0]= "<<f[t][0]<<endl;
                                cout <<"ireal="<<ireal<<" iimag="<<iimag<<" i["<<t<<"]="<<i.calc[t]<<endl;
                        }

                        for (int n=0;n<CubeNum;n++)
                        {
                                dotproduct=mux*cubes[n].x+muy*cubes[n].y+muz*cubes[n].z;
                                IrealCube+=cos(dotproduct)*cubes[n].density;
                                IimagCube+=sin(dotproduct)*cubes[n].density;
                        }

                        if (t==0)
                        {
                                cout <<"product= "<<product<<endl;
                                cout <<"IrealCube["<<t<<"]= "<<IrealCube<<endl;
                                cout <<"IimagCube["<<t<<"]= "<<IimagCube<<endl;
                        }
                        if (t==10)
                        {
                                cout <<"IrealCube["<<t<<"]= "<<IrealCube<<endl;
                                cout <<"IimagCube["<<t<<"]= "<<IimagCube<<endl;
                        }
                        ireal+=(IrealCube*product);
                        iimag+=(IimagCube*product);

                        if (t==0)
                        {
                                cout <<"ireal="<<ireal<<" iimag="<<iimag<<" i["<<t<<"]="<<i.calc[t]<<endl;
                                cout <<endl<<endl;
                        }

                        i.calc[t]=i.calc[t]+ireal*ireal+iimag*iimag;
                }

                theta += 2.0*pi/thetainc;
                if (theta>2.0*pi)
                {
                        theta=pi/thetainc;
                        cosphi += 2.0/nslices;
                        if (cosphi>0.0-0.9/nslices)
                        {
                                break;
                        }
                }
        }
        cout <<"NumVectors= "<<NumVectors<<endl;
        for (int t=0;t<points;t++) i.calc[t]=i.calc[t]*0.5/NumVectors;

        for (int n=0;n<points;n++)
        {
                delete [] f[n];
        }
        delete [] f;
}

void FindCavitiesByCutOff(Real contrast, vector<CubeStruct> &cubes)
{
        //Assumes that if a cube is not surrounded by a certain number of 
        //other cubes then it must be in a cavity and eliminates it.
        int CubeNum=cubes.size();
        int *WithinCutOff;
        Real r;
        Real dx, dy, dz;
        Real xn, yn, zn;

        SafeArrayAlloc(WithinCutOff, CubeNum, "WithinCutOff");
        for (int n=0;n<CubeNum;n++)
        {
                WithinCutOff[n]=0;
        }
        Real CutOff=5.0;
        Real CutOff2=CutOff*CutOff;
        for (int n=0;n<CubeNum-1;n++)
        {
                if (cubes[n].density!=-contrast)
                {
                        xn=cubes[n].x;
                        yn=cubes[n].y;
                        zn=cubes[n].z;
                        for (int m=n+1;m<CubeNum;m++)
                        {
                                if (cubes[m].density!=-contrast)
                                {
                                        dx=xn-cubes[m].x;
                                        dy=yn-cubes[m].y;
                                        dz=zn-cubes[m].z;
                                        r=dx*dx+dy*dy+dz*dz;
                                        if (r<CutOff2)
                                        {
                                                WithinCutOff[n]++;
                                                WithinCutOff[m]++;
                                        }
                                }
                        }
                }
        }

        for (int n=0;n<CubeNum;n++)
        {
                if (WithinCutOff[n]<500.0 && cubes[n].density!=-contrast)
                {
                        cout <<"WithinCutOff["<<n<<"]= "<<WithinCutOff[n]<<endl;
                }
                if (WithinCutOff[n]<500)
                {
                        cubes[n].density=-contrast;	
                }
        }

}

void AssignCubeDensities1(vector<CubeStruct> &cubes, Real contrast, vector<AtomStruct> Atoms, ParamStruct &params)
{
        //This is an old function for generating the hydration shell.  It
        //has been made obsolete by AssignCubeDensities2.
        int CubeNum, distance, Location, nearest;
        int MaxXBin, MaxYBin, MaxZBin;
        Real dist;
        Real dx, dy, dz;
        Real mindist;
        Real xgrid, ygrid, zgrid;
        Real xmax, ymax, zmax;
        Real xmin, ymin, zmin;
        vector< vector<Real> > AtomTypesDensityMatrix, ElementsDensityMatrix;
        CubeStruct cube;
        CubeNum=0;
        MinMax(xmin, ymin, zmin, xmax, ymax, zmax, Atoms);

        //MakeGrid(xmin, ymin, zmin, xmax, ymax, zmax, params.CubeSize, cubes);
        ExpandMinMax(xmin, ymin, zmin, xmax, ymax, zmax, params);

        MaxXBin=int(floor((xmax-xmin)/params.CubeSize+0.5));
        MaxYBin=int(floor((ymax-ymin)/params.CubeSize+0.5));
        MaxZBin=int(floor((zmax-zmin)/params.CubeSize+0.5));

        for (xgrid=xmin;xgrid<=xmax;xgrid+=params.CubeSize)
        {
                cout <<"xgrid= "<<xgrid<<" CubeNum= "<<CubeNum<<endl;
                for (ygrid=ymin;ygrid<=ymax;ygrid+=params.CubeSize)
                {
                        for (zgrid=zmin;zgrid<=zmax;zgrid+=params.CubeSize)
                        {
                                Location=location(xgrid, ygrid, zgrid, nearest, mindist, Atoms, params.maxhs);
                                if (Location==InProtein || Location==InHydrationShell)
                                {
                                        cube.x=xgrid;
                                        cube.y=ygrid;
                                        cube.z=zgrid;
                                        dx=xgrid-Atoms[nearest].x;
                                        dy=ygrid-Atoms[nearest].y;
                                        dz=zgrid-Atoms[nearest].z;

                                        dist=sqrt(dx*dx+dy*dy+dz*dz);
                                        distance=int(floor(dist/params.RecBin+0.5));
                                        if (!params.UniformHydrationShell)
                                        {
                                                GetDensity(cube.density, Atoms[nearest], distance, AtomTypesDensityMatrix, ElementsDensityMatrix, 0, params);
                                                cube.density-=contrast;
                                        }
                                        else
                                        {
                                                if (Location==InProtein) cube.density=-contrast;
                                                else if (Location==InHydrationShell) cube.density=params.hsdensity;
                                                else cout <<"ERROR determining hydration shell density"<<endl;
                                        }
                                        cube.AtomType=Atoms[nearest].AtomType;
                                        cube.IntDist=distance;
                                        SafePushBack(cubes, cube, "cubes");
                                        CubeNum++;
                                }
                        }
                }
        }
}

void SubtractBulkDensity(vector<CubeStruct> &cubes, Real density)
{
        //Subtracts of the bulk solvent density from every cube.  This is
        //done to account for the effect of the bulk solvent.  This is a 
        //result of Babinet's principle.
        int NumCubes=cubes.size();
        for (int j=0;j<NumCubes;j++) cubes[j].density-=density;
}

void SubtractBulkDensity(vector<CubeStruct> &cubes, Real density, Real PrdfBulkDensity)
{
        //Subtracts one density from within the protein and another from the
        //hydration shell.  This is because the pRDFs may have a slightly 
        //higher or lower bulk density.
        int NumCubes=cubes.size();
        for (int j=0;j<NumCubes;j++) 
        {
                if (cubes[j].LocatedIn==InHydrationShell) 
                {
                        cubes[j].density-=PrdfBulkDensity;
                }
                else if (cubes[j].LocatedIn==InProtein) 
                {
                        cubes[j].density-=density;
                }
                else
                {
                        cout <<"cube "<<j<<" should not exist"<<endl;
                        PrintCubeInfo(cubes[j]);
                        exit(EXIT_FAILURE);
                }
        }
}

void RemoveZeroDensityCubes(vector<CubeStruct> &cubes)
{
        //Eliminates cubes with zero electron density.
        int NumCubes=cubes.size();
        vector<CubeStruct> TempCubes;

        for (int j=0;j<NumCubes;j++)
        {
                if (cubes[j].density!=0) SafePushBack(TempCubes, cubes[j], "TempCubes");
        }
        cubes=TempCubes;
}

void SetNonZeroCubesToHydrationShell(vector<CubeStruct> &cubes, Real contrast)
{
        //Sets the location of all cubes with a density greater than negative
        //bulk density to the hydration shell.
        int ncube=cubes.size();

        for (int i=0;i<ncube;i++)
        {
                if (cubes[i].density!=-contrast) cubes[i].LocatedIn=InHydrationShell;
        }
}

void AssignCubeDensities(vector<CubeStruct> &cubes, Real contrast, Real hsdensity, vector<AtomStruct> &Atoms, ParamStruct &params)
{
        //Either calls HyPred to hydrate the protein or reads in hydration
        //density map from a pdb.
        Matrix AtomTypesDensityMatrix, ElementsDensityMatrix, NormalDensityMatrix, NormalElementsDensityMatrix;
        DeleteVector(cubes);
        RemoveHydrationShell(Atoms);
        if (params.AssignCubes!="FromPdb")
        {	
                ReadGofRFile(params.AtomTypesGofRFile, 0, AtomTypesDensityMatrix, params.RecBin);
                ReadGofRFile(params.ElementsGofRFile, 0, ElementsDensityMatrix, params.RecBin);
                ReadGofRFile(params.NormalGofRFile, 0, NormalDensityMatrix, params.RecBin);
                ReadGofRFile(params.NormalElementsGofRFile, 0, NormalElementsDensityMatrix, params.RecBin);
        }
        if (params.AssignCubes=="Location3") AssignCubeDensities1(cubes, contrast, Atoms, params);
        else if (params.AssignCubes=="Location4") 
        {
                BuildHydrationShell(cubes, Atoms, AtomTypesDensityMatrix, ElementsDensityMatrix, NormalDensityMatrix, NormalElementsDensityMatrix, contrast, params.UniformHydrationShell, hsdensity, params);
                //SubtractBulkDensity(cubes, contrast);
                SubtractBulkDensity(cubes, contrast, params.PrdfBulkDensity);
        }
        else if (params.AssignCubes=="FromPdb") 
        {
                ReadPdbCube(contrast, cubes, params);
                SetNonZeroCubesToHydrationShell(cubes, contrast);
        }
        else 
        {
                cout <<"Error. "<<params.AssignCubes<<" is not a valid option for AssignCubes."<<endl;
                exit(EXIT_FAILURE);
        }
        cout <<"hsdensity= "<<hsdensity<<endl;
        if (params.UniformHydrationShell)
        {
                int NumCubes=cubes.size();
                for (int j=0;j<NumCubes;j++)
                {
                        if (cubes[j].LocatedIn==InHydrationShell) cubes[j].density=hsdensity;
                }
        }
        if (contrast!=0) RemoveZeroDensityCubes(cubes);
        //cout <<"CubeNum= "<<CubeNum<<endl;

}

void Hydrate(Real contrast, Real hsdensity, vector<AtomStruct> &Atoms, Real atomr[], ParamStruct params)
{
        //Build the hydration shell using cubes, and then converts the cubes to
        //atoms.  This is useful for p(r) calculations.
        vector<CubeStruct> cubes;
        cout <<"In Hydrate contrast= "<<contrast<<" hsdensity= "<<hsdensity<<endl;
        AssignCubeDensities(cubes, contrast, hsdensity, Atoms, params);
        //DecreaseCubeResolution(cubes, params);
        ConvertCubesToAtoms(cubes, atomr, Atoms, params);
}

void ConvertGrid2(Real xmin, Real ymin, Real zmin, vector<CubeStruct> &cubes, Real ***CubeWeight3D, ParamStruct params)
{
        //Conversts a 1D vector of cubes to a 3D array with the densities.
        int CubeNum=cubes.size();
        int xbin, ybin, zbin;
        int MaxXBin, MaxYBin, MaxZBin;
        Real xmax, ymax, zmax;

        MinMax(xmin, ymin, zmin, xmax, ymax, zmax, cubes);

        MaxXBin=int(floor((xmax-xmin)/params.CubeSize+0.5))+1;
        MaxYBin=int(floor((ymax-ymin)/params.CubeSize+0.5))+1;
        MaxZBin=int(floor((zmax-zmin)/params.CubeSize+0.5))+1;
        cout <<"In ConvertGrid2"<<endl;

        //for (xbin=0; xbin<MaxXBin; xbin++)
        //{
        //        for (ybin=0; ybin<MaxYBin; ybin++)
        //        {
        //                for (zbin=0; zbin<MaxZBin; zbin++)
        //                {
        //                        CubeConversion[xbin][ybin][zbin]=0;
        //                }
        //        }
        //}
        cout <<"Reset CubeConversion"<<endl;
        cout <<"MaxXBin= "<<MaxXBin<<" MaxYBin= "<<MaxYBin<<" MaxZBin= "<<MaxZBin<<endl;
        cout <<"xmin= "<<xmin<<" ymin= "<<ymin<<" zmin= "<<zmin<<endl;
        for (int n=0;n<CubeNum;n++)
        {
                //cout <<"n= "<<n<<endl;
                //PrintCubeInfo(cubes[n]);
                //cout <<endl;
                xbin=int(floor((cubes[n].x-xmin)/params.CubeSize+0.5));
                ybin=int(floor((cubes[n].y-ymin)/params.CubeSize+0.5));
                zbin=int(floor((cubes[n].z-zmin)/params.CubeSize+0.5));
                //cout <<"xbin= "<<xbin<<" ybin= "<<ybin<<" zbin= "<<zbin<<endl;
                //CubeConversion[xbin][ybin][zbin]=n;
                //cout <<"CubeConversion= "<<CubeConversion[xbin][ybin][zbin]<<endl;
                CubeWeight3D[xbin][ybin][zbin]=cubes[n].density;
                //cout <<"CubeWeight3D= "<<CubeWeight3D[xbin][ybin][zbin]<<endl;
        }
        cout <<"Done with convert grid"<<endl;
}

void GeneralMinMax(int Num, vector<Real> Xcoor, vector<Real> Ycoor, vector<Real> Zcoor, Real &xmin, Real &ymin, Real &zmin, Real &xmax, Real &ymax, Real &zmax)
{
        //Returns the min and max of the first Num values of three vectors.

        xmin=Xcoor[0];
        ymin=Ycoor[0];
        zmin=Zcoor[0];
        xmax=Xcoor[0];
        ymax=Ycoor[0];
        zmax=Zcoor[0];
        for (int n=0;n<Num;n++)
        {
                if (Xcoor[n]>xmax) xmax=Xcoor[n];
                if (Xcoor[n]<xmin) xmin=Xcoor[n];
                if (Ycoor[n]>ymax) ymax=Ycoor[n];
                if (Ycoor[n]<ymin) ymin=Ycoor[n];
                if (Zcoor[n]>zmax) zmax=Zcoor[n];
                if (Zcoor[n]<zmin) zmin=Zcoor[n];
        }
}

void DensityMapToPdb(lattice &cubes, string PdbFile)
{
        //Prints out a density map in pdb format.  This should be moved to a
        //header file, since other programs use the same function.
        char CharPdbFile[1000];
        int MaxXBin, MaxYBin, MaxZBin;
        AddIndexToFile(PdbFile);
        MaxXBin=cubes.size();
        MaxYBin=cubes[0].size();
        MaxZBin=cubes[0][0].size();
        strcpy(CharPdbFile, PdbFile.c_str());
        for (int j=0;j<MaxXBin;j++)
        {
                for (int k=0;k<MaxYBin;k++)
                {
                        for (int l=0;l<MaxZBin;l++)
                        {
                                PrintPdbLine(PdbFile, true, cubes[j][k][l].IntDist, "O", "HOH", cubes[j][k][l].AtomType, cubes[j][k][l].x, cubes[j][k][l].y, cubes[j][k][l].z, 1.0, cubes[j][k][l].density);
                        }
                }
        }
        ofstream pdb(CharPdbFile, ios::app);
        pdb <<"END"<<endl;
}

void CubesToPdb(vector<CubeStruct> &cubes, string PdbFile)
{
        //Converts cubes to atoms and then prints the atoms to a pdb file.
        vector<AtomStruct> Atoms;

        ConvertCubesToAtoms(cubes, Atoms);
        PrintPdb(PdbFile, Atoms);
}

CubeStruct AverageCubes(lattice &cubes, int xStart, int yStart, int zStart, int xEnd, int yEnd, int zEnd)
{
        //Finds the average location and density of a portion of a grid of cubes
        Real NumAveraged;
        CubeStruct TempCube;
        //cout <<"In AverageCubes"<<endl;
        NumAveraged=Real((xEnd-xStart+1)*(yEnd-yStart+1)*(zEnd-zStart+1));
        TempCube=cubes[xStart][yStart][zStart];
        TempCube.x=0;
        TempCube.y=0;
        TempCube.z=0;
        TempCube.density=0;
        InitializeCube(TempCube);
        for (int j=xStart;j<=xEnd;j++)
        {
                for (int k=yStart;k<=yEnd;k++)
                {
                        for (int l=zStart;l<=zEnd;l++)
                        {
                                TempCube.x+=cubes[j][k][l].x;
                                TempCube.y+=cubes[j][k][l].y;
                                TempCube.z+=cubes[j][k][l].z;
                                TempCube.density+=cubes[j][k][l].density;	
                        }
                }
        }
        if (NumAveraged!=0)
        {
                TempCube.x/=NumAveraged;
                TempCube.y/=NumAveraged;
                TempCube.z/=NumAveraged;
                TempCube.density/=NumAveraged;
        }
        else 
        {
                cout <<"Warning zero cubes averaged"<<endl;
        }
        //cout <<"Leaving AverageCubes"<<endl;
        return TempCube;
}

void DecreaseCubeResolution(lattice &cubes, Real CubeSize)
{
        //Increases the sizes of cubes, by grouping cubes together.
        lattice TempLattice;
        int xbin, ybin, zbin;
        int xEnd, yEnd, zEnd;
        int MaxXBin, MaxYBin, MaxZBin;
        int MaxXBinTemp, MaxYBinTemp, MaxZBinTemp;
        cout <<"In DecreaseCubeResolution"<<endl;
        MaxXBin=cubes.size();
        MaxYBin=cubes[0].size();
        MaxZBin=cubes[0][0].size();

        MaxXBinTemp=int(Real(MaxXBin)*0.5+0.6);
        MaxYBinTemp=int(Real(MaxYBin)*0.5+0.6);
        MaxZBinTemp=int(Real(MaxZBin)*0.5+0.6);
        InitializeLattice(TempLattice, MaxXBinTemp, MaxYBinTemp, MaxZBinTemp);
        //SetCoordinates(TempLattice, xmin, ymin, zmin, CubeSize*2.0);
        for (int j=0;j<MaxXBin;j+=2)
        {
                for (int k=0;k<MaxYBin;k+=2)
                {
                        for (int l=0;l<MaxZBin;l+=2)
                        {
                                xbin=j/2;
                                ybin=k/2;
                                zbin=l/2;
                                xEnd=j+1;
                                yEnd=k+1;
                                zEnd=l+1;
                                if (xEnd>=MaxXBin) xEnd=MaxXBin-1;
                                if (yEnd>=MaxYBin) yEnd=MaxYBin-1;
                                if (zEnd>=MaxZBin) zEnd=MaxZBin-1;
                                TempLattice[xbin][ybin][zbin]=AverageCubes(cubes, j, k, l, xEnd, yEnd, zEnd);
                                TempLattice[xbin][ybin][zbin].x=cubes[j][k][l].x+CubeSize*0.5;
                                TempLattice[xbin][ybin][zbin].y=cubes[j][k][l].y+CubeSize*0.5;
                                TempLattice[xbin][ybin][zbin].z=cubes[j][k][l].z+CubeSize*0.5;
                        }
                }
        }
        cubes=TempLattice;
}

void DecreaseCubeResolution(vector<CubeStruct> &cubes, ParamStruct &params)
{
        //Increases the sizes of a vector of cubes.
        lattice CubeLattice;
        ConvertCubeVectorToLattice(cubes, CubeLattice, params.CubeSize);
        DecreaseCubeResolution(CubeLattice, params.CubeSize);
        ConvertCubeLatticeToCubeVector(CubeLattice, cubes);
        SafeResize(cubes, "cubes in DecreaseCubeResolution");
        params.CubeSize*=2.0;
}

void cube4d(int natom, Real atomr[], Real numatom[], Real contrast, Real hsdensity, IntensityStruct &i, vector<AtomStruct> &Atoms, ParamStruct params)
{
        //Sets up the calculation for the cube method.
        char CharInputCubeDensityPdbFile[1000];
        char CharOutputCubeDensityPdbFile[1000];
        int points=i.calc.size();
        int TotalParticles=Atoms.size();
        int CubeNum, NumCrossSections;
        int NumCubeLimit=8000000;
        vector< vector< vector<int> > > CubeConversion;
        vector<CubeStruct> cubes;
        Real CrossSection;
        Real xmax, ymax, zmax;
        Real xmin, ymin, zmin;
        string PdbOutputFile; //Initialize This
        //PdbOutputFile="/home/jouko/project/WAXS/pdb/ubq_NoIons.pdb";
        //PrintPdb(PdbOutputFile, Atoms);

        //center(atomm, natom, 0, atomr);

        if (params.MoveProtein) MoveAtoms(Atoms, -params.XOrigin, -params.YOrigin, -params.ZOrigin);
        if (params.CenterAtomsInBox) CenterAtomsInBox(Atoms);
        cout <<"natom= "<<natom<<" TotalParticles= "<<TotalParticles<<endl;

        PrintPdb(PdbOutputFile, Atoms);


        CubeNum=0;

        cout <<"Allocated memory"<<endl;

        strcpy(CharInputCubeDensityPdbFile, params.InputCubeDensityPdbFile.c_str());
        strcpy(CharOutputCubeDensityPdbFile, params.OutputCubeDensityPdbFile.c_str());
        cout <<"CharOutputCubeDensityPdbFile= "<<CharOutputCubeDensityPdbFile<<endl;
        for (int t=0;t<points;t++) i.calc[t]=0;

        AssignCubeDensities(cubes, contrast, hsdensity, Atoms, params);
        cout <<"CubeNum= "<<cubes.size()<<endl;
        if (params.AssignCubes!="FromPdb" && params.OutputCubeDensityPdbFile!="")
        {
                cout <<"In if (params.AssignCubes!='FromPdb')"<<endl;
                char CharCubeDensityPdbFilePre[1000];
                string CubeDensityPdbFilePre;

                CubeDensityPdbFilePre=params.OutputCubeDensityPdbFile+"Pre";
                cout <<"CubeDensityPdbFilePre= "<<CubeDensityPdbFilePre<<endl;
                strcpy(CharCubeDensityPdbFilePre, CubeDensityPdbFilePre.c_str());
                cout <<"After strcpy"<<endl;
                ofstream pdb1(CharCubeDensityPdbFilePre, ios::app);
                pdb1 <<"HEADER    Created by "<<Version<<endl;
                CubeNum=cubes.size();
                cout <<"About to print density map"<<endl;
                for (int n=0;n<CubeNum;n++)
                {
                        PrintPdbLine(pdb1, true, 1, "O", "HOH", 1, cubes[n].x, cubes[n].y, cubes[n].z, 1.0, cubes[n].density);
                }

                pdb1  << "END"<<endl;
                pdb1.close();
                cout <<"Printed density map"<<endl;
                MinMax(xmin, ymin, zmin, xmax, ymax, zmax, cubes);
                cout <<"After MinMax"<<endl;
                if (params.RemoveCavities=="ByConnectivity") EliminateCavities(xmin, ymin, zmin, contrast, cubes, params.CubeSize);
                cout <<"After EliminateCavities"<<endl;
                if (params.RemoveCavities=="ByCutOff") FindCavitiesByCutOff(contrast, cubes);
                //cout <<"CharOutputCubeDensityPdbFile= "<<CharOutputCubeDensityPdbFile<<endl;
                ofstream pdb2(CharOutputCubeDensityPdbFile, ios::app);
                pdb2 <<"HEADER    Created by "<<Version<<endl;
                cout <<"CubeNum= "<<cubes.size()<<endl;
                CubeNum=cubes.size();
                for (int n=0;n<CubeNum;n++)
                {
                        //cout <<"cubes["<<n<<"].IntDist= "<<cubes[n].IntDist<<endl;
                        //cout <<"cubes.AtomType= "<<cubes[n].AtomType<<endl;
                        //cout <<"cubes.density= "<<cubes[n].density<<endl;
                        //cout <<"cubes.x= "<<cubes[n].x<<endl;
                        //cout <<"cubes.y= "<<cubes[n].y<<endl;
                        //cout <<"cubes.z= "<<cubes[n].z<<endl;
                        PrintPdbLine(params.OutputCubeDensityPdbFile, true, cubes[n].IntDist, "O", "HOH", cubes[n].AtomType, cubes[n].x, cubes[n].y, cubes[n].z, 1.0, cubes[n].density);
                        //PrintPdbLine(params.OutputCubeDensityPdbFile, true, 1, "O", "HOH", 1, cubes[n].x, cubes[n].y, cubes[n].z, 1.0, cubes[n].density);
                }

                pdb2  << "END"<<endl;
                string CubeDensityPdbFileB=params.OutputCubeDensityPdbFile+"B";
                cout <<"CubeDensityPdbFileB= "<<CubeDensityPdbFileB<<endl;
                NumCrossSections=int((zmax-zmin)/params.CubeSize+0.5)+1;
                if (isOdd(NumCrossSections)) CrossSection=(zmin+zmax)*0.5;
                else CrossSection=(zmin+zmax)*0.5+0.5*params.CubeSize;
                for (int n=0;n<CubeNum;n++)
                {
                        Real ScaledDensity;
                        if (cubes[n].density>0.0) ScaledDensity=cubes[n].density*150.0+50.0;
                        if (cubes[n].density==0.0) ScaledDensity=50.0;
                        if (cubes[n].density<0.0) ScaledDensity=cubes[n].density*300.0+50.0;
                        //cout <<"zmin= "<<zmin<<" zmax= "<<zmax<<" CrossSection= "<<CrossSection
                        //<<" cubes["<<n<<"].z= "<<cubes[n].z<<endl;
                        Real d=0.1;
                        if (cubes[n].density!=-contrast && abs(cubes[n].z-CrossSection)<d)
                        {
                                PrintPdbLine(CubeDensityPdbFileB, true, 1, "O", "HOH", 1, cubes[n].x, cubes[n].y, cubes[n].z, 1.0, ScaledDensity);
                        }
                }
        }
        FindProteinElectrons(Atoms, cubes, params.CubeSize);

        if (int(cubes.size())>NumCubeLimit)
        {
                DecreaseCubeResolution(cubes, params);
        }

        if (params.BoolFFT) 
        {
                if (params.FFTType=="3D") CubeFFT3D(numatom, contrast, i, cubes, Atoms, params);
                if (params.FFTType=="1D") CubeFFT3(contrast, i, cubes, Atoms, params);
        }
        else CalcCubeIntensity(cubes, i, Atoms, params);
        cout <<"cubes.size= "<<cubes.size()<<endl;
        if (points>100) cout <<"i[100]="<<i.calc[100]<<endl;
        //DecreaseCubeResolution(cubes, params);
        cout <<"cubes.size= "<<cubes.size()<<endl;
        if (params.CalcPrFromStructure || params.CalcPrFromStructure2 || params.CalcRg)
        {
                ConvertCubesToAtoms(cubes, atomr, Atoms, params);
        }
        //PrintPdb("/home/jouko/project/WAXS/test/pdb/DecreasedResolution.pdb", Atoms);
}


void InterpolateFFT(vector< vector< vector<Real> > > &TempOutputReal, vector< vector< vector<Real> > > &TempOutputImag, Real IrealCube[], Real IimagCube[], int points, IntensityStruct &i, int xmax, int ymax, Real dx, Real dy, Real dz, Real inc, int Max, ParamStruct &params, Real cosphi, Real theta)
{
        //FFT calculates the scattering at certain discreate values of s.
        //To obtain the scattering at other values of s interpolation must be
        //done.  This function carries out the interpolation.
        //cosphi and theta are only passed for debugging.
        int u;
        Real st0, st1, dotproduct, InvBinSize;
        Real *TempIrealCube, *TempIimagCube;

        SafeArrayAlloc(TempIrealCube, points, "TempIrealCube");
        SafeArrayAlloc(TempIimagCube, points, "TempIimagCube");

        InvBinSize=-dz*Real(Max)*params.CubeSize/(2.0*pi);
        if (cosphi<-0.64 && cosphi>-0.66 && theta<0.16 && theta>0.14) cout <<"dz= "<<dz<<" Max= "<<Max<<" params.CubeSize= "<<params.CubeSize<<" InvBinSize= "<<InvBinSize<<endl;
        //cout <<"Step size= "<<-(2.0*pi)/(dz*Real(Max)*CubeSize)<<endl;

        for (int xbin=0; xbin<xmax ; xbin++)
        {
                //cout <<"xbin= "<<xbin<<endl;
                for (int ybin=0; ybin<ymax ; ybin++)
                {
                        for (int t=0;t<points;t++)
                        {
                                u=int(floor(i.s[t]*InvBinSize));
                                //if (cosphi<-0.64 && cosphi>-0.66 && theta<0.16 && theta>0.14) cout <<"u= "<<u<<" i.s["<<t<<"]= "<<i.s[t]<<" InvBinSize= "<<InvBinSize<<endl;
                                st0=-Real(u)*2.0*pi/(dz*Real(Max)*params.CubeSize);
                                st1=-Real(u+1)*2.0*pi/(dz*Real(Max)*params.CubeSize);
                                //if (cosphi<-0.64 && cosphi>-0.66 && theta<0.16 && theta>0.14) 
                                //{
                                //	cout <<"st0= "<<st0<<" st1= "<<st1<<" s["<<t<<"]= "<<i.s[t]<<" u= "<<u<<endl;
                                //	cout <<"xbin= "<<xbin<<" ybin= "<<ybin<<endl;
                                //}			

                                if (u+1<Max)
                                {
                                        TempIrealCube[t]=(TempOutputReal[xbin][ybin][u+1]-TempOutputReal[xbin][ybin][u])*(i.s[t]-st0)*InvBinSize+TempOutputReal[xbin][ybin][u];
                                        TempIimagCube[t]=(TempOutputImag[xbin][ybin][u+1]-TempOutputImag[xbin][ybin][u])*(i.s[t]-st0)*InvBinSize+TempOutputImag[xbin][ybin][u];
                                }
                                //if (xbin==35 && ybin==50) cout <<"TempIrealCube["<<t<<"]= "<<TempIrealCube[t]<<" TempOutputReal["<<u+1<<"]= "<<TempOutputReal[u+1]<<" TempOutputReal["<<u<<"]= "<<TempOutputReal[u]<<" s["<<t<<"]= "<<s[t]<<" st0= "<<st0<<" st1= "<<st1<<endl;
                        }

                        for (int t=0;t<points;t++)
                        {
                                dotproduct=(dx*Real(xbin)*params.CubeSize+dy*Real(ybin)*params.CubeSize)*i.s[t];
                                IrealCube[t]+=TempIrealCube[t]*cos(dotproduct)-TempIimagCube[t]*sin(dotproduct);
                                IimagCube[t]+=TempIimagCube[t]*cos(dotproduct)+TempIrealCube[t]*sin(dotproduct);
                        }
                }
        }
}

Real InterpolateFFT(int xbin, int ybin, int zbin, Real dx, Real dy, Real dz, Real ***TempOutputReal, Real CubeSize, Real Max)
{
        //This interpolates the scattering when 3D FFT is used.  This funcion
        //may not be working.  3D FFT does not work.
        vector<Real> variable;
        vector< vector<Real> > coefficient;
        SafeAlloc(variable, 4, "variable");
        CreateMatrix(coefficient, 5, 4);
        //cout <<"First equation"<<endl;	
        coefficient[0][0]=Real(xbin+1)*2.0*pi/(CubeSize*Real(Max));
        coefficient[1][0]=Real(ybin)*2.0*pi/(CubeSize*Real(Max));
        coefficient[2][0]=Real(zbin)*2.0*pi/(CubeSize*Real(Max));
        coefficient[3][0]=1.0;
        coefficient[4][0]=TempOutputReal[xbin+1][ybin][zbin];
        //cout <<"Second equation"<<endl;
        coefficient[0][1]=Real(xbin)*2.0*pi/(CubeSize*Real(Max));
        coefficient[1][1]=Real(ybin+1)*2.0*pi/(CubeSize*Real(Max));
        coefficient[2][1]=Real(zbin)*2.0*pi/(CubeSize*Real(Max));
        coefficient[3][1]=1.0;
        coefficient[4][1]=TempOutputReal[xbin][ybin+1][zbin];
        //cout <<"Third equation"<<endl;
        coefficient[0][2]=Real(xbin)*2.0*pi/(CubeSize*Real(Max));
        coefficient[1][2]=Real(ybin)*2.0*pi/(CubeSize*Real(Max));
        coefficient[2][2]=Real(zbin+1)*2.0*pi/(CubeSize*Real(Max));
        coefficient[3][2]=1.0;
        coefficient[4][2]=TempOutputReal[xbin][ybin][zbin+1];
        //cout <<"About to enter equation"<<endl;	
        coefficient[0][3]=Real(xbin)*2.0*pi/(CubeSize*Real(Max));
        coefficient[1][3]=Real(ybin)*2.0*pi/(CubeSize*Real(Max));
        coefficient[2][3]=Real(zbin)*2.0*pi/(CubeSize*Real(Max));
        coefficient[3][3]=1.0;
        coefficient[4][3]=TempOutputReal[xbin][ybin][zbin];
        //PrintMatrix(coefficient);
        Equation(coefficient, variable);
        //cout <<"Left equation"<<endl;
        //cout <<endl;
        //cout <<"variable[0]= "<<variable[0]<<" variable[1]= "<<variable[1]<<" variable[2]= "<<variable[2]<<endl;
        return variable[0]*dx+variable[1]*dy+variable[2]*dz+variable[3];
}

void InterpolateFFT3D(Real ***TempOutputReal, Real ***TempOutputImag, Real IrealCube[], Real IimagCube[], int points, IntensityStruct i, int xmax, int ymax, Real dx, Real dy, Real dz, Real inc, int Max, ParamStruct params)
{
        //This interpolates the scattering when 3D FFT is used.  This funcion
        //may not be working.  3D FFT does not work.
        int xbin, ybin, zbin;
        Real InvBinSizeX, InvBinSizeY, InvBinSizeZ;
        vector<Real> TempIrealCube, TempIimagCube;

        SafeAlloc(TempIrealCube, points, "TempIrealCube");
        SafeAlloc(TempIimagCube, points, "TempIimagCube");

        InvBinSizeX=abs(dx)*Real(Max)*params.CubeSize/(2.0*pi);
        InvBinSizeY=abs(dy)*Real(Max)*params.CubeSize/(2.0*pi);
        InvBinSizeZ=abs(dz)*Real(Max)*params.CubeSize/(2.0*pi);

        for (int t=0;t<points;t++)
        {
                xbin=int(floor(i.s[t]*InvBinSizeX));
                ybin=int(floor(i.s[t]*InvBinSizeY));
                zbin=int(floor(i.s[t]*InvBinSizeZ));
                if (dx<0 && xbin!=0) xbin=abs(Max-xbin);
                if (dy<0 && ybin!=0) ybin=abs(Max-ybin);
                if (dz<0 && zbin!=0) zbin=abs(Max-zbin);
                //cout <<"t= "<<t<<" xbin= "<<xbin<<" ybin= "<<ybin<<" zbin= "<<zbin<<endl;
                IrealCube[t]=TempOutputReal[xbin][ybin][zbin];
                IimagCube[t]=-TempOutputImag[xbin][ybin][zbin];
                if (t==10)
                {
                        cout <<"TempOutputReal["<<xbin<<"]["<<ybin<<"]["<<zbin<<"]= "<<TempOutputReal[xbin][ybin][zbin]<<endl;
                        cout <<"TempOutputReal["<<ybin<<"]["<<xbin<<"]["<<zbin<<"]= "<<TempOutputReal[ybin][xbin][zbin]<<endl;
                        cout <<"TempOutputReal["<<zbin<<"]["<<xbin<<"]["<<ybin<<"]= "<<TempOutputReal[zbin][xbin][ybin]<<endl;
                        cout <<"TempOutputReal["<<xbin<<"]["<<zbin<<"]["<<ybin<<"]= "<<TempOutputReal[xbin][zbin][ybin]<<endl;
                        cout <<"TempOutputReal["<<ybin<<"]["<<zbin<<"]["<<xbin<<"]= "<<TempOutputReal[ybin][zbin][xbin]<<endl;
                        cout <<"TempOutputReal["<<zbin<<"]["<<ybin<<"]["<<xbin<<"]= "<<TempOutputReal[zbin][ybin][xbin]<<endl;
                }
                //IrealCube[t]=InterpolateFFT(xbin, ybin, zbin, i.s[t]*dx, i.s[t]*dy, i.s[t]*dz, TempOutputReal, params.CubeSize, Real(Max));
                //IimagCube[t]=InterpolateFFT(xbin, ybin, zbin, i.s[t]*dx, i.s[t]*dy, i.s[t]*dz, TempOutputImag, params.CubeSize, Real(Max));
        }

}

bool AllocateMemoryForFFT(Array3D &TempOutputImagX, Array3D &TempOutputRealX, Array3D &TempOutputImagY, Array3D &TempOutputRealY, Array3D &TempOutputImagZ, Array3D &TempOutputRealZ, int xmax, int ymax, int zmax, int Max)
{
        //Attempts to allocate memory for FFT.  If memory allocation fails a
        //more memory efficient algorithm is tried.
        bool MemoryAllocated=true;
        cout <<"In AllocateMemoryForFFT"<<endl;
        cout <<"xmax= "<<xmax<<" ymax= "<<ymax<<" zmax= "<<zmax<<" Max= "<<Max<<endl;
        cout <<"In try"<<endl;
        Safe3DAlloc(TempOutputImagX, ymax, zmax, Max, "TempOutputImagX");	
        Safe3DAlloc(TempOutputRealX, ymax, zmax, Max, "TempOutputRealX");	
        cout <<"Allocated TempOutputRealX"<<endl;
        Safe3DAlloc(TempOutputImagY, xmax, zmax, Max, "TempOutputImagY");	
        Safe3DAlloc(TempOutputRealY, xmax, zmax, Max, "TempOutputRealY");	
        cout <<"Allocated TempOutputRealY"<<endl;
        Safe3DAlloc(TempOutputImagZ, xmax, ymax, Max, "TempOutputImagZ");	
        MemoryAllocated=Safe3DAlloc(TempOutputRealZ, xmax, ymax, Max);	
        cout <<"Allocated Z"<<endl;
        cout <<"xmax= "<<xmax<<" ymax= "<<ymax<<" zmax= "<<zmax<<" Max= "<<Max<<endl;
        Print3DVectorSize(TempOutputImagX, "TempOutputImagX");
        Print3DVectorSize(TempOutputRealX, "TempOutputRealX");
        Print3DVectorSize(TempOutputImagY, "TempOutputImagY");
        Print3DVectorSize(TempOutputRealY, "TempOutputRealY");
        Print3DVectorSize(TempOutputImagZ, "TempOutputImagZ");
        Print3DVectorSize(TempOutputRealZ, "TempOutputRealZ");
        cout <<"MemoryAllocated= "<<MemoryAllocated<<endl;
        return MemoryAllocated;
}

void CubeFFT(lattice &CubeWeight3D, Array3D &TempOutputImagX, Array3D &TempOutputRealX, Array3D &TempOutputImagY, Array3D &TempOutputRealY, Array3D &TempOutputImagZ, Array3D &TempOutputRealZ, int xmax, int ymax, int zmax, int Max, Real contrast, IntensityStruct &i, vector<CubeStruct> cubes, vector<AtomStruct> &Atoms, ParamStruct params)
{
        //Calculates scattering from the cubes using FFT.
        int TotalParticles=Atoms.size();
        int CubeNum=cubes.size();
        int xbin, ybin, zbin;
        int points=i.calc.size();
        Real dx, dy, dz;
        Real mux, muy, muz;
        Real *ireal, *iimag;
        Real *InputReal, *InputImag, *OutputReal, *OutputImag;
        Real *IrealCube, *IimagCube;
        Real dotproduct, product;
        Real **f;
        Real inc;
        Real NumVectors;
        Real volume;
        Real thetainc;
        Real theta;
        Real nslices, cosphi, sinphi;

        SafeArrayAlloc(InputReal, Max, "InputReal");
        SafeArrayAlloc(InputImag, Max, "InputImag");
        SafeArrayAlloc(OutputReal, Max, "OutputReal");
        SafeArrayAlloc(OutputImag, Max, "OutputImag");
        SafeArrayAlloc(ireal, points, "ireal");
        SafeArrayAlloc(iimag, points, "iimag");
        SafeArrayAlloc(IrealCube, points, "IrealCube");
        SafeArrayAlloc(IimagCube, points, "IimagCube");

        for (zbin=0;zbin<Max;zbin++)
        {
                InputReal[zbin]=0.0;
                InputImag[zbin]=0.0;
                OutputReal[zbin]=0.0;
                OutputImag[zbin]=0.0;
        }

        cout <<"TotalParticles= "<<TotalParticles<<endl;

        for (int t=0;t<points;t++) i.calc[t]=0;

        cout <<"CubeNum="<<CubeNum<<endl;

        volume=0;
        Safe2DArrayAlloc(f, points, NumAtomTypes, "f");

        for (int t=0;t<points;t++)
        {
                for (int n=0;n<NumAtomTypes;n++)
                {
                        f[t][n]=0;
                }
        }

        for (int n=0;n<NumAtomTypes;n++)
        {
                for (int t=0;t<points;t++)
                {
                        f[t][n]=SolventCorrectedScattering(n, 0, i.s[t], i.f[t][n], 0, 0, params);
                }
                cout <<"f[0]["<<n<<"]= "<<f[0][n]<<endl;
        }

        cout <<"f[0][0]="<<f[0][0]<<endl;
        NumVectors=0.0;
        nslices=params.VectorsPerInclination;
        cosphi=1.0/nslices-1.0;
        thetainc=params.VectorsPerInclination;
        theta=pi/thetainc;

        //cout <<"About to average"<<endl;

        for (ybin=0; ybin<ymax ; ybin++)
        {
                //cout <<"ybin= "<<ybin<<endl;
                for (zbin=0; zbin<zmax ; zbin++)
                {
                        for (xbin=0; xbin<xmax ; xbin++)
                        {
                                //cout <<"CubeWeight3D["<<xbin<<"]["<<ybin<<"]["<<zbin<<"]= "<<endl;
                                InputReal[xbin]=CubeWeight3D[xbin][ybin][zbin].density;
                                //if (ybin==35 && zbin==50 )cout <<"InputReal["<<xbin<<"]= "<<InputReal[xbin]<<endl;
                        }
                        for (xbin=xmax; xbin<Max; xbin++) InputReal[xbin]=0.0;
                        FFT(InputReal, InputImag, OutputReal, OutputImag, Max);
                        //cout <<"Finnished with FFT"<<endl;
                        for (int t=0;t<Max;t++)
                        {
                                //if (ybin==0) cout <<"OutputReal["<<t<<"]= "<<OutputReal[t]<<endl;
                                //if (ybin==0) cout <<"OutputImag["<<t<<"]= "<<OutputImag[t]<<endl;
                                //if (ybin==35 && zbin==50) cout <<"OutputReal["<<t<<"]= "<<OutputReal[t]<<" OutputImag["<<t<<"]= "<<OutputImag[t]<<endl;
                                TempOutputRealX[ybin][zbin][t]=OutputReal[t];
                                TempOutputImagX[ybin][zbin][t]=OutputImag[t];
                        }
                }
        }

        for (xbin=0; xbin<xmax ; xbin++)
        {
                for (zbin=0; zbin<zmax ; zbin++)
                {
                        for (ybin=0; ybin<ymax ; ybin++)
                        {
                                InputReal[ybin]=CubeWeight3D[xbin][ybin][zbin].density;
                        }
                        for (ybin=ymax; ybin<Max; ybin++) InputReal[ybin]=0.0;
                        FFT(InputReal, InputImag, OutputReal, OutputImag, Max);
                        for (int t=0;t<Max;t++)
                        {
                                TempOutputRealY[xbin][zbin][t]=OutputReal[t];
                                TempOutputImagY[xbin][zbin][t]=OutputImag[t];
                        }
                }
        }
        Real sum=0;
        for (xbin=0; xbin<xmax ; xbin++)
        {
                for (ybin=0; ybin<ymax ; ybin++)
                {
                        for (zbin=0; zbin<zmax ; zbin++)
                        {
                                InputReal[zbin]=CubeWeight3D[xbin][ybin][zbin].density;
                        }
                        for (zbin=zmax; zbin<Max; zbin++) InputReal[zbin]=0.0;
                        FFT(InputReal, InputImag, OutputReal, OutputImag, Max);
                        sum+=OutputReal[0];
                        cout <<"xbin= "<<xbin<<" ybin= "<<ybin<<" OutputReal[0]= "<<OutputReal[0]<<endl;
                        for (int t=0;t<Max;t++)
                        {
                                TempOutputRealZ[xbin][ybin][t]=OutputReal[t];
                                TempOutputImagZ[xbin][ybin][t]=OutputImag[t];
                        }
                }
        }
        cout <<"sum= "<<sum<<endl;
        cout <<"Finnished with FFTs"<<endl;

        while (true)
        {
                NumVectors+=1.0;
                cout <<"cosphi="<<cosphi<<" theta="<<theta<<endl;
                sinphi=sqrt(1-cosphi*cosphi);
                //cosphi=-0.65;
                //theta=0.15708;
                dx=sinphi*cos(theta);
                dy=sinphi*sin(theta);
                dz=cosphi;

                Real DX=abs(dx);
                Real DY=abs(dy);
                Real DZ=abs(dz);

                if (DX>DY && DX>DZ)
                {
                        if (dx>0)
                        {
                                dx=-dx;
                                dy=-dy;
                                dz=-dz;
                        }
                }
                else if (DY>DX && DY>DZ)
                {
                        if (dy>0)
                        {
                                dx=-dx;
                                dy=-dy;
                                dz=-dz;
                        }
                }
                else if (DZ>DX && DZ>DY)
                {
                        if (dz>0)
                        {
                                dx=-dx;
                                dy=-dy;
                                dz=-dz;
                        }
                }

                for (int t=0;t<points;t++)
                {
                        ireal[t]=0;
                        iimag[t]=0;
                        IrealCube[t]=0;
                        IimagCube[t]=0;
                        mux=i.s[t]*dx;
                        muy=i.s[t]*dy;
                        muz=i.s[t]*dz;
                        //cout <<"s["<<t<<"]= "<<s[t]<<" dx= "<<dx<<" dy= "<<dy<<" dz= "<<dz<<endl;
                        for (int n=0;n<TotalParticles;n++)
                        {
                                //cout <<"x["<<n<<"]= "<<Atoms[n].x<<" y= "<<y[n]<<" z= "<<z[n]<<endl;
                                dotproduct=mux*Atoms[n].x+muy*Atoms[n].y+muz*Atoms[n].z;
                                ireal[t]+=f[t][Atoms[n].atomid]*cos(dotproduct);
                                iimag[t]+=f[t][Atoms[n].atomid]*sin(dotproduct);
                                if (n<10)
                                {
                                        //cout <<"ireal["<<t<<"]= "<<ireal[t]<<" f["<<t<<"]["<<atomid[n]<<"]= "<<f[t][atomid[n]]<<endl;
                                        //cout <<"cos("<<dotproduct<<")= "<<cos(dotproduct)<<" sin("<<dotproduct<<")= "<<sin(dotproduct)<<endl;
                                        //cout <<"iimag["<<t<<"]= "<<iimag[t]<<endl;
                                }
                        }
                        if (t==10)
                        {
                                cout <<"mux= "<<mux<<" muy= "<<muy<<" muz= "<<muz<<endl;
                        }
                }
                cout <<"ireal[10]= "<<ireal[10]<<" iimag[10]= "<<iimag[10]<<endl;
                cout <<"dx= "<<dx<<" dy= "<<dy<<" dz= "<<dz<<endl;
                cout <<"f["<<10<<"][0]= "<<f[10][0]<<endl;
                if (dx<dy && dx<dz)
                {
                        InterpolateFFT(TempOutputRealX, TempOutputImagX, IrealCube, IimagCube, points, i, ymax, zmax, dy, dz, dx, inc, Max, params, cosphi, theta);
                }
                else if (dy<dx && dy<dz)
                {
                        InterpolateFFT(TempOutputRealY, TempOutputImagY, IrealCube, IimagCube, points, i, xmax, zmax, dx, dz, dy, inc, Max, params, cosphi, theta);
                }
                else if (dz<dx && dz<dy)
                {
                        InterpolateFFT(TempOutputRealZ, TempOutputImagZ, IrealCube, IimagCube, points, i, xmax, ymax, dx, dy, dz, inc, Max, params, cosphi, theta);
                }
                else
                {
                        cout <<"Error in CubeFFT3"<<endl;
                }


                for (int t=0;t<points;t++)
                {
                        mux=i.s[t]*dx;
                        muy=i.s[t]*dy;
                        muz=i.s[t]*dz;
                        product=8.0*sin(mux*params.CubeSize*0.5)*sin(muy*params.CubeSize*0.5)*sin(muz*params.CubeSize*0.5)/(mux*muy*muz);
                        if (t==10)
                        {
                                cout <<"IrealCube["<<t<<"]= "<<IrealCube[t]<<endl;
                                cout <<"IimagCube["<<t<<"]= "<<IimagCube[t]<<endl;
                        }
                        ireal[t]=ireal[t]+IrealCube[t]*product;
                        iimag[t]=iimag[t]+IimagCube[t]*product;
                        i.calc[t]=i.calc[t]+ireal[t]*ireal[t]+iimag[t]*iimag[t];
                }

                theta += 2.0*pi/thetainc;
                if (theta>2.0*pi)
                {
                        theta=pi/thetainc;
                        cosphi += 2.0/nslices;
                        if (cosphi>0.0-0.9/nslices)
                        {
                                break;
                        }
                }
        }

        if (points>100) cout <<"i[100]="<<i.calc[100]<<endl;

        for (int n=0;n<points;n++) i.calc[n]=i.calc[n]*0.5/NumVectors;
        cout <<"cubes.size()= "<<cubes.size()<<endl;
        //DensityMapToPdb(CubeWeight3D, "/home/jouko/CubeWeight3D.pdb");
        //ConvertCubeLatticeToCubeVector(CubeWeight3D, cubes); This should be uncommented if cubes is changed to &cubes
        cout <<"cubes.size()= "<<cubes.size()<<endl;
        cout <<"Done with ConverCubeLatticeToCubeVector"<<endl;

        delete [] InputReal;
        delete [] InputImag;
        delete [] OutputReal;
        delete [] OutputImag;
        delete [] ireal;
        delete [] iimag;
        delete [] IrealCube;
        delete [] IimagCube;

        for (int t=0;t<points;t++)
        {
                delete [] f[t];
        }
        delete [] f;
}

void InterpolateFFT(vector< vector<Real> > &ireal, vector< vector<Real> > &iimag, vector<VectorStruct> &v, IntensityStruct &i, vector< vector<Real> > &IrealCube, vector< vector<Real> > &IimagCube, int Max, Real CubeSize)
{
        //FFT only calculates the scattering at certain discreate points.  
        //Inorder to calculate scattering at other points interpolation is used.
        int u, NumVectors, points;
        Real dz, sk0, sk1, InvBinSize;
        NumVectors=v.size();
        points=i.s.size();
        for (int j=0;j<NumVectors;j++)
        {
                if (abs(v[j].z)>abs(v[j].x) && abs(v[j].z)>abs(v[j].y)) dz=v[j].z;
                if (abs(v[j].x)>abs(v[j].y) && abs(v[j].x)>abs(v[j].z)) dz=v[j].x;
                if (abs(v[j].y)>abs(v[j].x) && abs(v[j].y)>abs(v[j].z)) dz=v[j].y;
                if (dz>0) dz=-dz;
                InvBinSize=-dz*Real(Max)*CubeSize/(2.0*pi);
                for (int k=0;k<points;k++)
                {
                        u=int(floor(i.s[k]*InvBinSize));
                        sk0=-Real(u)*2.0*pi/(dz*Real(Max)*CubeSize);
                        sk1=-Real(u+1)*2.0*pi/(dz*Real(Max)*CubeSize);

                        if (u+1<Max)
                        {
                                if (false)
                                {
                                        cout <<"IrealCube["<<j<<"]["<<k<<"]= "<<IrealCube[j][k]<<endl;
                                        cout <<"ireal["<<j<<"]["<<u+1<<"]= "<<ireal[j][u+1]<<endl;
                                        cout <<"IimagCube["<<j<<"]["<<k<<"]= "<<IimagCube[j][k]<<endl;
                                        cout <<"iimag["<<j<<"]["<<u+1<<"]= "<<iimag[j][u+1]<<endl;
                                        cout <<"i.s["<<k<<"]= "<<i.s[k]<<" sk0= "<<sk0<<" InvBinSize= "<<InvBinSize<<" dz= "<<dz<<endl;
                                }
                                IrealCube[j][k]=(ireal[j][u+1]-ireal[j][u])*(i.s[k]-sk0)*InvBinSize+ireal[j][u];
                                IimagCube[j][k]=(iimag[j][u+1]-iimag[j][u])*(i.s[k]-sk0)*InvBinSize+iimag[j][u];
                        }
                }
        }
}

bool HasNonZero(Real a[], int Max)
{
        //Return true if one of the first Max elements of an array is not zero.
        //Returns false otherwise.
        for (int j=0;j<Max;j++)
        {
                if (a[j]!=0) return true;
        }
        return false;
}

void X1DFFT(int Max, int xmax, int ymax, int zmax, lattice &CubeWeight3D, vector<VectorStruct> &v, ParamStruct &params, vector< vector<Real> > &ireal, vector< vector<Real> > &iimag)
{
        //In order to calculate scattering FFT is performed on a series of 1D
        //vectors.  This function performs FFT on one dimensional vectors of
        //cubes alligned in the x direction.
        int NumVectors=v.size();
        Real s, ds, product, dotproduct, TempIReal, TempIImag;
        Real mux, muy, muz;
        Real *InputReal, *InputImag, *OutputReal, *OutputImag;
        cout <<"In X1DFFT"<<endl;
        SafeArrayAlloc(InputReal, Max, "InputReal");
        SafeArrayAlloc(InputImag, Max, "InputImag");
        SafeArrayAlloc(OutputReal, Max, "OutputReal");
        SafeArrayAlloc(OutputImag, Max, "OutputImag");
        cout <<"Allocated Memory"<<endl;
        for (int ybin=0;ybin<ymax;ybin++)
        {
                for (int zbin=0;zbin<zmax;zbin++)
                {
                        for (int xbin=0;xbin<xmax;xbin++)
                        {
                                InputReal[xbin]=CubeWeight3D[xbin][ybin][zbin].density;
                        }
                        for (int xbin=xmax;xbin<Max;xbin++) InputReal[xbin]=0.0;
                        if (HasNonZero(InputReal, xmax))
                        {
                                FFT(InputReal, InputImag, OutputReal, OutputImag, Max);
                                //cout <<"ybin= "<<ybin<<" zbin= "<<zbin<<endl;
                                for (int j=0;j<NumVectors;j++)
                                {
                                        if (abs(v[j].x)>abs(v[j].y) && abs(v[j].x)>abs(v[j].z))
                                        {
                                                //cout <<"Adding to intensity X v["<<j<<"].x= "<<v[j].x<<" v.y= "<<v[j].y<<" v.z= "<<v[j].z<<endl;
                                                if (v[j].x>0)
                                                {
                                                        v[j].x=-v[j].x;
                                                        v[j].y=-v[j].y;
                                                        v[j].z=-v[j].z;
                                                }
                                                for (int k=0;k<Max;k++)
                                                {
                                                        ds=-2.0*pi/(v[j].x*Real(Max)*params.CubeSize);
                                                        s=ds*Real(k);
                                                        mux=s*v[j].x;
                                                        muy=s*v[j].y;
                                                        muz=s*v[j].z;
                                                        if (s!=0)
                                                        {
                                                                product=8.0*sin(mux*params.CubeSize*0.5)*sin(muy*params.CubeSize*0.5)*sin(muz*params.CubeSize*0.5)/(mux*muy*muz);
                                                        }
                                                        else
                                                        {
                                                                product=params.CubeSize*params.CubeSize*params.CubeSize;
                                                        }
                                                        dotproduct=muy*Real(ybin)*params.CubeSize+muz*Real(zbin)*params.CubeSize;
                                                        TempIReal=OutputReal[k]*cos(dotproduct)-OutputImag[k]*sin(dotproduct);
                                                        TempIImag=OutputImag[k]*cos(dotproduct)+OutputReal[k]*sin(dotproduct);
                                                        /* 
                                                           if (k==0)
                                                           {
                                                           cout <<"OutputReal["<<k<<"]= "<<OutputReal[k]<<" OuputImag["<<k<<"]= "<<OutputImag[k]<<endl;
                                                           cout <<"dotproduct= "<<dotproduct<<" cos= "<<cos(dotproduct)<<" sin= "<<sin(dotproduct)<<endl;
                                                           cout <<"TempIReal= "<<TempIReal<<" TempIImag= "<<TempIImag<<endl;
                                                           cout <<"ireal["<<j<<"]["<<k<<"]= "<<ireal[j][k]<<" iimag["<<j<<"]["<<k<<"]= "<<iimag[j][k]<<endl;
                                                           }
                                                           */
                                                        ireal[j][k]+=TempIReal*product;
                                                        iimag[j][k]+=TempIImag*product;
                                                }
                                        }
                                }
                        }
                }
        }
        delete [] InputReal;
        delete [] InputImag;
        delete [] OutputReal;
        delete [] OutputImag;
}

void Y1DFFT(int Max, int xmax, int ymax, int zmax, lattice &CubeWeight3D, vector<VectorStruct> &v, ParamStruct &params, vector< vector<Real> > &ireal, vector< vector<Real> > &iimag)
{
        //In order to calculate scattering FFT is performed on a series of 1D
        //vectors.  This function performs FFT on one dimensional vectors of
        //cubes alligned in the y direction.
        int NumVectors=v.size();
        Real s, ds, product, dotproduct, TempIReal, TempIImag;
        Real mux, muy, muz;
        Real *InputReal, *InputImag, *OutputReal, *OutputImag;
        cout <<"In Y1DFFT"<<endl;
        SafeArrayAlloc(InputReal, Max, "InputReal");
        SafeArrayAlloc(InputImag, Max, "InputImag");
        SafeArrayAlloc(OutputReal, Max, "OutputReal");
        SafeArrayAlloc(OutputImag, Max, "OutputImag");
        cout <<"Allocated Memory"<<endl;
        for (int xbin=0;xbin<xmax;xbin++)
        {
                for (int zbin=0;zbin<zmax;zbin++)
                {
                        for (int ybin=0;ybin<ymax;ybin++)
                        {
                                InputReal[ybin]=CubeWeight3D[xbin][ybin][zbin].density;
                        }
                        for (int ybin=ymax;ybin<Max;ybin++) InputReal[ybin]=0.0;
                        if (HasNonZero(InputReal, ymax))
                        {
                                FFT(InputReal, InputImag, OutputReal, OutputImag, Max);
                                for (int j=0;j<NumVectors;j++)
                                {
                                        if (abs(v[j].y)>=abs(v[j].x) && abs(v[j].y)>=abs(v[j].z))
                                        {
                                                //cout <<"Adding to intensity Y v["<<j<<"].x= "<<v[j].x<<" v.y= "<<v[j].y<<" v.z= "<<v[j].z<<endl;
                                                if (v[j].y>0)
                                                {
                                                        v[j].x=-v[j].x;
                                                        v[j].y=-v[j].y;
                                                        v[j].z=-v[j].z;
                                                }
                                                for (int k=0;k<Max;k++)
                                                {
                                                        ds=-2.0*pi/(v[j].y*Real(Max)*params.CubeSize);
                                                        s=ds*Real(k);
                                                        mux=s*v[j].x;
                                                        muy=s*v[j].y;
                                                        muz=s*v[j].z;
                                                        if (s!=0)
                                                        {
                                                                product=8.0*sin(mux*params.CubeSize*0.5)*sin(muy*params.CubeSize*0.5)*sin(muz*params.CubeSize*0.5)/(mux*muy*muz);
                                                        }
                                                        else
                                                        {
                                                                product=params.CubeSize*params.CubeSize*params.CubeSize;
                                                        }
                                                        dotproduct=mux*Real(xbin)*params.CubeSize+muz*Real(zbin)*params.CubeSize;
                                                        TempIReal=OutputReal[k]*cos(dotproduct)-OutputImag[k]*sin(dotproduct);
                                                        TempIImag=OutputImag[k]*cos(dotproduct)+OutputReal[k]*sin(dotproduct);
                                                        if (k==0)
                                                        {
                                                                //cout <<"mux= "<<mux<<" muy= "<<muy<<endl;
                                                                //cout <<"TempIReal= "<<TempIReal<<" TempIImag= "<<TempIImag<<endl;
                                                        }
                                                        if (k==0)
                                                        {
                                                                //cout <<"OutputReal["<<k<<"]= "<<OutputReal[k]<<" OuputImag["<<k<<"]= "<<OutputImag[k]<<endl;
                                                                //cout <<"dotproduct= "<<dotproduct<<" cos= "<<cos(dotproduct)<<" sin= "<<sin(dotproduct)<<endl;
                                                                //cout <<"TempIReal= "<<TempIReal<<" TempIImag= "<<TempIImag<<endl;
                                                                //cout <<"ireal["<<j<<"]["<<k<<"]= "<<ireal[j][k]<<" iimag["<<j<<"]["<<k<<"]= "<<iimag[j][k]<<endl;
                                                        }
                                                        ireal[j][k]+=TempIReal*product;
                                                        iimag[j][k]+=TempIImag*product;
                                                }
                                        }
                                }
                        }
                }
        }
        delete [] InputReal;
        delete [] InputImag;
        delete [] OutputReal;
        delete [] OutputImag;
}

void Z1DFFT(int Max, int xmax, int ymax, int zmax, lattice &CubeWeight3D, vector<VectorStruct> &v, ParamStruct &params, vector< vector<Real> > &ireal, vector< vector<Real> > &iimag)
{
        //In order to calculate scattering FFT is performed on a series of 1D
        //vectors.  This function performs FFT on one dimensional vectors of
        //cubes alligned in the z direction.
        int NumVectors=v.size();
        Real s, ds, product, dotproduct, TempIReal, TempIImag;
        Real mux, muy, muz;
        Real *InputReal, *InputImag, *OutputReal, *OutputImag;
        cout <<"In Z1DFFT"<<endl;
        SafeArrayAlloc(InputReal, Max, "InputReal");
        SafeArrayAlloc(InputImag, Max, "InputImag");
        SafeArrayAlloc(OutputReal, Max, "OutputReal");
        SafeArrayAlloc(OutputImag, Max, "OutputImag");
        cout <<"Allocated Memory for Z1DFFT"<<endl;
        for (int xbin=0;xbin<xmax;xbin++)
        {
                cout <<"xbin= "<<xbin<<endl;
                for (int ybin=0;ybin<ymax;ybin++)
                {
                        //cout <<"ybin= "<<ybin<<endl;
                        for (int zbin=0;zbin<zmax;zbin++)
                        {
                                InputReal[zbin]=CubeWeight3D[xbin][ybin][zbin].density;
                        }
                        for (int zbin=zmax;zbin<Max;zbin++) 
                        {
                                InputReal[zbin]=0.0;
                        }
                        if (HasNonZero(InputReal, zmax))
                        {
                                FFT(InputReal, InputImag, OutputReal, OutputImag, Max);
                                for (int j=0;j<NumVectors;j++)
                                {
                                        if (abs(v[j].z)>abs(v[j].x) && abs(v[j].z)>abs(v[j].y))
                                        {
                                                //cout <<"Adding to intensity Z v["<<j<<"].x= "<<v[j].x<<" v.y= "<<v[j].y<<" v.z= "<<v[j].z<<endl;
                                                for (int k=0;k<Max;k++)
                                                {
                                                        ds=-2.0*pi/(v[j].z*Real(Max)*params.CubeSize);
                                                        s=ds*Real(k);
                                                        mux=s*v[j].x;
                                                        muy=s*v[j].y;
                                                        muz=s*v[j].z;
                                                        if (s!=0)
                                                        {
                                                                product=8.0*sin(mux*params.CubeSize*0.5)*sin(muy*params.CubeSize*0.5)*sin(muz*params.CubeSize*0.5)/(mux*muy*muz);
                                                        }
                                                        else
                                                        {
                                                                product=params.CubeSize*params.CubeSize*params.CubeSize;
                                                        }
                                                        dotproduct=mux*Real(xbin)*params.CubeSize+muy*Real(ybin)*params.CubeSize;
                                                        TempIReal=OutputReal[k]*cos(dotproduct)-OutputImag[k]*sin(dotproduct);
                                                        TempIImag=OutputImag[k]*cos(dotproduct)+OutputReal[k]*sin(dotproduct);
                                                        if (k==0)
                                                        {
                                                                //cout <<"mux= "<<mux<<" muy= "<<muy<<endl;
                                                                //cout <<"TempIReal= "<<TempIReal<<" TempIImag= "<<TempIImag<<endl;
                                                        }
                                                        ireal[j][k]+=TempIReal*product;
                                                        iimag[j][k]+=TempIImag*product;
                                                }
                                        }
                                }
                        }
                }
        }
        delete [] InputReal;
        delete [] InputImag;
        delete [] OutputReal;
        delete [] OutputImag;
}

void MemoryEfficientCubeFFT(lattice &CubeWeight3D, int xmax, int ymax, int zmax, int Max, Real contrast, int points, IntensityStruct &i, vector<CubeStruct> &cubes, vector<AtomStruct> &Atoms, ParamStruct &params)
{
        //Calculates the scattering from the cubes using FFT.  This function
        //is more memory efficient than the other FFT functions, because
        //it performs FFT on a 1D vector of cubes, calculates the contribution
        //to the scattering from the 1D vector of cubes and then moves onto 
        //the next 1D vector instead of calculating FFT on all 1D vectors at 
        //once or allocating memory for 3D FFT.
        int TotalParticles=Atoms.size();
        int CubeNum=cubes.size();
        int NumVectors;
        Real mux, muy, muz;
        Real dotproduct;
        Real **f;
        Real thetainc;
        Real nslices;
        Real RealPart, ImagPart;
        vector< vector<Real> > ireal, iimag;
        vector< vector<Real> > IrealCube, IimagCube;
        vector< vector<Real> > IrealAtom, IimagAtom;
        vector<VectorStruct> v;
        //Max=2048;
        Max=4096;


        cout <<"TotalParticles= "<<TotalParticles<<endl;

        for (int t=0;t<points;t++) i.calc[t]=0;

        cout <<"CubeNum="<<CubeNum<<endl;

        Safe2DArrayAlloc(f, points, NumAtomTypes, "f");

        for (int t=0;t<points;t++)
        {
                for (int n=0;n<NumAtomTypes;n++)
                {
                        f[t][n]=0;
                }
        }

        for (int n=0;n<NumAtomTypes;n++)
        {
                for (int t=0;t<points;t++)
                {
                        f[t][n]=SolventCorrectedScattering(n, 0, i.s[t], i.f[t][n], 0, 0, params);
                }
                cout <<"f[0]["<<n<<"]= "<<f[0][n]<<endl;
        }

        cout <<"f[0][0]="<<f[0][0]<<endl;
        nslices=params.VectorsPerInclination;
        thetainc=params.VectorsPerInclination;
        MakePointsOnAHemiSphere(thetainc, nslices, v);
        NumVectors=v.size();
        Safe2DAlloc(ireal, NumVectors, Max, "ireal in MemoryEfficientCubeFFT");
        Safe2DAlloc(iimag, NumVectors, Max, "iimag in MemoryEfficientCubeFFT");
        Safe2DAlloc(IrealCube, NumVectors, Max, "IrealCube in MemoryEfficientCubeFFT");
        Safe2DAlloc(IimagCube, NumVectors, Max, "IimagCube in MemoryEfficientCubeFFT");
        Safe2DAlloc(IrealAtom, NumVectors, Max, "IrealAtom in MemoryEfficientCubeFFT");
        Safe2DAlloc(IimagAtom, NumVectors, Max, "IimagAtom in MemoryEfficientCubeFFT");

        cout <<"About to enter X1DFFT"<<endl;
        X1DFFT(Max, xmax, ymax, zmax, CubeWeight3D, v, params, ireal, iimag);
        cout <<"About to enter Y1DFFT"<<endl;
        Y1DFFT(Max, xmax, ymax, zmax, CubeWeight3D, v, params, ireal, iimag);
        cout <<"About to enter Z1DFFT"<<endl;
        Z1DFFT(Max, xmax, ymax, zmax, CubeWeight3D, v, params, ireal, iimag);
        cout <<"Left Z1DFFT"<<endl;
        for (int j=0;j<NumVectors;j++)
        {
                //cout <<"ireal["<<j<<"][0]= "<<ireal[j][0]<<" iimag["<<j<<"][0]= "<<iimag[j][0]<<" v.x= "<<v[j].x<<" v.y= "<<v[j].y<<" v.z= "<<v[j].z<<endl;
        }
        InterpolateFFT(ireal, iimag, v, i, IrealCube, IimagCube, Max, params.CubeSize);
        for (int j=0;j<NumVectors;j++)
        {
                //cout <<"IrealCube["<<j<<"][0]= "<<IrealCube[j][0]<<" IimagCube["<<j<<"][0]= "<<IimagCube[j][0]<<endl;
        }

        for (int j=0;j<NumVectors;j++)
        {
                for (int k=0;k<points;k++)
                {
                        mux=v[j].x*i.s[k];
                        muy=v[j].y*i.s[k];
                        muz=v[j].z*i.s[k];
                        for (int n=0;n<TotalParticles;n++)
                        {
                                dotproduct=mux*Atoms[n].x+muy*Atoms[n].y+muz*Atoms[n].z;
                                IrealAtom[j][k]+=f[k][Atoms[n].atomid]*cos(dotproduct);
                                IimagAtom[j][k]+=f[k][Atoms[n].atomid]*sin(dotproduct);
                        }
                }
        }

        for (int j=0;j<NumVectors;j++)
        {
                for (int k=0;k<points;k++)
                {
                        RealPart=IrealAtom[j][k]+IrealCube[j][k];
                        ImagPart=IimagAtom[j][k]+IimagCube[j][k];
                        i.calc[k]+=RealPart*RealPart+ImagPart*ImagPart;
                }
        }

        if (points>100) cout <<"i[100]="<<i.calc[100]<<endl;

        for (int n=0;n<points;n++) i.calc[n]=i.calc[n]*4.0*pi/Real(NumVectors);
        cout <<"cubes.size()= "<<cubes.size()<<endl;
        //DensityMapToPdb(CubeWeight3D, "/home/jouko/CubeWeight3D.pdb");
        //ConvertCubeLatticeToCubeVector(CubeWeight3D, cubes); This should be uncommented if cubes is changed to &cubes
        cout <<"cubes.size()= "<<cubes.size()<<endl;
        cout <<"Done with ConverCubeLatticeToCubeVector"<<endl;


        for (int t=0;t<points;t++)
        {
                delete [] f[t];
        }
        delete [] f;
}
void CubeFFT3(Real contrast, IntensityStruct &i, vector<CubeStruct> &cubes, vector<AtomStruct> &Atoms, ParamStruct &params)
{
        //Calculates the scattering from cubes using FFT.  Don't get this
        //function mixed up with CubeFFT3D.  Should rename function.
        int CubeNum=cubes.size(), Max;
        int MaxXBin, MaxYBin, MaxZBin;
        int xbin, ybin, zbin;
        int xmax, ymax, zmax;
        int points=i.calc.size();
        Real CubeXMax, CubeYMax, CubeZMax;
        Real CubeXMin, CubeYMin, CubeZMin;

        Real numatom[NumAtomTypes];
        vector< vector< vector<Real> > > TempOutputRealX, TempOutputImagX;
        vector< vector< vector<Real> > > TempOutputRealY, TempOutputImagY;
        vector< vector< vector<Real> > > TempOutputRealZ, TempOutputImagZ;
        lattice CubeWeight3D;
        Max=2048;

        FindNumAtom(numatom, Atoms);
        cout <<"In CubeFFT"<<endl;
        MinMax(CubeXMin, CubeYMin, CubeZMin, CubeXMax, CubeYMax, CubeZMax, cubes);
        cout <<"CubeXMin= "<<CubeXMin<<" CubeYMin= "<<CubeYMin<<" CubeZMin= "<<CubeZMin<<endl;
        xmax=int((CubeXMax-CubeXMin)/params.CubeSize+0.5)+1;
        ymax=int((CubeYMax-CubeYMin)/params.CubeSize+0.5)+1;
        zmax=int((CubeZMax-CubeZMin)/params.CubeSize+0.5)+1;
        MaxXBin=xmax;
        MaxYBin=ymax;
        MaxZBin=zmax;
        cout <<"xmax= "<<xmax<<" ymax= "<<ymax<<" zmax= "<<zmax<<endl;
        FindProteinElectrons(Atoms, cubes, params.CubeSize);

        MoveAtoms(Atoms, -CubeXMin, -CubeYMin, -CubeZMin);

        for (int n=0;n<CubeNum;n++)
        {
                cubes[n].x-=CubeXMin;
                cubes[n].y-=CubeYMin;
                cubes[n].z-=CubeZMin;
        }

        ConvertCubeVectorToLattice(cubes, CubeWeight3D, params.CubeSize);
        Real TotalMass=0;
        for (xbin=0;xbin<xmax;xbin++)
        {
                for (ybin=0;ybin<ymax;ybin++)
                {
                        for (zbin=0;zbin<zmax;zbin++)
                        {
                                //cout <<"xbin= "<<xbin<<" ybin= "<<ybin<<" zbin= "<<zbin<<endl;
                                TotalMass+=CubeWeight3D[xbin][ybin][zbin].density*params.CubeSize*params.CubeSize*params.CubeSize;
                        }
                }
        }
        cout <<"Hydration+Excluded= "<<TotalMass<<endl;

        while (true)
        {
                bool Doubled=false;
                if (xmax>Max)
                {
                        Max*=2;
                        Doubled=true;
                }
                if (ymax>Max)
                {
                        Max*=2;
                        Doubled=true;
                }
                if (zmax>Max)
                {
                        Max*=2;
                        Doubled=true;
                }
                if (!Doubled) break;
        }
        cout  <<"Max= "<<Max<<endl;
        /*
           MemoryAllocated=AllocateMemoryForFFT(TempOutputImagX, TempOutputRealX, TempOutputImagY, TempOutputRealY, TempOutputImagZ, TempOutputRealZ, xmax, ymax, zmax, Max);
           if (MemoryAllocated) CubeFFT(CubeWeight3D, TempOutputImagX, TempOutputRealX, TempOutputImagY, TempOutputRealY, TempOutputImagZ, TempOutputRealZ, xmax, ymax, zmax, Max, contrast, points, i, cubes, Atoms, params);
           else 
           {
           TempOutputImagX.clear();
           TempOutputRealX.clear();
           TempOutputImagY.clear();
           TempOutputRealY.clear();
           TempOutputImagZ.clear();
           TempOutputRealZ.clear();
           ConvertCubeLatticeToCubeVector(CubeWeight3D, cubes);
           CalcCubeIntensity(cubes, i, points, Atoms, params);
           }
           */
        MemoryEfficientCubeFFT(CubeWeight3D, xmax, ymax, zmax, Max, contrast, points, i, cubes, Atoms, params);
}

void CubeFFT3D(Real numatom[], Real contrast, IntensityStruct &i, vector<CubeStruct> &cubes, vector<AtomStruct> &Atoms, ParamStruct params)
{
        //Calculates the scattering from the cubes using 3D FFT.  This function
        //does not appear to be working.  It should not be used.
        int TotalParticles=Atoms.size();
        int CubeNum=cubes.size(), Max;
        int MaxXBin, MaxYBin, MaxZBin;
        int xbin, ybin, zbin;
        int xmax, ymax, zmax;
        int points=i.calc.size();
        Real CubeXMax, CubeYMax, CubeZMax;
        Real CubeXMin, CubeYMin, CubeZMin;
        Real dx, dy, dz;
        Real mux, muy, muz;
        Real *ireal, *iimag;
        Real *IrealCube, *IimagCube;
        Real dotproduct, product;
        Real **f;
        Real inc;
        Real NumVectors;
        Real volume, TotalDensity;
        Real thetainc;
        Real theta;

        Real cosphi;
        Real nslices;
        Real sinphi;
        Real ***InputReal, ***InputImag, ***OutputReal, ***OutputImag;
        Real ***CubeWeight3D;
        Max=1024;
        vector<AtomStruct> AtomCubes;

        cout <<"In CubeFFT3D"<<endl;
        MinMax(CubeXMin, CubeYMin, CubeZMin, CubeXMax, CubeYMax, CubeZMax, cubes);
        cout <<"CubeXMin= "<<CubeXMin<<" CubeYMin= "<<CubeYMin<<" CubeZMin= "<<CubeZMin<<endl;
        cout <<"CubeXMax= "<<CubeXMax<<" CubeYMax= "<<CubeYMax<<" CubeZMax= "<<CubeZMax<<endl;
        xmax=int((CubeXMax-CubeXMin)/params.CubeSize+0.5)+1;
        ymax=int((CubeYMax-CubeYMin)/params.CubeSize+0.5)+1;
        zmax=int((CubeZMax-CubeZMin)/params.CubeSize+0.5)+1;
        MaxXBin=xmax;
        MaxYBin=ymax;
        MaxZBin=zmax;
        cout <<"xmax= "<<xmax<<" ymax= "<<ymax<<" zmax= "<<zmax<<endl;
        while (true)
        {
                bool Doubled=false;
                if (xmax>Max)
                {
                        Max*=2;
                        Doubled=true;
                }
                if (ymax>Max)
                {
                        Max*=2;
                        Doubled=true;
                }
                if (zmax>Max)
                {
                        Max*=2;
                        Doubled=true;
                }
                if (!Doubled) break;
        }
        cout  <<"Max= "<<Max<<endl;
        cout <<"Allocated TempOutputRealY"<<endl;

        SafeArrayAlloc(ireal, points, "ireal");
        SafeArrayAlloc(iimag, points, "iimag");
        SafeArrayAlloc(IrealCube, points, "IrealCube");
        SafeArrayAlloc(IimagCube, points, "IimagCube");

        FindProteinElectrons(Atoms, cubes, params.CubeSize);
        Safe3DArrayAlloc(InputReal, Max, Max, Max, "InputReal");
        Safe3DArrayAlloc(InputImag, Max, Max, Max, "InputImag");
        Safe3DArrayAlloc(OutputReal, Max, Max, Max, "OutputReal");
        Safe3DArrayAlloc(OutputImag, Max, Max, Max, "OutputImag");

        Safe3DArrayAlloc(CubeWeight3D, xmax, ymax, zmax, "CubeWeight3D");

        for (xbin=0;xbin<xmax;xbin++)
        {
                for (ybin=0;ybin<ymax;ybin++)
                {
                        for (zbin=0;zbin<zmax;zbin++)
                        {
                                CubeWeight3D[xbin][ybin][zbin]=0;
                        }
                }
        }

        for (xbin=0;xbin<Max;xbin++)
        {
                for (ybin=0;ybin<Max;ybin++)
                {
                        for (zbin=0;zbin<Max;zbin++)
                        {
                                InputReal[xbin][ybin][zbin]=0;
                                InputImag[xbin][ybin][zbin]=0;
                                OutputReal[xbin][ybin][zbin]=0;
                                OutputImag[xbin][ybin][zbin]=0;
                        }
                }
        }
        //AllocateCubeConversion(CubeXMin, CubeYMin, CubeZMin, CubeXMax, CubeYMax, CubeZMax);
        cout <<"CubeXMax= "<<CubeXMax<<" CubeYMax= "<<CubeYMax<<" CubeZMax= "<<CubeZMax<<endl;
        ConvertGrid2(CubeXMin, CubeYMin, CubeZMin, cubes, CubeWeight3D, params);
        Real TotalMass=0;
        for (xbin=0;xbin<xmax;xbin++)
        {
                for (ybin=0;ybin<ymax;ybin++)
                {
                        for (zbin=0;zbin<zmax;zbin++)
                        {
                                //cout <<"xbin= "<<xbin<<" ybin= "<<ybin<<" zbin= "<<zbin<<endl;
                                TotalMass+=CubeWeight3D[xbin][ybin][zbin]*params.CubeSize*params.CubeSize*params.CubeSize;
                                TotalDensity+=CubeWeight3D[xbin][ybin][zbin];
                        }
                }
        }
        Real NumElectrons[NumAtomTypes];
        SetNumElectrons(NumElectrons);
        for (int n=0;n<NumAtomTypes;n++)
        {
                cout <<"NumElectrons["<<n<<"]= "<<NumElectrons[n]<<" numatom["<<n<<"]= "<<numatom[n]<<endl;
                TotalMass+=NumElectrons[n]*numatom[n];
        }

        cout <<"TotalMass= "<<TotalMass<<endl;
        MoveAtoms(Atoms, -CubeXMin, -CubeYMin, -CubeZMin);

        for (int n=0;n<CubeNum;n++)
        {
                cubes[n].x-=CubeXMin;
                cubes[n].y-=CubeYMin;
                cubes[n].z-=CubeZMin;
        }
        //ConvertCubesToAtoms(cubes, AtomCubes);
        //AppendVector(AtomCubes, Atoms);
        //PrintPdb("/home/jouko/project/WAXS/pdb/TestCubes.pdb", Atoms);
        cout <<"TotalParticles= "<<TotalParticles<<endl;

        for (int t=0;t<points;t++) i.calc[t]=0;


        cout <<"CubeNum="<<CubeNum<<endl;

        volume=0;
        Safe2DArrayAlloc(f, points, NumAtomTypes, "f");

        for (int t=0;t<points;t++)
        {
                for (int n=0;n<NumAtomTypes;n++)
                {
                        f[t][n]=0;
                }
        }

        for (int n=0;n<NumAtomTypes;n++)
        {
                for (int t=0;t<points;t++)
                {
                        f[t][n]=SolventCorrectedScattering(n, 0, i.s[t], i.f[t][n], 0, 0, params);
                }
        }

        cout <<"f[0][0]="<<f[0][0]<<endl;
        NumVectors=0.0;
        nslices=20.0;
        cosphi=1.0/nslices-1.0;
        thetainc=20.0;
        theta=pi/thetainc;

        cout <<"About to average"<<endl;

        for (xbin=0; xbin<Max ; xbin++)
        {
                cout <<"xbin= "<<xbin<<endl;
                for (ybin=0; ybin<Max ; ybin++)
                {
                        for (zbin=0; zbin<Max ; zbin++)
                        {
                                //cout <<"CubeWeight3D["<<xbin<<"]["<<ybin<<"]["<<zbin<<"]= "<<endl;
                                if (xbin<xmax && ybin<ymax && zbin<zmax) InputReal[xbin][ybin][zbin]=CubeWeight3D[xbin][ybin][zbin];
                                else InputReal[xbin][ybin][zbin]=0.0;
                                //if (ybin==35 && zbin==50 )cout <<"InputReal["<<xbin<<"]= "<<InputReal[xbin]<<endl;
                        }
                }
        }
        cout <<"About to enter FFT3D"<<endl;
        FFT3D(InputReal, InputImag, OutputReal, OutputImag, Max);
        cout <<"OutputReal[0][0][0]= "<<OutputReal[0][0][0]<<" OutputImag[0][0][0]= "<<OutputImag[0][0][0]<<endl;
        cout <<"TotalDensity= "<<TotalDensity<<endl;
        cout <<"Finnished with FFTs"<<endl;

        while (true)
        {
                NumVectors+=1.0;
                cout <<"cosphi="<<cosphi<<" theta="<<theta<<endl;
                sinphi=sqrt(1-cosphi*cosphi);
                dx=sinphi*cos(theta);
                dy=sinphi*sin(theta);
                dz=cosphi;

                for (int t=0;t<points;t++)
                {
                        ireal[t]=0;
                        iimag[t]=0;
                        IrealCube[t]=0;
                        IimagCube[t]=0;
                        mux=i.s[t]*dx;
                        muy=i.s[t]*dy;
                        muz=i.s[t]*dz;
                        //cout <<"s["<<t<<"]= "<<s[t]<<" dx= "<<dx<<" dy= "<<dy<<" dz= "<<dz<<endl;
                        for (int n=0;n<TotalParticles;n++)
                        {
                                //cout <<"x["<<n<<"]= "<<Atoms[n].x<<" y= "<<y[n]<<" z= "<<z[n]<<endl;
                                dotproduct=mux*Atoms[n].x+muy*Atoms[n].y+muz*Atoms[n].z;
                                ireal[t]+=f[t][Atoms[n].atomid]*cos(dotproduct);
                                iimag[t]+=f[t][Atoms[n].atomid]*sin(dotproduct);
                                if (n<10)
                                {
                                        //cout <<"ireal["<<t<<"]= "<<ireal[t]<<" f["<<t<<"]["<<atomid[n]<<"]= "<<f[t][atomid[n]]<<endl;
                                        //cout <<"cos("<<dotproduct<<")= "<<cos(dotproduct)<<" sin("<<dotproduct<<")= "<<sin(dotproduct)<<endl;
                                        //cout <<"iimag["<<t<<"]= "<<iimag[t]<<endl;
                                }
                        }
                }
                cout <<"ireal[10]= "<<ireal[10]<<" iimag[10]= "<<iimag[10]<<endl;
                InterpolateFFT3D(OutputReal, OutputImag, IrealCube, IimagCube, points, i, xmax, ymax, dy, dz, dx, inc, Max, params);

                for (int t=0;t<points;t++)
                {
                        mux=i.s[t]*dx;
                        muy=i.s[t]*dy;
                        muz=i.s[t]*dz;
                        product=8.0*sin(mux*params.CubeSize*0.5)*sin(muy*params.CubeSize*0.5)*sin(muz*params.CubeSize*0.5)/(mux*muy*muz);
                        if (t==10)
                        {
                                cout <<"IrealCube["<<t<<"]= "<<IrealCube[t]<<endl;
                                cout <<"IimagCube["<<t<<"]= "<<IimagCube[t]<<endl;
                        }
                        ireal[t]=ireal[t]+IrealCube[t]*product;
                        iimag[t]=iimag[t]+IimagCube[t]*product;
                        i.calc[t]=i.calc[t]+ireal[t]*ireal[t]+iimag[t]*iimag[t];
                }

                theta += 2.0*pi/thetainc;
                if (theta>2.0*pi)
                {
                        theta=pi/thetainc;
                        cosphi += 2.0/nslices;
                        if (cosphi>0.0-0.9/nslices)
                        {
                                break;
                        }
                }
        }

        if (points>100) cout <<"i[100]="<<i.calc[100]<<endl;

        for (int n=0;n<points;n++) i.calc[n]=i.calc[n]*0.5/NumVectors;

        for (xbin=0;xbin<Max;xbin++)
        {
                for (ybin=0;ybin<Max;ybin++)
                {
                        delete [] InputReal[xbin][ybin];
                        delete [] InputImag[xbin][ybin];
                        delete [] OutputReal[xbin][ybin];
                        delete [] OutputImag[xbin][ybin];
                }
                delete [] InputReal[xbin];
                delete [] InputImag[xbin];
                delete [] OutputReal[xbin];
                delete [] OutputImag[xbin];
        }

        delete [] InputReal;
        delete [] InputImag;
        delete [] OutputReal;
        delete [] OutputImag;

        delete [] ireal;
        delete [] iimag;
        delete [] IrealCube;
        delete [] IimagCube;

        for (int t=0;t<points;t++)
        {
                delete [] f[t];
        }
        delete [] f;
        for (xbin=0;xbin<xmax;xbin++)
        {
                for (ybin=0;ybin<ymax;ybin++)
                {
                        delete [] CubeWeight3D[xbin][ybin];
                }
                delete [] CubeWeight3D[xbin];
        }
        delete [] CubeWeight3D;
}
int max(int a, int b)
{
        if (a>b) return a;
        else return b;
}

int min(int a, int b)
{

        if (a<b) return a;
        else return b;
}

Real max(Real a, Real b)
{

        if (a>b) return a;
        else return b;
}

int CrysolSolvent(Real numatom[], Real atomm[], int natom, Real contrast, Real atomr[], vector<AtomStruct> &Atoms, ParamStruct params)
{
        //Uses a method similar to that of Crysol to evaluate surface, and
        //add in dummy hydration shell atoms.
        AtomStruct Atom;
        int TotalParticles;
        Real nslices;
        Real cosphi, sinphi;
        Real thetainc, theta;
        Real current;
        Real ProbeRadius=1.4;
        Real Projection;
        Real solventatomr[5];
        Real VectX, VectY, VectZ;
        Real xn, yn, zn;

        int temp=Atoms.size();

        for (int n=temp;0<=n;n--)
        {
                if (Atoms[n].atomid==HydrationShell || Atoms[n].atomid==ExcludedVolume) Atoms.pop_back();
        }

        numatom[HydrationShell]=0;

        center(atomm, contrast, atomr, Atoms, params);

        solventatomr[HYDROGEN]=0.787609;
        solventatomr[CARBON]=1.540862;
        solventatomr[NITROGEN]=1.295872;
        solventatomr[OXYGEN]=1.205371;
        solventatomr[SULFUR]=1.582549;

        atomr[7]=1.5;

        nslices=50.0;
        cosphi=1.0/nslices-1.0;
        thetainc=50.0;
        theta=pi/thetainc;

        TotalParticles=natom;

        while (true)
        {
                sinphi=sqrt(1.0-cosphi*cosphi);
                VectX=sinphi*cos(theta);
                VectY=sinphi*sin(theta);
                VectZ=cosphi;
                theta += 2.0*pi/thetainc;
                if (theta>2.0*pi)
                {
                        theta=pi/thetainc;
                        cosphi+=2.0/nslices;
                        if (cosphi>1.0)
                        {
                                break;
                        }
                }

                current=0;
                for (int n=0;n<natom;n++)
                {
                        xn=Atoms[n].x;
                        yn=Atoms[n].y;
                        zn=Atoms[n].z;
                        Projection=VectX*xn+VectY*yn+VectZ*zn;
                        current=max(current, Projection+solventatomr[Atoms[n].atomid]*0.5+ProbeRadius);
                }

                if (current <= 0)
                {
                        cout <<"current="<<current<<endl;
                }
                Atom.x=VectX*current;
                Atom.y=VectY*current;
                Atom.z=VectZ*current;
                Atom.AtomNumber=TotalParticles;
                Atom.AtomName="So";
                Atom.ResidueName="SOL";
                Atom.ChainName="A";
                Atom.ResidueNum=1;
                Atom.Occupancy=1;
                Atom.BFactor=0;
                Atom.SegID="SOLV";
                Atom.ID="So";
                Atom.atomid=HydrationShell;
                SafePushBack(Atoms, Atom, "Atoms");
                numatom[HydrationShell]+=1.0;
                TotalParticles++;
        }

        PrintPdb("/home/jouko/WAXS/pdb/CrysolSolvent", Atoms);

        return TotalParticles;
}

int FindGreatestDistance(Real bin, vector<AtomStruct> &Atoms)
{
        //Obtains a value which is certain to be greater then Dmax.  Then
        //calculates the integer value that corresponds to that distance
        //in the p(r) array.  This is used to initialize p(r).
        Real dx, dy, dz, dist;
        Real MinX, MinY, MinZ;
        Real MaxX, MaxY, MaxZ;

        MinMax(MinX, MinY, MinZ, MaxX, MaxY, MaxZ, Atoms);

        dx=MaxX-MinX;
        dy=MaxY-MinY;
        dz=MaxZ-MinZ;

        dist=sqrt(dx*dx+dy*dy+dz*dz);

        return int(floor(dist/bin));

}

void FourierTransformIntensity(IntensityStruct &i, vector<Real> &PR, int GreatestDistance, Real bin)
{
        //Calculates p(r) from I(s) using Fourier transformation.
        int points=i.s.size();
        Real r=bin, sinc=i.s[1]-i.s[0];
        for (int distance=1;distance<GreatestDistance;distance++)
        {
                for (int n=0;n<points;n++)
                {
                        PR[distance]+=i.calc[n]*i.s[n]*sin(i.s[n]*r)*sinc;
                }
                PR[distance]=r*PR[distance]*bin/(2.0*pi);
                r+=bin;
        }
}

void FourierTransformPr(IntensityStruct &i, vector<Real> &PR, int GreatestDistance, Real bin)
{
        //Calculates I(s) from p(r) using Fourier transformation.
        int points=i.s.size();
        Real dist;
        for (int m=0;m<points;m++)
        {
                dist=bin;
                i.calc[m]=0;
                for (int n=1;n<GreatestDistance;n++)
                {
                        i.calc[m]+=PR[n]*sin(i.s[m]*dist)*4.0*pi*pi*bin/(i.s[m]*dist);
                        dist+=bin;
                }
        }
}

void EquateIntensity(IntensityStruct &i1, IntensityStruct &i2)
{
        i2.calc=i1.calc;
        i2.s=i1.s;
}

Real CheckS(Real atomr[], Real s, Real density, Real hsdensity, ParamStruct &params, Real f0[])
{
        //Checks to make sure that the scattering factor of an atom at I(0)
        //equals the number of electrons in the atom.
        Real f, max=0, quotient;
        for (int m=0;m<NumAtomTypes;m++)
        {
                f=SolventCorrectedScattering(m, atomr[m], s, density, hsdensity, params);
                quotient=f/f0[m];
                if (quotient>max) max=quotient;
        }
        return max;
}

Real DetermineSmax(Real atomr[], Real density, Real hsdensity, ParamStruct &params)
{
        //Finds a value of s for which f(s) is below a certain fraction of f(0).
        Real f0[NumAtomTypes];
        Real small_s=1e-5, test_s, criterion=1e-1, max_quotient;
        for (int j=0;j<NumAtomTypes;j++) f0[j]=0;
        for (int m=0;m<NumAtomTypes;m++)
        {
                f0[m]=SolventCorrectedScattering(m, atomr[m], small_s, density, hsdensity, params);
        }
        test_s=small_s;
        while (true)
        {
                max_quotient=CheckS(atomr, test_s, density, hsdensity, params, f0);
                if (max_quotient<criterion) return test_s;
                else test_s*=1.2;
        }
}
void CalcPrFromScattering(vector<Real> &PR, Real atomr[], Real bin, Real contrast, Real hsdensity, string PrFromStructureFile, vector<AtomStruct> &Atoms, ParamStruct params)
{
        //Calculates scattering intensity to very large s and then finds p(r)
        //using Fourier transformation.
        string PrDir, IntensityDir;
        int points;
        int GreatestDistance;
        int window;
        Real sinc, smax;
        Real mindistance[NumAtomTypes][NumAtomTypes];
        Real maxdistance[NumAtomTypes][NumAtomTypes];
        IntensityStruct i, i2;

        cout <<"In CalcPrFromScattering"<<endl;

        GreatestDistance=FindGreatestDistance(bin, Atoms);
        SafeAlloc(PR, GreatestDistance+1, "PR");

        for (int n=0;n<NumAtomTypes;n++)
        {
                for (int m=0;m<NumAtomTypes;m++)
                {
                        mindistance[n][m]=GreatestDistance;
                        maxdistance[n][m]=0;
                }
        }

        sinc=0.001;
        //smax=DetermineSmax(atomr, contrast, hsdensity, params);
        //cout <<"smax= "<<smax<<endl;
        smax=2.0*pi/bin;
        cout <<"smax= "<<smax<<endl;
        points=int(smax/sinc);

        InitializeIntensity(i, points);
        i.s[0]=sinc;
        for (int n=1;n<points;n++) i.s[n]=i.s[n-1]+sinc;
        AtomScatteringFactor(i, params.ScatteringType);

        window=int(5.0/bin);
        /*
           for (int n=0;n<NumAtomTypes;n++)
           {
           for (int m=0;m<NumAtomTypes;m++)
           {
           fn=SolventCorrectedScattering(n, atomr[n], 40.0, contrast, hsdensity, params);
           fm=SolventCorrectedScattering(m, atomr[m], 40.0, contrast, hsdensity, params);
           cout <<fn*fm<<"\t";
           }
           cout <<endl;
           }
           */
        cout <<"About to enter findpr"<<endl;
        if (pr.size()==0) findpr(bin, maxdistance, mindistance, Atoms, params.UniformHydrationShell);
        MinMaxPr(maxdistance, mindistance);
        InitializeSinc(bin, maxdistance, i);
        EquateIntensity(i, i2);
        //PrintPdb("/home/jouko/project/WAXS/test/log/InCalcPrFromScattering.pdb", Atoms);
        cout <<"About to enter intensity"<<endl;
        Debye(atomr, i, contrast, hsdensity, bin, points, maxdistance, mindistance, Atoms, params);
        cout <<"Left intensity"<<endl;

        cout <<"Calculated PR"<<endl;

        FourierTransformIntensity(i, PR, GreatestDistance, bin);
        //FourierTransformPr(i2, PR, GreatestDistance, bin);
        NormalizeArea(PR, bin);

        cout << endl;
        for (int m=0;m<2000;m++)
        {
                cout <<i.s[m]<<"\t"<<i.calc[m]<<"\t"<<i2.calc[m]<<endl;
        }
}

void TestScatteringFactorPr(Real atomr, ParamStruct &params)
{
        //Calculates p(r) for a single atom, by Fourier transforming its
        //scattering factor.
        int GreatestDistance, points;
        Real smax, sinc=0.001, bin=0.05, r=0;
        vector<Real> PR;
        IntensityStruct i;

        cout <<"In TestScatteringFactorPr"<<endl;
        GreatestDistance=200;
        SafeAlloc(PR, GreatestDistance, "PR");
        smax=2.0*pi/bin;
        points=int(smax/sinc);
        cout <<"points= "<<points<<endl;
        InitializeIntensity(i, points);
        cout <<"After InitializeIntensity"<<endl;
        i.s[0]=sinc;
        for (int n=1;n<points;n++) i.s[n]=i.s[n-1]+sinc;
        AtomScatteringFactor(i, params.ScatteringType);
        cout <<"After AtomScatteringFactor"<<endl;
        for (int j=0;j<points;j++)
        {
                i.calc[j]=SolventCorrectedScattering(ExcludedVolume, atomr, i.s[j], 1.0, 1.0, params);
        }
        cout <<"After i.calc"<<endl;
        FourierTransformIntensity(i, PR, GreatestDistance, bin);
        cout <<"After FourierTransformIntensity"<<endl;
        ofstream file;
        file.open("/home/jouko/project/WAXS/pr/cube_pr", ios::app);

        for (int j=0;j<GreatestDistance;j++)
        {
                file <<r<<"\t"<<PR[j]<<endl;
                r+=bin;
        }
}

void CalcPrFromScattering(Real atomr[], Real bin, Real contrast, Real hsdensity, string PrFromStructureFile, vector<AtomStruct> &Atoms, ParamStruct params)
{
        //Calculates scattering intensity for different solvent effects and
        //calculates the p(r) from the scattering intensity.
        char CharPrDir[1000];
        int GreatestDistance;
        vector<Real> PR, PR_NoHydration, PR_vacuum;
        Real dist;
        CalcPrFromScattering(PR, atomr, bin, contrast, hsdensity, PrFromStructureFile, Atoms, params);
        params.UniformHydrationShell=true;
        CalcPrFromScattering(PR_NoHydration, atomr, bin, contrast, 0, PrFromStructureFile, Atoms, params);
        CalcPrFromScattering(PR_vacuum, atomr, bin, 0, 0, PrFromStructureFile, Atoms, params);
        GreatestDistance=PR.size();
        dist=0;
        cout <<"PrFromStructureFile= "<<PrFromStructureFile<<endl;
        strcpy(CharPrDir, PrFromStructureFile.c_str());
        ofstream pofr2(CharPrDir, ios::app);
        pofr2 <<"r\tUser P(r)\tNo hydration shell P(r)\tVacuum P(r)"<<endl;
        for (int n=0;n<GreatestDistance;n++)
        {
                pofr2 << dist<<"\t"<<PR[n]
                        <<"\t"<<PR_NoHydration[n]
                        <<"\t"<<PR_vacuum[n]<<endl;
                dist+=bin;
        }
        pofr2 << endl;
}

void NormalizeArea(vector<Real> &PR, Real bin)
{
        //Set the area under a curve to 1.
        int GreatestDistance=PR.size();
        Real TotalPr=0;
        cout <<"In NormalizeArea"<<endl;
        for (int n=0;n<GreatestDistance;n++) TotalPr+=PR[n];
        cout <<"In NormalizeArea TotalPr= "<<TotalPr<<" bin= "<<bin<<endl;	
        for (int n=0;n<GreatestDistance;n++) PR[n]=PR[n]/(TotalPr*bin);
}

void Pr3DToPr1D(vector< vector< vector<Real> > > &PR, vector<Real> &PR2, vector<Real> &PR_NoHydration, vector<Real> &PR_Vacuum, vector<Real> &PR_Excluded, Real NumElectrons[], Real atomr[], ParamStruct &params)
{
        //Creates a 1D p(r) whrich only has information about electron pair
        //correlations from a 3D p(r, AtomType1, AtomType2), which is a 
        //collection of p(r)s for different atom type pairs.
        int GreatestDistance=PR[0][0].size();
        Real MassProduct, VacuumMassProduct;
        Real SolventCorrectedElectrons[NumAtomTypes], ExcludedElectrons[NumAtomTypes];
        CalculateSolventCorrectedElectrons(SolventCorrectedElectrons, params.contrast, atomr, params);
        CalculateExcludedElectrons(ExcludedElectrons, params.contrast, atomr, params); 
        cout <<"GreatestDistance= "<<GreatestDistance<<endl;
        for (int n=0;n<NumAtomTypes;n++)
        {
                cout <<"SolventCorrectedElectrons["<<n<<"]= "<<SolventCorrectedElectrons[n]<<" NumElectrons= "<<NumElectrons[n]<<endl;
                for (int m=n;m<NumAtomTypes;m++)
                {
                        MassProduct=SolventCorrectedElectrons[n]*SolventCorrectedElectrons[m];
                        cout <<"SolventCorrectedElectrons["<<n<<"]= "<<SolventCorrectedElectrons[n]<<endl;
                        cout <<"SolventCorrectedElectrons["<<m<<"]= "<<SolventCorrectedElectrons[m]<<endl;
                        cout <<"MassProduct= "<<MassProduct<<endl;
                        VacuumMassProduct=NumElectrons[n]*NumElectrons[m];
                        for (int distance=0;distance<GreatestDistance;distance++)
                        {
                                PR2[distance]+=MassProduct*PR[n][m][distance];
                                //cout <<"m= "<<m<<" n= "<<n<<" PR2["<<distance<<"]= "<<PR2[distance]<<endl;
                                if (n!=HydrationShell && m!=HydrationShell) 
                                {
                                        PR_NoHydration[distance]+=MassProduct*PR[n][m][distance];
                                        cout <<"MassProduct= "<<MassProduct<<endl;
                                        cout <<"PR["<<n<<"]["<<m<<"]["<<distance<<"]= "<<PR[n][m][distance]<<endl;
                                        cout <<"distance= "<<distance<<endl;
                                        cout <<"PR_NoHydration.size()= "<<PR_NoHydration.size()<<endl;
                                        cout <<"PR_NoHydration= "<<PR_NoHydration[distance]<<endl;
                                }
                                if (n!=HydrationShell && m!=HydrationShell) PR_Vacuum[distance]+=VacuumMassProduct*PR[n][m][distance];
                                if (n==ExcludedVolume && m==ExcludedVolume)
                                {
                                        PR_Excluded[distance]+=ExcludedElectrons[n]*ExcludedElectrons[m]*PR[n][m][distance];
                                }
                        }
                }
        }
}

Real SumOverRange(vector<Real> &v, int start, int end)
{
        //Returns the sum of a part of a vector.
        int Size=v.size();
        Real sum=0;

        if (start<0) start=0;
        if (end>=Size) end=Size-1;
        for (int j=start;j<end;j++) sum+=v[j];
        return sum;
}

void ChangeBinSize(vector<Real> &pr1, vector<Real> &pr2, Real bin1, Real bin2)
{
        //Increases the bin size of a 1D p(r).
        int end, start, BinRatio=int(bin2/bin1+0.5);
        int GreatestDistance1=pr1.size();
        int GreatestDistance2=GreatestDistance1/BinRatio;
        int GreatestDistanceb=pr2.size();
        cout <<"BinRatio= "<<BinRatio<<endl;
        cout <<"GreatestDistance2= "<<GreatestDistance2<<" GreatestDistanceb= "<<GreatestDistanceb<<endl;
        for (int j=0;j<GreatestDistanceb;j++)
        {
                start=j*BinRatio;
                end=(j+1)*BinRatio;
                cout <<"j= "<<j<<" start= "<<start<<" end= "<<end<<endl;
                pr2[j]=SumOverRange(pr1, start, end);
        }
}

void ChangeBinSize(vector< vector< vector<Real> > > &pr1, vector< vector< vector<Real> > > &pr2, Real bin1, Real bin2)
{
        //Changes the bin size of a 3D p(r).
        cout <<"In ChangeBinSize"<<endl;
        for (int j=0;j<NumAtomTypes;j++)
        {
                for (int k=0;k<NumAtomTypes;k++)
                {
                        ChangeBinSize(pr1[j][k], pr2[j][k], bin1, bin2);
                }
        }
}

void ConvoluteWithGaussian(vector<Real> &pr, Real width, Real bin)
{
        //Covolutes a 1D p(r) or any other vector with a Gaussian.
        int size=pr.size();
        Real dr;
        vector<Real> TempPr;

        SafeAlloc(TempPr, size, "TempPr in CovoluteWithGaussian");

        for (int j=0;j<size;j++)
        {
                for (int k=0;k<size;k++)
                {
                        dr=Real(k-j)*bin;
                        TempPr[k]+=pr[j]*exp(-dr*dr/width);
                }
        }
        pr=TempPr;
}

void Convolute3DPrWithGaussian(vector< vector< vector<Real> > > &pr, Real bin)
{
        //Convolutes a 3D p(r) with a Gaussian.  The width of the Gaussian
        //should depend on the atom types, but this was not implemented.
        int GreatestDistance;

        if (pr.size()>0)
        {
                if (pr[0].size()>0)
                {
                        GreatestDistance=pr[0][0].size();
                }
        }
        else
        {
                cout <<"ERROR: Pr size is 0"<<endl;
                exit(EXIT_FAILURE);
        }

        for (int j=0;j<NumAtomTypes;j++)
        {
                for (int k=0;k<NumAtomTypes;k++)
                {
                        ConvoluteWithGaussian(pr[j][k], 2.0, bin);
                }
        }

}

void CalcPr2(Real atomr[], Real bin, Real contrast, Real hsdensity, string PrFromStructureFile2, vector<AtomStruct> &Atoms, ParamStruct params, vector< vector< vector<Real> > > &PR, vector<Real> &PR2, vector<Real> &PR_NoHydration, vector<Real> &PR_Vacuum, vector<Real> &PR_Excluded)
{
        //Calculates 1D p(r) from the structure of the protein.
        //Simultaneously calculates p(r) for the hydrated protein, unhydrated
        //protein, protein in vacuum, and the p(r) of the dummy excluded volume
        //atoms.
        int TotalParticles=Atoms.size();
        int GreatestDistance;
        Real dist, invbin;
        Real dx, dy, dz;
        Real xmax, ymax, zmax;
        Real xmin, ymin, zmin;
        Real NumElectrons[NumAtomTypes], ExcludedElectrons[NumAtomTypes];
        Real maxdistance[NumAtomTypes][NumAtomTypes], mindistance[NumAtomTypes][NumAtomTypes];
        cout << "In CalcPr2"<<endl;

        SetNumElectrons(NumElectrons);
        invbin=1.0/bin;

        CalculateExcludedElectrons(ExcludedElectrons, contrast, atomr, params);
        MinMax(xmin, ymin, zmin, xmax, ymax, zmax, Atoms);
        dx=xmax-xmin;
        dy=ymax-ymin;
        dz=zmax-zmin;
        dist=sqrt(dx*dx+dy*dy+dz*dz);
        GreatestDistance=int(floor(dist*invbin+0.5));

        SafeAlloc(PR2, GreatestDistance, "PR2");
        SafeAlloc(PR_NoHydration, GreatestDistance, "PR_NoHydration");
        SafeAlloc(PR_Vacuum, GreatestDistance, "PR_Vacuum");
        SafeAlloc(PR_Excluded, GreatestDistance, "PR_Excluded");
        Safe3DAlloc(PR, NumAtomTypes, NumAtomTypes, GreatestDistance, "PR");
        cout <<"atomr[ExcludedVolume]= "<<atomr[ExcludedVolume]<<endl;
        for (int n=0;n<NumAtomTypes;n++)
        {
                Real f, g;
                g=ScatteringFactor(n, 0.00001, params.ScatteringType);
                f=SolventCorrectedScattering(n, atomr[n], 0.00001, g, params.contrast, params.hsdensity, params);
                //cout <<"SolventCorrectedElectrons["<<n<<"]= "
                //        <<SolventCorrectedElectrons[n]
                //        <<" f= "<<f<<" ExcludedElectrons["<<n<<"]= "
                //        <<ExcludedElectrons[n]<<endl;
        }
        cout <<"CubeSize= "<<params.CubeSize<<endl;
        cout <<"TotalParticles= "<<TotalParticles<<endl;
        if (pr.size()==0) findpr(params.bin, maxdistance, mindistance, Atoms, params.UniformHydrationShell);
        if (params.bin!=bin) 
        {
                cout <<"About to ChangeBinSize"<<endl;
                ChangeBinSize(pr, PR, params.bin, bin);	
                cout <<"Changed bin size"<<endl;
        }
        else 
        {
                cout <<"About to PR=pr"<<endl;
                PR=pr;
                cout <<"Set PR=pr"<<endl;
        }
        cout <<"About to enter Pr3DToPr1D"<<endl;
        Convolute3DPrWithGaussian(PR, bin);
        Pr3DToPr1D(PR, PR2, PR_NoHydration, PR_Vacuum, PR_Excluded, NumElectrons, atomr, params);
        NormalizeArea(PR2, bin);
        NormalizeArea(PR_NoHydration, bin);
        NormalizeArea(PR_Vacuum, bin);
}

void CalcPr2(Real atomr[], Real bin, Real contrast, Real hsdensity, string PrFromStructureFile2, vector<AtomStruct> &Atoms, ParamStruct params)
{
        //Calculates p(r) for various solvation conditions and print the 
        //output.
        char CharPrFile[1000];
        int GreatestDistance;
        Real dist;
        vector<Real> PR2, PR_NoHydration, PR_Vacuum, PR_Excluded;
        vector< vector< vector<Real> > > PR;

        CalcPr2(atomr, bin, contrast, hsdensity, PrFromStructureFile2, Atoms, params, PR, PR2, PR_NoHydration, PR_Vacuum, PR_Excluded);
        GreatestDistance=PR[0][0].size();
        AddIndexToFile(PrFromStructureFile2);
        strcpy(CharPrFile, PrFromStructureFile2.c_str());
        cout <<"GreatestDistance= "<<GreatestDistance<< " CharPrFile= "<<CharPrFile<<endl;
        cout <<"PrFromStructureFile2= "<<PrFromStructureFile2<<endl;
        dist=0;
        ofstream PrOut(CharPrFile, ios::app);
        //PrOut <<"Pr for protein "<<params.xyz<<endl;
        PrOut <<"r\tUser P(r)\tNo hydration shell P(r)\tVacuum P(r)\tExcluded P(r)"<<endl;
        for (int distance=0;distance<GreatestDistance;distance++)
        {
                ofstream PrOut(CharPrFile, ios::app);
                PrOut << dist << "\t" 
                        << PR2[distance] << "\t" 
                        << PR_NoHydration[distance] << "\t" 
                        << PR_Vacuum[distance] << "\t"
                        << PR_Excluded[distance] << endl;
                dist+=bin;
        }
}

void ConvertCubesToAtoms(vector<CubeStruct> &cubes, Real atomr[], vector<AtomStruct> &Atoms, ParamStruct &params)
{
        //Copies some of the information in a vector of cube to a vector of
        //atoms.
        AtomStruct Atom;
        int natom=Atoms.size(), CubeNum=cubes.size();
        cout <<"In ConvertCubesToAtoms"<<endl;
        Atoms.resize(natom+CubeNum);
        cout <<"Resized Atoms"<<endl;
        Atom.ResidueName="HOH";
        Atom.AtomName="O";
        for (int n=0;n<CubeNum;n++)
        {
                Atom.x=cubes[n].x;
                Atom.y=cubes[n].y;
                Atom.z=cubes[n].z;
                if (cubes[n].density==-params.contrast)
                {
                        Atom.atomid=ExcludedVolume;
                        Atom.weight=1.0;
                }
                else
                {
                        Atom.atomid=HydrationShell;
                        Atom.weight=cubes[n].density;
                }
                Atoms[natom+n]=Atom;
        }
        cout <<"Atoms.size()= "<<Atoms.size()<<endl;
        for (int n=0;n<NumAtomTypes;n++) atomr[n]=0;
        atomr[ExcludedVolume]=params.CubeSize*CubeToSphere;
        if (params.ExcludedVolumeSphereType=="GaussianSphere") atomr[ExcludedVolume]*=HardSphereToGaussianSphere;
        atomr[HydrationShell]=atomr[ExcludedVolume];
}

void PrintMaxdistance(Real maxdistance[][NumAtomTypes])
{
        //Print out the Dmax array of the 3D p(r).  I.e the maximum distance
        //between each pair of atom types.
        cout <<"maxdistance= "<<endl;
        for (int n=0;n<NumAtomTypes;n++)
        {
                for (int m=0;m<NumAtomTypes;m++)
                {
                        cout <<"maxdistance["<<n<<"]["<<m<<"]= "<<maxdistance[n][m]<<endl;
                }
        }
}

void MinMaxPr(Real maxdistance[][NumAtomTypes], Real mindistance[][NumAtomTypes])
{
        //Finds the Dmin and Dmax arrays of the the 3D p(r).  I.e. the 
        //minimum and maximum distance between each pair of atom types.
        bool first;
        int GreatestDistance;

        GreatestDistance=pr[0][0].size();
        for (int n=0;n<NumAtomTypes;n++)
        {
                for (int m=n;m<NumAtomTypes;m++)
                {
                        first=true;
                        for (int distance=0;distance<GreatestDistance;distance++)
                        {
                                //cout <<"pr["<<n<<"]["<<m<<"]["<<distance<<"]= "<<pr[n][m][distance]<<endl;
                                if (pr[n][m][distance]!=0)
                                {
                                        maxdistance[n][m]=distance;
                                }

                                if (pr[n][m][distance]!=0 && first)
                                {
                                        mindistance[n][m]=distance;
                                        first=false;
                                }
                        }
                }
        }
}

void findpr(Real bin, Real maxdistance[][NumAtomTypes], Real mindistance[][NumAtomTypes], vector<AtomStruct> &Atoms, bool UniformHydrationShell)
{
        //Calculates the 3D p(r) from the protein structure.
        //I.e p(r, AtomType1, AtomType2)
        int distance;
        int atomidm, atomidn;
        int GreatestDistance;
        int TotalParticles=Atoms.size();
        Real dist;
        Real invbin;
        Real weightn;
        Real xn, yn, zn;
        Real dx, dy, dz;
        Real xmax, ymax, zmax;
        Real xmin, ymin, zmin;
        cout <<"TotalParticles= "<<TotalParticles<<endl;
        invbin=1.0/bin;
        MinMax(xmin, ymin, zmin, xmax, ymax, zmax, Atoms);	
        cout <<"xmin= "<<xmin<<" ymin= "<<ymin<<" zmin= "<<zmin<<" xmax= "<<xmax<<" ymax= "<<ymax<<" zmax= "<<zmax<<endl;
        dx=xmax-xmin;
        dy=ymax-ymin;
        dz=zmax-zmin;
        dist=sqrt(dx*dx+dy*dy+dz*dz);
        GreatestDistance=int(dist*invbin+0.5)+1;
        cout <<"GreatestDistance= "<<GreatestDistance<<endl;
        Safe3DAlloc(pr, NumAtomTypes, NumAtomTypes, GreatestDistance, "pr");
        cout <<"Assigned pr"<<endl;
        for (int n=0;n<NumAtomTypes;n++)
        {
                for (int m=0;m<NumAtomTypes;m++)
                {
                        maxdistance[n][m]=0;
                        mindistance[n][m]=0;
                }
        }
        cout <<"Initialized maxdistance"<<endl;
        for (int n=0;n<NumAtomTypes;n++)
        {
                for (int m=0;m<NumAtomTypes;m++)
                {
                        mindistance[n][m]=1000;
                }
        }
        cout <<"Initialized mindistance"<<endl;
        if (!UniformHydrationShell)
        {
                for (int n=0;n<TotalParticles-1;n++)
                {
                        xn=Atoms[n].x;
                        yn=Atoms[n].y;
                        zn=Atoms[n].z;
                        atomidn=Atoms[n].atomid;
                        weightn=Atoms[n].weight;
                        //PrintAtomInfo(Atoms[n]);
                        for (int m=n+1;m<TotalParticles;m++)
                        {
                                dx=xn-Atoms[m].x;
                                dy=yn-Atoms[m].y;
                                dz=zn-Atoms[m].z;
                                dist=sqrt(dx*dx+dy*dy+dz*dz);
                                distance=int(floor(dist*invbin+0.5));
                                atomidm=Atoms[m].atomid;
                                //cout <<"n= "<<n<<" m= "<<m<<" atomidm= "<<atomidm<<" atomidn= "<<atomidn<<" distance= "<<distance<<endl;
                                //PrintAtomInfo(Atoms[m]);
                                if (atomidn+1>atomidm) pr[atomidm][atomidn][distance]+=weightn*Atoms[m].weight;
                                else pr[atomidn][atomidm][distance]+=weightn*Atoms[m].weight;
                        }
                }
        }
        else
        {
                for (int n=0;n<TotalParticles-1;n++)
                {
                        xn=Atoms[n].x;
                        yn=Atoms[n].y;
                        zn=Atoms[n].z;
                        atomidn=Atoms[n].atomid;
                        if (n%1000==0) cout <<"n= "<<n<<" atomid= "<<Atoms[n].atomid<<endl;
                        for (int m=n+1;m<TotalParticles;m++)
                        {
                                dx=xn-Atoms[m].x;
                                dy=yn-Atoms[m].y;
                                dz=zn-Atoms[m].z;
                                dist=sqrt(dx*dx+dy*dy+dz*dz);
                                distance=int(floor(dist*invbin+0.5));
                                atomidm=Atoms[m].atomid;
                                if (n==2568)
                                {
                                        //cout <<"atomidm= "<<atomidm<<" atomidn= "<<atomidn<<" distance= "<<distance<<endl;
                                        //PrintAtomInfo(Atoms[n]);
                                        //PrintAtomInfo(Atoms[m]);
                                        //cout <<endl;
                                }
                                if (atomidn+1>atomidm) pr[atomidm][atomidn][distance]+=1.0;
                                else pr[atomidn][atomidm][distance]+=1.0;
                        }
                }
        }
        cout <<"Found pr"<<endl;
        //Print3DVectorSize(pr, "pr");

        MinMaxPr(maxdistance, mindistance);
        if (GreatestDistance>250)
        {
                for (int n=0;n<NumAtomTypes;n++)
                {
                        for (int m=0;m<NumAtomTypes;m++)
                        {
                                cout <<pr[n][m][250]<<"\t";
                        }
                        cout <<endl;
                }
                cout <<endl;
        }

        if (GreatestDistance>300)
        {
                for (int n=0;n<NumAtomTypes;n++)
                {
                        for (int m=0;m<NumAtomTypes;m++)
                        {
                                cout <<pr[n][m][300]<<"\t";
                        }
                        cout <<endl;
                }
                cout <<endl;
        }
}

void InitializeSinc(Real bin, Real maxdistance[][NumAtomTypes], IntensityStruct &i)
{
        //Creates a look up table for the sinc function.  sinc(x)=sin(x)/x.
        int points1=i.calc.size();
        int points2=i.expi.size();
        int points;
        int max, distance;
        Real dist, product;
        vector<Real> temp;
        cout <<"In InitializeSinc"<<endl;
        max=0;
        points=getMax(points1, points2);
        i.sinc.clear();
        for (int m=0;m<NumAtomTypes;m++)
        {
                for (int n=0;n<NumAtomTypes;n++)
                {
                        cout <<"m= "<<m<<" n= "<<n<<endl;
                        if (maxdistance[m][n]>max)
                        {
                                cout <<"max= "<<max<<" maxdistance["<<m<<"]["<<n<<"]= "<<maxdistance[m][n]<<endl;
                                max=int(maxdistance[m][n]);
                        }
                }
        }
        cout <<"max= "<<max<<endl;
        for (int t=0;t<points;t++)
        {
                cout <<"t= "<<t<<endl;
                dist=0;
                for (distance=0;distance<=max+1;distance++)
                {
                        //cout <<"distance= "<<distance<<endl;
                        product=dist*i.s[t];
                        if (product!=0) SafePushBack(temp, sin(product)/product, "temp");
                        else SafePushBack(temp, Real(1.0), "temp");
                        dist+=bin;
                }
                cout <<"Before SafePushBackTemp"<<endl;
                SafePushBack(i.sinc, temp, "i.sinc");
                cout <<"After SafePushBackTemp"<<endl;
                temp.clear();
        }
}

void ApplyScatteringFactors(Real atomr[], IntensityStruct &i, Real contrast, Real hsdensity, ParamStruct params, int points)
{
        //Multiplies the scattering calculated using point particles by the 
        //scattering factors of the the atoms to obtain the scattering 
        //intensity.
        Real f[NumAtomTypes], fn;
        for (int n=0;n<NumAtomTypes;n++) f[n]=0;
        cout <<"In ApplyScatteringFactors"<<endl;
        cout <<"contrast= "<<contrast<<endl;
        for (int n=0;n<NumAtomTypes;n++) cout <<"atomr["<<n<<"]= "<<atomr[n]<<endl;
        cout <<"SphereType= "<<params.ExcludedVolumeSphereType<<endl;
        cout <<"HydrationShell atom Scattering"<<endl;
        cout <<"hsdensity= "<<hsdensity<<" Uniform= "<<params.UniformHydrationShell<<endl;
        for (int t=0;t<points;t++)
        {
                for (int n=0;n<NumAtomTypes;n++)
                {
                        f[n]=SolventCorrectedScattering(n, atomr[n], i.s[t], i.f[t][n], contrast, hsdensity, params);
                }
                cout <<i.s[t]<<"\t"<<i.f[t][HydrationShell]<<"\t"<<f[HydrationShell]<<endl;
                for (int n=0;n<NumAtomTypes;n++)
                {
                        fn=f[n];
                        for (int m=n;m<NumAtomTypes;m++)
                        {
                                i.CrossTerm[t][m][n]=i.point[t][m][n]*fn*f[m]*8.0*pi;
                                i.calc[t]+=i.CrossTerm[t][m][n];
                                //if (t%100==0)
                                //{
                                //	cout <<"i.CrossTerm["<<t<<"]["<<m<<"]["<<n<<"]= "<<i.CrossTerm[t][m][n]<<"\t"
                                //		<<"i.point["<<t<<"]["<<m<<"]["<<n<<"]= "<<i.point[t][m][n]<<"\t"
                                //		<<"f["<<n<<"]= "<<f[n]<<"\t"
                                //		<<"i.calc["<<t<<"]= "<<i.calc[t]<<endl;
                                //}
                        }
                }
        }
}

void ZeroIntensity(IntensityStruct &i)
{
        //Sets the scattering intensities to zero.
        int points=i.calc.size();
        for (int t=0;t<points;t++)
        {
                i.calc[t]=0;
                for (int n=0;n<NumAtomTypes;n++)
                {
                        for (int m=0;m<NumAtomTypes;m++)
                        {
                                i.point[t][n][m]=0;
                                i.CrossTerm[t][n][m]=0;
                        }
                }
        }
}

void Debye(Real atomr[], IntensityStruct &i, Real contrast, Real hsdensity, Real bin, int points, Real maxdistance[NumAtomTypes][NumAtomTypes], Real mindistance[NumAtomTypes][NumAtomTypes], vector<AtomStruct> &Atoms, ParamStruct params)
{
        //Calculates the scattering intensity using the Debye formula.
        int TotalParticles=Atoms.size();
        int distance;
        Real dist;
        Real product;
        Real numatom[NumAtomTypes];

        cout <<"Allocated memory"<<endl;
        FindNumAtom(numatom, Atoms);
        ZeroIntensity(i);	

        cout <<"Initialized CrossTerm"<<endl;
        cout <<"TotalParticles= "<<TotalParticles<<endl;
        Print3DVectorSize(i.point, "i.point");
        Print3DVectorSize(pr, "pr");
        Print2DVectorSize(i.sinc, "i.sinc");
        if (i.sinc.size() == 0)
        {
                for (int t=0;t<points;t++)
                {
                        for (int n=0;n<NumAtomTypes;n++)
                        {
                                for (int m=n;m<NumAtomTypes;m++)
                                {
                                        dist=mindistance[n][m]*bin;
                                        for (distance=int(mindistance[n][m]);distance<=maxdistance[n][m];distance++)
                                        {
                                                product=i.s[t]*dist;
                                                i.point[t][m][n]+=pr[n][m][distance]*sin(product)/(product);
                                                dist+=bin;
                                        }
                                }
                                i.point[t][n][n]+=numatom[n]*0.5;
                        }
                }
        }
        else
        {
                for (int t=0;t<points;t++)
                {
                        for (int n=0;n<NumAtomTypes;n++)
                        {
                                for (int m=n;m<NumAtomTypes;m++)
                                {
                                        for (distance=int(mindistance[n][m]);distance<=maxdistance[n][m];distance++)
                                        {
                                                //if (n==0 && m==0 && distance>9990)
                                                //{
                                                //	cout <<"i.point["<<t<<"]["<<m<<"]["<<n<<"]= "<<i.point[t][m][n]<<"\t"
                                                //	<<"pr["<<n<<"]["<<m<<"]["<<distance<<"]= "<<pr[n][m][distance]<<"\t"
                                                //	<<"i.sinc["<<t<<"]["<<distance<<"]= "<<i.sinc[t][distance]<<endl;
                                                //}
                                                //cout <<"distance= "<<distance<<" t= "<<t<<" m= "<<m<<" n= "<<n<<endl;
                                                i.point[t][m][n]+=pr[n][m][distance]*i.sinc[t][distance];
                                        }
                                }
                                i.point[t][n][n]+=numatom[n]*0.5;
                        }
                }
        }
        ApplyScatteringFactors(atomr, i, contrast, hsdensity, params, points);
}

void nohistogram(Real atomr[], IntensityStruct &i, Real contrast, Real hsdensity, Real bin, int points, Real numatom[], vector<AtomStruct> &Atoms, ParamStruct params )
{
        //Calculates the scattering intensity using the Debue formula without 
        //histograming distances.
        int TotalParticles=Atoms.size();
        int atomidm, atomidn;
        Real dist;
        Real fn;
        Real xn, yn, zn;
        Real dx, dy, dz;
        Real dists;
        Real f[1000][NumAtomTypes];


        for (int t=0;t<points;t++)
        {
                i.calc[t]=0;
                for (int n=0;n<NumAtomTypes;n++)
                {
                        f[t][n]=0;
                }
        }

        for (int t=0;t<points;t++)
        {
                for (int n=0;n<NumAtomTypes;n++)
                {
                        f[t][n]=SolventCorrectedScattering(n, atomr[n], i.s[t], i.f[t][n], contrast, hsdensity, params);
                }
        }

        for (int n=0;n<TotalParticles-1;n++)
        {
                atomidn=Atoms[n].atomid;
                xn=Atoms[n].x;
                yn=Atoms[n].y;
                zn=Atoms[n].z;
                for (int m=n+1;m<TotalParticles;m++)
                {
                        atomidm=Atoms[m].atomid;
                        dx=xn-Atoms[m].x;
                        dy=yn-Atoms[m].y;
                        dz=zn-Atoms[m].z;
                        dist=sqrt(dx*dx+dy*dy+dz*dz);
                        for (int t=0;t<points;t++)
                        {
                                dists=dist*i.s[t];
                                i.calc[t]+=f[t][atomidn]*f[t][atomidm]*sin(dists)/(dists);
                        }
                }
        }

        for (int t=0;t<points;t++)
        {
                for (int n=0;n<NumAtomTypes;n++)
                {
                        fn=f[t][n];
                        i.calc[t]+=fn*fn*numatom[n]*0.5;
                }
        }
}

void nohistogram2(Real atomr[], IntensityStruct &i, Real contrast, Real hsdensity, Real bin, int points, Real numatom[], vector<AtomStruct> &Atoms, ParamStruct params)
{
        //Calculates the scattering intensity using the Debue formula without 
        //histograming distances.  Each atom has a different atomic radius.
        int TotalParticles=Atoms.size();
        int natom;
        Real mindist;
        Real dist;
        Real xn, yn, zn;
        Real dx, dy, dz;
        Real dists;
        Real iatomrn;
        Real volume;
        Real totalvolume;
        Real scale;
        Real *iatomr;
        Real **f;

        Safe2DArrayAlloc(f, 1000, TotalParticles+1, "f");
        SafeArrayAlloc(iatomr, TotalParticles, "iatomr");

        natom=0;
        for (int n=0;n<6;n++)
        {
                natom+=int(numatom[n]);
        }

        for (int n=0;n<natom;n++)
        {
                iatomr[n]=atomr[Atoms[n].atomid];
        }

        for (int k=0;k<100;k++)
        {
                for (int n=0;n<natom;n++)
                {
                        xn=Atoms[n].x;
                        yn=Atoms[n].y;
                        zn=Atoms[n].z;
                        mindist=1000;
                        for (int m=0;m<natom;m++)
                        {
                                if (m!=n)
                                {
                                        dx=xn-Atoms[m].x;
                                        dy=yn-Atoms[m].y;
                                        dz=zn-Atoms[m].z;
                                        dist=dx*dx+dy*dy+dz*dz;
                                        if (dist<100)
                                        {
                                                dist=sqrt(dist)-iatomr[m];
                                                if (dist<mindist)
                                                {
                                                        mindist=dist;
                                                }
                                        }
                                }
                        }
                        iatomr[n]=(mindist+iatomr[n])*0.5;
                }
        }

        volume=findvolume(natom, atomr, numatom, Atoms);

        totalvolume=0;

        for (int n=0;n<natom;n++)
        {
                iatomrn=iatomr[n];
                totalvolume+=5.568328*iatomrn*iatomrn*iatomrn;
        }

        scale=exp(log(volume/totalvolume)/3.0);
        for (int n=0;n<natom;n++)
        {
                iatomr[n]=iatomr[n]*scale;
        }

        for (int t=0;t<points;t++)
        {
                i.calc[t]=0;
                for (int n=0;n<TotalParticles;n++)
                {
                        f[t][n]=0;
                }
        }

        for (int t=0;t<points;t++)
        {
                for (int n=0;n<TotalParticles;n++)
                {
                        f[t][n]=SolventCorrectedScattering(n, iatomr[n], i.s[t], i.f[t][n], contrast, hsdensity, params);
                }
        }

        for (int n=0;n<TotalParticles-1;n++)
        {
                xn=Atoms[n].x;
                yn=Atoms[n].y;
                zn=Atoms[n].z;
                for (int m=n+1;m<TotalParticles;m++)
                {
                        dx=xn-Atoms[m].x;
                        dy=yn-Atoms[m].y;
                        dz=zn-Atoms[m].z;
                        dist=sqrt(dx*dx+dy*dy+dz*dz);
                        for (int t=0;t<points;t++)
                        {
                                dists=dist*i.s[t];
                                i.calc[t]+=f[t][n]*f[t][m]*sin(dists)/(dists);
                        }
                }
        }

        for (int n=0;n<TotalParticles;n++)
        {
                for (int t=0;t<points;t++)
                {
                        i.calc[t]+=f[t][n]*f[t][n]*0.5;
                }
        }
}

inline Real Distance(Real x1, Real y1, Real z1, Real x2, Real y2, Real z2)
{
        //Returns the distance between two points.
        Real dx, dy, dz;
        dx=x2-x1;
        dy=y2-y1;
        dz=z2-z1;
        return sqrt(dx*dx+dy*dy+dz*dz);
}

void charge(Real atomr[], IntensityStruct &i, Real contrast, Real hsdensity, Real bin, int points, Real numatom[], vector<AtomStruct> &Atoms, ParamStruct params)
{
        //Takes into account the partial charges of the atoms when calculating
        //scattering intensity.  This function probably does not work.
        int TotalParticles=Atoms.size();
        int atomidm, atomidn;
        Real dist;
        Real dists;
        Real xn, yn, zn;
        Real partialchargen, partialchargem;
        Real f[1000][NumAtomTypes], g[1000][NumAtomTypes];

        for (int t=0;t<points;t++)
        {
                i.calc[t]=0;
                for (int n=0;n<NumAtomTypes;n++)
                {
                        f[t][n]=0;
                        g[t][n]=0;
                }
        }

        for (int t=0;t<points;t++)
        {
                for (int n=0;n<NumAtomTypes;n++)
                {
                        f[t][n]=SolventCorrectedScattering(n, atomr[n], i.s[t], i.f[t][n], 0, hsdensity, params);
                        //g[t][n]=-SphericalScattering(atomr[n], i.s[t], i.f[t][n], contrast, params);
                }
        }

        for (int n=0;n<TotalParticles-1;n++)
        {
                xn=Atoms[n].x;
                yn=Atoms[n].y;
                zn=Atoms[n].z;
                atomidn=Atoms[n].atomid;
                partialchargen=Atoms[n].charge;
                for (int m=n+1;m<TotalParticles;m++)
                {
                        dist=Distance(xn, yn, zn, Atoms[m].x, Atoms[m].y, Atoms[m].z); 
                        atomidm=Atoms[m].atomid;
                        partialchargem=Atoms[m].charge;
                        for (int t=0;t<points;t++)
                        {
                                dists=dist*i.s[t];
                                i.calc[t]+=(f[t][atomidn]*(partialchargen*(3.0-i.s[t])+i.s[t])/3.0+g[t][atomidn])*(f[t][atomidm]*(partialchargem*(3.0-i.s[t])+i.s[t])/3.0+g[t][atomidm])*sin(dists)/(dists);
                        }
                }
        }

        for (int t=0;t<points;t++)
        {
                for (int n=0;n<TotalParticles;n++)
                {
                        atomidn=Atoms[n].atomid;
                        partialchargen=Atoms[n].charge;
                        i.calc[t]+=(f[t][atomidn]*(partialchargen*(3.0-i.s[t])+i.s[t])/3.0+g[t][atomidn])*(f[t][atomidn]*(partialchargen*(3.0-i.s[t])+i.s[t])/3.0+g[t][atomidn])*0.5;
                }
        }
}

void SphericalHarmonic(Real phi, Real theta, int lvalue, int mvalue, Real &RealOut, Real &ImagOut)
{
        //Evaluates a spherical harmonic.  Untested.

        int l;
        Real cthetan;
        vector<Real> r;
        Real phin;
        Real sthetan;
        Real stheta2;
        Real ylm;
        Real wlm[100][100];
        Real pow[100];
        Real sign;
        long double longl;
        long double longm;
        long double sqrtm[100];
        long double invsqrtm[100];

        for (int n=0;n<100;n++)
        {
                pow[n]=0;
                sqrtm[n]=0;
                invsqrtm[n]=0;
                for (int m=0;m<100;m++)
                {
                        wlm[n][m]=0;
                }
        }

        longl=1.0;
        wlm[0][0]=1.0;
        for (l=1;l<lvalue+1;l++)
        {
                wlm[l][l]=wlm[l-1][l-1]*(1.0-1.0/(2.0*longl));
                longl=longl+1.0;
        }

        longl=1.0;
        sign=-1.0;

        for (l=1;l<lvalue+1;l++)
        {
                wlm[l][l]=sign*sqrt((2*longl+1.0)*wlm[l][l]/(4.0*pi));
                longl=longl+1.0;
                sign=-sign;
        }

        wlm[0][0]=wlm[0][0]/sqrt(4.0*pi);

        pow[0]=1.0;

        longl=0;

        longm=0;

        for (int m=0;m<lvalue-1;m++)
        {
                sqrtm[m]=sqrt((longl+longm+2.0)*(longl-longm-1.0));
                invsqrtm[m]=1.0/sqrt((longl-longm)*(longl+longm+1.0));
                longm=longm+1.0;
        }

        stheta2=sthetan*sthetan;

        if (l>0)
        {
                wlm[l][l-1]=-sqrt(2.0*longl)*cthetan*wlm[l][l];
        }

        for (int t=1;t<=lvalue;t++)
        {
                pow[t]=pow[t-1]*sthetan;
        }

        longm=2.0*longl-2.0;
        for (int m=lvalue-2;m>-lvalue;m-=1)
        {
                wlm[lvalue][m]=-(sqrtm[m]*stheta2*wlm[l][m+2]+longm*cthetan*wlm[l][m+1])*invsqrtm[m];
                longm=longm-2.0;
        }


        ylm=wlm[lvalue][mvalue]*pow[mvalue];
        longm=phin*Real(mvalue);
        RealOut=ylm*cos(longm);
        ImagOut=ylm*sin(longm);
}

void FitSphericalHarmonics(Real numatom[], int natom, Real atomr[], IntensityStruct &i, Real contrast, Real hsdensity, int lmax, int points, 
                vector<AtomStruct> &Atoms, ParamStruct params)
{
        //Finds the coefficients of spherical harmonics that form a surface.

        int l, TotalParticles=Atoms.size();
        Real rn, dr, besn;
        Real strn;
        Real ylm;
        Real reallmn, imaglmn;
        Real NumElectrons[NumAtomTypes];
        Real *r, *phi, *stheta, *ctheta, *F;
        Real ***real, ***imag, ***breal, ***bimag;
        Real wlm[100][100];
        long double cthetan, sthetan, stheta2, phin;
        long double pow[100];
        long double sign;
        long double longl;
        long double longm;
        long double sqrtm[100];
        long double invsqrtm[100];
        long double sterm[150];

        lmax=3;
        dr=0.1;

        for (int n=0;n<100;n++)
        {
                pow[n]=0;
                sqrtm[n]=0;
                invsqrtm[n]=0;
                for (int m=0;m<100;m++)
                {
                        wlm[n][m]=0;
                }
        }
        SetNumElectrons(NumElectrons);
        TotalParticles=CrysolSolvent(numatom, NumElectrons, natom, contrast, atomr, Atoms, params);

        Safe3DArrayAlloc(real, lmax+2, lmax+2, TotalParticles+2, "real");
        Safe3DArrayAlloc(imag, lmax+2, lmax+2, TotalParticles+2, "imag");
        Safe3DArrayAlloc(breal, lmax+2, lmax+2, points+2, "breal");
        Safe3DArrayAlloc(bimag, lmax+2, lmax+2, points+2, "bimag");

        SafeArrayAlloc(r, TotalParticles, "r");
        SafeArrayAlloc(phi, TotalParticles, "phi");
        SafeArrayAlloc(stheta, TotalParticles, "stheta");
        SafeArrayAlloc(ctheta, TotalParticles, "ctheta");
        SafeArrayAlloc(F, TotalParticles, "F");

        for (int n=0;n<TotalParticles;n++)
        {
                //F[n]=0;
        }

        for (l=0;l<=lmax;l++)
        {
                for (int m=0;m<=lmax;m++)
                {
                        for (int n=0;n<TotalParticles;n++)
                        {
                                real[l][m][n]=0;
                                imag[l][m][n]=0;
                        }

                        for (int n=0;n<points;n++)
                        {
                                breal[l][m][n]=0;
                                bimag[l][m][n]=0;
                        }
                }
        }

        cout <<"Initialized real"<<endl;

        cart2sphere(r, phi, stheta, ctheta, Atoms);

        longl=1.0;
        wlm[0][0]=1.0;
        for (l=1;l<lmax;l++)
        {
                wlm[l][l]=wlm[l-1][l-1]*(1.0-1.0/(2.0*longl));
                longl=longl+1.0;
        }

        longl=1.0;
        sign=-1.0;

        for (l=1;l<lmax;l++)
        {
                wlm[l][l]=sign*sqrt((2*longl+1.0)*wlm[l][l]/(4.0*LongPi));
                cout <<"wlm["<<l<<"]["<<l<<"]= "<<wlm[l][l]<<endl;
                longl=longl+1.0;
                sign=-sign;
        }

        wlm[0][0]=wlm[0][0]/sqrt(4.0*LongPi);

        pow[0]=1.0;
        longl=0;
        for (l=0;l<lmax;l++)
        {
                longm=0;

                for (int m=0;m<l-1;m++)
                {
                        sqrtm[m]=sqrt((longl+longm+2.0)*(longl-longm-1.0));
                        invsqrtm[m]=1.0/sqrt((longl-longm)*(longl+longm+1.0));
                        longm=longm+1.0;
                }

                for (int n=natom;n<TotalParticles;n++)
                {
                        cthetan=ctheta[n];
                        sthetan=stheta[n];
                        phin=phi[n];

                        stheta2=sthetan*sthetan;

                        if (l>0)
                        {
                                wlm[l][l-1]=-sqrt(2.0*longl)*cthetan*wlm[l][l];
                        }

                        for (int t=1;t<l+1;t++)
                        {
                                pow[t]=pow[t-1]*sthetan;
                        }

                        longm=2.0*longl-2.0;
                        for (int m=l-2;m>-1;m=m-1)
                        {
                                wlm[l][m]=-(sqrtm[m]*stheta2*wlm[l][m+2]+longm*cthetan*wlm[l][m+1])*invsqrtm[m];
                                longm-=2.0;
                        }

                        longm=0.0;
                        for (int m=0;m<l+1;m++)
                        {
                                ylm=wlm[l][m]*pow[m];
                                real[l][m][n]=ylm*cos(longm);
                                real[l][m][n]=ylm*sin(longm);
                                longm+=phin;
                        }
                }
                longl+=1.0;
        }

        for (l=0;l<lmax;l++)
        {
                for (int m=0;m<=l;m++)
                {
                        for (int n=natom;n<TotalParticles;n++)
                        {
                                reallmn=real[l][m][n];
                                for (rn=r[n];rn<r[n]+3.0;rn+=dr)
                                {
                                        for (int t=0;t<points;t++)
                                        {
                                                strn=i.s[t]*rn;
                                                if (strn<20) besn=bessel(strn, l, sterm);
                                                else besn=bessel2(strn, l);
                                                besn=besn*rn*rn;
                                                breal[l][m][t]+=besn;
                                        }
                                }

                                reallmn=real[l][m][n];
                                imaglmn=imag[l][m][n];
                                for (int t=0;t<points;t++)
                                {
                                        breal[l][m][t]+=breal[l][m][t]*reallmn;
                                        bimag[l][m][t]+=breal[l][m][t]*imaglmn;
                                }
                        }
                }
        }

        for (l=0;l<lmax;l++)
        {
                for (int m=0;m<=l;m++)
                {
                        for (int t=0;t<points;t++)
                        {
                                breal[l][m][t]=breal[l][m][t]*dr/Real(TotalParticles-natom);
                                bimag[l][m][t]=breal[l][m][t]*dr/Real(TotalParticles-natom);
                                i.calc[t]+=(breal[l][m][t]*breal[l][m][t]+bimag[l][m][t]*bimag[l][m][t]);
                        }
                }
        }

        /****************************************************************************************
          for (int n=natom;n<TotalParticles;n++)
          {
        //cout <<"r["<<n<<"]= "<<r[n]<<" F["<<n<<"]= "<<F[n]<<endl;
        }



        for (int n=natom;n<TotalParticles;n++)
        {
        PrintPdbLine("/home/jouko/WAXS/SphericalHarmonics.pdb", n, SH, SPH, 1, F[n]*stheta[n]*cos(phi[n]), F[n]*stheta[n]*sin(phi[n]), F[n]*ctheta[n], 1.0, 0.0);
        }
        ofstream pdb2("/home/jouko/WAXS/SphericalHarmonics.pdb", ios::app);
        pdb2  << "END";



        for (int n=natom;n<TotalParticles;n++)
        {
        l=29;
        m=10;

        Index=(l+1)*(l+2)/2-l+m-1;
        PrintPdbLine("/home/jouko/WAXS/SphericalHarmonics_29_10.pdb", n, SH, SPH, 1, 10.0*real[Index][n]*stheta[n]*cos(phi[n]), 10.0*real[Index][n]*stheta[n]*sin(phi[n]), 10.0*real[Inedex][n]*ctheta[n], 1.0, 0.0);
        }
        ofstream pdb3("/home/jouko/WAXS/SphericalHarmonics_29_10.pdb", ios::app);
        pdb3  << "END";
         ********************************************************************************************/
}

void multipole(Real atomr[], IntensityStruct &i, Real contrast, Real hsdensity, int points, vector<AtomStruct> Atoms, ParamStruct params)
{
        //Finds intensity using multipole expansion this is the same method as
        //CRYSOL.
        //Spherical harmonics calculated using algorithm form:
        //Libbrecht, K.G. 1985. Practical Considerations for the Generation of Large-Order Spherical Harmonics
        //This needs to be tested.
        int l, TotalParticles=Atoms.size();
        int atomidn, lmax, natom;
        Real cthetan, sthetan, stheta2, phin;
        Real rn, dr, reallmn, imaglmn;
        Real strn;
        Real f[2000][NumAtomTypes];
        Real *r, *ctheta, *stheta, *phi;
        Real **bess;
        Real ***real, ***imag, ***areal, ***aimag, ***breal, ***bimag;
        Real besn;
        Real ylm;
        Real wlm[100][100];
        Real pow[100];
        Real sign;
        long double longl;
        long double longm;
        long double j[100];
        long double sqrtm[100];
        long double invsqrtm[100];
        long double sterm[150];
        time_t seconds;
        time_t seconds2;

        natom=FindNumProteinAtoms(Atoms);
        lmax=params.lmax;

        cout <<"natom= "<<natom<<" TotalParticles= "<<TotalParticles<<endl;

        SafeArrayAlloc(r, TotalParticles, "r");
        SafeArrayAlloc(ctheta, TotalParticles, "ctheta");
        SafeArrayAlloc(stheta, TotalParticles, "stheta");
        SafeArrayAlloc(phi, TotalParticles, "phi");

        Safe2DArrayAlloc(bess, lmax+2, points+2, "bess");

        cout <<"Initialized 1"<<endl;
        Safe3DArrayAlloc(real, lmax+2, lmax+2, TotalParticles+2, "real");
        Safe3DArrayAlloc(imag, lmax+2, lmax+2, TotalParticles+2, "imag");
        Safe3DArrayAlloc(areal, lmax+2, lmax+2, points+2, "areal");
        Safe3DArrayAlloc(aimag, lmax+2, lmax+2, points+2, "aimag");
        Safe3DArrayAlloc(breal, lmax+2, lmax+2, points+2, "breal");
        Safe3DArrayAlloc(bimag, lmax+2, lmax+2, points+2, "bimag");


        cout <<"Initialized 2"<<endl;


        for (l=0;l<=lmax;l++)
        {
                for (int m=0;m<=lmax;m++)
                {
                        for (int n=0;n<TotalParticles;n++)
                        {
                                real[l][m][n]=0;
                                imag[l][m][n]=0;
                        }

                        for (int n=0;n<points;n++)
                        {
                                areal[l][m][n]=0;
                                aimag[l][m][n]=0;
                                breal[l][m][n]=0;
                                bimag[l][m][n]=0;
                        }
                }
        }

        cout <<"Initialized 4"<<endl;

        for (l=0;l<=lmax;l++)
        {
                for (int t=0;t<points;t++)
                {
                        bess[l][t]=0;
                }
        }

        for (int n=0;n<100;n++)
        {
                pow[n]=0;
                sqrtm[n]=0;
                invsqrtm[n]=0;
                for (int m=0;m<100;m++)
                {
                        wlm[n][m]=0;
                }
        }

        cout <<"Initialized 5"<<endl;

        l=0;

        for (int n=0;n<100;n++) j[n]=0;

        for (int n=0;n<150;n++) sterm[n]=0;

        for (int t=0;t<points;t++)
        {
                i.calc[t]=0;
                for (int n=0;n<NumAtomTypes;n++)
                {
                        f[t][n]=0;
                }
        }

        cout <<"Initialized"<<endl;

        center(Atoms);

        cart2sphere(r, phi, stheta, ctheta, Atoms);

        cout <<"Centered and converted"<<endl;
        dr=0.1;
        if (hsdensity!=0)
        {
                cout <<"Normization= "<<dr*sqrt(2.0/pi)*hsdensity/Real(TotalParticles-natom)<<endl;
        }
        for (int t=0;t<points;t++)
        {
                for (int n=0;n<NumAtomTypes;n++)
                {
                        f[t][n]=SolventCorrectedScattering(n, atomr[n], i.s[t], i.f[t][n], contrast, hsdensity, params);
                }
        }

        longl=1.0;
        wlm[0][0]=1.0;
        for (l=1;l<=lmax;l++)
        {
                wlm[l][l]=wlm[l-1][l-1]*(1.0-1.0/(2.0*longl));
                longl=longl+1.0;
        }

        longl=1.0;
        sign=-1.0;
        for (l=1;l<=lmax;l++)
        {
                wlm[l][l]=sign*sqrt((2*longl+1.0)*wlm[l][l]/(4.0*LongPi));
                longl=longl+1.0;
                sign=-sign;
        }
        wlm[0][0]=wlm[0][0]/sqrt(4.0*LongPi);

        seconds=time(NULL);

        pow[0]=1.0;
        longl=0;
        for (l=0;l<=lmax;l++)
        {
                longm=1;

                for (int m=1;m<150;m++)
                {
                        sterm[m]=(longm+longl)/(2.0*(2.0*longm+2.0*longl+1.0)*(longm+longl)*longm);
                        longm=longm+1.0;
                }

                longm=0;

                for (int m=0;m<l-1;m++)
                {
                        sqrtm[m]=sqrt((longl+longm+2.0)*(longl-longm-1.0));
                        invsqrtm[m]=1.0/sqrt((longl-longm)*(longl+longm+1.0));
                        longm=longm+1.0;
                }

                for (int n=0;n<TotalParticles;n++)
                {
                        cthetan=ctheta[n];
                        phin=phi[n];
                        sthetan=stheta[n];
                        stheta2=sthetan*sthetan;
                        rn=r[n];
                        if (l>0)
                        {
                                wlm[l][l-1]=-sqrt(2.0*longl)*cthetan*wlm[l][l];
                        }

                        for (int t=1;t<l+1;t++)
                        {
                                pow[t]=pow[t-1]*sthetan;
                        }//This maybe were the error is

                        longm=2.0*longl-2.0;
                        for (int m=l-2;m>-1;m=m-1)
                        {
                                wlm[l][m]=-(sqrtm[m]*stheta2*wlm[l][m+2]+longm*cthetan*wlm[l][m+1])*invsqrtm[m];
                                longm=longm-2.0;
                        }

                        longm=phin;
                        for (int m=1;m<l+1;m++)
                        {
                                ylm=wlm[l][m]*pow[m];
                                if (n==0)
                                {
                                        cout <<"l= "<<l<<" m= "<<m<<" cthetan= "<<cthetan<<" phin= "<<phin<<" ylm= "<<ylm<<endl;
                                }
                                real[l][m][n]=ylm*cos(longm);
                                imag[l][m][n]=ylm*sin(longm);
                                longm+=phin;
                        }

                        real[l][0][n]=wlm[l][0];
                        imag[l][0][n]=0;
                }
                longl+=1.0;
        }

        cout <<"About to calculate areal"<<endl;
        longl=0;
        for (l=0;l<=lmax;l++)
        {
                longm=1;
                for (int m=1;m<150;m++)
                {
                        sterm[m]=(longm+longl)/(2.0*(2.0*longm+2.0*longl+1.0)*(longm+longl)*longm);
                        longm+=1.0;
                }
                longl+=1.0;

                for (int n=0;n<natom;n++)
                {
                        atomidn=Atoms[n].atomid;
                        rn=r[n];
                        for (int t=0;t<points;t++)
                        {
                                strn=i.s[t]*rn;
                                if (strn<20) besn=bessel(strn, l, sterm);
                                else besn=bessel2(strn, l);
                                if (n==0)
                                {
                                        cout <<"strn= "<<strn<<" l= "<<l<<" besn= "<<besn<<endl;
                                }
                                besn=besn*f[t][atomidn];
                                for (int m=0;m<=l;m++)
                                {
                                        areal[l][m][t]=areal[l][m][t]+besn*real[l][m][n];
                                        aimag[l][m][t]=aimag[l][m][t]+besn*imag[l][m][n];
                                }
                        }
                }
        }

        Real volume;
        volume=0;

        for (int n=natom;n<TotalParticles;n++)
        {
                for (rn=r[n];rn<r[n]+3.0;rn+=dr)
                {
                        volume+=rn*rn;
                }
        }

        if (hsdensity!=0)
        {
                cout <<"HydrationShell volume= "<<volume*4*pi*dr/Real(TotalParticles-natom)<<endl;
        }
        cout <<"About to calculate breal"<<endl;
        dr=0.1;
        longl=0;
        for (l=0;l<=lmax;l++)
        {
                longm=1;
                for (int m=1;m<150;m++)
                {
                        sterm[m]=(longm+longl)/(2.0*(2.0*longm+2.0*longl+1.0)*(longm+longl)*longm);
                        longm+=1.0;
                }
                longl+=1.0;
                cout <<"l= "<<l<<endl;
                for (int n=natom;n<TotalParticles;n++)
                {
                        for (rn=r[n];rn<r[n]+3.0;rn+=dr)
                        {
                                for (int t=0;t<points;t++)
                                {
                                        strn=i.s[t]*rn;
                                        if (strn<20) besn=bessel(strn, l, sterm);
                                        else besn=bessel2(strn, l);

                                        besn=besn*rn*rn;
                                        bess[l][t]+=besn;
                                }
                        }
                        for (int m=0;m<=l;m++)
                        {
                                reallmn=real[l][m][n];
                                imaglmn=imag[l][m][n];
                                for (int t=0;t<points;t++)
                                {
                                        breal[l][m][t]+=bess[l][t]*reallmn;
                                        bimag[l][m][t]+=bess[l][t]*imaglmn;
                                }
                        }
                }
        }

        for (l=0;l<=lmax;l++)
        {
                for (int t=0;t<points;t++)
                {
                        areal[l][0][t]=areal[l][0][t]*0.5;
                        breal[l][0][t]=breal[l][0][t]*0.5;
                }
        }

        for (l=0;l<=lmax;l++)
        {
                for (int m=0;m<=l;m++)
                {
                        for (int t=0;t<points;t++)
                        {
                                areal[l][m][t]=areal[l][m][t]*4.0*pi;
                                aimag[l][m][t]=aimag[l][m][t]*4.0*pi;
                                if (hsdensity!=0)
                                {
                                        breal[l][m][t]=breal[l][m][t]*dr*sqrt(2.0/pi)*hsdensity/Real(TotalParticles-natom);
                                        bimag[l][m][t]=breal[l][m][t]*dr*sqrt(2.0/pi)*hsdensity/Real(TotalParticles-natom);
                                }					
                                i.calc[t]+=( (areal[l][m][t]+breal[l][m][t])*(areal[l][m][t]+breal[l][m][t])+(aimag[l][m][t]+bimag[l][m][t])*(aimag[l][m][t]+bimag[l][m][t]) );
                                //i[t]+=(areal[l][m][t]*areal[l][m][t] + aimag[l][m][t]*aimag[l][m][t]);
                                //i[t]+=(breal[l][m][t]*breal[l][m][t] + bimag[l][m][t]*bimag[l][m][t]);
                        }
                }
        }

        for (int m=0;m<lmax+2;m++)
        {
                delete [] bess[m];
        }
        delete [] bess;

        for (int m=0;m<lmax+2;m++)
        {
                for (int n=0;n<lmax+2;n++)
                {
                        delete [] real[m][n];
                        delete [] imag[m][n];
                        delete [] areal[m][n];
                        delete [] aimag[m][n];
                        delete [] breal[m][n];
                        delete [] bimag[m][n];
                }

                delete [] real[m];
                delete [] imag[m];
                delete [] areal[m];
                delete [] aimag[m];
                delete [] breal[m];
                delete [] bimag[m];
        }

        delete [] real;
        delete [] imag;
        delete [] areal;
        delete [] aimag;
        delete [] breal;
        delete [] bimag;

        delete [] r;
        delete [] ctheta;
        delete [] stheta;
        delete [] phi;


        seconds2=time(NULL);
}

void implicitresidues()
{
        Real intensity;

        intensity=0;
}

void InitializeDistance(Real maxdistance[][NumAtomTypes])
{
        //Initializes the 2D Dmax array.  That is the maximum distances
        //between the different atom types.
        for (int m=0;m<NumAtomTypes;m++)
        {
                for (int n=0;n<NumAtomTypes;n++) maxdistance[m][n]=0;
        }
}

void ApplyPeriodicBoundaryConditions(vector<AtomStruct> &Atoms, ParamStruct &params)
{
        ApplyPeriodicBoundaryConditions(Atoms, params.XBoxLength, params.YBoxLength, params.ZBoxLength);
}

void watersphere(Real atomr[], IntensityStruct &i, Real contrast, Real hsdensity, Real bin, int points, Real numatom[], vector<AtomStruct> &Atoms, ParamStruct params)
{
        //This is not the primary function for using water sphere.
        //Elliminates atoms outside of a sphere and calculates the scattering
        //from the remaining atoms.  The difference between this and the other
        //water sphere function is that the excluded volume dummy atom is not
        //histogrammed.
        int m, TotalParticles=Atoms.size();
        Real *i2;
        Real fst;
        Real st;
        Real satomr;
        Real xn, yn, zn;
        Real r;
        Real rs;
        Real sradius, sradius2;
        Real maxdistance[NumAtomTypes][NumAtomTypes], mindistance[NumAtomTypes][NumAtomTypes];
        Real f[NumAtomTypes][2000];
        Real fs[2000];

        InitializeDistance(maxdistance);
        InitializeDistance(mindistance);

        SafeArrayAlloc(i2, points, "i2");

        for (int t=0;t<points;t++)
        {
                for (int n=0;n<NumAtomTypes;n++)
                {
                        f[n][t]=0;
                }
                fs[t]=0;
                i2[t]=0;
        }

        for (int t=0;t<points;t++)
        {
                for (int n=0;n<NumAtomTypes;n++)
                {
                        f[n][t]=SolventCorrectedScattering(n, atomr[n], i.s[t], i.f[t][n], 0, 0, params);
                }
        }

        if (!CommandLine)
        {
                cout <<"Enter desired resolution: ";
                cin >> bin;
                cout <<"Enter water ball radius: ";
                cin >> sradius;
        }
        else
        {
                sradius=params.SphereRadius;
        }

        m=0;
        sradius2=sradius*sradius;
        center(Atoms);
        ApplyPeriodicBoundaryConditions(Atoms, params);
        for (int n=0;n<TotalParticles;n++)
        {
                xn=Atoms[n].x;
                yn=Atoms[n].y;
                zn=Atoms[n].z;
                if ( (xn*xn+yn*yn+zn*zn) <sradius2 )
                {
                        Atoms[m].x=xn;
                        Atoms[m].y=yn;
                        Atoms[m].z=zn;
                        Atoms[m].atomid=Atoms[n].atomid;
                        m++;
                }
                else
                {
                        numatom[Atoms[n].atomid]-=1.0;
                }
        }

        for (int n=0;n<NumAtomTypes;n++)
        {
                cout <<"natom["<<n<<"]="<<numatom[n]<<endl;
        }

        findpr(bin, maxdistance, mindistance, Atoms, params.UniformHydrationShell);
        InitializeSinc(bin, maxdistance, i);
        Debye(atomr, i, 0, 0, bin, points, maxdistance, mindistance, Atoms, params);

        cout <<"i[0]="<<i.calc[0]<<endl;

        for (int t=0;t<points;t++)
        {
                st=i.s[t];
                satomr=st*sradius;
                fs[t]=-contrast*4.0*pi*(sin(satomr)-satomr*cos(satomr))/(st*st*st);
        }

        for (int t=0;t<points;t++)
        {
                fst=fs[t];
                i.calc[t]+=fst*fst*0.5;
        }
        cout <<"i[0]="<<i.calc[0]<<endl;

        for (int n=0;n<m;n++)
        {
                xn=Atoms[n].x;
                yn=Atoms[n].y;
                zn=Atoms[n].z;
                r=sqrt(xn*xn+yn*yn+zn*zn);
                for (int t=0;t<points;t++)
                {
                        rs=r*i.s[t];
                        i2[t]+=fs[t]*f[Atoms[n].atomid][t]*sin(rs)/rs;
                }
        }

        for (int t=0;t<points;t++) i.calc[t]+=i2[t];

        cout <<"i[0]="<<i.calc[0]<<endl;
}

Real ContinuumScattering(Real sradius, Real s, Real density[], Real contrast)
{
        //Calculates the scattering from an object with spherical symmetry
        //with an arbitrary density profile.
        Real f;
        Real distance;

        f=0;

        f=SphereScattering(sradius-2.0, contrast, s);
        distance=0;
        for (int n=19;n>-1;n--)
        {
                //cout <<"density["<<n<<"]="<<density[n]<<endl;
                f+=SphereScattering(sradius-2.0+distance+0.1, contrast-density[n], s);
                f-=SphereScattering(sradius-2.0+distance, contrast-density[n], s);
                distance+=0.1;
        }

        distance=0;

        for (int n=1;n<20;n++)
        {
                //cout <<"density["<<n<<"]="<<density[n]<<endl;
                f+=SphereScattering(sradius+distance+0.1, density[n], s);
                f-=SphereScattering(sradius+distance, density[n], s);
                distance+=0.1;
        }

        return f;

}

void FindNumAtom(Real numatom[], vector<AtomStruct> &Atoms)
{
        //Counts the number of atoms belonging to each atom type.
        int TotalParticles=Atoms.size();
        for (int n=0;n<NumAtomTypes;n++) numatom[n]=0;
        for (int n=0;n<TotalParticles;n++) numatom[Atoms[n].atomid]+=Atoms[n].weight*Atoms[n].weight;
}

void MakeWaterSphere(Real atomr[], Real numatom[], vector<AtomStruct> &Atoms, ParamStruct params)
{
        //Elliminates waters outside of a specified radius.  This is for
        //calculating scattering from an explicit solvent MD simulation.
        AtomStruct Atom;
        bool FindCenterOfMass=true;
        char CharTotalElectrons[1000];
        int TotalParticles=Atoms.size();
        int test=4000;
        Real farthest;
        Real sradius, sradius2;
        Real xn, yn, zn;
        Real NumElectrons[NumAtomTypes];
        vector<AtomStruct> TempAtoms;
        cout <<"In MakeWaterSphere"<<endl;
        SetNumElectrons(NumElectrons);
        if (test<TotalParticles) 
        {
                PrintAtomInfo(Atoms[test]);
                cout <<"Atoms["<<test<<"].ChainName= "<<Atoms[test].ChainName<<endl;
        }
        for (int n=0;n<NumAtomTypes;n++) atomr[n]=0.0;
        GetResidueIDs(Atoms);
        atomr[WaterSphere]=params.SphereRadius;
        if (FindCenterOfMass) center(Atoms);
        cout <<"verbose= "<<verbose<<endl;
        if (verbose) cout <<"FoundCenterOfMass"<<endl;
        PrintAtomInfo(Atoms[0]);
        if (!FindCenterOfMass) MoveAtoms(Atoms, -params.XOrigin, -params.YOrigin, -params.ZOrigin);
        ApplyPeriodicBoundaryConditions(Atoms, params.XBoxLength, params.YBoxLength, params.ZBoxLength);
        if (params.SphereRadius*2.0>params.XBoxLength)
        {
                GetAtomsWithinRadius(Atoms, 0, 0, 0, params.SphereRadius);
                cout <<"params.SphereRadius= "<<params.SphereRadius<<endl;
                cout <<"XBoxLength= "<<params.XBoxLength<<" YBoxLength= "<<params.YBoxLength<<" ZBoxLength= "<<params.ZBoxLength<<endl<<endl;
                cout <<"Warning: Sphere diameter is larger than the box dimensions.  Using images."<<endl;
                Image(Atoms, X, params.XBoxLength, 3);
        }
        if (params.SphereRadius*2.0>params.YBoxLength)
        {
                GetAtomsWithinRadius(Atoms, 0, 0, 0, params.SphereRadius);
                cout <<"params.SphereRadius= "<<params.SphereRadius<<endl;
                cout <<"XBoxLength= "<<params.XBoxLength<<" YBoxLength= "<<params.YBoxLength<<" ZBoxLength= "<<params.ZBoxLength<<endl<<endl;
                cout <<"Warning: Sphere diameter is larger than the box dimensions.  Using images."<<endl;
                Image(Atoms, Y, params.YBoxLength, 3);
        }
        if (params.SphereRadius*2.0>params.ZBoxLength)
        {
                GetAtomsWithinRadius(Atoms, 0, 0, 0, params.SphereRadius);
                cout <<"params.SphereRadius= "<<params.SphereRadius<<endl;
                cout <<"XBoxLength= "<<params.XBoxLength<<" YBoxLength= "<<params.YBoxLength<<" ZBoxLength= "<<params.ZBoxLength<<endl<<endl;
                cout <<"Warning: Sphere diameter is larger than the box dimensions.  Using images."<<endl;
                Image(Atoms, Z, params.ZBoxLength, 3);
        }
        GetAtomsWithinRadius(Atoms, 0, 0, 0, params.SphereRadius);
        PrintAtomInfo(Atoms[0]);
        string OutputFile="/home/jouko/project/WAXS/pdb/SphereTest_";
        OutputFile+=RealToStr(params.SphereRadius)+".pdb";
        //PrintPdb(OutputFile, Atoms);
        if (verbose) cout <<"Applied periodic boundary conditions"<<endl;
        Atom.x=0;
        Atom.y=0;
        Atom.z=0;
        Atom.atomid=WaterSphere;
        Atom.weight=1.0;
        Atom.ResidueName="SPH";
        Atom.residueid=0;
        Atom.AtomName="SPH";
        Atoms.resize(Atoms.size()+1);
        Atoms[Atoms.size()-1]=Atom;
        //SafePushBack(Atoms, Atom, "Atoms in MakeWaterSphere");
        sradius=params.SphereRadius;
        sradius2=sradius*sradius;
        farthest=0.0;
        Real TotalElectrons=0.0;
        TotalParticles=Atoms.size();
        cout <<"TotalParticles= "<<TotalParticles<<endl;
        for (int n=0;n<TotalParticles;n++)
        {
                //if (verbose) cout <<"n= "<<n<<endl;
                xn=Atoms[n].x;
                yn=Atoms[n].y;
                zn=Atoms[n].z;

                if ( (xn*xn+yn*yn+zn*zn) < sradius2 )
                {
                        if (Atoms[n].residueid>20) TotalElectrons+=NumElectrons[Atoms[n].atomid];
                        if (Atoms[n].residueid<21)
                        {
                                if (sqrt(xn*xn+yn*yn+zn*zn)>farthest)
                                {
                                        farthest=sqrt(xn*xn+yn*yn+zn*zn);
                                }
                        }
                        SafePushBack(TempAtoms, Atoms[n], "TempAtoms in MakeWaterSphere");
                }
                else if (Atoms[n].residueid<21) 
                {
                        PrintAtomInfo(Atoms[n]);
                        cout <<"Error protein cut"<<endl;
                }
        }
        CopyVector(TempAtoms, Atoms);
        cout <<"farthest= "<<farthest<<endl;

        strcpy(CharTotalElectrons, params.TotalElectronsFile.c_str());
        ofstream File(CharTotalElectrons, ios::app);

        File << TotalElectrons<<endl;

        cout <<"TotalParticles= "<<TotalParticles<<endl;

        FindNumAtom(numatom, Atoms);
        for (int n=0;n<NumAtomTypes;n++) cout <<"numatom["<<n<<"]= "<<numatom[n]<<endl;
        /*
           string PdbOut, pdb, RadiusStr;
           char CharPdbOut[1000];
           size_t position;
           position=params.xyz.rfind("/");
           pdb=params.xyz.substr(position+1, 100);
           stringstream out;
           out << params.SphereRadius;
           RadiusStr=out.str();
           PdbOut="/home/jouko/WAXS/Hydration_" + RadiusStr;
           PdbOut+=pdb;
           cout <<"PdbOut= "<<PdbOut<<endl;
           strcpy(CharPdbOut, PdbOut.c_str());
           PrintPdb(CharPdbOut, Atoms);
           */
}

void Image(vector<AtomStruct> &Atoms, ParamStruct &params, int NumPerEdge)
{
        //Takes a system and makes copies of it around the original system.
        AtomStruct Atom;
        int TotalParticles2;
        int TotalParticles=Atoms.size();
        int a, b, c;
        Real Xlength, Ylength, Zlength;
        vector<AtomStruct> TempAtoms;

        TotalParticles2=NumPerEdge*NumPerEdge*NumPerEdge*TotalParticles;
        Xlength=Real(NumPerEdge)*params.XBoxLength;
        Ylength=Real(NumPerEdge)*params.YBoxLength;
        Zlength=Real(NumPerEdge)*params.ZBoxLength;

        for (a=1;a<=NumPerEdge;a++)
        {
                for (b=1;b<=NumPerEdge;b++)
                {
                        for (c=1;c<=NumPerEdge;c++)
                        {
                                cout <<"a= "<<a<<" b= "<<b<<" c= "<<c<<endl;
                                for (int n=0;n<TotalParticles;n++)
                                {
                                        CopyAtom(Atoms[n], Atom);
                                        Atom.x=Atoms[n].x+Real(a)*params.XBoxLength;
                                        Atom.y=Atoms[n].y+Real(b)*params.YBoxLength;
                                        Atom.z=Atoms[n].z+Real(c)*params.ZBoxLength;
                                        SafePushBack(TempAtoms, Atom, "TempAtoms");
                                }
                        }
                }
        }
        CopyVector(TempAtoms, Atoms);

        cout <<"XBoxLength="<<params.XBoxLength<<" YBoxLength="<<params.YBoxLength<<" ZBoxLength="<<params.ZBoxLength<<endl;
        MoveAtoms(Atoms, -Real(NumPerEdge-1)*params.XBoxLength*0.5, -Real(NumPerEdge-1)*params.YBoxLength*0.5, -Real(NumPerEdge-1)*params.ZBoxLength*0.5);
        ApplyPeriodicBoundaryConditions(Atoms, Xlength, Ylength, Zlength);
        params.XBoxLength=Xlength;
        params.YBoxLength=Ylength;
        params.ZBoxLength=Zlength;
        cout <<"XBoxLength="<<params.XBoxLength<<" YBoxLength="<<params.YBoxLength<<" ZBoxLength="<<params.ZBoxLength<<endl;
        cout <<"TotalParticles= "<<Atoms.size()<<endl;
}

void OutputIntensity(IntensityStruct &i, string IntensityFile)
{
        char CharIntensityFile[1000];
        int points=i.calc.size();
        AddIndexToFile(IntensityFile);
        strcpy(CharIntensityFile, IntensityFile.c_str());
        ofstream file(CharIntensityFile, ios::app);
        file <<"s\tcalc\tCalcError\tExperiment\tError\tvacuum\tNoHydration"<<endl;
        for (int t=0;t<points;t++)
        {
                file <<i.s[t]<<"\t"<<i.calc[t]<<"\t"<<i.CalcError[t]
                        <<"\t"<<i.expi[t]<<"\t"<<i.error[t]<<"\t"
                        <<i.vacuum[t]<<"\t"<<i.NoHydration[t]<<endl;
        }
        file.close();
}

void trajectory(Real atomr[], Real solventatomr[], IntensityStruct &i, Real contrast, Real hsdensity, Real bin, int points, Real scale, ParamStruct params)
{
        //Calulates a WAXS pattern for an ensemble of structures from a 
        //trajectory.  Should add in ability to read tinker arc files.
        bool IsStructure;
        AtomStruct Atom;
        vector<AtomStruct> Atoms, PSFAtoms;
        string base, DCDFile="", dir, IntensityOutputFile;
        int nthstruct=1, nInDcd, nthDcd=1, NumStructures=0;
        int natom, NumCalculated=0;
        Real StdDev;
        Real maxdistance[NumAtomTypes][NumAtomTypes];
        Real mindistance[NumAtomTypes][NumAtomTypes];
        Real *TotalI, *TotalI2;
        Real numatom[NumAtomTypes];

        SafeArrayAlloc(TotalI, points, "TotalI");
        SafeArrayAlloc(TotalI2, points, "TotalI2");

        for (int t=0;t<points;t++)
        {
                TotalI[t]=0;
                TotalI2[t]=0;
                i.calc[t]=0;
        }

        if (verbose) cout <<"DCDFiles= "<<params.DCDFiles<<endl;
        if (params.PsfFile!="") ReadPsf(params.PsfFile, PSFAtoms);
        if (verbose) cout <<"params.nDCDFiles= "<<params.nDCDFiles<<endl;
        if (verbose) cout <<"nstructures="<<params.nstructures<<endl;
        dir=RemoveExtension(params.IntensityOutputFile) + "/";
        base=GetBase(params.IntensityOutputFile);
        IntensityOutputFile=dir + base;
        cout <<"dir= "<<dir<<" base= "<<base<<endl;
        cout <<"IntensityOutputFile= "<<IntensityOutputFile<<endl;
        while (true)
        {
                if (verbose) cout <<"NumStructures= "<<NumStructures<<" nstructures= "<<params.nstructures<<endl;
                if (params.PsfFile!="") CopyVector(PSFAtoms, Atoms);
                cout <<"Copied PSFAtoms"<<endl;
                IsStructure=GetNextFrame(Atoms, nthstruct, nInDcd, nthDcd, DCDFile, params.DCDFiles, params.PdbList, params);
                //for (int i=0;i<Atoms.size();i++) PrintAtomInfo(Atoms[i]);
                NumStructures++;
                if (!IsStructure) break;
                if (NumStructures>params.nstructures) break;
                if (verbose) cout <<"After IsStructure"<<endl;
                cout <<"NumStructures= "<<NumStructures<<" params.skip= "<<params.skip<<" mod= "<<NumStructures%params.skip<<endl;
                if (NumStructures%params.skip==0)
                {
                        preprocess(params, natom, numatom, atomr, Atoms, maxdistance, mindistance, i);
                        ChooseIntensityCalculation(natom, atomr, solventatomr, i, contrast, hsdensity, bin, points, numatom, maxdistance, mindistance, Atoms, params);
                        OutputIntensity(i, IntensityOutputFile);
                        for (int t=0;t<points;t++) 
                        {
                                TotalI[t]+=i.calc[t];
                                TotalI2[t]+=i.calc[t]*i.calc[t];
                        }
                        NumCalculated++;
                }
        }
        NumStructures--;
        for (int t=0;t<points;t++) 
        {
                i.calc[t]=TotalI[t]/Real(NumCalculated);
                StdDev=sqrt(TotalI2[t]/Real(NumCalculated)-TotalI[t]*TotalI[t]/(Real(NumCalculated)*Real(NumCalculated)));
                i.CalcError[t]=StdDev/sqrt(Real(NumCalculated));
        }
        delete [] TotalI;
        delete [] TotalI2;
}

void SetWeightsForMerge(vector<AtomStruct> &Atoms, Real NumStructures)
{
        //When superimposing structures on top of each other the weight
        //of each structure has to be one over the number of structures.
        int TotalParticles=Atoms.size();
        for (int j=0;j<TotalParticles;j++)
        {
                if (IsProtein(Atoms[j]) || isNA(Atoms[j])) Atoms[j].weight=1.0;
                else Atoms[j].weight=1.0/NumStructures;
        }
}

void Merge(vector<AtomStruct> &Atoms, ParamStruct params)
{
        //Superimposes snapshots from a trajectory on top of each other.
        //This is not useful.  I just used this as a test.
        bool IsStructure;
        AtomStruct Atom;
        vector<AtomStruct> PSFAtoms, TempAtoms;
        string base, DCDFile="", dir, IntensityOutputFile;
        int nthstruct=1, nInDcd, nthDcd=1, NumStructures=0;

        if (params.PsfFile!="") ReadPsf(params.PsfFile, PSFAtoms);
        while (true)
        {
                if (verbose) cout <<"NumStructures= "<<NumStructures<<" nstructures= "<<params.nstructures<<endl;
                if (params.PsfFile!="") CopyVector(PSFAtoms, TempAtoms);
                IsStructure=GetNextFrame(TempAtoms, nthstruct, nInDcd, nthDcd, DCDFile, params.DCDFiles, params.PdbList, params);
                //for (int i=0;i<Atoms.size();i++) PrintAtomInfo(Atoms[i]);
                if (NumStructures!=0) RemoveProteinAtoms(TempAtoms);
                NumStructures++;
                if (!IsStructure) break;
                if (NumStructures>params.nstructures) break;
                if (verbose) cout <<"After IsStructure"<<endl;
                GetResidueIDs(TempAtoms);
                AppendVector(TempAtoms, Atoms);
        }
        NumStructures--;
        SetWeightsForMerge(Atoms, Real(NumStructures));
}

Real normalize(vector<Real> &calc, vector<Real> expi, vector<Real> error, int beginfit, int endfit)
{
        //Calculates the fit between two scattering patterns.
        int points=calc.size();
        Real chisqr;
        Real isqr;
        Real crossterm;
        Real scale;
        Real sumsqr;
        Real * inverror2;

        SafeArrayAlloc(inverror2, endfit+1, "inverror2");

        for (int t=0;t<endfit+1;t++) inverror2[t]=0;

        //cout <<"error[0]="<<error[0]<<" error[1]="<<error[1]<<endl;

        for (int t=0;t<endfit+1;t++)
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

        for (int t=beginfit;t<endfit+1;t++)
        {
                //cout <<"In normalize.  i["<<t<<"]="<<i[t]<<endl;
                isqr+=calc[t]*calc[t]*inverror2[t];
                crossterm+=calc[t]*expi[t]*inverror2[t];
        }

        if (isqr!=0) scale=crossterm/isqr;
        else
        {
                scale=1.0;
                cout <<"Error scalling\n";
        }

        for (int t=0;t<points;t++) calc[t]=calc[t]*scale;

        for (int t=beginfit;t<endfit+1;t++)
        {
                sumsqr=sumsqr+(calc[t]-expi[t])*(calc[t]-expi[t])*inverror2[t];
        }

        chisqr=sumsqr/float(endfit-beginfit+1);

        delete [] inverror2;

        return chisqr;
}

Real radiusgyration(Real atomr[], Real contrast, Real hsdensity, vector<AtomStruct> &Atoms, ParamStruct params)
{
        //Calculates the solvent corrected radius of gyration from the
        //structure.
        int TotalParticles=Atoms.size();
        Real avex, avey, avez;
        Real mass;
        Real dist2;
        Real distsqr;
        Real rg;
        Real totalmass;
        Real dx, dy, dz;
        Real SolventCorrectedElectrons[NumAtomTypes], NumElectrons[NumAtomTypes];
        avex=0;
        avey=0;
        avez=0;
        mass=0;
        dist2=0;
        distsqr=0;
        totalmass=0;
        SetNumElectrons(NumElectrons);
        CalculateSolventCorrectedElectrons(SolventCorrectedElectrons, contrast, atomr, params);
        for (int i=0;i<NumAtomTypes;i++)
        {
                cout <<"SolventCorrectedElectrons["<<i<<"]= "<<SolventCorrectedElectrons[i]<<endl;
        }
        FindProteinCenter(NumElectrons, contrast, atomr, Atoms, params, avex, avey, avez);
        cout <<"avex= "<<avex<<" avey= "<<avey<<" avez= "<<avez<<endl;
        for (int n=0;n<TotalParticles;n++)
        {
                dx=avex-Atoms[n].x;
                dy=avey-Atoms[n].y;
                dz=avez-Atoms[n].z;
                dist2=dx*dx+dy*dy+dz*dz;
                totalmass+=SolventCorrectedElectrons[Atoms[n].atomid]*Atoms[n].weight;
                distsqr+=dist2*SolventCorrectedElectrons[Atoms[n].atomid]*Atoms[n].weight;
        }
        cout <<"avex= "<<avex<<" avey= "<<avey<<" avez= "<<avez<<endl;
        cout <<"distsqr= "<<distsqr<<" totalmass= "<<totalmass<<endl;
        if (distsqr/totalmass>0) rg=sqrt(distsqr/totalmass);
        else rg=0;

        return rg;
}

Real calcVacuumRg(vector<AtomStruct> &Atoms)
{
        Real density=0, hsdensity=0;
        Real atomr[NumAtomTypes];
        ParamStruct params;

        SetDefaultParameters(params);
        SetDefaultAtomicRadii(atomr, params);
        params.contrast=0;
        params.hsdensity=0;
        return radiusgyration(atomr, density, hsdensity, Atoms, params);
        
}

Real calcExcludedVolumeRg(Real atomr[], Real density, vector<AtomStruct> &Atoms, ParamStruct params)
{
        Real hsdensity=0;
        params.hsdensity=0.0;
        return radiusgyration(atomr, density, hsdensity, Atoms, params);
}

Real gyrationfromscattering(vector<Real> i, vector<Real> s)
{
        //Calculates the radius of gyration from the first two scattering
        //points.
        Real point1;
        Real point2;
        Real rg;

        point1=log(i[0]);
        point2=log(i[1]);

        rg=(point1-point2)/(s[1]*s[1]-s[0]*s[0]);

        if (rg>0) rg=sqrt(rg*3.0);
        else rg=0;

        return rg;
}

void FitPr4(int NumVariables, int points, vector<Real> s, vector<Real> expi, vector<Real> error, ParamStruct params)
{
        //Calculates p(r) from I(s) using indirect Fourier transformation.
        //I am not sure how well this is working.
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
        Real *PR, *BestPR;
        vector<Real> i;

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

        SafeArrayAlloc(coefficient, NumVariables, "coefficient");
        SafeArrayAlloc(BestCoefficient, NumVariables, "BestCoefficient");
        SafeArrayAlloc(PR, NumPrPoints, "PR");
        SafeArrayAlloc(BestPR, NumPrPoints, "BestPR");
        Safe2DArrayAlloc(imn, NumVariables, points, "imn");

        cout <<"Allocated memory"<<endl;
        for (int n=0;n<NumPrPoints;n++)
        {
                PR[n]=0;
                BestPR[n]=0;
        }

        for (int n=0;n<NumVariables;n++)
        {
                coefficient[n]=0;
                BestCoefficient[n]=0;
        }

        coefficient[0]=1.0;
        cout <<"Initialized arrays"<<endl;

        for (int m=0;m<NumVariables;m++)
        {
                for (int n=0;n<points;n++)
                {
                        dist=inc;
                        for (int t=1;t<NumPrPoints;t++)
                        {
                                imn[m][n]+=sin(pi*Real((m+1)*t)/Real(NumPrPoints-1))*sin(s[n]*dist)/(s[n]*dist);
                                dist+=inc;
                        }
                }
        }

        for (int m=0;m<NumVariables;m++)
        {
                for (int n=0;n<params.EndFit;n++)
                {
                        i[n]+=coefficient[m]*imn[m][n];
                }
        }	

        cout <<"beginfit= "<<params.BeginFit<<" endfit= "<<params.EndFit<<endl;
        chisqr_old=normalize(i, expi, error, params.BeginFit, params.EndFit);
        /**
          for (int n=0;n<NumVariables;n++);
          {
          cout <<"coefficeint["<<n<<"]= "<<coefficient[n]<<endl;
          }
          cout <<endl;
         **/
        NumIterations=0;
        while (NumIterations<MaxIterations)
        {
                RandNum=rand();
                pick=int(Real(RandNum*NumVariables)/Real(RAND_MAX));
                move=(Real(rand())/Real(RAND_MAX)-0.5)*MaxMove;
                NumIterations++;

                coefficient[pick]+=move;

                for (int n=0;n<NumPrPoints;n++)
                {
                        PR[n]=0;
                }

                for (int m=0;m<NumVariables;m++)
                {
                        for (int n=params.BeginFit;n<params.EndFit;n++)
                        {
                                i[n]+=coefficient[m]*imn[m][n];
                        }
                }

                chisqr=normalize(i, expi, error, params.BeginFit, params.EndFit);

                for (int m=0;m<NumVariables;m++)
                {
                        for (int n=0;n<NumPrPoints;n++)
                        {
                                PR[n]+=coefficient[m]*sin(pi*Real((m+1)*n)/Real(NumPrPoints-1));
                        }
                }
                TotalPr=0;
                for (int n=0;n<NumPrPoints;n++)
                {
                        TotalPr+=PR[n];
                }

                for (int n=0;n<NumPrPoints;n++)
                {
                        PR[n]=PR[n]/TotalPr;
                }

                TotalPenalty=0;
                for (int n=1;n<NumPrPoints-1;n++)
                {
                        derivative2=(PR[n-1]-2.0*PR[n]+PR[n+1])/(inc*inc);
                        TotalPenalty+=derivative2*derivative2;
                }
                Real calc=0;
                Real product, product2;
                for (int n=0;n<NumVariables;n++)
                {
                        product=pi*Real(n+1)/Dmax;
                        product2=product*product;
                        calc+=0.5*coefficient[n]*coefficient[n]*product2*product2*Dmax;
                }

                cout <<"chisqr= "<<chisqr<<" TotalPenalty= "<<TotalPenalty*params.PenaltyCoefficient<<" calc= "<<calc*params.PenaltyCoefficient<<endl;
                chisqr+=TotalPenalty*params.PenaltyCoefficient;

                TotalChisqr+=chisqr;
                TotalChisqrSqr+=chisqr*chisqr;


                if (chisqr<BestChisqr)
                {
                        cout <<"New best chisqr= "<<chisqr<<endl;
                        BestChisqr=chisqr;
                        for (int n=0;n<NumVariables;n++)
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

        for (int n=0;n<points;n++) i[n]=0;

        for (int m=0;m<NumVariables;m++)
        {
                for (int n=0;n<points;n++)
                {
                        i[n]+=BestCoefficient[m]*imn[m][n];
                }
        }

        normalize(i, expi, error, params.BeginFit, params.EndFit);

        for (int m=0;m<NumVariables;m++)
        {
                for (int n=0;n<NumPrPoints;n++)
                {
                        PR[n]+=BestCoefficient[m]*sin(pi*Real((m+1)*n)/Real(NumPrPoints-1));
                }
        }

        dist=0;
        for (int n=0;n<NumPrPoints;n++)
        {
                ofstream pofr("/home/jouko/WAXS/pofr_ubiquitin.txt", ios::app);
                pofr << dist<<"\t"<<PR[n]<<endl;
                dist+=inc;
        }

        for (int n=0;n<points;n++)
        {
                ofstream intensity("/home/jouko/WAXS/pofr_intensity.txt", ios::app);
                intensity <<s[n]<<"\t"<<i[n]<<"\t"<<expi[n]<<endl;
        }
}

Real FitPr5(int NumVariables, int points, int NumPrPoints, Real Dmax, vector<Real> s, vector<Real> i, vector<Real> expi, Real PR[], ParamStruct params)
{
        //Calculates p(r) from I(s) using indirect Fourier transformation.
        //I am not sure how well this is working.
        cout <<"In FitPr5"<<endl;
        Real chisqr;
        Real dist, inc;
        cout <<"**imn"<<endl;
        Real **imn;
        cout <<"*variable, Coefficient"<<endl;
        vector<long double> variable;
        vector< vector<long double> > Coefficient;
        cout <<"In FitPr5"<<endl;	
        inc=Dmax/Real(NumPrPoints-1);

        SafeAlloc(variable, NumVariables, "variable");
        CreateMatrix(Coefficient, NumVariables+1, NumVariables);
        Safe2DArrayAlloc(imn, NumVariables, points, "imn");
        cout <<"Allocated memory"<<endl;
        for (int m=0;m<NumVariables;m++)
        {
                for (int n=0;n<points;n++)
                {
                        imn[m][n]=0;
                }
        }
        cout <<"Initialized imn"<<endl;
        for (int m=0;m<500;m++)
        {
                for (int n=0;n<500;n++)
                {
                        Coefficient[m][n]=0;
                }
        }
        cout <<"Initialized Coefficient"<<endl;
        cout <<"Allocated memory"<<endl;
        for (int n=0;n<NumPrPoints;n++)
        {
                PR[n]=0;
        }
        cout <<"Initialized PR"<<endl;
        for (int n=0;n<NumVariables;n++)
        {
                variable[n]=0;
        }

        cout <<"Initialized arrays"<<endl;

        for (int m=0;m<NumVariables;m++)
        {
                for (int n=0;n<points;n++)
                {
                        dist=inc;
                        for (int t=1;t<NumPrPoints;t++)
                        {
                                imn[m][n]+=sin(pi*Real((m+1)*t)/Real(NumPrPoints-1))*sin(s[n]*dist)/(s[n]*dist);
                                dist+=inc;
                        }
                }
        }

        for (int m=0;m<NumVariables;m++)
        {
                for (int n=0;n<points;n++)
                {
                        imn[m][n]*=inc;
                }
        }

        for (int m=0;m<NumVariables;m++)
        {
                for (int n=0;n<NumVariables;n++)
                {
                        for (int t=params.BeginFit;t<params.EndFit;t++)
                        {
                                Coefficient[m][n]+=imn[m][t]*imn[n][t];
                        }	
                }
        }

        //for (int n=0;n<NumVariables;n++)
        //{
        //	cout <<"Coefficient["<<n<<"]["<<n<<"]= "<<Coefficient[n][n]<<endl;
        //}
        //cout <<endl;

        for (int n=0;n<NumVariables;n++)
        {
                for (int t=params.BeginFit;t<params.EndFit;t++)
                {
                        Coefficient[NumVariables][n]+=expi[t]*imn[n][t];
                }
        }

        for (int n=0;n<NumVariables;n++)
        {
                //Coefficient[n][n]+=Real((n+1)*(n+1)*(n+1)*(n+1))*PenaltyCoefficient;
                Coefficient[n][n]+=Real((n+1)*(n+1))*params.PenaltyCoefficient/Dmax;
        }

        //for (int n=0;n<NumVariables;n++)
        //{
        //        cout <<"Coefficient["<<n<<"]["<<n<<"]= "<<Coefficient[n][n]<<endl;
        //}
        //cout <<endl;


        Equation(Coefficient, variable);

        for (int n=0;n<points;n++) i[n]=0;

        for (int m=0;m<NumVariables;m++)
        {
                for (int n=0;n<points;n++)
                {
                        i[n]+=variable[m]*imn[m][n];
                }
        }

        chisqr=0;
        for (int n=params.BeginFit;n<params.EndFit;n++)
        {
                chisqr+=(expi[n]-i[n])*(expi[n]-i[n]);
        }
        chisqr/=(params.EndFit-params.BeginFit+1);
        //chisqr=normalize(i, expi, beginfit, endfit, points, error);

        Real calc=0;
        for (int n=0;n<NumVariables;n++)
        {
                calc+=variable[n]*variable[n]*Real((n+1)*(n+1))*params.PenaltyCoefficient/Dmax;
        }

        cout <<"chisqr= "<<chisqr<<" penalty= "<<calc<<" total= "<<chisqr+calc<<endl;
        chisqr+=calc;

        for (int m=0;m<NumVariables;m++)
        {
                for (int n=0;n<NumPrPoints;n++)
                {
                        PR[n]+=variable[m]*sin(pi*Real((m+1)*n)/Real(NumPrPoints-1));
                }
        }

        return chisqr;

}

void PrFromScattering(int points, vector<Real> s, vector<Real> expi, vector<Real> error, ParamStruct params)
{
        //Finds P(r) which optimizes the fit to experimental data and is smooth.
        //Also can optimize Dmax.
        int NumPrPoints, NumFunctions, NumIterations;
        bool OptimizeDmax=true;
        char CharIntensityFile[1000], CharPrFile[1000];
        string IntensityFile;
        Real BestChisqr, BestDmax;
        Real chisqr, chisqr_old, dChisqr;
        Real dist, Dmax, inc;
        Real MaxMove, move, Temp, step;
        Real Rand;
        Real *PR;
        vector<Real> i;

        NumPrPoints=5000;
        NumFunctions=60;
        NumIterations=0;
        BestChisqr=10000000;

        SafeArrayAlloc(PR, NumPrPoints, "PR");
        SafeAlloc(i, points, "i");

        for (int n=0;n<NumPrPoints;n++)
        {
                PR[n]=0;
        }
        step=0.5;
        MaxMove=10.0;
        Dmax=40.0;
        Temp=10000.0;
        cout <<"About to enter FitPr5"<<endl;
        cout <<"NumFunctions= "<<NumFunctions<<" points= "<<points<<" NumPrPoints= "<<NumPrPoints<<endl;
        cout <<" Dmax= "<<Dmax<<" s[100]= "<<s[100]<<" i[100]= "<<i[100]<<" expi[100]= "<<expi[100]<<" PR[100]= "<<PR[100]<<endl;

        FitPr5(NumFunctions, points, NumPrPoints, Dmax, s, i, expi, PR, params);
        cout <<"Finished with FirPr5"<<endl;
        while (NumIterations<100)
        {
                move=MaxMove*(Real(rand()/Real(RAND_MAX)-0.5));
                Dmax+=move;

                chisqr=FitPr5(NumFunctions, points, NumPrPoints, Dmax, s, i, expi, PR, params);
                dChisqr=chisqr-chisqr_old;

                Rand=Real(rand())/Real(RAND_MAX);

                if (Rand < exp(-dChisqr/Temp))
                {
                        chisqr_old=chisqr;
                }
                else
                {
                        Dmax-=move;
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

        while (step*step>0.00001 && OptimizeDmax)
        {
                Dmax+=step;

                chisqr_old=chisqr;
                chisqr=FitPr5(NumFunctions, points, NumPrPoints, Dmax, s, i, expi, PR, params);
                cout <<"chisqr= "<<chisqr<<" chisqr_old= "<<chisqr_old<<" Dmax= "<<Dmax<<endl;

                if (chisqr>chisqr_old)
                {
                        Dmax-=step;
                        step=-step*0.5;
                        chisqr=chisqr_old;
                }
                else step=step*1.2;
        }

        IntensityFile=params.PrFile+"Intensity.txt";

        strcpy(CharPrFile, params.PrFile.c_str());
        strcpy(CharIntensityFile, IntensityFile.c_str());

        inc=Dmax/Real(NumPrPoints-1);

        dist=0;
        for (int n=0;n<NumPrPoints;n++)
        {
                ofstream pofr(CharPrFile, ios::app);
                pofr << dist<<"\t"<<PR[n]<<endl;
                dist+=inc;
        }

        for (int n=0;n<points;n++)
        {
                ofstream intensity(CharIntensityFile, ios::app);
                intensity <<s[n]<<"\t"<<i[n]<<"\t"<<expi[n]<<endl;
        }
}

void FitPolynomial(int degree, int FirstPoint, int LastPoint, Real InputX[], Real InputY[], Real output[])
{
        //Finds a polynomial which best fits a set of points.
        //This puts the y-axis intercept at zero.
        int i, p;
        int NumFitted;
        Real **Coefficient;
        Real Coefficientnn;
        Real InputXP, OutputP;
        Real scale;
        Real Sum;
        Real *variable;
        Real **xi;


        NumFitted=LastPoint;

        if (NumFitted>degree+1)
        {
                Safe2DArrayAlloc(Coefficient, degree+2, degree+2, "Coefficient");
                Safe2DArrayAlloc(xi, NumFitted, degree+2, "xi");

                for (int n=0;n<NumFitted;n++)
                {
                        for (int m=0;m<degree+1;m++)
                        {
                                xi[n][m]=0;
                        }
                }

                SafeArrayAlloc(variable, degree+1, "variable");
                for (int n=0;n<degree+1;n++)
                {
                        variable[n]=0;
                        for (int m=0;m<degree+1;m++)
                        {
                                Coefficient[m][n]=0;
                        }
                }

                for (int n=0;n<NumFitted;n++)
                {
                        xi[n][0]=1.0;
                        for (int m=1;m<degree+1;m++)
                        {
                                xi[n][m]=xi[n][m-1]*InputX[n];
                        }
                }

                for (int n=0;n<degree+1;n++)
                {
                        for (int m=0;m<degree+1;m++)
                        {
                                for (p=FirstPoint;p<LastPoint;p++)
                                {
                                        Coefficient[n][m]+=xi[p][n]*xi[p][m];
                                }
                        }

                        Sum=0;
                        for (p=FirstPoint;p<LastPoint;p++)
                        {
                                Sum+=xi[p][n]*InputY[p];
                        }
                        Coefficient[n][degree+1]=Sum;
                }

                for (int n=0;n<degree;n++)
                {
                        Coefficientnn=Coefficient[n][n];
                        if (Coefficientnn!=0)
                        {
                                for (int m=n+1;m<degree+1;m++)
                                {
                                        scale=Coefficient[m][n]/Coefficientnn;
                                        for (i=n;i<degree+2;i++)
                                        {
                                                Coefficient[m][i]-=Coefficient[n][i]*scale;
                                        }
                                }
                        }
                }

                for (int n=degree;n>-1;n=n-1)
                {
                        Coefficientnn=Coefficient[n][n];
                        if (Coefficientnn!=0)
                        {
                                variable[n]=Coefficient[n][degree+1]/Coefficient[n][n];
                                for (int m=n-1;m>-1;m=m-1)
                                {
                                        Coefficient[m][degree+1]=Coefficient[m][degree+1]-variable[n]*Coefficient[m][n];
                                }
                        }
                }

                for (p=0;p<NumFitted;p++)
                {
                        InputXP=InputX[p];
                        OutputP=0;
                        for (int m=degree;m>-1;m--)
                        {
                                OutputP=OutputP*InputXP;
                                OutputP+=variable[m];
                        }
                        output[p]=OutputP;
                }
        }
        else
        {
                cout <<"Error the degree of the polynomial exceeds the number of points."<<endl;
        }
}

void Smooth(int degree, int SegmentSize, int Step, int StartFrom, int EndAt, Real InputX[], Real InputY[], Real output[])
{       
        //Smooths a curve by fitting polynomials to different sections of it.
        int FirstPoint, LastPoint;
        Real *NumPoints;
        Real *Total;

        SafeArrayAlloc(NumPoints, EndAt, "NumPoints");
        SafeArrayAlloc(Total, EndAt, "Total");

        for (int n=0;n<EndAt;n++)
        {
                NumPoints[n]=0;
                Total[n]=0;
        }

        FirstPoint=StartFrom;
        LastPoint=FirstPoint+SegmentSize;

        while(LastPoint<EndAt)
        {
                FitPolynomial(degree, FirstPoint, LastPoint, InputX, InputY, output);
                for (int n=FirstPoint;n<LastPoint;n++)
                {
                        Total[n]+=output[n];
                        NumPoints[n]+=1.0;
                }
                FirstPoint+=Step;
                LastPoint+=Step;
        }

        for (int n=StartFrom;n<EndAt;n++)
        {
                if (NumPoints[n]!=0)
                {
                        output[n]=Total[n]/NumPoints[n];
                }
                else
                {
                        output[n]=InputY[n];
                }
        }

}

void ExtrapolateIntensity(vector<Real> s, vector<Real> expi, int points, ParamStruct params)
{
        //Extrapolates intensity to s=0;
        int point1, point2;
        Real expi1, expi2, expi3;
        Real i0;
        Real logi0, logi1, logi2;
        Real maxi;
        Real maxpoint, minpoint;
        Real rg2;


        for (int n=0;n<points;n++)
        {
                if (expi[n]>maxi)
                {
                        maxi=expi[n];
                        maxpoint=n;
                }
        }
        cout <<"maxi= "<<maxi<<" maxpoint= "<<maxpoint<<" expi[100]= "<<expi[100]<<endl;
        for (int n=int(maxpoint);n<points;n++)
        {
                expi3=expi2;
                expi2=expi1;
                expi1=expi[n];

                if ( (expi2 < expi1) && (expi2 < expi3) )
                {
                        minpoint=n-1;
                        break;
                }

                n++;
        }

        point1=int(floor(maxpoint+(minpoint-maxpoint)/3.0));
        point2=int(floor(minpoint-(minpoint-maxpoint)/3.0));

        if (point1==point2)
        {
                point1=27;
                point2=35;
        }
        cout <<"point1= "<<point1<<" point2= "<<point2<<endl;
        logi1=log(expi[point1]);
        logi2=log(expi[point2]);

        rg2=(logi1-logi2)/(s[point2]*s[point2]-s[point1]*s[point1]);
        logi0=logi1+s[point1]*s[point1]*rg2;
        i0=exp(logi0);

        for (int n=0;n<params.BeginFit;n++)
        {
                expi[n]=i0/exp(s[n]*s[n]*rg2);
                //cout <<"In extrapolatedi["<<n<<"]="<<extrapolatedi[n]<<endl;
        }
}

void fourierpr(vector<Real> s, vector<Real> expi, vector<Real> error, int points, char prfile[1000], ParamStruct params)
{
        //Finds p(r) by fourier transforming extrapolated experimental 
        //scattering data
        const long int MAX=1000000;
        int distance;
        int NumFitted;
        int Peak;
        int point1, point2;
        int SegmentSize, Step;
        bool CalcRinc=true;
        Real ds, dist;
        Real expi1, expi2, expi3;
        Real i0;
        Real logi0, logi1, logi2;
        Real maxi;
        Real maxpoint, minpoint;
        Real MaxPR;
        Real r, rinc;
        Real rg, rg2;
        Real rsphere;
        Real rsn;
        Real scale, scale2, scale4;
        Real sinc, Smax;
        Real sn, sn2, sn4;
        Real TotalPR2;
        Real PR[1000], DistanceArray[1000];
        bool BelowZero;

        //Real angle[MAX];
        //Real extrapolatedi[MAX];
        //Real is[MAX];

        Real *angle;
        vector<Real> extrapolatedi;
        Real *is;

        SafeArrayAlloc(angle, MAX, "angle");
        SafeAlloc(extrapolatedi, MAX, "extrapolatedi");
        SafeArrayAlloc(is, MAX, "is");

        for (int n=0;n<MAX;n++)
        {
                is[n]=0;
                angle[n]=0;
        }

        for (int n=0;n<params.EndFit;n++)
        {
                angle[n]=s[n];
        }

        for (int n=0;n<1000;n++)
        {
                PR[n]=0;
        }

        expi1=0;
        expi2=0;
        expi3=0;
        maxi=0;
        sinc=0.001;

        if (CalcRinc)
        {
                rinc=2.0*pi/s[params.EndFit];
        }
        else rinc=0.1;

        Smax=2.0*pi*10.0/rinc;

        for (int n=0;n<points;n++)
        {
                if (expi[n]>maxi)
                {
                        maxi=expi[n];
                        maxpoint=n;
                }
        }
        cout <<"maxi= "<<maxi<<" maxpoint= "<<maxpoint<<" expi[100]= "<<expi[100]<<endl;
        for (int n=int(maxpoint);n<points;n++)
        {
                expi3=expi2;
                expi2=expi1;
                expi1=expi[n];

                if ( (expi2 < expi1) && (expi2 < expi3) )
                {
                        minpoint=n-1;
                        break;
                }

                n++;
        }

        point1=int(floor(maxpoint+(minpoint-maxpoint)/3.0));
        point2=int(floor(minpoint-(minpoint-maxpoint)/3.0));

        if (point1==point2)
        {
                point1=27;
                point2=35;
        }
        cout <<"point1= "<<point1<<" point2= "<<point2<<endl;
        logi1=log(expi[point1]);
        logi2=log(expi[point2]);

        rg2=(logi1-logi2)/(s[point2]*s[point2]-s[point1]*s[point1]);
        logi0=logi1+s[point1]*s[point1]*rg2;
        i0=exp(logi0);

        if (rg2>0)
        {
                rg=sqrt(rg2*3.0);
                rsphere=rg*sqrt(5.0/3.0);
        }
        else
        {
                rg=0;
                cout <<"Radius of gyration error\n";
        }
        cout <<"rg= "<<rg<<endl;
        for (int n=0;n<maxpoint;n++)
        {
                extrapolatedi[n]=i0/exp(s[n]*s[n]*rg2);
                //cout <<"In extrapolatedi["<<n<<"]="<<extrapolatedi[n]<<endl;
        }

        for (int n=int(maxpoint);n<params.EndFit;n++)
        {
                extrapolatedi[n]=expi[n];
                //cout <<"extrapolatedi["<<n<<"]="<<extrapolatedi[n]<<endl;
        }

        scale=angle[params.EndFit-1];
        scale2=scale*scale;
        scale4=scale2*scale2*expi[params.EndFit-1];

        for (int n=params.EndFit; n<MAX ;n++)
        {
                angle[n]=angle[n-1]+0.001;
                sn=angle[n];
                sn2=sn*sn;
                sn4=sn2*sn2;
                //fsphere=3*(sin(s[n]*rsphere)-s[n]*rsphere*cos(s[n]*rsphere))/(s[n]*s[n]*s[n]*rsphere*rsphere*rsphere);
                //extrapolatedi[n]=i0*fsphere*fsphere;
                extrapolatedi[n]=scale4/sn4;
                //cout <<"extrapolatedi2["<<n<<"n= "<<extrapolatedi[n]<<endl;
                //extrapolatedi[n]=expi[points-1]*exp(s[points-1]*s[points-1])/exp(s[n]*s[n]);

        }
        for (int n=0;n<MAX;n++)
        {
                ofstream check("/home/jouko/WAXS/Extrapolatedi_initial.txt", ios::app);
                check << extrapolatedi[n] <<endl;
        }

        cout <<"extrapolatedi[1000]= "<<extrapolatedi[1000];
        dist=0;
        for (int n=0;n<1000;n++)
        {
                DistanceArray[n]=dist;
                dist+=0.1;
        }

        for (int m=0; m<1; m++)
        {
                for (int n=0; n<MAX; n++)
                {
                        is[n]=extrapolatedi[n]*angle[n];
                }

                r=0.0;
                for (distance=0;distance<1000;distance++)
                {
                        //PR[distance]=is[0]*sin(s[0]*r)*s[1]*0.5;
                        for (int n=1;n<MAX-2;n++)
                        {
                                if (angle[n]<Smax)
                                {
                                        ds=(angle[n+1]-angle[n-1])*0.5;
                                        PR[distance]+=is[n]*sin(angle[n]*r)*ds;
                                }
                                //cout <<"PR["<<distance<<"]="<<PR[distance]<<" is["<<n<<"]="<<is[n]<<" s["<<n<<"]="<<s[n]<<" r="<<r<<" ds="<<ds<<endl;
                        }
                        PR[distance]=r*PR[distance]/(2.0*pi*pi);
                        //cout <<"PR["<<distance<<"]= "<<PR[distance]<<" r= "<<r<<endl;
                        r+=rinc;
                }

                MaxPR=PR[0];
                for (int n=0;n<1000;n++)
                {
                        if (PR[n]>MaxPR)
                        {
                                MaxPR=PR[n];
                                Peak=n;
                        }
                }

                BelowZero=false;
                for (int n=Peak;n<1000;n++)
                {
                        if (PR[n]<0)
                        {
                                BelowZero=true;
                        }

                        if (PR[n]>0 && BelowZero)
                        {
                                NumFitted=n;
                        }
                }

                if (BelowZero==false)
                {
                        NumFitted=1000;
                }

                for (distance=NumFitted;distance<1000;distance++)
                {
                        //PR[distance]=0;
                }

                //Smooth(degree, SegmentSize, Step, StartFrom, EndAt, InputX[], InputY[], output[])
                SegmentSize=NumFitted/10;
                Step=NumFitted/100;
                //Smooth(2, SegmentSize, Step, 0, NumFitted, DistanceArray, PR, PR2);

                for (distance=0;distance<NumFitted;distance++)
                {
                        //PR[distance]=PR2[distance];
                }

                for (int n=0;n<MAX;n++)
                {
                        extrapolatedi[n]=0;
                        r=rinc;
                        for (distance=1;distance<1000;distance++)
                        {
                                rsn=r*angle[n];
                                extrapolatedi[n]+=PR[distance]*sin(rsn)/rsn;
                                r+=rinc;
                        }
                        extrapolatedi[n]=extrapolatedi[n]*4.0*pi*0.1;
                }

                normalize(extrapolatedi, expi, error, int(maxpoint), params.EndFit);

                for (int n=int(maxpoint);n<params.EndFit;n++)
                {
                        //extrapolatedi[n]=expi[n];
                }
        }

        MaxPR=PR[0];
        for (int n=0;n<1000;n++)
        {
                if (PR[n]>MaxPR)
                {
                        MaxPR=PR[n];
                }
        }

        if (MaxPR!=0)
        {
                for (int n=0;n<1000;n++)
                {
                        PR[n]=PR[n]/MaxPR;
                }
        }

        TotalPR2=0;
        for (int n=0;n<1000;n++)
        {
                TotalPR2+=PR[n]*PR[n]*0.1;
        }

        //cout <<"TotalElectrons="<<TotalElectrons<<endl;
        Real TotalElectrons; //Make sure TotalElectrons is assigned a value
        if (TotalPR2!=0 && TotalElectrons!=0)
        {
                for (int n=0;n<1000;n++)
                {
                        PR[n]=PR[n]*TotalElectrons*TotalElectrons/(TotalPR2*4.0*pi*pi);
                }
        }

        for (int n=0;n<MAX;n++)
        {
                ofstream check("/home/jouko/WAXS/Extrapolatedi.txt", ios::app);
                check << angle[n]<<" "<<extrapolatedi[n] <<endl;
        }
        r=0;
        for (int n=0;n<1000;n++)
        {
                ofstream checkpr(prfile, ios::app);
                checkpr <<r<<" "<< PR[n] <<endl;
                r+=rinc;
        }

        //delete[] extrapolatedi;
        //delete[] angle;
        //delete[] PR;
        //delete[] PR2;
}

void fourierpr2(vector<Real> s, vector<Real> i, int points, char prfile[1000], Real bin)
{
        //Finds p(r) by fourier transforming calculated intensity
        int distance;
        Real ds;
        Real r;
        Real PR[2000];
        Real *is;

        SafeArrayAlloc(is, points+10, "is");

        for (int n=0;n<2000;n++)
        {
                PR[n]=0;
        }

        for (int n=0; n<points; n++) is[n]=0;

        for (int n=0; n<points; n++) is[n]=i[n]*s[n];

        r=0;
        for (distance=0;distance<2000;distance++)
        {
                //PR[distance]=is[0]*sin(s[0]*r)*s[1]*0.5;
                for (int n=1;n<points-2;n++)
                {
                        ds=(s[n+1]-s[n-1])*0.5;
                        PR[distance]+=is[n]*sin(s[n]*r)*ds;
                        //cout <<"PR["<<distance<<"]="<<PR[distance]<<" is["<<n<<"]="<<is[n]<<" s["<<n<<"]="<<s[n]<<" r="<<r<<" ds="<<ds<<endl;
                }
                PR[distance]=r*PR[distance]/(2.0*pi*pi);
                r+=bin;
        }

        for (int n=0;n<2000;n++)
        {
                ofstream checkpr(prfile, ios::app);
                checkpr << PR[n] <<endl;
        }
}

void GetFitOptions(vector<bool> fit, ParamStruct params)
{
        //Interogates the used about which parameters to optimize.
        bool temp;
        int optimizeradii;

        cout <<"At what point do you want the fit to begin? ";
        cin >> params.BeginFit;
        cout <<"At what point do you want the fit to end? ";
        cin >> params.EndFit;

        cout <<"What parameters do you want to optimize?\n";
        cout <<"Optimize excluded volume scale factor: ";
        cin >> temp;
        fit[0]=temp;
        cout <<"Optimize hydration shell concentration: ";
        cin >> temp;
        fit[1]=temp;
        cout <<"Optimize solvent radii separetly: ";
        cin >> temp;
        fit[2]=temp;
        cout <<"Optimize atomic group excluded volumes individually: ";
        cin >> optimizeradii;

        if (optimizeradii==1)
        {
                cout <<"Optimize hydrogen excluded volume: ";
                cin >> temp;
                fit[3]=temp;
                cout <<"Optimize carbon excluded volume: ";
                cin >> temp;
                fit[4]=temp;
                cout <<"Optimize nitrogen excluded volume: ";
                cin >> temp;
                fit[5]=temp;
                cout <<"Optimize oxygen excluded volume: ";
                cin >> temp;
                fit[6]=temp;
                cout <<"Optimize sulfur excluded volume: ";
                cin >> temp;
                fit[7]=temp;
                cout <<"Optimize iron excluded volume: ";
                cin >> temp;
                fit[8]=temp;
        }
}

void solvent(Real atomr[], vector<AtomStruct> &Atoms, Real contrast, ParamStruct &params)
{
        solvent(atomr, Atoms, contrast, params.AtomTypesGofRFile, params.ElementsGofRFile, params.RecBin, params.maxhs, params.UniformHydrationShell, params.ExcludedVolumeSphereType, params.excluded, params);
}

Real fit(int natom, Real atomr[], IntensityStruct &i, Real contrast, Real &hsdensity, Real bin, int points, Real numatom[], vector<AtomStruct> &Atoms, ParamStruct params)
{
        //Fix this
        string experiment;
        string line;
        int TotalParticles;
        vector<bool> fit;
        Real mindistance[NumAtomTypes][NumAtomTypes];
        Real maxdistance[NumAtomTypes][NumAtomTypes];
        Real firstchisqr;
        Real chisqr1, chisqr2, chisqr3;
        Real total;
        Real volume;
        Real fitatomr[100];
        Real solventatomr[100];
        Real parameter[100];
        Real step[100];
        Real chisqr[100];

        for (int t=0;t<NumAtomTypes;t++)
        {
                fitatomr[t]=0;
                solventatomr[t]=0;
                for (int n=0;n<NumAtomTypes;n++)
                {
                        maxdistance[n][t]=0;
                        mindistance[n][t]=1000;
                }
        }

        for (int t=0;t<50;t++)
        {
                chisqr[t]=0;
                fit[t]=false;
                parameter[t]=0;
                step[t]=0;
        }

        for (int t=0;t<50;t++)
        {
                step[t]=0.0001;
        }


        ReadIntensityFile(params.ExperimentFile, i.s, i.expi, i.error);

        GetFitOptions(fit, params);

        if (params.VolumeOption=="FindRadii")
        {
                findradii(natom, numatom, atomr, Atoms);
                volume=findvolume(natom, atomr, numatom, Atoms);
        }

        for (int n=0;n<6;n++)
        {
                solventatomr[n]=atomr[n];
        }

        solvent(atomr, Atoms, contrast, params);
        TotalParticles=Atoms.size();
        findpr(bin, maxdistance, mindistance, Atoms, params.UniformHydrationShell);
        InitializeSinc(bin, maxdistance, i);

        for (int n=0;n<NumAtomTypes;n++) fitatomr[n]=atomr[n];

        parameter[0]=1.0;
        parameter[1]=hsdensity;
        parameter[2]=1.0;
        parameter[3]=1.0;

        for (int u=0;u<10;u++)
        {
                solvent(atomr, Atoms, contrast, params);
                findpr(bin, maxdistance, mindistance, Atoms, params.UniformHydrationShell);
                InitializeSinc(bin, maxdistance, i);

                ChooseIntensityCalculation(natom, atomr, solventatomr, i, contrast, hsdensity, bin, points, numatom, maxdistance, mindistance, Atoms, params);

                firstchisqr=normalize(i.calc, i.expi, i.error, params.BeginFit, params.EndFit);

                cout <<firstchisqr<<endl;

                for (int t=0;t<50;t++)
                {
                        step[t]=0.00001;
                }

                step[2]=0.01;

                for (int t=0;t<3;t++)
                {
                        if (fit[t]==1)
                        {
                                parameter[t]=parameter[t]+step[t];

                                for (int n=0;n<6;n++)
                                {
                                        atomr[n]=fitatomr[n]*parameter[0];
                                        solventatomr[n]=fitatomr[n]*parameter[2];
                                }

                                hsdensity=parameter[1];

                                if (t==2)
                                {
                                        solvent(atomr, Atoms, contrast, params);
                                        findpr(bin, maxdistance, mindistance, Atoms, params.UniformHydrationShell);
                                        InitializeSinc(bin, maxdistance, i);
                                }

                                ChooseIntensityCalculation(natom, atomr, solventatomr, i, contrast, hsdensity, bin, points, numatom, maxdistance, mindistance, Atoms, params);

                                chisqr[t]=normalize(i.calc, i.expi, i.error, params.BeginFit, params.EndFit);

                                parameter[t]=parameter[t]-step[t];

                                cout <<chisqr[t]<<"\t"<<t<<endl;
                        }
                }

                cout <<endl;

                total=0;

                for (int t=0;t<3;t++)
                {
                        if (fit[t]==1)
                        {
                                total += (firstchisqr-chisqr[t])*(firstchisqr-chisqr[t])/step[t];
                        }
                }

                for (int t=0;t<3;t++)
                {
                        if (fit[t]==1)
                        {
                                step[t]=(firstchisqr-chisqr[t])*0.00001/(sqrt(total*step[t]));
                        }
                }

                chisqr1=firstchisqr;
                chisqr2=0;
                chisqr3=0;

                while (true)
                {
                        for (int t=0;t<3;t++)
                        {
                                if (fit[t]==1)
                                {
                                        parameter[t] += step[t];
                                }
                        }

                        for (int n=0;n<6;n++)
                        {
                                atomr[n]=fitatomr[n]*parameter[0];
                        }

                        hsdensity=parameter[1];

                        //if (t==2)
                        //{
                        solvent(atomr, Atoms, contrast, params);
                        findpr(bin, maxdistance, mindistance, Atoms, params.UniformHydrationShell);
                        InitializeSinc(bin, maxdistance, i);
                        //}

                        Debye(atomr, i, contrast, hsdensity, bin, points, maxdistance, mindistance, Atoms, params);

                        chisqr2=chisqr1;
                        chisqr1=normalize(i.calc, i.expi, i.error, params.BeginFit, params.EndFit);

                        cout <<chisqr1<<"\t"<<parameter[0]<<"\t"<<parameter[1]<<"\t"<<parameter[2]<<endl;

                        if (chisqr2<chisqr1)
                        {
                                for (int t=0;t<3;t++)
                                {
                                        parameter[t]=parameter[t]-step[t];
                                        step[t]=-step[t]*0.5;
                                }
                        }

                        if (chisqr1<chisqr2)
                        {
                                for (int t=0;t<3;t++)
                                {
                                        step[t]=step[t]*1.2;
                                }
                        }

                        if ((chisqr2-chisqr1)*(chisqr2-chisqr1)<0.0000001)
                        {
                                break;
                        }
                }
        }

        return chisqr1;
}

Real OptimizeParameter(IntensityStruct &i, vector<AtomStruct> &Atoms, Real contrast, Real bin, Real maxdistance[][NumAtomTypes], Real mindistance[][NumAtomTypes], Real hsdensity, Real atomr[], vector<Real> &parameter, ParamStruct &params, int index, int natom)
{
        int NumParameters=parameter.size();
        Real chisqr, chisqr_old, rg, rg2, step, criterion=1e-10;
        Real BestParameter, BestChisqr;;
        Real numatom[NumAtomTypes];
        vector<Real> fitatomr, fitatomr2;

        chisqr_old=normalize(i.calc, i.expi, i.error, params.BeginFit, params.EndFit);
        BestChisqr=chisqr_old;
        SafeAlloc(fitatomr, NumParameters, "fitatomr");
        FindNumAtom(numatom, Atoms);

        step=0.001;
        BestParameter=parameter[index];
        while (true)
        {
                parameter[index]=parameter[index]+step;
                fitatomr[0]=parameter[3];
                fitatomr[1]=parameter[4];
                fitatomr[2]=parameter[5];
                fitatomr[3]=parameter[6];
                fitatomr[4]=parameter[7];

                for (int n=0;n<6;n++)
                {
                        atomr[n]=fitatomr[n]*parameter[0];
                        //solventatomr[n]=fitatomr2[n]*parameter[2];
                }

                hsdensity=parameter[1];

                if (index==2)
                {
                        solvent(atomr, Atoms, contrast, params);
                        if (params.histogram=="yes")
                        {
                                findpr(bin, maxdistance, mindistance, Atoms, params.UniformHydrationShell);
                                InitializeSinc(bin, maxdistance, i);
                        }
                }

                if (!params.DoTrajectory)
                {
                        ChooseIntensityCalculation(natom, atomr, solventatomr, i, contrast, hsdensity, bin, params.EndFit+1, numatom, maxdistance, mindistance, Atoms, params);
                }
                else
                {
                        for (int n=0;n<params.EndFit+1;n++) i.calc[n]=0;

                        trajectory(atomr, solventatomr, i, contrast, hsdensity, bin, params.EndFit+1, parameter[0], params);
                }

                chisqr=normalize(i.calc, i.expi, i.error, params.BeginFit, params.EndFit);
                if (!params.DoTrajectory)
                {
                        rg=radiusgyration(atomr, contrast, hsdensity, Atoms, params);
                        rg2=gyrationfromscattering(i.calc, i.s);
                }

                cout <<setw(9)<<left<<"Chisqr"<<setw(9)<<chisqr<<setw(9)<<left<<index<<setw(12)<<parameter[index]<<setw(9)<<"Rg"<<setw(9)<<rg<<setw(9)<<"Rg"<<setw(9)<<rg2<<endl;
                cout <<"chisqr= "<<chisqr<<" chisqr_old= "<<chisqr_old<<" dChisqr= "<<chisqr-chisqr_old<<" step= "<<step<<endl;
                if (chisqr_old>chisqr)
                {
                        step=step*1.2;
                        parameter[index]=parameter[index]+step;
                        chisqr_old=chisqr;
                        BestParameter=parameter[index];
                        BestChisqr=chisqr;
                }
                else
                {
                        parameter[index]-=step;
                        step=-step*0.5;
                }
                if (step*step<criterion) break;
        }
        parameter[index]=BestParameter;
        return BestChisqr;
}

Real ItterateOverParameters(vector<bool> fit, IntensityStruct &i, vector<AtomStruct> &Atoms, Real contrast, Real bin, Real maxdistance[][NumAtomTypes], Real mindistance[][NumAtomTypes], Real hsdensity, Real atomr[], vector<Real> &parameter, ParamStruct &params, int natom)
{
        int NumParameters=parameter.size();
        Real chisqr;
        for (int t=0;t<NumParameters;t++)
        {
                if (fit[t] && !(t==0 && params.excluded=="Lattice") ) 
                {
                        chisqr=OptimizeParameter(i, Atoms, contrast, bin, maxdistance, mindistance, hsdensity, atomr, parameter, params, t, natom);
                }
        }
        return chisqr;
}

Real fit2(int natom, Real atomr[], IntensityStruct &i, Real contrast, Real &hsdensity, Real bin, int points, Real numatom[], vector<AtomStruct> &Atoms, ParamStruct params)
{
        vector<bool> fit;
        int TotalParticles, NumParameters=8;
        Real mindistance[NumAtomTypes][NumAtomTypes];
        Real maxdistance[NumAtomTypes][NumAtomTypes];
        Real chisqr, chisqr_old, criterion=0.5;
        Real bestchisqr;
        Real fitatomr[100];
        Real fitatomr2[NumAtomTypes];
        vector<Real> parameter;

        SafeAlloc(fit, NumParameters, "fit");
        SafeAlloc(parameter, NumParameters, "parameter");

        for (int t=0;t<NumAtomTypes;t++)
        {
                for (int n=0;n<NumAtomTypes;n++)
                {
                        maxdistance[n][t]=0;
                        mindistance[n][t]=1000;
                }
        }


        ReadIntensityFile(params.ExperimentFile, i.s, i.expi, i.error);

        if (params.excluded=="Superimpose" || params.excluded=="Lattice")
        {
                if (params.VolumeOption=="FindRadii")
                {
                        FindRadii3(natom, atomr, numatom, Atoms);
                        //findradii(natom, numatom, atomr);
                        findvolume(natom, atomr, numatom, Atoms);
                }


                solvent(atomr, Atoms, contrast, params);
                //TotalParticles=CrysolSolvent(numatom, atomm, natom, contrast, atomr);
        }
        else
        {
                hsdensity=0;
                TotalParticles=natom;
        }
        //GetFitOptions(fit, params); //This should be uncommented;

        for (int j=0;j<NumParameters;j++)
        {
                fit[j]=false;
        }

        fit[HYDRATION_OPT]=true;
        fit[EXCLUDED_VOLUME_SCALE_OPT]=true;

        for (int n=0;n<NumAtomTypes;n++)
        {
                cout <<atomr[n]<<endl;
        }

        if (params.histogram=="yes")
        {
                findpr(bin, maxdistance, mindistance, Atoms, params.UniformHydrationShell);
                InitializeSinc(bin, maxdistance, i);
        }

        for (int n=0;n<6;n++)
        {
                fitatomr[n]=atomr[n];
                fitatomr2[n]=atomr[n];
        }

        parameter[0]=1.0;
        parameter[1]=hsdensity;
        parameter[2]=1.0;
        parameter[3]=atomr[0];
        parameter[4]=atomr[1];
        parameter[5]=atomr[2];
        parameter[6]=atomr[3];
        parameter[7]=atomr[4];

        chisqr_old=normalize(i.calc, i.expi, i.error, params.BeginFit, params.EndFit);
        bestchisqr=chisqr;
        while (true)
        {
                chisqr=ItterateOverParameters(fit, i, Atoms, contrast, bin, maxdistance, mindistance, hsdensity, atomr, parameter, params, natom);
                cout <<"chisqr_old= "<<chisqr_old<<" chisqr= "<<chisqr<<" dChisqr= "<<chisqr_old-chisqr<<" criterion= "<<criterion<<endl;
                if (chisqr_old-chisqr<criterion) break;
                chisqr_old=chisqr;
        }

        if (!params.DoTrajectory)
        {
                ChooseIntensityCalculation(natom, atomr, solventatomr, i, contrast, hsdensity, bin, points, numatom, maxdistance, mindistance, Atoms, params);
        }
        else
        {
                for (int n=0;n<params.EndFit+1;n++) i.calc[n]=0;

                trajectory(atomr, solventatomr, i, contrast, hsdensity, bin, points, parameter[0], params);
        }

        chisqr=normalize(i.calc, i.expi, i.error, params.BeginFit, params.EndFit);

        //printpr(atomr, bin, contrast, hsdensity, TotalPartilces);
        //fourierpr(s, i, points, "/home/jouko/WAXS/WAXS2/fouriercalc.txt");

        return chisqr;
}

Real fit3(int natom, Real atomr[], IntensityStruct &i, Real contrast, Real &hsdensity, Real bin, int points, Real numatom[], vector<AtomStruct> &Atoms, ParamStruct params)
{
        //This needs to be fixed
        string experiment;
        string line;
        Real mindistance[NumAtomTypes][NumAtomTypes];
        Real maxdistance[NumAtomTypes][NumAtomTypes];
        Real accept;
        Real size;
        Real temp;
        Real chisqr1;
        Real chisqr2;
        Real total;
        Real rg, rg2;
        Real fitatomr[100];
        Real fitatomr2[NumAtomTypes];
        Real solventatomr[100];
        Real parameter[100];
        Real parameter2[100];
        Real parameter3[100];
        Real fitparameter[100];
        Real step[100];
        Real chisqr[100];

        temp=10;

        for (int t=0;t<NumAtomTypes;t++)
        {
                fitatomr[t]=0;
                fitatomr2[t]=0;
                solventatomr[t]=0;
                for (int n=0;n<NumAtomTypes;n++)
                {
                        maxdistance[n][t]=0;
                        mindistance[n][t]=1000;
                }
        }

        for (int t=0;t<50;t++)
        {
                chisqr[t]=0;
                parameter[t]=0;
                parameter2[t]=0;
                parameter3[t]=0;
                fitparameter[t]=0;
                step[t]=0;
        }

        cout <<"At what point do you want the fit to begin? ";
        cin >> params.BeginFit;
        cout <<"At what point do you want the fit to end? ";
        cin >> params.EndFit;

        if (params.VolumeOption=="FindRadii")
        {
                findradii(natom, numatom, atomr, Atoms);
                findvolume(natom, atomr, numatom, Atoms);
        }

        if (params.VolumeOption=="FindRadii")
        {
                findvolume(natom, atomr, numatom, Atoms);
        }

        solvent(atomr, Atoms, contrast, params);
        if (params.histogram=="yes")
        {
                findpr(bin, maxdistance, mindistance, Atoms, params.UniformHydrationShell);
                InitializeSinc(bin, maxdistance, i);
        }

        for (int n=0;n<6;n++)
        {
                fitatomr[n]=atomr[n];
                fitatomr2[n]=atomr[n];
        }

        parameter[0]=1.0;
        parameter[1]=hsdensity;
        parameter[2]=1.0;
        parameter[3]=atomr[0];
        parameter[4]=atomr[1];
        parameter[5]=atomr[2];
        parameter[6]=atomr[3];
        parameter[7]=atomr[4];

        for (int u=0;u<100;u++)
        {
                for (int t=0;t<8;t++)
                {
                        step[t]=rand()-RAND_MAX*0.5;
                        total=total+step[t]*step[t];
                }

                total=sqrt(total);

                size=rand()/RAND_MAX;

                for (int t=0;t<8;t++)
                {
                        step[t]=step[t]*size/total;
                        parameter[t]=parameter[t]+step[t];
                }

                for (int n=0;n<6;n++)
                {
                        atomr[n]=fitatomr[n]*parameter[0];
                        solventatomr[n]=fitatomr2[n]*parameter[2];
                }

                hsdensity=parameter[1];

                solvent(atomr, Atoms, contrast, params);
                if (params.histogram=="yes")
                {
                        findpr(bin, maxdistance, mindistance, Atoms, params.UniformHydrationShell);
                        InitializeSinc(bin, maxdistance, i);
                }

                ChooseIntensityCalculation(natom, atomr, solventatomr, i, contrast, hsdensity, bin, points, numatom, maxdistance, mindistance, Atoms, params);

                chisqr2=chisqr1;
                chisqr1=normalize(i.calc, i.expi, i.error, params.BeginFit, params.EndFit);

                rg=radiusgyration(atomr, contrast, hsdensity, Atoms, params);
                rg2=gyrationfromscattering(i.calc, i.s);

                cout <<chisqr1<<"\t"<<rg<<"\t"<<rg2<<endl;

                chisqr1=chisqr2-chisqr1;

                if (accept*RAND_MAX<rand()*temp)
                {
                        chisqr1=chisqr2;
                        for (int t=0;t<8;t++)
                        {
                                parameter[t]=parameter[t]-step[t];
                        }
                }
        }

        return chisqr1;
}

void ChooseIntensityCalculation(int natom, Real atomr[], Real solventatomr[], IntensityStruct &i, Real contrast, Real hsdensity, Real bin, int points, Real numatom[], Real maxdistance[][NumAtomTypes], Real mindistance[][NumAtomTypes], vector<AtomStruct> &Atoms, ParamStruct params)
{

        cout <<"In choose intensity"<<endl;
        if (params.excluded=="Superimpose")
        {
                cout <<"excluded= "<<params.excluded<<endl;
        }
        if (params.excluded=="WaterSphere" || params.excluded=="Superimpose" || params.excluded=="Lattice")
        {
                if (params.histogram=="no")
                {
                        if (!params.UseCharges) nohistogram(atomr, i, contrast, hsdensity, bin, points, numatom, Atoms, params);
                        else charge(atomr, i, contrast, hsdensity, bin, points, numatom, Atoms, params);
                }
                else if (params.histogram=="yes") Debye(atomr, i, contrast, hsdensity, bin, points, maxdistance, mindistance, Atoms, params);
                else if (params.histogram=="multipole") multipole(atomr, i, contrast, hsdensity, points, Atoms, params);
                else
                {
                        cout <<"Error unrecognized value for histogram "<<params.histogram<<endl;
                        exit(EXIT_FAILURE);
                }
        }
        else if (params.excluded=="VectorAverage") VectorAverage(natom, atomr, contrast, hsdensity, points, i, params, Atoms);
        else if (params.excluded=="Cube") cube4d(natom, atomr, numatom, contrast, hsdensity, i, Atoms, params);
        //else if (excluded==2) cube4c(natom, TotalParticles, atomr, numatom, contrast, hsdensity, points, s, i);
        //else if (params.excluded=="WaterSphere") watersphere3(atomr, i, contrast, bin, points, numatom, Atoms, params);
        //else if (excluded==5) Image2(TotalParticles, atomr, s, i, contrast, hsdensity, bin, points, numatom, maxdistance, mindistance);
        //else if (excluded=="CubeFFT") CubeFFT3(natom, TotalParticles, atomr, numatom, contrast, hsdensity, points, s, i);
        else 
        {
                cout <<"Error "<<params.excluded<<" is not a valid option for calculating density"<<endl;
                exit(EXIT_FAILURE);
        }
        cout <<"Leaving ChooseIntensityCalculation"<<endl;
}

void ReadParameterFile(string ParameterFile, Real atomr[])
{
        //Complete later
        int NumLines;
        string par, parameter;
        vector<string> lines, str;

        ReadLines(ParameterFile, lines, "parameter file");
        ProcessLinesForReadingParameters(lines);
        NumLines=lines.size();

        for (int i=0;i<NumLines;i++)
        {
                Tokenize2(lines[i], " ", str);
                par=str[0];
                parameter=str[1];
                if (par=="RadiusHydrogen") atomr[HYDROGEN]=StrToFloat(parameter);
                else if (par=="RadiusCarbon") atomr[CARBON]=StrToFloat(parameter);
                else if (par=="RadiusNitrogen") atomr[NITROGEN]=StrToFloat(parameter);
                else if (par=="RadiusOxygen") atomr[OXYGEN]=StrToFloat(parameter);
                else if (par=="RadiusSulfur") atomr[SULFUR]=StrToFloat(parameter);
                else if (par=="RadiusPhosphorus") atomr[PHOSPHORUS]=StrToFloat(parameter);
                else if (par=="RadiusIron") atomr[IRON]=StrToFloat(parameter);
                else if (par=="RadiusSodium") atomr[NA_Plus]=StrToFloat(parameter);
                else if (par=="RadiusHydrationShell") atomr[HydrationShell]=StrToFloat(parameter);
        }
}

void SetDefaultParameters(ParamStruct &params)
{
        params.AngleBin=10;
        params.AngularCutOff=120.0;
        params.AngularDependence=false;
        params.AngularDependence2=false;
        params.ConvexConcave=false;
        params.ConvexConcave2=false;
        params.ElectrostaticPotential=false;
        params.ElectrostaticPotentialBinSize=0.5;
        params.MaxDistForZeroDensity=2.0;
        params.MaxElectrostaticPotential=9.0;
        params.MinElectrostaticPotential=-9.0;
        params.PhobicPhilic=false;
        params.AngleFile="";
        params.ChargeRadiiScale=0.0;
        params.AssignCubes="Location4";
        params.AtomTypesGofRFile="";
        params.BeginFit=30;
        params.bin=0.1;
        params.BoolFFT=true;
        params.BoolMerge=false;
        params.BoxSize=85.0;
        params.CalcNoHydrationIntensity=false;
        params.CalcPrFromExperiment=false;
        params.CalcPrFromStructure2=false;
        params.CalcPrFromStructure=false;
        params.CalcRg=false;
        params.CalcVacuumIntensity=false;
        params.CenterAtomsInBox=false;
        params.contrast=0.334269;
        params.RadiusHydrogen=0.90225;
        params.RadiusCarbon=1.43693;
        params.RadiusNitrogen=1.24527;
        params.RadiusOxygen=1.22099;
        params.RadiusSulfur=2.19596;
        params.RadiusPhosphorus=1.11;
        params.RadiusIron=0;
        params.RadiusSodium=0;
        params.RadiusHydrationShell=0.75;
        params.RoundGrid=true;
        params.CubeSize=0.5;
        params.DCDFiles="";
        params.DoFit=false;
        params.DoTrajectory="false";
        params.ElementsGofRFile="";
        params.EndFit=850;
        params.excluded="Superimpose";
        params.ExcludedVolumeSphereType="HardSphere";
        params.ExperimentFile="";
        params.FFTType="1D";
        params.GofRFile="";
        params.histogram="yes";
        params.hsdensity=0.03;
        params.HydrationRadiusOption="ByAtomType";
        params.HydrationRadiusScale=0.53;
        params.InputCubeDensityPdbFile="";
        params.IntensityOutputFile="";
        params.OutputCubeDensityPdbFile="";
        params.maxhs=3.0;
        params.MaxS=2.65;
        params.MoveProtein=false;
        params.nDCDFiles=1;
        params.NormalGofRFile="";
        params.NormalElementsGofRFile="";
        params.nstructures=1;
        params.NumPerEdge=3;
        params.PenaltyCoefficient=0;
        params.points=900;
        params.PrdfBulkDensity=0.334;
        params.PrFromStructureFile="";
        params.PrFromStructureFile2="";
        params.PrmFile="";
        params.PsfFile="";
        params.RecBin=0.1;
        params.ReadAngleFile=true;
        params.ReadUnknownResidues=true;
        params.ReadWaters=true;
        params.RemoveCavities="ByConnectivity";
        params.RtfFile="";
        params.ScatteringType="X-ray";
        params.SecondNearestNeighbor=false;
        params.skip=1;
        params.SolventOption="AtomTypes";
        params.SphereRadius=50.0;
        params.TotalElectronsFile="";
        params.UniformHydrationShell=false;
        params.UseCharges=false;
        params.UseChargeRadii=false;
        params.UseWaterSphere=false;
        params.VectorsPerInclination=20;
        params.VDWOffSet=0;
        params.VDWScale=1.0;
        params.VDWRadiusOption="ByAtomType";
        params.verbose=false;
        params.VolumeOption="Default";
        params.XBoxLength=85.0;
        params.XOrigin=0.0;
        params.xyz="";
        params.YBoxLength=85.0;
        params.YBoxLength=85.0;
        params.YOrigin=0.0;
        params.ZOrigin=0.0;
        params.ZBoxLength=0.0;
}

void PrintParams(ParamStruct &params)
{
        cout <<"Trajectory= "<<BoolToStr(params.DoTrajectory);
}

void PrintHelp()
{
        cout    
                <<"This program calculates a SWAXS pattern given an input "
                <<"protein structure.  The possible input parameters are listed "
                <<"bellow."<<endl<<endl

                <<"AngularDependence: Uses one of two sets of pRDFs depending "
                <<"upon whether the angle is less than or greater than a cut "
                <<"off angle"<<endl<<endl

                <<"AngularDependecee2: Uses a pRDF depending upon the angle "
                <<"formed by the cube, closest atom, and atom to which the "
                <<"closest atom is bonded"<<endl<<endl

                <<"BoolFFT: When the cube option is used this determines if "
                <<"FFT is used to calculate the scattering"<<endl<<endl

                <<"BoolMerge: Can be used to supperimpose many structures "
                <<"from a trajectory file on top of each other.  "
                <<"There is no reason to use this option"<<endl<<endl

                <<"CalcNoHydrationIntensity: In addition to the intensity "
                <<"from the given parameters this calculates the scattering "
                <<"from a protein with no hydration shell.  This option "
                <<"may not be compatible with all other options."<<endl<<endl

                <<"CalcPrFromExperiment: Given an inputed SWAXS pattern, "
                <<"calculates a p(r) using indirect Fourier transform.  "
                <<"This is a little buggy and not as good as GNOM."<<endl<<endl

                <<"ConvexConcave:  Uses one of two sets of pRDFs depending "
                <<"upon whether the cube is in a convex or concave region."

                <<"ConvexConcave2: Same as ConcaveConvex, but uses a different "
                <<"method to determine the concave and convex regions."
                <<endl<<endl

                <<"CalcPrFromStructure: Calculates the SWAXS pattern of a "
                <<"protein and calculates the p(r) from the SWAXS pattern."
                <<endl<<endl

                <<"CalcPrFromStructure2: Calculates the p(r) directly from "
                <<"the protein structure."<<endl<<endl

                <<"CalcRg: Calculates Rg."<<endl<<endl

                <<"CalcVacuumIntensity:  In addition to the intensity "
                <<"from the given parameters this calculated the scattering "
                <<"from a protein in vacuum."<<endl<<endl

                <<"CenterAtomsInBox:  Can be used with the cube method."
                <<"Does not effect the scattering but is useful when "
                <<"comparing protein hydration shells."<<endl<<endl

                <<"DoFit:  Adjusts parameters to fit to given experimental "
                <<"data."<<endl<<endl

                <<"DoTrajectory: Calculates a SWAXS pattern for every "
                <<"structure in a trajectory."<<endl<<endl

                <<"ElectrostaticPotential: Makes HyPred use the electrostatic "
                <<"potential when calculating the hydration shell density map."
                <<endl<<endl

                <<"PhobicPhilic: Makes HyPred use the hydrophobicity of the "
                <<"second nearest neighbor."<<endl<<endl

                <<"MoveProtein: Moves the system.  This does not effect the "
                <<"scattering but it can be usefull when comparing hydration "
                <<"shells.  This is just a boolean with acceptable values of "
                <<"yes and no.  To actually move the protein XOrigin... have "
                <<"to be set."<<endl<<endl

                <<"ReadAngelFile: Reads in the set of values of q for which "
                <<"the scattering intensity will be calculated."<<endl<<endl

                <<"ReadUnknownResidues:  Determines if residues other than the "
                <<"standard set of 20 residues, water, and HEM are used in the "
                <<"calculation."<<endl<<endl

                <<"ReadWaters: Determines if waters in the input pdb file are "
                <<"used.  This should be set to yes when the WaterSphere "
                <<"option is used and to no otherwise."<<endl<<endl

                <<"SecondNearestNeighbor: Determines if HyPred uses the "
                <<"identity of the second nearest solute atom to predict "
                <<"the hydration shell density map."<<endl<<endl

                <<"UniformHydrationShell: Determines if HyPred is used to "
                <<"predict the hydration shell density map or if the "
                <<"hydration shell density is set to a constant."

                <<"UseCharges: Determines if charges are used in determining "
                <<"scattering factors.  Untested."<<endl<<endl

                <<"UseWaterSphere: Determines if atoms outside of a certain "
                <<"radius should be elliminated.  Should be set to yes if "
                <<"calculating scattering from an all atom MD simulation."
                <<endl<<endl

                <<"verbose: Determines if debugging statements are printed out."
                <<endl<<endl<<endl

                <<"The following parameters should be string:"<<endl<<endl

                <<"AngleFile: The path to the file containing the values of q "
                <<"for which scattering is calculated."<<endl<<endl

                <<"AssignCubes: Determines how the densities of the cubes in "
                <<"the hydration shell are set.  Acceptable values are "
                <<"FromPdb Location3 and Location4."<<endl<<endl

                <<"AtomTypesGofRFile: The path to the file containing the "
                <<"detailed set of pRDFs used by HyPred."<<endl<<endl

                <<"DCDFiles: The path to file containing the paths to the DCD "
                <<"files to be used when DoTrajector is set to true."
                <<endl<<endl

                <<"ElementsGofRFile: The path to the file containing the "
                <<"coarse set of pRDFs used by HyPred."<<endl<<endl

                <<"ExperimentFile: The path to a file containing experimental "
                <<"data.  If DoFit is set to yes the calculated scattering "
                <<"will be fit to this data.  If not the chisqr of the fit "
                <<"will still be printed even if though no parameters are "
                <<"optimized."<<endl<<endl

                <<"excluded: The way in which the excluded volume is dealt with.  "
                <<"Acceptable values are Superimpose, WaterSphere, Lattice, "
                <<"VectorAverage, and Cube."<<endl
                <<"Superimpose: Places dummy excluded atoms on top of real atoms."
                <<endl
                <<"WaterSphere: This should be used with all atom MD simulations."
                <<endl
                <<"Lattice:  Places dummy excluded atoms on a lattice on top the "
                <<"protein."
                <<endl
                <<"VectorAverage: Instead of using the Debye formula, scattering "
                <<"is calculated by numerical averaging.  This option should not"
                <<"be in excluded."
                <<endl
                <<"Cube:  Places cubes on top of the protein and hydration shell."
                <<endl<<endl

                <<"ExcludedVolumeSphereType: This is the type of dummy atom to be "
                <<"used when excluded is set to superimpose.  The acceptable "
                <<"values are HardSphere and Gaussian."<<endl<<endl 

                <<"FFTType: This should be set 1D as 3D does not work."<<endl<<endl

                <<"GofRFile: Old parameter."<<endl<<endl

                <<"histogram: If set to yes the interatomic distances are binned.  "
                <<"If set to no the interatomic distances are not binned.  "
                <<"Setting this to yes speeds up the calculation by orders of "
                <<"magnitude.  If set to multipole, multipole expansion is used, "
                <<"instead of the Debye formula."<<endl<<endl

                <<"HydrationRadiusOption: Determines how the radii used by HyPred "
                <<"are set.  The two acceptable values are ByElement and "
                <<"ByAtomType."<<endl<<endl


                <<"InputCubeDensityPdbFile: The path to a  pdb file containing "
                <<"the hydration shell densities in the b-factor coloumn.  This "
                <<"needs to be set when AssingCubes is set to FromPdb."
                <<endl<<endl

                <<"IntensityOutputFile: The file to which the intensities are "
                <<"outputed to.  If none is specified the intensity file will "
                <<"be named after the pdb file."<<endl<<endl

                <<"NormalGofRFile: The path to the file containing the pRDFs "
                <<"with only distance and atom type information."<<endl<<endl

                <<"PdbList: The path to the file containing a list of pdbs "
                <<"for which SWAXS scattering patterns will be be calculated if "
                <<"DoTrajector is set to yes."

                <<"PrFile: The file to which the p(r) calculated from SWAXS data "
                <<"is outputted."<<endl<<endl

                <<"PrFromStructureFile: The file to which the p(r) calculated "
                <<"from the Fourier transform of a SWAXS scattering pattern."
                <<endl<<endl

                <<"PrFromStructureFile2: The file to which the p(r) calculated "
                <<"from the protein structure is outputted."<<endl<<endl

                <<"PrmFile: A Charmm .prm file.  This contains the vdw radii "
                <<"used to determine the hydration shell."<<endl<<endl

                <<"PsfFile: A Charmm .psf file.  This is used by HyPred to get "
                <<"the atom connectivity needed for angular dependence."
                <<endl<<endl

                <<"RemoveCavities: The method by which hydration shell density in "
                <<"the protein is elliminated.  The two acceptable values are "
                <<"ByConnectivity and ByCutOff.  ByConnectivity is reccomended."
                <<endl<<endl

                <<"RtfFile: A Charmm .rtf file."<<endl<<endl

                <<"TotalElectronsFile: When the WaterSphere option is used "
                <<"this is the file to which the number of solvent electrons "
                <<"within the sphere is outputted."<<endl<<endl

                <<"ScatteringType: The type of scattering to be calculated. "
                <<"The acceptable values are X-ray and Neutron."<<endl<<endl

                <<"SolventOption: Does this do anything?  I will have to check."
                <<endl<<endl

                <<"VolumeOption: Determines what set of radii are used for the "
                <<"excluded volumes of the dummy atoms.  The acceptable values are "
                <<"Default, FraserSuzukiMacRea, FindRadii, and UserDefined."
                <<endl<<endl

                <<"xyz: The cartessian coordinate file.  Usually left blank.  "
                <<"Instead xyz is the second parameter passed to WAXS."
                <<endl<<endl<<endl

                <<"The following parameters should be integers."<<endl<<endl

                <<"BeginFit: When a fit is performed this specifies the number of "
                <<"the point at which the fit begins."<<endl<<endl

                <<"EndFit: When a fit is performed this specifies the number of "
                <<"the point at which the fit ends."<<endl<<endl

                <<"lmax: When excluded is set to multipole this determines the "
                <<"number of spherical harmonics used."<<endl<<endl

                <<"nDCDFiles: When DoTrajectory is used this determines the "
                <<"number of .dcd files to be read in."<<endl<<endl

                <<"NumPerEdge: When all atom explicit water MD simulations are "
                <<"used this can be used to create images of the original "
                <<"simulation."<<endl<<endl

                <<"nstructures: When DoTrajectory is used this determines the "
                <<"number of structures to be read in."<<endl<<endl

                <<"points: This is the number of values of q for which intensities "
                <<"will be calculated."<<endl<<endl

                <<"skip: When DoTrajectory is set to true.  Every skipth structure "
                <<"will be used."<<endl<<endl

                <<"VectorsPerInclination: Determines the number of vectors used "
                <<"when excluded is set to Cube or VectorAverage is set to yes.  "
                <<"The actual number of vectors is 0.5*VectorsPerInclination^2."
                <<endl<<endl

                <<"AngleBin: When AngularDependence2 is set to yes this determines "
                <<"the angular bin size used by HyPred."<<endl<<endl

                <<"AngularCutOff: When AngularDependence is set to yes this sets "
                <<"the angular cut off used by HyPred."<<endl<<endl

                <<"bin: The bin size when histogram is set to yes."<<endl<<endl

                <<"BoxSize: This has been replaced by XBoxLength..."<<endl<<endl

                <<"contrast: The density of bulk solvent."<<endl<<endl

                <<"CubeSize: The size of the cubes used when excluded is set to "
                <<"Cubes."<<endl<<endl

                <<"ElectrostaticPotentialBinSize: When HyPred uses the "
                <<"electrostatic potential it uses this bin size"<<endl<<endl

                <<"hsdensity: The hydration shell density."<<endl<<endl

                <<"HydrationSpacing: Not sure what this is."<<endl<<endl

                <<"HydrationScale: Scales the vdw radii used by HyPred."<<endl<<endl

                <<"MaxElectrostaticPotential: Electrostatic potentials larger "
                <<"than this are considered to be the same."<<endl<<endl

                <<"maxhs: The thickness of the hydration shell layer."<<endl<<endl

                <<"MaxS: The maximum value of q for which I(q) is calculated"
                <<endl<<endl

                <<"MinElectrostaticPotential: Electrostatic potentials smaller "
                <<"than this are considered to be the same."<<endl<<endl

                <<"PenaltyCoefficient: When p(r) is calculated by the indirect "
                <<"Fourier transform this penalizes high second derivatives."
                <<endl<<endl

                <<"PrdfBulkDensity: This is the density subtracted from the "
                <<"hydration shell when HyPred is used.  Ideally it should be "
                <<"the same as contrast."<<endl<<endl

                <<"RecBin: The bin size used by HyPred."<<endl<<endl

                <<"RadiusHydrogen...: The excluded volume radii."<<endl<<endl

                <<"SphereRadius: This is the size of the sphere when the option "
                <<"UseWaterSphere is set to true."<<endl<<endl

                <<"VDWOffSet: When setting VDW radii this adds a constant to all "
                <<"the VDW radii."<<endl<<endl

                <<"VDWScale: When setting VDW radii this scales the VDW radii."
                <<endl<<endl 

                <<"XBoxLength: The length of the x dimension of the simulation "
                <<"box.  Used when calculating WAXS from an MD simulation. "
                <<"Can also be used to set the displacement of the grid when "
                <<"using the cube method."<<endl<<endl

                <<"YBoxLength: Same as XBoxLength."<<endl<<endl

                <<"ZBoxLength: Same as XBoxLength."<<endl<<endl

                <<"XOrigin: Moves protein by XOrigin angstroms."<<endl<<endl

                <<"YOrigin: Same as XOrigin."<<endl<<endl

                <<"ZOrigin: Same as XOrigin."<<endl<<endl;
}

void GetAtoms(vector<AtomStruct> &Atoms, Real numatom[], ParamStruct params)
{
        string str;
        int natom, position;
        if (!params.BoolMerge)
        {
                position=params.xyz.rfind(".");
                str=params.xyz.substr(position+1,4);
                if (str=="pdb") natom=ReadPdb(params.xyz, Atoms);
                else ReadXYZ(params.xyz, Atoms);
                if (!params.ReadWaters) RemoveWaters(Atoms);
                if (!params.ReadUnknownResidues) RemoveUnknownResidues(Atoms);
                AssignAtomIDs(Atoms);
                GetAtomType(Atoms);
                GetAtomType2(Atoms);
        }
        else Merge(Atoms, params);
}

void SetDefaultAtomicRadii(Real atomr[], ParamStruct params)
{
        cout <<"In SetDefaultAtomicRadii"<<endl;
        for (int n=0;n<NumAtomTypes;n++) atomr[n]=0;
        atomr[HYDROGEN]=0.902255;
        atomr[CARBON]=1.43693;
        atomr[NITROGEN]=1.24527;
        atomr[OXYGEN]=1.22099;
        atomr[SULFUR]=2.19596;
        atomr[PHOSPHORUS]=1.11;
        cout <<"ExcludedVolumeSphereType= "<<params.ExcludedVolumeSphereType<<endl;
        if (params.ExcludedVolumeSphereType=="HardSphere")
        {
                for (int n=0;n<NumAtomTypes;n++) atomr[n]*=GaussianToHardSphere;
        }
        atomr[ExcludedVolume]=0.75*packing;
        atomr[HydrationShell]=0.75*packing;
        cout <<"atomr[HydrationShell]= "<<atomr[HydrationShell]<<endl;
        if (params.ExcludedVolumeSphereType=="GaussianSphere")
        {
                atomr[ExcludedVolume]/=GaussianToHardSphere;
                atomr[HydrationShell]/=GaussianToHardSphere;
        }
        cout <<"atomr[HydrationShell]= "<<atomr[HydrationShell]<<endl;
}

void SetFraserSuzukiMacReaAtomicRadii(Real atomr[])
{
        atomr[HYDROGEN]=1.07;
        atomr[CARBON]=1.58;
        atomr[NITROGEN]=0.84;
        atomr[OXYGEN]=1.30;
        atomr[SULFUR]=1.68;
        atomr[PHOSPHORUS]=1.11;
        atomr[ExcludedVolume]=0.75;
        atomr[HydrationShell]=0.75;
}

void GetAngles(IntensityStruct &i, ParamStruct params)
{
        string line;
        if (params.ReadAngleFile)
        {

                fstream theta;
                OpenFile(params.AngleFile, theta, "angle file");

                for (int t=0;t<params.points;t++)
                {
                        getline(theta,line);
                        i.s[t]=StrToFloat(line);
                }
                theta.close();
        }
        else
        {
                for (int t=0;t<params.points;t++)
                {
                        i.s[t]=Real(t+1)*params.MaxS/Real(params.points);
                }
        }
}

void InitializeIntensity(IntensityStruct &i, int points)
{
        vector<Real> array;
        vector< vector<Real> > matrix;
        cout <<"In InitializeIntensity"<<endl;
        SafeAlloc(i.calc, points, "i.calc");
        SafeAlloc(i.CalcError, points, "i.CalcError");
        SafeAlloc(i.error, points, "i.error");
        SafeAlloc(i.expi, points, "i.expi");
        SafeAlloc(i.s, points, "i.s");
        SafeAlloc(i.vacuum, points, "i.vacuum");
        SafeAlloc(i.NoHydration, points, "i.NoHydration");
        SafeAlloc(array, NumAtomTypes, "array");
        CreateMatrix(matrix, NumAtomTypes, NumAtomTypes);
        cout <<"Created Matrix"<<endl;
        for (int n=0;n<points;n++) 
        {
                SafePushBack(i.f, array, "i.f");
                SafePushBack(i.point, matrix, "i.point");
                SafePushBack(i.CrossTerm, matrix, "i.CrossTerm");
        }
        cout <<"Allocated others"<<endl;
}

void interogate(ParamStruct params, Real atomr[])
{
        string StrUniformHydrationShell;
        cout <<"In if (!CommandLine)"<<endl;
        cout << "Enter the xyz file path: ";
        getline(cin, params.xyz);
        cout <<"Enter average solvent electron concentration: ";
        cin >> params.contrast;
        cout <<"Use uniform hydration shell density?";
        cin.ignore();
        getline(cin, StrUniformHydrationShell);
        params.UniformHydrationShell=StrToBool(StrUniformHydrationShell);
        cout <<"UniformHydrationShell= "<<params.UniformHydrationShell<<endl;

        if (params.UniformHydrationShell)
        {
                cout <<"Enter hydration shell contrast: ";
                cin >> params.hsdensity;
        }
        else
        {
                cin.ignore();
                cout <<"Enter file containinig g(r): ";
                getline(cin, params.GofRFile);
        }

        cout <<"0: Superimpose dummy solvent atoms on top of real atoms."<<endl;
        cout <<"1: Place dummy solvent atoms on grid."<<endl;
        cout <<"2: Use cube method."<<endl;
        cout <<"3: Use a water sphere."<<endl;
        cout <<"4: Use a water box."<<endl;
        cout <<"5: Use an image."<<endl;
        cout <<"Enter option: ";
        cin >> params.excluded;

        if (params.excluded=="Superimpose" || params.excluded=="Lattice")
        {
                cout <<"0: Use Debye formula without histogram.\n";
                cout <<"1: Use Debye formula with histogram.\n";
                cout <<"2: Use mulipole expansion.\n";
                cout <<"Enter option: ";
                cin >> params.histogram;
        }
        else
        {
                params.histogram="multipole";
        }

        if (params.histogram=="yes")
        {
                cout <<"Enter your desired resolution: ";
                cin >> params.bin;
        }

        if (params.excluded=="Superimpose")
        {
                cout <<"0: Calculate excluded volumes of atoms.\n";
                cout <<"1: Use default values.\n";
                cout <<"2: Enter your own values.\n";
                cout <<"3: Use Fraser-MacRae-Suzuki values.\n";
                cout <<"Enter option: ";
                cin >> params.VolumeOption;

                if (params.VolumeOption=="UserDefined")
                {
                        cout <<"Enter hydrogen radius: ";
                        cin >> atomr[HYDROGEN];
                        cout <<"Enter carbon radius: ";
                        cin >> atomr[CARBON];
                        cout <<"Enter nitrogen radius: ";
                        cin >> atomr[NITROGEN];
                        cout <<"Enter oxygen radius: ";
                        cin >> atomr[OXYGEN];
                        cout <<"Enter sulfur radius: ";
                        cin >> atomr[SULFUR];
                        cout <<"Enter iron radius: ";
                        cin >>atomr[IRON];
                        cout <<"Enter phosphorus radius: ";
                        cin >>atomr[PHOSPHORUS];
                }

                cout <<"0: Use gaussian density spheres."<<endl;
                cout <<"1: Use hard spheres."<<endl;
                cout <<"2: Use trapazoidal density spheres."<<endl;
                cin >> params.ExcludedVolumeSphereType;
        }
        if (params.histogram=="no")
        {
                cout <<"Account for partial charges: ";
                cin >>params.UseCharges;
        }
        cout <<"Enter the number of points to be calculated: ";
        cin >> params.points;

        cout <<"Enter file path containing angles: ";
        cin.ignore();
        getline(cin, params.AngleFile);
        cout <<"0: Calculate scattering from one structure.\n";
        cout <<"1: Calculate scattering from trajectory.\n";
        cout <<"Enter option: ";
        cin >> params.DoTrajectory;

        cout <<"0: Don't optimize scattering to fit to experimental data.\n";
        cout <<"1: Optimize scattering to fit to experimental data.\n";
        cout <<"Enter option: ";
        cin >> params.DoFit;

        if (params.DoFit)
        {
                cin.ignore();
                cout <<"Enter experimental data file path: ";
                getline(cin, params.ExperimentFile);
        }

        if (params.DoTrajectory)
        {
                cout <<"Enter the number of structures in the trajectory: ";
                cin >> params.nstructures;
                cout <<"Use every nth structure: ";
                cin >> params.skip;
        }
} 

void SetHydrationRadius(vector<AtomStruct> &Atoms, Real solventatomr[])
{
        int natom=Atoms.size();
        for (int n=0;n<natom;n++) Atoms[n].HydrationRadius=solventatomr[Atoms[n].atomid];
}

void SetHydrationRadius(vector<AtomStruct> &Atoms)
{
        int natom=Atoms.size();
        for (int n=0;n<natom;n++) Atoms[n].HydrationRadius=Atoms[n].vdw;
}

void SetDefaultVDW(vector<AtomStruct> &Atoms, Real solventatomr[])
{
        cout <<"In SetDefaultVDW"<<endl;
        int natom=Atoms.size();
        for (int n=0;n<natom;n++) Atoms[n].vdw=solventatomr[Atoms[n].atomid];
}

void SetDefaultSolventAtomR(Real solventatomr[])
{
        solventatomr[HYDROGEN]=0.504367;
        solventatomr[CARBON]=1.36141;
        solventatomr[NITROGEN]=1.38984;
        solventatomr[OXYGEN]=0.887856;
        solventatomr[SULFUR]=1.22271;
}

void ScaleHydrationRadius(Real HydrationRadiusScale, vector<AtomStruct> &Atoms)
{
        int natom=Atoms.size();
        for (int j=0;j<natom;j++) Atoms[j].HydrationRadius*=HydrationRadiusScale;
}

void SetVDWRadius(vector<AtomStruct> &Atoms, ParamStruct params, Real solventatomr[])
{
        cout <<"In SetVDWRadius"<<endl;
        if (params.PrmFile!="" && params.RtfFile!="") 
        {
                SetVDW(params.PrmFile, params.RtfFile, Atoms);
        }
        else SetDefaultVDW(Atoms, solventatomr);
        if (params.HydrationRadiusOption=="ByElement") SetHydrationRadius(Atoms, solventatomr);
        else if (params.HydrationRadiusOption=="ByAtomType") SetHydrationRadius(Atoms);
        else 
        {
                cout <<"ERROR: Unrecognized HydrationRadiusOption "<<params.HydrationRadiusOption<<endl;
                exit(EXIT_FAILURE);
        }
        if (params.VDWRadiusOption=="ByElement") SetDefaultVDW(Atoms, solventatomr);
        else if (params.HydrationRadiusOption=="ByAtomType") SetVDW(params.PrmFile, params.RtfFile, Atoms);
        else 
        {
                cout <<"ERROR: Unrecognized VDWRadiusOption "<<params.VDWRadiusOption<<endl;
                exit(EXIT_FAILURE);
        }
        ScaleVDW(params.VDWScale, params.VDWOffSet, Atoms);
        ScaleHydrationRadius(params.HydrationRadiusScale, Atoms);
        if (params.UseChargeRadii) AddChargeRadii(Atoms, params.ChargeRadiiScale);
}

void preprocess(ParamStruct params, int natom, Real numatom[], Real atomr[], vector<AtomStruct> &Atoms, Real maxdistance[][NumAtomTypes], Real mindistance[][NumAtomTypes], IntensityStruct &i)
{
        cout <<"In preprocess"<<endl;
        if (params.excluded=="Superimpose" || params.excluded=="Lattice")
        {
                cout <<"volumeoption= "<<params.VolumeOption<<endl;
                if (params.VolumeOption=="FindRadii")
                {
                        //findradii(natom, numatom, atomr);
                        FindRadii2(natom, numatom, atomr, Atoms);
                        findvolume(natom, atomr, numatom, Atoms);
                }

                if (params.hsdensity!=0 || params.excluded=="Lattice")
                {
                        cout <<"About to add solvent"<<endl;
                        solvent(atomr, Atoms, params.contrast, params);
                        //PrintPdb("/home/jouko/WAXS/pdb/WAXS145Hydration.pdb", Atoms);
                        cout <<"Left solvent"<<endl;
                        //CrysolSolvent(numatom, atomm, natom, contrast, atomr);
                }
        }
        else params.hsdensity=0;
        if (params.UseWaterSphere) MakeWaterSphere(atomr, numatom, Atoms, params);
        if (params.histogram=="yes" && (params.UseWaterSphere || params.excluded=="Superimpose" || params.excluded=="Lattice"))
        {
                cout <<"About to enter findpr"<<endl;
                findpr(params.bin, maxdistance, mindistance, Atoms, params.UniformHydrationShell);
                cout <<"Found pr"<<endl;
                InitializeSinc(params.bin, maxdistance, i);
                cout <<"Initialized sinc"<<endl;
        }
}

void waxs(ParamStruct params, int natom, Real atomr[], IntensityStruct &i, IntensityStruct &IntensityVacuum, IntensityStruct &IntensityNoHydration, vector<AtomStruct> &Atoms, Real numatom[])
{
        Real chisqr, scale;
        Real maxdistance[NumAtomTypes][NumAtomTypes], mindistance[NumAtomTypes][NumAtomTypes];
        cout <<"In waxs"<<endl;
        //PrintPdb("/home/jouko/project/WAXS/test/log/InWaxs.pdb", Atoms);
        for (int m=0;m<NumAtomTypes;m++)
        {
                for (int n=0;n<NumAtomTypes;n++)
                {
                        mindistance[n][m]=1000.0;
                        maxdistance[n][m]=0;
                }
        }
        if (params.DoFit)
        {
                //chisqr=fit(natom, atomr, s, i, contrast, hsdensity, bin, points, numatom);
                chisqr=fit2(natom, atomr, i, params.contrast, params.hsdensity, params.bin, params.points, numatom, Atoms, params);
        }
        else
        {
                if (!params.DoTrajectory)
                {
                        cout <<"preparing"<<endl;
                        preprocess(params, natom, numatom, atomr, Atoms, maxdistance, mindistance, i);
                        cout <<"About to enter ChooseIntensityCalculation"<<endl;
                        ChooseIntensityCalculation(natom, atomr, solventatomr, i, params.contrast, params.hsdensity, params.bin, params.points, numatom, maxdistance, mindistance, Atoms, params);
                }

                if (params.DoTrajectory)
                {
                        trajectory(atomr, solventatomr, i, params.contrast, params.hsdensity, params.bin, params.points, scale, params);
                }
        }

        if (params.CalcVacuumIntensity) 
        {
                bool temp=params.UniformHydrationShell;
                params.UniformHydrationShell=true;
                ChooseIntensityCalculation(natom, atomr, solventatomr, IntensityVacuum, 0, 0, params.bin, params.points, numatom, maxdistance, mindistance, Atoms, params);
                params.UniformHydrationShell=temp;
        }

        if (params.CalcNoHydrationIntensity) 
        {
                bool temp=params.UniformHydrationShell;
                params.UniformHydrationShell=true;
                ChooseIntensityCalculation(natom, atomr, solventatomr, IntensityNoHydration, params.contrast, 0, params.bin, params.points, numatom, maxdistance, mindistance, Atoms, params);
                params.UniformHydrationShell=temp;
        }
        //if (params.excluded=="Cube") Hydrate(params.contrast, params.hsdensity, Atoms, atomr, params);
}

void SetIntensityOutputFileName(string xyz, string &intensityfile, string &logfile)
{
        bool valid;
        char charintensity[1000], charlogfile[1000];
        string str, StrFileIndex;
        int FileIndex;
        size_t position;

        position=xyz.find(".");
        str=xyz.substr(0,position);
        intensityfile=str+"Scattering.txt";
        logfile=str+"Log.txt";
        strcpy(charintensity, intensityfile.c_str());
        strcpy(charlogfile, logfile.c_str());

        valid=exist(charlogfile);
        FileIndex=2;

        while (valid)
        {
                stringstream out;
                out << FileIndex;
                StrFileIndex=out.str();
                logfile=str+"Log"+StrFileIndex+".txt";
                strcpy(charlogfile, logfile.c_str());
                intensityfile=str+"Scattering"+StrFileIndex+".txt";
                strcpy(charintensity, intensityfile.c_str());
                valid=exist(charlogfile);
                FileIndex++;
        }
}

void HydrationEffects(IntensityStruct &i, IntensityStruct &IntensityVacuum, IntensityStruct &IntensityNoHydration, ParamStruct params)
{
        int points=i.calc.size();

        for (int t=0;t<points;t++)
        {
                if (params.CalcNoHydrationIntensity) 
                {
                        i.NoHydration[t]=IntensityNoHydration.calc[t];
                }
                if (params.CalcVacuumIntensity)
                {
                        i.vacuum[t]=IntensityVacuum.calc[t];
                }
        }
}

void ParamsToAtomicRadii(ParamStruct &params, Real atomr[])
{
        atomr[HYDROGEN]=params.RadiusHydrogen;
        atomr[CARBON]=params.RadiusCarbon;
        atomr[NITROGEN]=params.RadiusNitrogen;
        atomr[OXYGEN]=params.RadiusOxygen;
        atomr[SULFUR]=params.RadiusSulfur;
        atomr[IRON]=params.RadiusIron;
        atomr[PHOSPHORUS]=params.RadiusPhosphorus;
        atomr[NA_Plus]=params.RadiusSodium;
        atomr[HydrationShell]=params.RadiusHydrationShell*packing;
        if (params.ExcludedVolumeSphereType=="GaussianSphere") atomr[HydrationShell]*=HardSphereToGaussianSphere;
}

int main (int argc, char *argv[])
{
        vector<AtomStruct> Atoms;
        ParamStruct params;
        bool BoolUseOptionFile, valid;
        char charintensity[1000], charlogfile[1000], paramfile[1000];
        string intensityfile;
        string str, StrFileIndex;
        string param, ParameterFile, PrFromStructureFile2;
        string UseOptionFile;
        string line, logfile;
        int TotalParticles, natom;
        Real chisqr;
        Real VacuumRg, ExcludedRg, Rg;
        Real rfactor, sum;
        Real numatom[NumAtomTypes], atomr[NumAtomTypes], solventatomr[NumAtomTypes];
        Real NumElectrons[NumAtomTypes];
        IntensityStruct i, IntensityVacuum, IntensityNoHydration;
        time_t seconds1, seconds2;

        natom=0;
        TotalParticles=0;
        cout <<"Starting"<<endl;

        for (int n=0;n<NumAtomTypes;n++)
        {
                numatom[n]=0;
                atomr[n]=0;
                solventatomr[n]=0;
        }

        SetDefaultParameters(params);
        SetDefaultAtomicRadii(atomr, params);
        cout <<"After SetDefaultAtomicRadii atomr[HydrationShell]= "<<atomr[HydrationShell]<<endl;
        SetDefaultSolventAtomR(solventatomr);
        SetNumElectrons(NumElectrons);
        if (argc==1)
        {
                CommandLine=false;
                cout <<"Use option file: ";
                getline(cin, UseOptionFile);
                BoolUseOptionFile=StrToBool(UseOptionFile);
        }


        if (argc>1 || BoolUseOptionFile)
        {
                CommandLine=true;
                if (argc>1)
                {
                        params.xyz=argv[1];
                        ParameterFile=argv[2];
                }
                else
                {
                        cout <<"Enter option file: ";
                        getline(cin, ParameterFile);
                        cin.ignore();
                }
                ReadParameterFile(ParameterFile, params);
                cout <<"About to enter ReadParameterFile(ParameterFile, atomr);"<<endl;
                ReadParameterFile(ParameterFile, atomr);
        }
        cout <<"UniformHydrationShell= "<<params.UniformHydrationShell<<endl;
        cout <<"atomr[HydrationShell]= "<<atomr[HydrationShell]<<endl;
        if (argc>1)
        {
                cout <<params.xyz<<endl;
        }
        cout <<"Before if (!CommandLine)"<<endl;
        cout <<"CommandLine= "<<CommandLine<<endl;
        if (!CommandLine) interogate(params, atomr);
        verbose=params.verbose;
        if (params.UniformHydrationShell) params.maxhs=3.0;
        else params.hsdensity=1.0;


        valid=exist(params.xyz);
        //if (!valid) FileNotValid(params.xyz, "cartesian coordinate file");
        if (!valid)
        {
                cout <<"ERROR: Unable to open cartesian coordinate file"<<endl;
                exit(EXIT_FAILURE);
        }

        cout <<"checked for cartesian coordinate file"<<endl;

        cout <<"skipped over manual input"<<endl;

        if (!params.UniformHydrationShell) params.maxhs=8.0;

        if (params.VolumeOption=="Default") SetDefaultAtomicRadii(atomr, params);
        else if (params.VolumeOption=="FraserSuzukiMacRea") SetFraserSuzukiMacReaAtomicRadii(atomr);
        else
        {
                cout <<"ERROR: Unrecognized VolumeOption "<<params.VolumeOption<<endl;
                cout <<"There are other VolumeOptions.  This should be fixed."<<endl;
        }
        cout <<"After VolumeOption==FraseSuzukiMacRea"<<endl;
        cout <<"atomr[HydrationShell]= "<<atomr[HydrationShell]<<endl;

        if (params.UseCharges)
        {
                param="/home/jouko/WAXS/charmm.prm";
                cout <<"Enter parameter file path: ";
                cin.ignore();
                getline(cin, param);
                strcpy(paramfile, param.c_str());
                valid=exist(paramfile);
                //if (!valid) FileNotValid(paramfile, "tinker parameter file");
                if (!valid) 
                {
                        cout <<"ERROR: Unable to open tinker parameter file"<<endl;
                        exit(EXIT_FAILURE);
                }
                GetCharges(paramfile, Atoms);
        }
        InitializeIntensity(i, params.points);
        GetAngles(i, params);
        AtomScatteringFactor(i, params.ScatteringType);
        if (params.CalcVacuumIntensity)
        {
                InitializeIntensity(IntensityVacuum, params.points);
                GetAngles(IntensityVacuum, params);
                AtomScatteringFactor(IntensityVacuum, params.ScatteringType);
        }
        if (params.CalcNoHydrationIntensity)
        {
                InitializeIntensity(IntensityNoHydration, params.points);
                GetAngles(IntensityNoHydration, params);
                AtomScatteringFactor(IntensityNoHydration, params.ScatteringType);
        }
        cout <<"After ReadAngleFile"<<endl;

        if (params.ExperimentFile!="") 
        {
                ReadIntensityFile(params.ExperimentFile, i.s, i.expi, i.error);
                params.points=i.s.size();
        }
        cout <<"params.points= "<<params.points<<" i.s.size()= "<<i.s.size()<<endl;
        seconds1=time(NULL);

        cout <<"about to read pdb"<<endl;
        GetAtoms(Atoms, numatom, params);
        ExcludedRg=calcExcludedVolumeRg(atomr, params.contrast, Atoms, params);
        if (params.CalcPrFromExperiment) PrFromScattering(params.points, i.s, i.expi, i.error, params);
        cout <<"TotalParticles= "<<TotalParticles<<endl;
        cout <<"About to select option for intensity\n";

        SetVDWRadius(Atoms, params, solventatomr);
        if (params.excluded=="Lattice")
        {
                for (int j=0;j<NumAtomTypes;j++)
                {
                        if (j!=ExcludedVolume && j!=HydrationShell) atomr[j]=0;
                }
        }
        cout <<"Before waxs. atomr[HydrationShell]= "<<atomr[HydrationShell]<<endl;
        waxs(params, natom, atomr, i, IntensityVacuum, IntensityNoHydration, Atoms, numatom);
        HydrationEffects(i, IntensityVacuum, IntensityNoHydration, params);	
        FindNumAtom(numatom, Atoms);

        cout <<endl;
        if (params.IntensityOutputFile=="")
        {
                SetIntensityOutputFileName(params.xyz, intensityfile, logfile);
        }
        else
        {
                intensityfile=params.IntensityOutputFile;
                AddIndexToFile(intensityfile);
                logfile=RemoveExtension(intensityfile)+"Log.txt";
        }
        if (params.PrFromStructureFile=="") params.PrFromStructureFile=RemoveExtension(intensityfile)+"PrFT.txt";
        if (params.PrFromStructureFile2=="") PrFromStructureFile2=RemoveExtension(intensityfile)+"Pr.txt";
        else PrFromStructureFile2=params.PrFromStructureFile2;
        if (params.CalcPrFromStructure) CalcPr2(atomr, 0.5, params.contrast, params.hsdensity, PrFromStructureFile2, Atoms, params);
        if (params.CalcPrFromStructure2) CalcPrFromScattering(atomr, params.bin, params.contrast, params.hsdensity, params.PrFromStructureFile, Atoms, params);
        //TestScatteringFactorPr(atomr[ExcludedVolume], params);

        if (params.CalcRg)
        {
                Real Rg2;

                VacuumRg=calcVacuumRg(Atoms);
                ExcludedRg=calcExcludedVolumeRg(atomr, params.contrast, Atoms, params);
                Rg=radiusgyration(atomr, params.contrast, params.hsdensity, Atoms, params);
                Rg2=gyrationfromscattering(i.calc, i.s);
                cout <<"VacuumRg= "<<VacuumRg<<endl;
                cout <<"ExcludedRg= "<<ExcludedRg<<endl;
                cout <<"Rg= "<<Rg<<endl;
                cout <<"Rg from scattering= "<<Rg2<<endl;
        }
        DeleteVector(pr);
        if (argc>1)
        {
                if (params.ExperimentFile!="")
                {
                        chisqr=log_normalize(i.calc, i.expi, params.BeginFit, params.EndFit, i.error);
                        rfactor=0;
                        sum=0;

                        for (int t=params.BeginFit;t<params.EndFit;t++)
                        {
                                rfactor += abs(i.calc[t]-i.expi[t]);
                                sum += i.expi[t];
                        }
                        rfactor=rfactor/sum;
                        cout <<"chisqr= "<<chisqr<<endl;
                }
        }

        strcpy(charintensity, intensityfile.c_str());
        strcpy(charlogfile, logfile.c_str());
        cout <<"Writing intensities to "<<charintensity<<endl;	
        OutputIntensity(i, intensityfile);

        seconds2=time(NULL);

        time_t now=time(0);

        tm* localtm=localtime(&now);

        ofstream waxslog(charlogfile, ios::app);
        //waxslog <<"Protein file: "<< params.xyz <<endl;
        waxslog <<"Version: "<<Version<<endl;
        waxslog <<"Date: "<<asctime(localtm);
        //waxslog <<"Optionfile: "<<ParameterFile<<endl;
        waxslog <<"Solvent density: "<<params.contrast<<endl;
        waxslog <<"Excluded volume sphere type: "<<params.ExcludedVolumeSphereType<<endl;
        waxslog <<"Uniform Hydration: "<<BoolToStr(params.UniformHydrationShell)<<endl;
        if (params.CalcRg)
        {
                waxslog <<"Rg in vacuum: "<<VacuumRg<<endl;
                waxslog <<"Rg of protein without hydration shell: "<<ExcludedRg<<endl;
                waxslog <<"Rg with user defined parameters: "<<Rg<<endl;
        }
        else
        {
                waxslog <<"Hydration shell density: "<<params.hsdensity<<endl;
        }
        waxslog <<"Solvent method: "<<params.excluded<<endl
                <<"Histogram: "<<params.histogram<<endl;
        if (params.excluded=="Superimpose" || params.excluded=="Lattice")
        {
                waxslog <<"Hydrogen radius: "<<atomr[HYDROGEN]<<endl
                        <<"Carbon radius: "<<atomr[CARBON]<<endl
                        <<"Nitrogen radius: "<<atomr[NITROGEN]<<endl
                        <<"Oxygen radius: "<<atomr[OXYGEN]<<endl
                        <<"Sulfur radius: "<<atomr[SULFUR]<<endl
                        <<"Iron radius: "<<atomr[IRON]<<endl
                        <<"Phosphorus: "<<atomr[PHOSPHORUS]<<endl;
        }
        if (params.AngleFile!="" && params.ReadAngleFile) waxslog <<"Angle file: "<<params.AngleFile<<endl;
        if (params.ExperimentFile!="") waxslog <<"Experimental data file: "<<params.ExperimentFile<<endl;
        waxslog <<"Trajectory: "<<BoolToStr(params.DoTrajectory)<<endl;
        waxslog <<"Optimization: "<<BoolToStr(params.DoFit)<<endl;
        if (params.DoFit || (argc>1 && params.ExperimentFile!=""))
        {
                waxslog <<"Begin Fit: "<<params.BeginFit<<endl
                        <<"End Fit: "<<params.EndFit<<endl
                        <<"Chi sqr: "<<chisqr<<endl;
        }
        waxslog <<"Run time: "<<seconds2-seconds1<<" seconds"<<endl;

        cout << "done"<<endl;

        if (argc==1)
        {
                cin >> TotalParticles;
        }
        else cout <<"Return"<<endl;

        return 0;
}
