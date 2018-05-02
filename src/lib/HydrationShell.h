#ifndef _HydrationShell_included_
#define _HydrationShell_included_

# include "AtomCode.h"
# include "concavity.h"
# include "Constants.h"
# include "ConvertStructures.h"
# include "GetAtomType.h"
# include "LinkedList.h"
# include "MinMax.h"
# include "ReadGofRFile.h"
# include "ResidueCode.h"
# include "StringUtils.h"
# include "Structures.h"
# include "TypeDef.h"
# include "VectorManip.h"
# include "WritePdb.h"

int CountNumCubesInHydrationShell(vector<CubeStruct> &Cubes)
{
        int ncube=Cubes.size();
        int HydrationShellCubes=0;
        for (int i=0;i<ncube;i++)
        {
                if (Cubes[i].LocatedIn==InHydrationShell)
                {
                        HydrationShellCubes++;
                }
        }
        cout <<"HydrationShellCubes= "<<HydrationShellCubes<<endl;
        return HydrationShellCubes;
}

Real FindLargestVDWRadius(vector<AtomStruct> &Atoms)
{
        int natom=FindNumProteinAtoms(Atoms);
        Real LargestRadius=0;

        for (int i=0;i<natom;i++)
        {
                //PrintAtomInfo(Atoms[i]);
                //cout <<"Atoms["<<i<<"].vdw= "<<Atoms[i].vdw<<endl;
                if (Atoms[i].vdw>LargestRadius) LargestRadius=Atoms[i].vdw;
        }
        return LargestRadius;
}

Real IsOverlap(vector<AtomStruct> Atoms, Real xgrid, Real ygrid, Real zgrid, int &nearest, int natom)
{
        int n;
        Real CutOff, dist, LargestRadius, mindist=1000.0;
        Real dx, dy, dz;
        CutOff=mindist*mindist;

        LargestRadius=FindLargestVDWRadius(Atoms);
        for (n=0;n<natom;n++)
        {
                dx=Atoms[n].x-xgrid;
                dy=Atoms[n].y-ygrid;
                dz=Atoms[n].z-zgrid;
                dist=dx*dx+dy*dy+dz*dz;
                //cout <<"dist= "<<dist<<endl;
                if (dist<CutOff)
                {
                        //if (Print) cout <<"n= "<<n<<" dist= "<<dist<<" CutOff= "<<CutOff<<endl;
                        //cout <<"dist= "<<dist<<" Atoms["<<n<<"].vdw= "<<Atoms[n].vdw<<" mindist= "<<mindist<<endl;
                        dist=sqrt(dist)-Atoms[n].vdw;
                        CutOff=(dist+LargestRadius)*(dist+LargestRadius);

                        if (dist<mindist)
                        {
                                //if (Print) cout <<"dist= "<<dist<<" mindist= "<<mindist<<endl;
                                mindist=dist;
                                nearest=n;
                                if (mindist<0) break;
                        }
                }
        }
        return mindist;
}

int location(Real xgrid, Real ygrid, Real zgrid, int &nearest, Real &mindist, vector<AtomStruct> &Atoms, Real maxhs)
{
        //This needs to be compared against MD simulations to validate.
        //Detemines if a point is inside the protein, in the hydration shell, or outside of the protein
        //bool Print=false;
        bool accept=false;
        int located, NearestTest, natom;
        //Real buffer=0;   //The hydration shell atom radius minus the distance a hydration shell atom can stick out of the hydtation shell
        Real ProteinHydrationBoundary, HydrationSolutionBoundary;
        Real MinDistTest;
        Real dx, dy, dz;
        Real xtest, ytest, ztest;
        Real ProbeRadius=1.4;
        Real r;
        Real radii;
        //cout <<"In location"<<endl;
        located=UNK;
        Real HydrationShellR=0.75;
        ProteinHydrationBoundary=ProbeRadius+HydrationShellR*0.5;
        HydrationSolutionBoundary=maxhs;	
        natom=FindNumProteinAtoms(Atoms);
        //cout <<"natom= "<<natom<<endl;
        mindist=IsOverlap(Atoms, xgrid, ygrid, zgrid, nearest, natom);
        //cout <<"mindist= "<<mindist<<endl;
        if (mindist<ProteinHydrationBoundary) located=InProtein;
        if (mindist>ProteinHydrationBoundary && mindist<HydrationSolutionBoundary) located=InHydrationShell;
        if (mindist>HydrationSolutionBoundary) located=InSolution;
        //if (Print) cout <<"mindist= "<<mindist<<" ProbeRadius= "<<ProbeRadius<<" maxhs= "<<maxhs<<" buffer= "<<buffer<<endl;
        if (mindist>ProteinHydrationBoundary && mindist<ProteinHydrationBoundary+ProbeRadius)
        {
                dx=xgrid-Atoms[nearest].x;
                dy=ygrid-Atoms[nearest].y;
                dz=zgrid-Atoms[nearest].z;
                r=1.0/sqrt(dx*dx+dy*dy+dz*dz);
                radii=Atoms[nearest].vdw+ProbeRadius;
                xtest=Atoms[nearest].x+radii*dx*r;
                ytest=Atoms[nearest].y+radii*dy*r;
                ztest=Atoms[nearest].z+radii*dz*r;
                MinDistTest=IsOverlap(Atoms, xtest, ytest, ztest, NearestTest, natom);

                if (MinDistTest>ProbeRadius) accept=true;

                if (!accept)
                {
                        xtest=xgrid+ProbeRadius*dx*r;
                        ytest=ygrid+ProbeRadius*dy*r;
                        ztest=zgrid+ProbeRadius*dz*r;
                        MinDistTest=IsOverlap(Atoms, xtest, ytest, ztest, NearestTest, natom);
                        if (MinDistTest>ProbeRadius) accept=true;
                }
                if (accept) located=InHydrationShell;
                else located=InProtein;
        }
        //if (Print) cout <<"location= "<<location<<endl;
        if (located==UNK)
        {
                cout <<"Error.  Unable to find location."<<endl;
                located=InProtein;
        }
        return located;
}

bool IsOverlapPairlist(vector<AtomStruct> &Atoms, int PairList[], int PairListNum, Real CutOff, Real CutOff2, Real ProbeRadius, Real xgrid, Real ygrid, Real zgrid)
{
        Real DistToSphere, r2;
        Real dx, dy, dz;
        for (int n=0;n<PairListNum;n++)
        {
                dz=Atoms[PairList[n]].z-zgrid;
                if ( dz<CutOff && dz>-CutOff )
                {
                        dx=Atoms[PairList[n]].x-xgrid;
                        dy=Atoms[PairList[n]].y-ygrid;
                        r2=dx*dx+dy*dy+dz*dz;
                        if ( r2<CutOff2 )
                        {
                                DistToSphere=sqrt(r2)-Atoms[PairList[n]].vdw;
                                if (DistToSphere<ProbeRadius) return true;
                        }
                }
        }
        return false;
}

void FindNearestVDWSurface(CubeStruct &Cube, vector<AtomStruct> &Atoms, int natom, Real LargestRadius)
{
        Real CutNearest=10000.0;
        Real r, dx, dy, dz;
        Cube.MinDist=10000.0;
        for (int m=0;m<natom;m++)
        {
                dx=Cube.x-Atoms[m].x;
                if (dx<CutNearest && dx>-CutNearest)
                {
                        dy=Cube.y-Atoms[m].y;
                        if (dy<CutNearest && dy>-CutNearest)
                        {
                                dz=Cube.z-Atoms[m].z;
                                if (dx<CutNearest && dz>-CutNearest)
                                {
                                        r=sqrt(dx*dx+dy*dy+dz*dz)-Atoms[m].vdw;
                                        if (r<Cube.MinDist)
                                        {
                                                Cube.MinDist=r;
                                                Cube.nearest=m;
                                                CutNearest=Cube.MinDist+LargestRadius;
                                        }
                                }
                        }
                }
        }
        Cube.MinDist+=Atoms[Cube.nearest].vdw;
}

void FindNearestHydrationRadius(CubeStruct &Cube, vector<AtomStruct> &Atoms, int natom, Real LargestRadius)
{
        bool BoolPrint=false;
        Real Maxdist=1000.0;
        Real r, dx, dy, dz;
        Cube.MinDist=1000.0;
        Cube.SecondMinDist=1000.0;
        Real d=0.1;
        if (abs(Cube.x-22.394)<d && abs(Cube.y-2.541)<d && abs(Cube.z-12.105)<d)
        {
                BoolPrint=true;
        }
        for (int n=0;n<natom;n++)
        {
                dx=Cube.x-Atoms[n].x;

                if (dx<Maxdist && dx>-Maxdist)
                {
                        dy=Cube.y-Atoms[n].y;
                        if (dy<Maxdist && dy>-Maxdist)
                        {
                                dz=Cube.z-Atoms[n].z;
                                if (dz<Maxdist && dz>-Maxdist)
                                {
                                        r=sqrt(dx*dx+dy*dy+dz*dz)-Atoms[n].HydrationRadius;
                                        if (r<Cube.MinDist)
                                        {
                                                Cube.SecondMinDist=Cube.MinDist;
                                                Cube.SecondNearest=Cube.nearest;
                                                Cube.MinDist=r;
                                                Cube.nearest=n;
                                                Maxdist=Cube.SecondMinDist+LargestRadius;
                                        }
                                        if (r>Cube.MinDist && r<Cube.SecondMinDist)
                                        {
                                                Cube.SecondMinDist=r;
                                                Cube.SecondNearest=n;
                                        }
                                        if (BoolPrint)
                                        {
                                                cout <<"n= "<<n<<endl;
                                                //Cubes[m].x+=0.000001;
                                                cout <<setprecision(10)<<"xgrid= "<<Cube.x<<" ygrid= "<<Cube.y<<" zgrid= "<<Cube.z<<endl;
                                                PrintAtomInfo(Atoms[n]);			
                                                cout <<"dx= "<<dx<<" dy= "<<dy<<" dz= "<<dz<<endl;
                                                cout <<"r= "<<r<<" SolventAtomR["<<Atoms[n].atomid<<"]= "<<Atoms[n].HydrationRadius<<endl;
                                                cout <<"Maxdist= "<<Maxdist<<" mindist= "<<Cube.MinDist<<" nearest= "<<Cube.nearest<<endl;
                                                cout <<endl;
                                        }
                                }
                        }
                }
        }
        if (BoolPrint)
        {
                cout <<"At end Cube.nearest= "<<Cube.nearest<<endl;
                cout <<"Atoms.HydrationRadius= "<<Atoms[Cube.nearest].HydrationRadius<<endl;
        }
        Cube.MinDist+=Atoms[Cube.nearest].HydrationRadius;
}

Real GetLargestRadius(vector<AtomStruct> &Atoms)
{
        int natom=Atoms.size();
        Real LargestRadius=0;

        for (int i=0;i<natom;i++)
        {
                if (Atoms[i].HydrationRadius>LargestRadius)
                {
                        LargestRadius=Atoms[i].HydrationRadius;
                }
        }
        return LargestRadius;
}

void FindNearestHydrationRadius(CubeStruct &Cube, vector<AtomStruct> &Atoms)
{
        int natom=Atoms.size();
        Real LargestRadius;

        natom=FindNumProteinAtoms(Atoms);
        LargestRadius=GetLargestRadius(Atoms);
        FindNearestHydrationRadius(Cube, Atoms, natom, LargestRadius);
}

void DetermineNearestForAllCubes(vector<CubeStruct> &Cubes, vector<AtomStruct> &Atoms, int natom, Real maxhs)
{
        bool BoolPrint=false;
        Real LargestRadius, Maxdist;
        int GridPoints=Cubes.size();

        LargestRadius=0;
        cout <<"natom= "<<natom<<endl;
        for (int n=0;n<natom;n++)
        {
                cout <<"n= "<<n<<" natom= "<<natom<<endl;
                if (Atoms[n].HydrationRadius>LargestRadius) 
                {
                        LargestRadius=Atoms[n].HydrationRadius;
                        cout <<"LargestRadius= "<<LargestRadius<<endl;
                        PrintAtomInfo(Atoms[n]);
                }
        }
        //Real d=0.1;
        for (int m=0;m<GridPoints;m++)
        {
                Maxdist=1000.0;
                Cubes[m].MinDist=1000.0;
                cout <<"m= "<<m<<" GridPoints= "<<GridPoints<<endl;
                if (Cubes[m].LocatedIn==InHydrationShell)
                {
                        FindNearestHydrationRadius(Cubes[m], Atoms, natom, LargestRadius);
                }
                BoolPrint=false;
        }
        for (int n=0;n<GridPoints;n++)
        {
                cout <<"n= "<<n<<" GridPoints= "<<GridPoints<<endl;
                if (Cubes[n].LocatedIn==InHydrationShell)
                {
                        if (Cubes[n].MinDist>maxhs) Cubes[n].LocatedIn=InSolution;
                }
                //cout <<"Cubes["<<n<<"].MinDist= "<<Cubes[n].MinDist
                //<<" nearest= "<<Cubes[n].nearest
                //<<" HydrationRadius= "<<Atoms[Cubes[n].nearest].HydrationRadius<<endl;
        }
}

void location4(int natom, Real xmin, Real xmax, Real ymin, Real ymax, Real zmin, Real zmax, vector<CubeStruct> &Cubes, vector<AtomStruct> Atoms, Real maxhs)
{
        //!!!This needs to be fixed before using
        //Determines which set of points in an array are in the protein, hydration shell, or bulk solvent.
        //Also determines the atom surface which is closest to each of the points.
        //Moves a sphere with radius=1.4 A on lattice.  If no overlaps with sphere then every grid point in
        //the sphere is in hydration shell or bulk solvent.
        cout <<"In location 4"<<endl;
        int GridPoints=Cubes.size();
        int PairListNum, PairListNum2, PairListNum3;
        int *PairList, *PairList2, *PairList3;
        bool SphereOverlap;
        Real CutOff, CutOff2;
        Real LargestRadius;
        Real inc;
        Real ProbeRadius, ProbeRadius2;
        Real xgrid, ygrid, zgrid;
        Real XPairMax, XPairMin;
        Real YPairMax, YPairMin;
        Real XPairMax2, XPairMin2;
        Real YPairMax2, YPairMin2;
        Real dx, dy, dz;
        clock_t end, start, TotalInIsOverlap=0, TotalInFindNearest=0;
        //clock_t TotalInFindAllNearest=0;
        clock_t TotalMakingPairList=0;
        clock_t TotalAfterSphereOverlap=0;
        time_t seconds1, seconds2;

        seconds1=time(NULL);
        PairList=new int[natom];
        PairList2=new int[GridPoints];
        PairList3=new int[GridPoints];

        for (int n=0;n<natom;n++) PairList[n]=0;

        for (int n=0;n<GridPoints;n++) 
        {
                PairList2[n]=0;
                PairList3[n]=0;
        }

        for (int n=0;n<GridPoints;n++) Cubes[n].LocatedIn=InProtein;
        LargestRadius=0.0;
        cout <<"natom= "<<Atoms.size()<<endl;
        cout <<"GridPoints= "<<GridPoints<<endl;
        for (int n=0;n<natom;n++)
        {
                if (Atoms[n].vdw>LargestRadius) LargestRadius=Atoms[n].vdw;
        }
        seconds2=time(NULL);
        cout <<"Initialization time was "<<seconds2-seconds1<<" seconds."<<endl;


        //ProbeRadius=0.75;
        ProbeRadius=1.4;
        inc=0.5;
        CutOff=ProbeRadius+LargestRadius;
        CutOff2=CutOff*CutOff;
        ProbeRadius2=ProbeRadius*ProbeRadius;

        for (xgrid=xmin;xgrid<xmax;xgrid+=inc)
        {
                //cout <<"xgrid= "<<xgrid<<endl;
                PairListNum3=0;
                XPairMax2=xgrid+ProbeRadius;
                XPairMin2=xgrid-ProbeRadius;
                for (int n=0;n<GridPoints;n++)
                {
                        if (Cubes[n].x>XPairMin2 && Cubes[n].x<XPairMax2)
                        {
                                PairList3[PairListNum3]=n;
                                PairListNum3++;
                        }
                }
                for (ygrid=ymin;ygrid<ymax;ygrid+=inc)
                {
                        start=clock();
                        PairListNum=0;
                        XPairMax=xgrid+CutOff;
                        XPairMin=xgrid-CutOff;
                        YPairMax=ygrid+CutOff;
                        YPairMin=ygrid-CutOff;

                        for (int n=0;n<natom;n++)
                        {
                                if (Atoms[n].x>XPairMin && Atoms[n].x<XPairMax && Atoms[n].y>YPairMin && Atoms[n].y<YPairMax)
                                {
                                        PairList[PairListNum]=n;
                                        PairListNum++;
                                }
                        }

                        PairListNum2=0;
                        YPairMax2=ygrid+ProbeRadius;
                        YPairMin2=ygrid-ProbeRadius;

                        for (int n=0;n<PairListNum3;n++)
                        {
                                if (Cubes[n].y>YPairMin2 && Cubes[n].y<YPairMax2)
                                {
                                        PairList2[PairListNum2]=PairList3[n];
                                        PairListNum2++;
                                }
                        }

                        end=clock();
                        TotalMakingPairList+=(end-start);
                        for (zgrid=zmin;zgrid<zmax;zgrid+=inc)
                        {
                                start=clock();
                                SphereOverlap=IsOverlapPairlist(Atoms, PairList, PairListNum, CutOff, CutOff2, ProbeRadius, xgrid, ygrid, zgrid);
                                end=clock();
                                TotalInIsOverlap+=(end-start);
                                start=clock();
                                if (!SphereOverlap)
                                {
                                        for (int n=0;n<PairListNum2;n++)
                                        {
                                                if (Cubes[PairList2[n]].LocatedIn==InProtein)
                                                {
                                                        dz=Cubes[PairList2[n]].z-zgrid;
                                                        if (dz<ProbeRadius && dz>-ProbeRadius)
                                                        {
                                                                dx=Cubes[PairList2[n]].x-xgrid;
                                                                dy=Cubes[PairList2[n]].y-ygrid;
                                                                if ( (dx*dx+dy*dy+dz*dz)<ProbeRadius2)
                                                                {
                                                                        start=clock();
                                                                        FindNearestVDWSurface(Cubes[PairList2[n]], Atoms, natom, LargestRadius);
                                                                        end=clock();
                                                                        TotalInFindNearest+=(end-start);
                                                                        if (Cubes[PairList2[n]].MinDist<maxhs) Cubes[PairList2[n]].LocatedIn=InHydrationShell;
                                                                        else Cubes[PairList2[n]].LocatedIn=InSolution;
                                                                }
                                                        }
                                                }
                                        }	
                                }
                                end=clock();
                                TotalAfterSphereOverlap+=(end-start);
                        }
                }
        }

        for (int n=0;n<GridPoints;n++)
        {
                //cout <<"Cubes["<<n<<"].MinDist= "<<Cubes[n].MinDist<<" Cubes.nearest= "<<Cubes[n].nearest<<" Atoms.vdw= "<<Atoms[Cubes[n].nearest].vdw<<endl;
                Cubes[n].MinDist+=Atoms[Cubes[n].nearest].vdw;
        }
        cout <<"Time in IsOverlap was "<<Real(TotalInIsOverlap)/Real(CLOCKS_PER_SEC)<<" seconds."<<endl;
        cout <<"Time in FindNearest was "<<Real(TotalInFindNearest)/Real(CLOCKS_PER_SEC)<<" seconds."<<endl;
        cout <<"Time in MakingPairList was "<<Real(TotalMakingPairList)/Real(CLOCKS_PER_SEC)<<" seconds."<<endl;
        cout <<"Time in AfterOverlap was "<<Real(TotalAfterSphereOverlap)/Real(CLOCKS_PER_SEC)<<" seconds."<<endl;
        seconds1=time(NULL);
        start=clock();
        DetermineNearestForAllCubes(Cubes, Atoms, natom, maxhs);
        end=clock();
        seconds2=time(NULL);
        cout <<"Time in FindAllNearest was "<<Real(end-start)/Real(CLOCKS_PER_SEC)<<" seconds."<<endl;
        cout <<"WallClock time in FindAllNearest was "<<seconds2-seconds1<<" seconds."<<endl;

        delete [] PairList;
        delete [] PairList2;
        delete [] PairList3;
        cout <<"Exiting location4"<<endl;
}

Real ASAOfAtomSharkeRupley(vector<AtomStruct> &Atoms, int index)
{
        int natom, accessible, points, nearest;
        Real ProbeRadius, d, r, sinphi, cosphi, theta, thetainc, nslices, asa;
        Real VectX, VectY, VectZ;
        Real xgrid, ygrid, zgrid;
        Real overlap;
        accessible=0;
        points=0;
        ProbeRadius=1.4;
        d=0.0001;                               //This is a hack, because the points should not overlap with the atoms they belong to.
        r=Atoms[index].vdw+ProbeRadius+d;
        nslices=10;
        thetainc=nslices;
        cosphi=1.0/nslices-1.0;
        theta=pi/thetainc;
        natom=FindNumProteinAtoms(Atoms);
        while (true)
        {
                sinphi=sqrt(1.0-cosphi*cosphi);
                VectX=sinphi*cos(theta);
                VectY=sinphi*sin(theta);
                VectZ=cosphi;
                theta+=2.0*pi/thetainc;
                if (theta>2.0*pi)
                {
                        theta=pi/thetainc;
                        cosphi+=2.0/nslices;
                        if (cosphi>1.0) break;
                }
                //cout <<"theta= "<<theta<<" cosphi= "<<cosphi<<endl;
                xgrid=Atoms[index].x+VectX*r;
                ygrid=Atoms[index].y+VectY*r;
                zgrid=Atoms[index].z+VectZ*r;
                overlap=IsOverlap(Atoms, xgrid, ygrid, zgrid, nearest, natom);
                //cout <<"overlap= "<<overlap<<endl;
                if (overlap>ProbeRadius) accessible++;
                points++;
        }
        //cout <<"accessible= "<<accessible<<" points= "<<points<<endl;
        asa=Real(accessible)*4.0*pi*r*r/Real(points);
        return asa;
}

Real ASASharkeRupley(vector<AtomStruct> &Atoms)
{
        int natom;
        Real asa=0;

        natom=FindNumProteinAtoms(Atoms);
        for (int i=0;i<natom;i++)
        {
                asa+=ASAOfAtomSharkeRupley(Atoms, i);
        }
        return asa;
}

Real FindSAOverlap(AtomStruct &Atom1, AtomStruct &Atom2)
{
        Real r1, r2, x, costheta, overlap, ProbeRadius;
        ProbeRadius=1.4;
        r1=Atom1.vdw+ProbeRadius;
        r2=Atom2.vdw+ProbeRadius;
        x=AtomDistance(Atom1, Atom2);
        if (r1+r2<x) return 0;
        costheta=(r2*r2+x*x-r1*r1)/(2.0*x*r2);
        overlap=2.0*pi*r2*r2*(1.0-costheta);
        //cout <<"r1= "<<r1<<" r2= "<<r2<<" x= "<<x<<" costheta= "<<costheta<<" overlap= "<<overlap<<endl;
        return overlap;
}

Real ASAOfAtomLCPO(vector<AtomStruct> &Atoms, int index)
{
        //Need to read actual paper
        int natom;
        Real ProbeRadius=1.4;
        Real overlap, r, asa, area;
        natom=FindNumProteinAtoms(Atoms);
        overlap=0;
        r=Atoms[index].vdw+ProbeRadius;
        area=4.0*pi*r*r;
        asa=area;
        for (int i=0;i<natom;i++)
        {
                overlap=FindSAOverlap(Atoms[i], Atoms[index]);
                if (i!=index) asa*=(area-overlap)/area;
        }
        return asa;
}

Real ASALCPO(vector<AtomStruct> &Atoms)
{
        int natom;
        Real asa=0;

        natom=FindNumProteinAtoms(Atoms);
        for (int i=0;i<natom;i++)
        {
                asa+=ASAOfAtomLCPO(Atoms, i);
        }
        return asa;
}

void QuickSurfaceASA(vector<AtomStruct> &Atoms, vector<AtomStruct> &SurfaceAtoms)
{
        int natom;
        Real threshold, asa;
        natom=FindNumProteinAtoms(Atoms);
        threshold=1.0;
        for (int i=0;i<natom;i++)
        {
                asa=ASAOfAtomLCPO(Atoms, i);
                cout <<"asa= "<<asa<<endl;
                if (asa>threshold) SafePushBack(SurfaceAtoms, Atoms[i], "SurfaceAtoms");
        }
}

void QuickSurface(vector<AtomStruct> Atoms, vector<AtomStruct> &Surface)
{
        vector<bool> BoolSurface;
        int natom, CurrentAtom;
        Real sinphi, cosphi, theta;
        Real VectX, VectY, VectZ;
        Real thetainc, nslices=100;
        Real current, projection;
        Real xn, yn, zn;
        center(Atoms);
        
        thetainc=nslices;
        cosphi=1.0/nslices-1.0;
        theta=pi/thetainc;
        natom=FindNumProteinAtoms(Atoms);
        SafeAlloc(BoolSurface, false, natom, "BoolSurface");

        while (true)
        {
                sinphi=sqrt(1.0-cosphi*cosphi);
                VectX=sinphi*cos(theta);
                VectY=sinphi*sin(theta);
                VectZ=cosphi;
                theta+=2.0*pi/thetainc;
                if (theta>2.0*pi)
                {
                        theta=pi/thetainc;
                        cosphi+=2.0/nslices;
                        if (cosphi>1.0) break;
                }
                cout <<"theta= "<<theta<<" cosphi= "<<cosphi<<endl;
                current=0;
                for (int n=0;n<natom;n++)
                {
                        xn=Atoms[n].x;
                        yn=Atoms[n].y;
                        zn=Atoms[n].z;
                        projection=VectX*xn+VectY*yn+VectZ*zn;
                        if (projection+Atoms[n].HydrationRadius>current)
                        {
                                current=projection+Atoms[n].HydrationRadius;
                                CurrentAtom=n;
                        }
                }
                BoolSurface[CurrentAtom]=true;
        }
        for (int i=0;i<natom;i++)
        {
                if (BoolSurface[i]) SafePushBack(Surface, Atoms[i], "Surface");
        }
}

bool DetermineLocalCurvature(vector<PlaneStruct> &planes, Real x, Real y, Real z)
{
        int nplane=planes.size();
        Real z2;
        //cout <<"In DetermineLocalCurvature nplane= "<<nplane<<endl;
        for (int i=0;i<nplane;i++)
        {
                z2=planes[i].coef[0]*x+planes[i].coef[1]*y+planes[i].coef[2];
                /*
                cout <<"z2= "<<z2<<" z= "<<z<<" planes["<<i<<"].location= "<<planes[i].location<<endl;
                cout <<"x= "<<x<<" y= "<<" z= "<<z<<endl;
                cout <<"planes["<<i<<"].coef[0]= "<<planes[i].coef[0]<<endl;
                cout <<"planes["<<i<<"].coef[1]= "<<planes[i].coef[1]<<endl;
                cout <<"planes["<<i<<"].coef[2]= "<<planes[i].coef[2]<<endl;
                */
                if (z2>z && planes[i].location=="AllAbove") return false;
                if (z2<z && planes[i].location=="AllBelow") return false;
        }
        return true;
}

bool DetermineLocalCurvature(vector<PlaneStruct> &planes, CubeStruct &Cube)
{
        return DetermineLocalCurvature(planes, Cube.x, Cube.y, Cube.z);
}

void DetermineLocalCurvature(vector<AtomStruct> &Atoms, lattice &Cubes)
{
        int MaxXBin, MaxYBin, MaxZBin;
        vector<AtomStruct> SurfaceAtoms;
        vector<PlaneStruct> planes;

        MaxXBin=Cubes.size();
        MaxYBin=Cubes[0].size();
        MaxZBin=Cubes[0][0].size();
        cout <<"In DetermineLocalCurvature"<<endl;
        QuickSurfaceASA(Atoms, SurfaceAtoms);
        //PrintPdb("/home2/jouko/project/WAXS/pdb/SurfaceAtoms.pdb", SurfaceAtoms);
        cout <<"Determined surface"<<endl;
        CalcPlanes(SurfaceAtoms, planes);
        cout <<"planes[0].coef[0]= "<<planes[0].coef[0]<<endl;
        cout <<"planes[0].coef[1]= "<<planes[0].coef[1]<<endl;
        cout <<"planes[0].coef[2]= "<<planes[0].coef[2]<<endl;
        for (int i=0;i<MaxXBin;i++)
        {
                for (int j=0;j<MaxYBin;j++)
                {
                        for (int k=0;k<MaxZBin;k++)
                        {
                                Cubes[i][j][k].concave=DetermineLocalCurvature(planes, Cubes[i][j][k]);
                        }
                }
        }
}

void DetermineLocalCurvature(vector<AtomStruct> &Atoms, vector<CubeStruct> &Cubes)
{
        int ncube;
        vector<AtomStruct> SurfaceAtoms;
        vector<PlaneStruct> planes;

        cout <<"In DetermineLocalCurvature"<<endl;
        QuickSurfaceASA(Atoms, SurfaceAtoms);
        PrintPdb("/home2/jouko/project/WAXS/pdb/SurfaceAtoms.pdb", SurfaceAtoms);
        cout <<"Determined surface"<<endl;
        CalcPlanes(SurfaceAtoms, planes);

        for (int i=0;i<ncube;i++)
        {
                Cubes[i].concave=DetermineLocalCurvature(planes, Cubes[i]);
        }
}

void Cavity(lattice &Cubes, ParamStruct params)
{
        bool XMinusBlocked, YMinusBlocked, ZMinusBlocked;
        bool XPlusBlocked, YPlusBlocked, ZPlusBlocked;
        bool BoolPrint=false;
        int BlockedSides;
        int nearest;
        int MaxXBin, MaxYBin, MaxZBin;
        int xbin, ybin, zbin;
        int xbin2, ybin2, zbin2;
        int xstart, ystart, zstart;
        int xend, yend, zend;
        int range=int(9.0/params.CubeSize);
        Real d=0.1;
        MaxXBin=Cubes.size();
        MaxYBin=Cubes[0].size();
        MaxZBin=Cubes[0][0].size();
        cout <<"In Cavity"<<endl;
        for (xbin=0; xbin<MaxXBin; xbin++)
        {
                cout <<"xbin= "<<xbin<<endl;
                for (ybin=0; ybin<MaxYBin; ybin++)
                {
                        for (zbin=0; zbin<MaxZBin; zbin++)
                        {
                                if (Cubes[xbin][ybin][zbin].LocatedIn!=InProtein)
                                {
                                        
                                        if (abs(Cubes[xbin][ybin][zbin].x-9.126)<d && abs(Cubes[xbin][ybin][zbin].y-5.440)<d && abs(Cubes[xbin][ybin][zbin].z+6.653)<d) BoolPrint=true;
                                        xstart=xbin-range;
                                        ystart=ybin-range;
                                        zstart=zbin-range;
                                        if (xstart<0) xstart=0;
                                        if (ystart<0) ystart=0;
                                        if (zstart<0) zstart=0;

                                        xend=xbin+range;
                                        yend=ybin+range;
                                        zend=zbin+range;
                                        if (xend>MaxXBin-1) xend=MaxXBin-1;
                                        if (yend>MaxYBin-1) yend=MaxYBin-1;
                                        if (zend>MaxZBin-1) zend=MaxZBin-1;

                                        if (BoolPrint)
                                        {
                                                cout <<"xstart= "<<xstart<<" ystart= "<<ystart<<" zstart= "<<zstart<<endl;
                                                cout <<"xend= "<<xend<<" yend= "<<yend<<" zend= "<<zend<<endl;
                                        }
                                        //if (xbin==29) cout <<"xend= "<<xend<<" yend= "<<yend<<" zend= "<<zend<<endl;
                                        //if (xbin==29) cout <<"xstart= "<<xstart<<" ystart= "<<ystart<<" zstart= "<<zstart<<endl;
                                        BlockedSides=0;
                                        XMinusBlocked=false;
                                        YMinusBlocked=false;
                                        ZMinusBlocked=false;
                                        XPlusBlocked=false;
                                        YPlusBlocked=false;
                                        ZPlusBlocked=false;
                                        for (xbin2=xstart; xbin2<xbin; xbin2++)
                                        {
                                                //if (xbin==29) cout <<"xbin2= "<<xbin2<<" GridValue["<<xbin2<<"]["<<ybin<<"]["<<zbin<<"]= "<<GridValue[xbin2][ybin][zbin]<<endl;
                                                if (Cubes[xbin2][ybin][zbin].LocatedIn==InProtein) 
                                                {
                                                        XMinusBlocked=true;
                                                        if (BoolPrint) cout <<"XMinus is Blocked"<<endl;
                                                        break;
                                                }
                                        }
                                        for (xbin2=xbin; xbin2<xend; xbin2++)
                                        {

                                                if (Cubes[xbin2][ybin][zbin].LocatedIn==InProtein)
                                                {
                                                        XPlusBlocked=true;
                                                        if (BoolPrint) cout <<"XPluss is Blocked"<<endl;
                                                        break;
                                                }
                                        }
                                        for (ybin2=ystart; ybin2<ybin; ybin2++)
                                        {
                                                if (BoolPrint)
                                                {
                                                        PrintCubeInfo(Cubes[xbin][ybin2][zbin]);
                                                }
                                                if (Cubes[xbin][ybin2][zbin].LocatedIn==InProtein)
                                                {
                                                        YMinusBlocked=true;
                                                        if (BoolPrint) cout <<"YMinus is Blocked"<<endl;
                                                        break;
                                                }
                                        }
                                        for (ybin2=ybin; ybin2<yend; ybin2++)
                                        {
                                                if (Cubes[xbin][ybin2][zbin].LocatedIn==InProtein)
                                                {
                                                        YPlusBlocked=true;
                                                        if (BoolPrint) cout <<"YPlus is Blocked"<<endl;
                                                        break;
                                                }
                                        }
                                        for (zbin2=zstart; zbin2<zbin; zbin2++)
                                        {
                                                if (Cubes[xbin][ybin][zbin2].LocatedIn==InProtein) 
                                                {
                                                        ZMinusBlocked=true;
                                                        if (BoolPrint) cout <<"ZMinus is Blocked"<<endl;
                                                        break;
                                                }
                                        }
                                        for (zbin2=zbin; zbin2<zend; zbin2++)
                                        {
                                                if (Cubes[xbin][ybin][zbin2].LocatedIn==InProtein)
                                                {
                                                        ZPlusBlocked=true;
                                                        if (BoolPrint) cout <<"ZPlus is Blocked"<<endl;
                                                        break;
                                                }
                                        }
                                        if ( (XMinusBlocked && XPlusBlocked) || (YMinusBlocked && YPlusBlocked) || (ZMinusBlocked && ZPlusBlocked)  ) Cubes[xbin][ybin][zbin].concave=true;
                                        else Cubes[xbin][ybin][zbin].concave=false;

                                        nearest=Cubes[xbin][ybin][zbin].nearest;
                                        BoolPrint=false;
                                }
                        }
                }
        }				
        cout <<"Done with cavity"<<endl;
}

void Cavity(vector<CubeStruct> &Cubes, ParamStruct params)
{
        lattice CubeLattice;
        for (int i=0;i<5;i++) PrintCubeInfo(Cubes[i]);
        cout <<endl;
        cout <<"Converting to lattice"<<endl;
        cout <<"Cubes.size()= "<<Cubes.size()<<endl;
        CountNumCubesInHydrationShell(Cubes);
        ConvertCubeVectorToLattice(Cubes, CubeLattice, params.CubeSize);
        Cavity(CubeLattice, params);
        ConvertCubeLatticeToCubeVector(CubeLattice, Cubes);
        cout <<"Done converting to Vector"<<endl;
        CountNumCubesInHydrationShell(Cubes);
        cout <<"Cubes.size()= "<<Cubes.size()<<endl;
        for (int i=0;i<5;i++) PrintCubeInfo(Cubes[i]);
}

void SafeMakeGrid(Real xmin, Real ymin, Real zmin, Real xmax, Real ymax, Real zmax, vector<AtomStruct> &Atoms, Real maxhs, Real CubeSize, vector<CubeStruct> &Cubes)
{
        bool MakeGridSafely=false;
        int LocatedIn, nearest;
        Real mindist;
        Real xgrid, ygrid, zgrid;
        CubeStruct Cube;
        DeleteVector(Cubes);
        InitializeCube(Cube);
        for (xgrid=xmin;xgrid<=xmax;xgrid+=CubeSize)
        {
                cout <<"xgrid= "<<xgrid<<" xmax= "<<xmax<<endl;
                for (ygrid=ymin;ygrid<=ymax;ygrid+=CubeSize)
                {
                        for (zgrid=zmin;zgrid<=zmax;zgrid+=CubeSize)
                        {
                                Cube.x=xgrid;
                                Cube.y=ygrid;
                                Cube.z=zgrid;
                                LocatedIn=location(xgrid, ygrid, zgrid, nearest, mindist, Atoms, maxhs);
                                Cube.LocatedIn=LocatedIn;
                                if (LocatedIn!=InSolution) SafePushBack(Cubes, Cube, "Cubes");
                                //cout <<"CubeNum= "<<Cubes.size();
                        }
                }
        }
}

bool MakeGrid(Real xmin, Real ymin, Real zmin, Real xmax, Real ymax, Real zmax, Real CubeSize, vector<CubeStruct> &Cubes)
{
        int ncube, nthcube=0;
        Real xgrid, ygrid, zgrid;
        CubeStruct Cube;
        DeleteVector(Cubes);
        InitializeCube(Cube);

        ncube=int((xmax-xmin)/CubeSize+1.5)*int((ymax-ymin)/CubeSize+1.5)*int((zmax-zmin)/CubeSize+1.5);
        if (!SafeAlloc(Cubes, Cube, ncube))
        {
                cout <<"Warning: Unable to allocate memory for cubes.  Using more memory efficient method"<<endl;
                return false;
        }
        for (xgrid=xmin;xgrid<=xmax;xgrid+=CubeSize)
        {
                for (ygrid=ymin;ygrid<=ymax;ygrid+=CubeSize)
                {
                        for (zgrid=zmin;zgrid<=zmax;zgrid+=CubeSize)
                        {
                                Cubes[nthcube].x=xgrid;
                                Cubes[nthcube].y=ygrid;
                                Cubes[nthcube].z=zgrid;
                                nthcube++;
                        }
                }
        }
        return true;
}

void Round(Real &x)
{
        x=Real(int(x*1000.0))/1000.0;
}

void ExpandMinMax(Real &xmin, Real &ymin, Real &zmin, Real &xmax, Real &ymax, Real &zmax, ParamStruct params)
{
        int n;
        n=int(floor((xmax+params.maxhs+params.XBoxLength*0.5)/params.CubeSize+1.0));
        xmax=-params.XBoxLength*0.5-params.CubeSize*0.5+params.CubeSize*Real(n);
        n=int(floor((xmin-params.maxhs+params.XBoxLength*0.5)/params.CubeSize+1.0));
        xmin=-params.XBoxLength*0.5-params.CubeSize*0.5+params.CubeSize*Real(n);
        n=int(floor((ymax+params.maxhs+params.YBoxLength*0.5)/params.CubeSize+1.0));
        ymax=-params.YBoxLength*0.5-params.CubeSize*0.5+params.CubeSize*Real(n);
        n=int(floor((ymin-params.maxhs+params.YBoxLength*0.5)/params.CubeSize+1.0));
        ymin=-params.YBoxLength*0.5-params.CubeSize*0.5+params.CubeSize*Real(n);
        n=int(floor((zmax+params.maxhs+params.ZBoxLength*0.5)/params.CubeSize+1.0));
        zmax=-params.ZBoxLength*0.5-params.CubeSize*0.5+params.CubeSize*Real(n);
        n=int(floor((zmin-params.maxhs+params.ZBoxLength*0.5)/params.CubeSize+1.0));
        zmin=-params.ZBoxLength*0.5-params.CubeSize*0.5+params.CubeSize*Real(n);
        //xmin+=0.0005;
        //zmin-=0.0005;
        //Round(xmin);
        //Round(ymin);
        //Round(zmin);
}

void EliminateCavities(vector<AtomStruct> &Atoms, Real HydrationShellRadius)
{
        bool changed=true;
        int m, n, TotalParticles=Atoms.size();
        int natom;
        int * group;
        Real dist, dsqr;
        Real dx, dy, dz;
        Real xn, yn, zn;
        group=new int[TotalParticles];
        natom=FindNumProteinAtoms(Atoms);
        for (n=0;n<TotalParticles;n++) group[n]=0;

        for (n=natom;n<TotalParticles;n++)
        {
                if (Atoms[n].atomid==HydrationShell) group[n]=n;
        }

        dsqr=4.0*HydrationShellRadius*HydrationShellRadius;
        while (changed)
        {
                changed=false;
                for (n=TotalParticles-1;n>natom-1;n=n-1)
                {
                        xn=Atoms[n].x;
                        yn=Atoms[n].y;
                        zn=Atoms[n].z;
                        if (Atoms[n].atomid==HydrationShell)
                        {
                                for (m=natom;m<TotalParticles;m++)
                                {
                                        if (Atoms[m].atomid==HydrationShell)
                                        {
                                                dx=xn-Atoms[m].x;
                                                dy=yn-Atoms[m].y;
                                                dz=zn-Atoms[m].z;
                                                dist=dx*dx+dy*dy+dz*dz;
                                                if (dist<3.0)
                                                {
                                                        if (group[n]>group[m])
                                                        {
                                                                group[m]=group[n];
                                                                changed=true;
                                                        }
                                                        else
                                                        {
                                                                group[n]=group[m];
                                                        }
                                                }
                                        }
                                }
                        }
                }
        }

        //tmax=0;
        /*************************************************************************************************************************************

          for (n=natom;n<w;n++)
          {
          if (atomid[n]==HydrationShell)
          {

          while(true)
          {
          for (t=0;t<tmax;t++)
          {
          if (group[n]==cavity[t])
          {
          cavityid[n]=t;
          break;
          }
          }

          cavity[tmax]=group[n];
          cavityid[n]=tmax;
          tmax++;
          }
          }
          }

         *****************************************************************************************************************************/

        for (n=natom;n<TotalParticles;n++)
        {
                if (Atoms[n].atomid==HydrationShell)
                {
                        if (group[n]<TotalParticles-1)
                        {
                                Atoms[n].atomid=ExcludedVolume;
                        }
                }
        }
        delete [] group;
}

void ConvertGrid(Real xmin, Real ymin, Real zmin, vector<CubeStruct> &cubes, Real CubeSize, vector< vector< vector<int> > > &CubeConversion)
{
        int m, n, CubeNum=cubes.size();
        int xbin, ybin, zbin;
        int MaxXBin, MaxYBin, MaxZBin;
        MaxXBin=CubeConversion.size();
        MaxYBin=CubeConversion[0].size();
        MaxZBin=CubeConversion[0][0].size();
        cout <<"In ConvertGrid"<<endl;
        cout <<"MaxXBin= "<<MaxXBin<<" MaxYBin= "<<MaxYBin<<" MaxZBin= "<<MaxZBin<<endl;
        for (xbin=0; xbin<MaxXBin; xbin++)
        {
                for (ybin=0; ybin<MaxYBin; ybin++)
                {
                        for (zbin=0; zbin<MaxZBin; zbin++)
                        {
                                CubeConversion[xbin][ybin][zbin]=-1;
                        }
                }
        }
        cout <<"Reset CubeConversion"<<endl;
        m=0;
        cout <<"xmin= "<<xmin<<" ymin= "<<ymin<<" zmin= "<<zmin<<endl;
        cout <<"CubeSize= "<<CubeSize<<endl;
        for (n=0;n<CubeNum;n++)
        {
                //cout <<"n= "<<n<<" CubeNum= "<<CubeNum<<endl;
                xbin=int(floor((cubes[n].x-xmin)/CubeSize+0.5));
                ybin=int(floor((cubes[n].y-ymin)/CubeSize+0.5));
                zbin=int(floor((cubes[n].z-zmin)/CubeSize+0.5));
                if (n==0)
                {
                        cout <<"xbin= "<<xbin<<" ybin= "<<ybin<<" zbin= "<<zbin<<endl;
                        cout <<"x= "<<cubes[n].x<<" y= "<<cubes[n].y<<" z= "<<cubes[n].z<<endl;
                }
                //cout <<"CubeConversion= "<<CubeConversion[xbin][ybin][zbin]<<endl;
                CubeConversion[xbin][ybin][zbin]=n;
        }
        cout <<"Done with convert grid"<<endl;
}

void AllocateCubeConversion(int MaxXBin, int MaxYBin, int MaxZBin, vector< vector< vector<int> > > &CubeConversion)
{
        int i;
        vector<int> v;
        vector< vector<int> > m;
        cout <<"MaxXBin= "<<MaxXBin<<endl;
        Create3DArray(CubeConversion, MaxXBin, MaxYBin, MaxZBin, "CubeConversion");
        cout <<"CubeConversion.size()= "<<CubeConversion.size()<<endl;
}

bool UpdateCubeGroups(vector<CubeStruct> &cubes, Real InvCubeSize, int GroupedWith[], vector< vector< vector<int> > > &CubeConversion)
{
        bool changed;
        int m, n, CubeNum=cubes.size();
        int MaxXBin, MaxYBin, MaxZBin;
        int xbin, ybin, zbin;
        int xbin2, ybin2, zbin2;
        int xstart, ystart, zstart;
        int xend, yend, zend;
        Real xmax, ymax, zmax;
        Real xmin, ymin, zmin;
        MaxXBin=CubeConversion.size();
        MaxYBin=CubeConversion[0].size();
        MaxZBin=CubeConversion[0][0].size();
        MinMax(xmin, ymin, zmin, xmax, ymax, zmax, cubes);
        changed=false;
        for (m=0;m<CubeNum;m++)
        {
                //cout <<"CubeLocation["<<m<<"]= "<<CubeLocation[m]<<endl;
                if (cubes[m].LocatedIn!=InProtein)
                {	
                        //cout <<"CubeX= "<<CubeX[m]<<" CubeY= "<<CubeY[m]<<" CubeZ= "<<CubeZ[m]<<endl;
                        xbin=int((cubes[m].x-xmin)*InvCubeSize+0.5);
                        ybin=int((cubes[m].y-ymin)*InvCubeSize+0.5);
                        zbin=int((cubes[m].z-zmin)*InvCubeSize+0.5);
                        //cout <<"xbin= "<<xbin<<" ybin= "<<ybin<<" zbin= "<<zbin<<endl;
                        //cout <<"MaxXBin= "<<MaxXBin<<" MaxYBin= "<<MaxYBin<<" MaxZBin= "<<MaxZBin<<endl;
                        xstart=xbin-1;
                        ystart=ybin-1;
                        zstart=zbin-1;

                        xend=xbin+1;
                        yend=ybin+1;
                        zend=zbin+1;

                        if (xstart<0) xstart=0;
                        if (ystart<0) ystart=0;
                        if (zstart<0) zstart=0;

                        if (xend>=MaxXBin) xend=MaxXBin-1;
                        if (yend>=MaxYBin) yend=MaxYBin-1;
                        if (zend>=MaxZBin) zend=MaxZBin-1;

                        for (xbin2=xstart;xbin2<=xend;xbin2++)
                        {
                                n=CubeConversion[xbin2][ybin][zbin];
                                if (n>=0 && cubes[n].LocatedIn!=InProtein)
                                {
                                        if (GroupedWith[m]!=GroupedWith[n]) 
                                        {
                                                changed=true;
                                                if (GroupedWith[m]<GroupedWith[n]) GroupedWith[n]=GroupedWith[m];
                                                else GroupedWith[m]=GroupedWith[n];
                                        }
                                }
                        }


                        for (ybin2=ystart;ybin2<=yend;ybin2++)
                        {
                                n=CubeConversion[xbin][ybin2][zbin];
                                if (n>=0 && cubes[n].LocatedIn!=InProtein)
                                {
                                        if (GroupedWith[m]!=GroupedWith[n]) 
                                        {
                                                changed=true;
                                                if (GroupedWith[m]<GroupedWith[n]) GroupedWith[n]=GroupedWith[m];
                                                else GroupedWith[m]=GroupedWith[n];
                                        }
                                }
                        }


                        for (zbin2=zstart;zbin2<=zend;zbin2++)
                        {
                                n=CubeConversion[xbin][ybin][zbin2];
                                if (n>=0 && cubes[n].LocatedIn!=InProtein)
                                {
                                        if (GroupedWith[m]!=GroupedWith[n]) 
                                        {
                                                changed=true;
                                                if (GroupedWith[m]<GroupedWith[n]) GroupedWith[n]=GroupedWith[m];
                                                else GroupedWith[m]=GroupedWith[n];
                                        }
                                }
                        }
                }
        }
        return changed;
}

void EliminateCavities(Real xmin, Real ymin, Real zmin, Real contrast, vector<CubeStruct> &cubes, Real CubeSize)
{
        bool changed;
        int m, n, CubeNum=cubes.size();
        int MaxXBin, MaxYBin, MaxZBin;
        int MaxPrint=1000, NumEliminated=0;
        int NumIterations;
        int *GroupedWith;
        vector< vector< vector<int> > > CubeConversion;
        Real InvCubeSize=1.0/CubeSize;
        Real xmax, ymax, zmax;
        time_t seconds1;
        time_t seconds2;


        seconds1=time(NULL);
        MinMax(xmin, ymin, zmin, xmax, ymax, zmax, cubes);
        MaxXBin=int((xmax-xmin)*InvCubeSize+0.5)+1;
        MaxYBin=int((ymax-ymin)*InvCubeSize+0.5)+1;
        MaxZBin=int((zmax-zmin)*InvCubeSize+0.5)+1;
        AllocateCubeConversion(MaxXBin, MaxYBin, MaxZBin, CubeConversion);
        ConvertGrid(xmin, ymin, zmin, cubes, CubeSize, CubeConversion);
        cout <<"In EliminateCavities"<<endl;
        cout <<"CubeNum= "<<CubeNum<<endl;
        GroupedWith=new int[CubeNum];
        cout <<"Allocated GroupedWith"<<endl;
        for (n=0;n<CubeNum;n++) GroupedWith[n]=0;
        for (n=0;n<CubeNum;n++) GroupedWith[n]=n;
        cout <<"Initialized GroupedWith"<<endl;

        NumIterations=0;
        while (true)
        {
                changed=false;
                NumIterations++;
                cout <<"NumIterations= "<<NumIterations<<endl;
                changed=UpdateCubeGroups(cubes, InvCubeSize, GroupedWith, CubeConversion);
                if (!changed) break;
        }
        cout <<"Found Groups"<<endl;
        cout <<"NumIterations= "<<NumIterations<<endl;
        for (m=0;m<CubeNum;m++)
        {
                if (GroupedWith[m]!=0 && cubes[m].LocatedIn!=InProtein)
                {
                        if (NumEliminated<MaxPrint)
                        {
                                cout <<"GroupedWith["<<m<<"]= "<<GroupedWith[m]<<endl;
                        }
                        NumEliminated++;
                        cubes[m].density=-contrast;
                        cubes[m].LocatedIn=InProtein;
                }
        }
        delete [] GroupedWith;
        seconds2=time(NULL);
        cout <<"NumEliminated= "<<NumEliminated<<endl;
        cout <<"EliminatedCavities"<<endl;
        cout <<"Eliminate Cavities took "<<seconds2-seconds1<<" seconds."<<endl;
}

Real InterpolateDensity(Matrix &DensityMatrix, int SolventType, Real mindist, Real RecBin)
{
        int distance=int(mindist/RecBin);
        int MaxDistance=DensityMatrix[0].size();
        Real dr, slope, density;

        if (DensityMatrix[SolventType][distance]!=UNK_DENSITY)
        {
                if (distance+1>=MaxDistance || DensityMatrix[SolventType][distance+1]==UNK_DENSITY)
                {
                        return DensityMatrix[SolventType][distance];
                }
                else 
                {
                        dr=mindist-Real(distance)*RecBin;
                        slope=(DensityMatrix[SolventType][distance+1]-DensityMatrix[SolventType][distance])/RecBin;
                        density=DensityMatrix[SolventType][distance]+dr*slope;
                        if (distance==30)
                        {
                        //cout <<"mindist= "<<mindist<<" RecBin= "<<RecBin<<endl;
                        //cout <<"DensityMatrix["<<SolventType<<"]["<<distance<<"]= "<<DensityMatrix[SolventType][distance]<<endl;
                        //cout <<"DensityMatrix["<<SolventType<<"]["<<distance+1<<"]= "<<DensityMatrix[SolventType][distance+1]<<endl;
                        //cout <<"dr= "<<dr<<" slope= "<<slope<<"density= "<<density<<endl;
                        //cout <<endl;
                        }
                        return density;
                }
        }
        return UNK_DENSITY;
}

inline bool GetDensity(Real &density, AtomStruct Atom, int distance, Matrix &AtomTypesDensityMatrix, Matrix &ElementsDensityMatrix, Real BulkDensity, Real RecBin)
{
        int SolventType;
        SolventType=Atom.AtomType;
        density=AtomTypesDensityMatrix[SolventType][distance];
        if (density==UNK_DENSITY) 
        {
                SolventType=Atom.AtomType2;
                cout <<"Warning unable to assign cube density from AtomTypes pRDF. SolventType="<< Atom.AtomType<<" distance= "<<distance<<" ResidueName= "<<Atom.ResidueName<<" AtomName= "<<Atom.AtomName<<endl;
                cout <<"AtomType2= "<<Atom.AtomType2<<endl;
                //cout <<"ElementsDensityMatrix= "<<ElementsDensityMatrix<<endl;
                density=ElementsDensityMatrix[SolventType][distance];
        }
        if (density==UNK_DENSITY)
        {
                cout <<"Warning unable to assign cube density.  Setting cube density to bulk solvent"<<endl;
                PrintAtomInfo(Atom);
                cout <<"distance= "<<distance<<endl;
                density=BulkDensity;
                return false;
        }
        return true;
}

inline bool GetDensity(Real &density, CubeStruct &Cube, Matrix &AtomTypesDensityMatrix, Matrix &ElementsDensityMatrix, Real BulkDensity, bool PhobicPhilic, Real RecBin)
{
        int SolventType;
        SolventType=Cube.AtomType;
        density=AtomTypesDensityMatrix[SolventType][Cube.IntDist];
        //density=InterpolateDensity(AtomTypesDensityMatrix, SolventType, Cube.MinDist, RecBin);
        cout <<"In GetDensity. density= "<<density<<" SolventType= "<<SolventType<<" Cube.IntDist= "<<Cube.IntDist<<endl;
        if (density==UNK_DENSITY && PhobicPhilic)
        {
                if (isEven(SolventType)) SolventType++;
                else SolventType--;
                density=AtomTypesDensityMatrix[SolventType][Cube.IntDist];
        }
        if (density==UNK_DENSITY) 
        {
                SolventType=Cube.ElementType;
                cout <<"Warning unable to assign cube density from AtomTypes pRDF. SolventType="<<Cube.AtomType<<" distance= "<<Cube.IntDist<<endl;
                PrintCubeInfo(Cube);
                //cout <<"ElementsDensityMatrix= "<<ElementsDensityMatrix<<endl;
                density=ElementsDensityMatrix[SolventType][Cube.IntDist];
                //density=InterpolateDensity(ElementsDensityMatrix, SolventType, Cube.MinDist, RecBin);
        }
        if (density==UNK_DENSITY && PhobicPhilic)
        {
                if (isEven(SolventType)) SolventType++;
                else SolventType--;
                density=AtomTypesDensityMatrix[SolventType][Cube.IntDist];
        }
        if (density==UNK_DENSITY)
        {
                cout <<"Warning unable to assign cube density.  Setting cube density to bulk solvent"<<endl;
                //PrintAtomInfo(Atom);
                cout <<"distance= "<<Cube.IntDist<<endl;
                density=BulkDensity;
                return false;
        }
        return true;
}

void RemoveHydrationShell(vector<AtomStruct> &Atoms)
{
        int TotalParticles=Atoms.size();
        vector<AtomStruct> TempAtoms;

        RemoveWaters(Atoms);

        for (int i=0;i<TotalParticles;i++)
        {
                if (Atoms[i].atomid!=HydrationShell && Atoms[i].atomid!=ExcludedVolume)
                {
                        SafePushBack(TempAtoms, Atoms[i], "TempAtoms in RemoveHydrationShell");
                }
        }
        Atoms=TempAtoms;
}

int solvent(Real atomr[], vector<AtomStruct> &Atoms, Real contrast, ParamStruct &params)
{
        AtomStruct Atom;
        int TotalParticles=Atoms.size();
        int distance;
        int LocatedIn, nearest;
        int natom;
        Real dx, dy, dz;
        Real inc;
        Real HydrationAtomSpacing;
        Real dist;
        Real xgrid, ygrid, zgrid;
        Real xmax, ymax, zmax;
        Real xmin, ymin, zmin;
        Real xstart[5], ystart[5], zstart[5];
        vector< vector<Real> > AtomTypesDensityMatrix, ElementsDensityMatrix;
        int temp=Atoms.size();

        Atom.ResidueName="HOH";
        Atom.AtomName="O";

        natom=FindNumProteinAtoms(Atoms);
        if (!params.UniformHydrationShell)
        {
                cout <<"About to readgofrfile"<<endl;
                ReadGofRFile3(params.AtomTypesGofRFile, contrast, AtomTypesDensityMatrix, params.RecBin);
                ReadGofRFile3(params.ElementsGofRFile, contrast, ElementsDensityMatrix, params.RecBin);
                cout <<"Finished reading gofrfile"<<endl;
                params.maxhs=8.0;
        }
        else params.maxhs=3.0;

        for (int n=temp-1;0<=n;n--)
        {
                if (Atoms[n].atomid==HydrationShell || Atoms[n].atomid==ExcludedVolume)
                {
                        Atoms.pop_back();
                }
        }


        for (int t=0;t<5;t++)
        {
                xstart[t]=0;
                ystart[t]=0;
                zstart[t]=0;
        }
        MinMax(xmin, ymin, zmin, xmax, ymax, zmax, Atoms);
        PrintAtomInfo(Atoms[0]);
        cout <<"xmin= "<<xmin<<" ymin= "<<ymin<<" zmin= "<<zmin<<" xmax= "<<xmax<<" ymax= "<<ymax<<" zmax= "<<zmax<<endl;
        HydrationAtomSpacing=atomr[HydrationShell]/packing;
        if (params.ExcludedVolumeSphereType=="GaussianSphere") inc/=GaussianToHardSphere;
        inc=atomr[HydrationShell]*2.0*sqrt(2.0)/packing;
        xstart[0]=0;
        ystart[0]=0;
        zstart[0]=0;

        xstart[1]=0;
        ystart[1]=-inc*0.5;
        zstart[1]=-inc*0.5;

        xstart[2]=-inc*0.5;
        ystart[2]=0;
        zstart[2]=-inc*0.5;

        xstart[3]=-inc*0.5;
        ystart[3]=-inc*0.5;
        zstart[3]=0;

        TotalParticles=natom;

        cout <<"atomr[HydrationShell]= "<<atomr[HydrationShell]<<" inc= "<<inc<<endl;
        for (int t=0;t<4;t++)
        {
                for (xgrid=xmin+xstart[t]-4.0;xgrid<xmax+4.0+inc;xgrid+=inc)
                {
                        cout <<"xgrid= "<<xgrid<<endl;
                        for (ygrid=ymin+ystart[t]-4.0;ygrid<ymax+4.0+inc;ygrid+=inc)
                        {
                                //cout <<"ygrid= "<<ygrid<<endl;
                                for (zgrid=zmin+zstart[t]-4.0;zgrid<zmax+4.0+inc;zgrid+=inc)
                                {
                                        //cout <<"zgrid= "<<zgrid<<endl;
                                        //cout <<"zmin= "<<zmin<<" zmax= "<<zmax<<endl;
                                        //cout <<"About to enter location"<<endl;
                                        LocatedIn=location(xgrid, ygrid, zgrid, nearest, dist, Atoms, params.maxhs);
                                        //cout <<"after LocatedIn"<<endl;
                                        //cout <<"zgrid= "<<zgrid<<endl;
                                        if (LocatedIn!=InProtein && LocatedIn!=InHydrationShell && LocatedIn!=InSolution) cout <<"Error location not identified"<<endl;
                                        if (LocatedIn==InHydrationShell)
                                        {
                                                //cout <<"value=InHydrationShell"<<endl;
                                        }
                                        //cout <<"xgrid= "<<xgrid<<" ygrid= "<<ygrid<<" zgrid= "<<zgrid<<" nearest= "<<nearest<<endl;
                                        //cout <<"Left location"<<endl;
                                        if (params.excluded=="Lattice")
                                        {
                                                if (LocatedIn==InProtein)
                                                {
                                                        Atom.x=xgrid;
                                                        Atom.y=ygrid;
                                                        Atom.z=zgrid;
                                                        Atom.atomid=ExcludedVolume;
                                                        Atom.weight=1.0;
                                                        SafePushBack(Atoms, Atom, "Atoms in solvent");
                                                }
                                        }
                                        //cout <<"After lattice zgrid= "<<zgrid<<endl;
                                        if (LocatedIn==InHydrationShell)
                                        {
                                                Atom.x=xgrid;
                                                Atom.y=ygrid;
                                                Atom.z=zgrid;
                                                if (!params.UniformHydrationShell)
                                                {
                                                        //cout <<"In if(!Uniform)"<<endl;
                                                        distance=int(floor(dist/params.RecBin+0.5));
                                                        //cout <<"nearest= "<<nearest<<" distance= "<<distance<<endl;
                                                        //cout <<"solventtype= "<<solventtype[nearest]<<endl;
                                                        //cout <<"densitym= "<<densitym[solventtype[nearest]][distance]<<endl;
                                                        GetDensity(Atom.weight, Atoms[nearest], distance, AtomTypesDensityMatrix, ElementsDensityMatrix, 0, params.RecBin);

                                                        //if (weight[TotalParticles]<-0.01 || weight[TotalParticles]>0.01)
                                                        //{
                                                        Atom.atomid=HydrationShell;
                                                        SafePushBack(Atoms, Atom, "Atoms in solvent");
                                                        //}
                                                }
                                                else
                                                {
                                                        dx=Atoms[nearest].x-xgrid;
                                                        dy=Atoms[nearest].y-ygrid;
                                                        dz=Atoms[nearest].z-zgrid;
                                                        dist=sqrt(dx*dx+dy*dy+dz*dz);
                                                        Atom.weight=1.0;
                                                        if (dist > atomr[Atoms[nearest].atomid] && dist < Atoms[nearest].vdw)
                                                        {
                                                                Atom.atomid=ExcludedVolume;
                                                                //cout <<"Made excluded volume dummy atom."<<endl;
                                                        }
                                                        else
                                                        {
                                                                Atom.atomid=HydrationShell;
                                                        }
                                                        SafePushBack(Atoms, Atom, "Atoms in solvent");
                                                }
                                        }
                                        //cout <<"After LocatedIn==InHydrationShell zgrid= "<<zgrid<<endl;
                                }
                        }
                }
        }

        EliminateCavities(Atoms, HydrationAtomSpacing);
        cout <<"HydrationAtomSpacing= "<<HydrationAtomSpacing<<endl;
        //cout <<"TotalParticles-natom="<<TotalParticles-natom<<endl;

        cout <<"About to enter PrintPDB"<<endl;

        //PrintPDB(TotalParticles, atomid);

        return TotalParticles;
}

void AssignCubeDensities(vector<CubeStruct> &Cubes, vector<AtomStruct> Atoms, Matrix &AtomTypesDensityMatrix, Matrix &ElementsDensityMatrix, Real BulkDensity, bool UniformHydrationShellDensity, Real HydrationShellDensity, bool PhobicPhilic, Real RecBin)
{
        bool SetDensity=true;
        int n, GridPoints=Cubes.size();
        Real density;
        vector<CubeStruct> TempCubes;
        clock_t start, end, total=0;
        cout <<"In AssignCubeDensities"<<endl;
        CountNumCubesInHydrationShell(Cubes);
        for (n=0;n<GridPoints;n++)
        {
                //if (n%100==0) cout <<"n= "<<n<<endl;

                if (Cubes[n].LocatedIn==InHydrationShell)
                {
                        Cubes[n].IntDist=int(Cubes[n].MinDist/RecBin+0.5);
                        start=clock();
                        if (UniformHydrationShellDensity) density=HydrationShellDensity;
                        else SetDensity=GetDensity(density, Cubes[n], AtomTypesDensityMatrix, ElementsDensityMatrix, BulkDensity, PhobicPhilic, RecBin);
                        if (!SetDensity)
                        {
                                PrintCubeInfo(Cubes[n]);
                                cout <<"AtomTypesDensityMatrix["<<Atoms[Cubes[n].nearest].AtomType<<"]["<<Cubes[n].IntDist<<"]= "<<AtomTypesDensityMatrix[Atoms[Cubes[n].nearest].AtomType][Cubes[n].IntDist]<<endl;
                        }
                        //Cubes[n].AtomType=Atoms[Cubes[n].nearest].AtomType;
                        end=clock();
                        total+=(end-start);
                        Cubes[n].density=density;
                        SafePushBack(TempCubes, Cubes[n], "TempCubes in AssignCubeDensities");
                }

                if (Cubes[n].LocatedIn==InProtein)
                {
                        Cubes[n].density=0;
                        SafePushBack(TempCubes, Cubes[n], "TempCubes in AssignCubeDensities");
                }
        }
        Cubes=TempCubes;
        cout <<"Time spent in GetDensity was "<<Real(total)/Real(CLOCKS_PER_SEC)<<" seconds."<<endl;
}

void SetCubeType(CubeStruct &Cube, vector<AtomStruct> &Atoms, vector<int> Hydrophobic, ParamStruct &params)
{
        int nearest=Cube.nearest;
        int SecondNearest;

        Cube.AllAtomType=nearest;
        Cube.AtomType=Atoms[nearest].AtomType;
        Cube.ElementType=Atoms[nearest].AtomType2;

        if (params.PhobicPhilic)
        {
                SecondNearest=Cube.SecondNearest;
                Cube.AllAtomType=Cube.AllAtomType*2+Hydrophobic[Atoms[SecondNearest].AtomType2];
                Cube.AtomType=Cube.AtomType*2+Hydrophobic[Atoms[SecondNearest].AtomType2];
                Cube.ElementType=Cube.ElementType*2+Hydrophobic[Atoms[SecondNearest].AtomType2];
        }
        
        if (params.ConvexConcave)
        {
                if (Cube.concave)
                {
                        Cube.AllAtomType=Cube.AllAtomType*2;
                        Cube.AtomType=Cube.AtomType*2;
                        Cube.ElementType=Cube.ElementType*2;
                }
                else
                {
                        Cube.AllAtomType=Cube.AllAtomType*2+1;
                        Cube.AtomType=Cube.AtomType*2+1;
                        Cube.ElementType=Cube.ElementType*2+1;
                }
        }

        if (params.AngularDependence)
        {
                if (Cube.angular)
                {
                        Cube.AllAtomType=Cube.AllAtomType*2;
                        Cube.AtomType=Cube.AtomType*2;
                        Cube.ElementType=Cube.ElementType*2;
                }
                else
                {
                        Cube.AllAtomType=Cube.AllAtomType*2+1;
                        Cube.AtomType=Cube.AtomType*2+1;
                        Cube.ElementType=Cube.ElementType*2+1;
                }
        }
}

void SetHydrophobic(vector<int> &Hydrophobic)
{
        Hydrophobic.clear();
        SafeAlloc(Hydrophobic, NumAtomTypes2, "Hydrophobic");
        Hydrophobic[HC]=HYDROPHOBIC;
        Hydrophobic[HN]=HYDROPHILIC;
        Hydrophobic[HO]=HYDROPHILIC;
        Hydrophobic[HS]=HYDROPHILIC;
        Hydrophobic[C]=HYDROPHOBIC;
        Hydrophobic[N]=HYDROPHOBIC;
        Hydrophobic[O]=HYDROPHILIC;
        Hydrophobic[S]=HYDROPHILIC;
}

Real VectorAngle(Real dx1, Real dy1, Real dz1, Real dx2, Real dy2, Real dz2)
{
        Real dotproduct, length1, length2, theta;
        dotproduct=dx1*dx2+dy1*dy2+dz1*dz2;
        length1=sqrt(dx1*dx1+dy1*dy1+dz1*dz1);
        length2=sqrt(dx2*dx2+dy2*dy2+dz2*dz2);
        theta=acos(dotproduct/(length1*length2));
        return theta;
}

bool Angular(vector<AtomStruct> &Atoms, CubeStruct &Cube, ParamStruct params)
{
        bool BoolPrint=false;
        Real d=0.1;
        int nearest, bonded;
        int xbin, ybin, zbin;
        Real theta, cutoff;
        Real dx1, dy1, dz1;
        Real dx2, dy2, dz2;
        cutoff=(params.AngularCutOff/360.0)*2.0*pi;

        if (abs(Cube.x+27.447)<d && abs(Cube.y+4.792)<d && abs(Cube.z+9.668)<d) BoolPrint=true;

        nearest=Cube.nearest;
        if (BoolPrint) 
        {
                cout <<"About to check bond"<<endl;
                PrintAtomInfo(Atoms[nearest]);
                cout <<"nearest= "<<nearest<<endl;
                cout <<"Atoms.connectivity.size()="<<Atoms[nearest].connectivity.size()<<endl;
        }
        bonded=Atoms[nearest].connectivity[0];
        if (BoolPrint) 
        {
                PrintAtomInfo(Atoms[bonded]);
                PrintCubeInfo(Cube);
        }
        dx1=Atoms[nearest].x-Atoms[bonded].x;
        dy1=Atoms[nearest].y-Atoms[bonded].y;
        dz1=Atoms[nearest].z-Atoms[bonded].z;
        dx2=Cube.x-Atoms[nearest].x;
        dy2=Cube.y-Atoms[nearest].y;
        dz2=Cube.z-Atoms[nearest].z;
        if (BoolPrint)
        {
                cout <<"dx1= "<<dx1<<" dy1= "<<dy1<<" dz1= "<<dz1<<endl;
                cout <<"dx2= "<<dx2<<" dy2= "<<dy2<<" dz2= "<<dz2<<endl;
        }
        theta=VectorAngle(dx1, dy1, dz1, dx2, dy2, dz2);
        if (BoolPrint) cout <<"theta= "<<theta<<" cutoff= "<<cutoff<<endl;
        if (theta<cutoff) return true;
        else return false;
}

void Angular(vector<AtomStruct> &Atoms, vector<CubeStruct> &Cubes, ParamStruct params)
{
        int ncube=Cubes.size();

        for (int i=0;i<ncube;i++) Cubes[i].angular=Angular(Atoms, Cubes[i], params);
}

void SetCubeTypes(vector<CubeStruct> &Cubes, vector<AtomStruct> &Atoms, ParamStruct &params)
{
        int ncube=Cubes.size();
        vector<int> Hydrophobic;
        if (params.ConvexConcave) Cavity(Cubes, params);
        if (params.PhobicPhilic) SetHydrophobic(Hydrophobic);
        if (params.AngularDependence) Angular(Atoms, Cubes, params);
        for (int i=0;i<ncube;i++) SetCubeType(Cubes[i], Atoms, Hydrophobic, params);
}

void BuildHydrationShell(vector<CubeStruct> &Cubes, vector<AtomStruct> &Atoms, Matrix &AtomTypesDensityMatrix, Matrix &ElementsDensityMatrix, Real BulkDensity, bool UniformHydrationShell, Real hsdensity, ParamStruct params)
{
        bool MadeGrid;
        DeleteVector(Cubes);
        int natom=Atoms.size();
        Real xmax, ymax, zmax;
        Real xmin, ymin, zmin;
        time_t seconds1, seconds2;
        MinMaxProtein(xmin, ymin, zmin, xmax, ymax, zmax, Atoms);
        cout <<"xmin= "<<xmin<<" ymin= "<<ymin<<" zmin= "<<zmin<<endl;
        ExpandMinMax(xmin, ymin, zmin, xmax, ymax, zmax, params);
        MadeGrid=MakeGrid(xmin, ymin, zmin, xmax, ymax, zmax, params.CubeSize, Cubes);	
        if (!MadeGrid) SafeMakeGrid(xmin, ymin, zmin, xmax, ymax, zmax, Atoms, params.maxhs, params.CubeSize, Cubes);
        cout <<"xmin= "<<xmin<<" ymin= "<<ymin<<" zmin= "<<zmin<<endl;
        seconds1=time(NULL);
        location4(natom, xmin, xmax, ymin, ymax, zmin, zmax, Cubes, Atoms, params.maxhs);
        CountNumCubesInHydrationShell(Cubes);
        seconds2=time(NULL);
        cout <<"location4 took "<<seconds2-seconds1<<" seconds."<<endl;
        seconds1=time(NULL);
        SetCubeTypes(Cubes, Atoms, params);
        AssignCubeDensities(Cubes, Atoms, AtomTypesDensityMatrix, ElementsDensityMatrix, BulkDensity, UniformHydrationShell, hsdensity, params.PhobicPhilic, params.RecBin);
        seconds2=time(NULL);
        cout <<"AssignCubeDensities took "<<seconds2-seconds1<<" seconds."<<endl;
        //if (!params.EliminateCavities) EliminateCavities(xmin, ymin, zmin, 0, Cubes, CubeSize);
}
#endif
