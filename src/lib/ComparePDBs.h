#ifndef _ComparePDBs_included_
#define _ComparePDBs_included_
 
Real ComparePDBs(Real contrast, vector<CubeStruct> &Cubes1, vector<CubeStruct> &Cubes2, vector<AtomStruct> &Atoms, Real bin, bool PrintDifferencePdb=false);
Real ComparePDBs(Real contrast, vector<AtomStruct> &Atoms1, vector<AtomStruct> &Atoms2, vector<AtomStruct> &DiffAtoms, Real CubeSize, bool PrintDifferencPdb=false);

# include "Structures.h"
# include "ConvertStructures.h"
# include "StringUtils.h"
# include "LinkedList.h"
# include "ReadPdb.h"
# include "TypeDef.h"
# include "VectorManip.h"
# include "WritePdb.h"

int FindMaxAtomType(vector<CubeStruct> &Cubes1, vector<CubeStruct> &Cubes2)
{
        int ncube1=Cubes1.size(), ncube2=Cubes2.size();
        int MaxAtomType=0;
        
        for (int i=0;i<ncube1;i++)
        {
                if (Cubes1[i].AtomType>MaxAtomType) MaxAtomType=Cubes1[i].AtomType;
        }
        for (int i=0;i<ncube2;i++)
        {
                if (Cubes2[i].AtomType>MaxAtomType) MaxAtomType=Cubes2[i].AtomType;
        }
        return MaxAtomType+1;
}

Real ComparePDBs(Real contrast, vector<CubeStruct> &Cubes1, vector<CubeStruct> &Cubes2, vector<AtomStruct> &Atoms, Real bin, bool PrintDifferencePdb)
{
        bool verbose=false;
	string DistFile, TypeFile;
	int m, n;
	int CubeNum, CubeNum2;
	vector< vector< vector<bool> > > Exists;
	int Xbins, Ybins, Zbins;
	int xgrid, ygrid, zgrid;
	int xtest, ytest, ztest;
        int MaxAtomType=0;
	int WorstCube, WorstMDCube;
	Real SumDiff, Sum2;
	Real MinX, MinY, MinZ;
	Real MaxX, MaxY, MaxZ;
	vector< vector< vector<Real> > > CubeDensity;
	Real CubeDensityXYZ;
	Real dWeight, dWeight2;
	Real MinCubeWeight, MinCubeWeight2;
        Real DiffSolvent, DiffProtein;
	Real SqrDiffSolvent, SqrDiffProtein;
	Real SqrDiff, NonZeroGridPoints;
        Real NonZeroPoints1, NonZeroPoints2;
	Real SumWeight, Sum;
	Real RFactor2, MisidentifiedPoint;
	Real SolventMisidentifiedAsProtein;
	Real ProteinMisidentifiedAsSolvent;
	Real *SumWeightDist, *SumDiffDist;
	Real *SumDist, *Sum2Dist;
	Real *SumWeightType, *SumDiffType;
	Real *SumType, *Sum2Type;
	AtomStruct Atom;
        Atoms.clear();
	Atom.AtomNumber=1;
	Atom.ResidueNum=1;
	Atom.weight=1.0;
	Atom.Occupancy=1.0;
	Atom.AtomName="O";
	Atom.ResidueName="HOH";
	Atom.ChainName="A";
        cout <<"In the real ComparePDBs"<<endl;
        MaxAtomType=FindMaxAtomType(Cubes1, Cubes2);
        cout <<"About to allocate memory"<<endl;
        cout <<"MaxAtomType= "<<MaxAtomType<<endl;
	SumWeightDist=new Real[1000];
	SumDiffDist=new Real[1000];
	SumDist=new Real[1000];
	Sum2Dist=new Real[1000];
        cout <<"MaxAtomType= "<<MaxAtomType<<endl;
        SafeArrayAlloc(SumWeightType, MaxAtomType, "SumWeightType");
        SafeArrayAlloc(SumDiffType, MaxAtomType, "SumDiffType");
        SafeArrayAlloc(SumType, MaxAtomType, "SumType");
        SafeArrayAlloc(Sum2Type, MaxAtomType, "Sum2Type");
	for (n=0;n<1000;n++)
	{
		SumWeightDist[n]=0;
		SumDiffDist[n]=0;
		SumDist[n]=0;
		Sum2Dist[n]=0;
	}
	MinX=Cubes1[0].x;
	MinY=Cubes1[0].y;
	MinZ=Cubes1[0].z;
	MaxX=Cubes1[0].x;
	MaxY=Cubes1[0].y;
	MaxZ=Cubes1[0].z;
	MisidentifiedPoint=0.0;
	SolventMisidentifiedAsProtein=0.0;
	ProteinMisidentifiedAsSolvent=0.0;
	SqrDiffSolvent=0;
	SqrDiffProtein=0;
	cout <<"In ComparePDBs"<<endl;
	
        CubeNum=Cubes1.size();
	CubeNum2=Cubes2.size();
	cout <<"CubeNum= "<<CubeNum<<" CubeNum2= "<<CubeNum2<<endl;
        cout <<"Cubes1.capacity()= "<<Cubes1.capacity()<<" Cubes2.capacity()= "<<Cubes2.capacity()<<endl;
	for (n=0;n<CubeNum;n++)
	{
		if (Cubes1[n].x<MinX) MinX=Cubes1[n].x;
		if (Cubes1[n].y<MinY) MinY=Cubes1[n].y;
		if (Cubes1[n].z<MinZ) MinZ=Cubes1[n].z;
		if (Cubes1[n].x>MaxX) MaxX=Cubes1[n].x;
		if (Cubes1[n].y>MaxY) MaxY=Cubes1[n].y;
		if (Cubes1[n].z>MaxZ) MaxZ=Cubes1[n].z;
	}
	cout <<"MinX= "<<MinX<<endl;
	for (n=0;n<CubeNum2;n++)
	{
		if (Cubes2[n].x<MinX) MinX=Cubes2[n].x;
		if (Cubes2[n].y<MinY) MinY=Cubes2[n].y;
		if (Cubes2[n].z<MinZ) MinZ=Cubes2[n].z;
		if (Cubes2[n].x>MaxX) MaxX=Cubes2[n].x;
		if (Cubes2[n].y>MaxY) MaxY=Cubes2[n].y;
		if (Cubes2[n].z>MaxZ) MaxZ=Cubes2[n].z;
	}
	cout <<"MinX= "<<MinX<<endl;
	Xbins=int(floor((MaxX-MinX)/bin+0.5))+1;
	Ybins=int(floor((MaxY-MinY)/bin+0.5))+1;
	Zbins=int(floor((MaxZ-MinZ)/bin+0.5))+1;
	cout <<"Xbins= "<<Xbins<<endl;
	Safe3DAlloc(CubeDensity, Xbins, Ybins, Zbins, "CubeDensity");
	Safe3DAlloc(Exists, Xbins, Ybins, Zbins, "Exists");

	MinCubeWeight=Cubes1[0].density;
	MinCubeWeight2=Cubes2[0].density;
	for (n=0;n<CubeNum;n++)
	{
		if (Cubes1[n].density<MinCubeWeight) MinCubeWeight=Cubes1[n].density;
	}

	for (n=0;n<CubeNum2;n++)
	{
		if (Cubes2[n].density<MinCubeWeight2) MinCubeWeight2=Cubes2[n].density;
	}
        cout <<"MinCubeWeight= "<<MinCubeWeight<<" MinCubeWeight2= "<<MinCubeWeight2<<endl;
	for (n=0;n<CubeNum;n++)
	{
		//Cubes1[n].density-=MinCubeWeight;
	}

	for (n=0;n<CubeNum2;n++)
	{
		//Cubes2[n].density-=MinCubeWeight2;
	}
	cout <<"Initialized Exist"<<endl;
	xtest=int((28.098-MinX)/bin+0.5);
	ytest=int((0.619-MinY)/bin+0.5);
	ztest=int((10.075-MinZ)/bin+0.5);
        int TestCube;
	for (n=0;n<CubeNum2;n++)
	{
		xgrid=int((Cubes2[n].x-MinX)/bin+0.5);
		ygrid=int((Cubes2[n].y-MinY)/bin+0.5);
		zgrid=int((Cubes2[n].z-MinZ)/bin+0.5);
		CubeDensity[xgrid][ygrid][zgrid]=Cubes2[n].density;
		//cout <<"CubeDensity["<<xgrid<<"]["<<ygrid<<"]["<<zgrid<<"]= "<<CubeDensity[xgrid][ygrid][zgrid]<<endl;
	        if (xgrid==xtest && ygrid==ytest && zgrid==ztest)
                {
                        TestCube=n;
                        cout <<"TestCube= "<<TestCube<<endl;
                        PrintCubeInfo(Cubes2[n]);
                }
                
                Exists[xgrid][ygrid][zgrid]=true;
	}
	//cout <<"CubeDensity["<<xtest<<"]["<<ytest<<"]["<<ztest<<"]= "<<CubeDensity[xtest][ytest][ztest]<<endl;
	cout <<"xgrid_-9.870= "<<int((-9.870-MinX)/bin+0.5)<<endl;
	cout <<"xgrid_-8.001= "<<int((-8.001-MinY)/bin+0.5)<<endl;
	cout <<"xgrid_-6.631= "<<int((-6.631-MinZ)/bin+0.5)<<endl;
	//cout <<"xgrid_16.631= "<<int(floor((16.631-MinX+0.5)/bin));
	//cout <<"xgrid_16.630= "<<int(floor((16.630-MinX+0.5)/bin));
	SumDiff=0;
	Sum2=0;
	cout <<"Assigned Exists"<<endl;
	n=0;
	//cout <<"Exists["<<xtest<<"]["<<ztest<<"]["<<ytest<<"]= "<<Exists[xtest][ytest][ztest]<<endl;

	SqrDiff=0;
	SolventMisidentifiedAsProtein=0;
	SqrDiffSolvent=0;
	ProteinMisidentifiedAsSolvent=0;
	NonZeroGridPoints=0;
        NonZeroPoints1=0;
        NonZeroPoints2=0;
	SqrDiffProtein=0;
        DiffProtein=0;
        DiffSolvent=0;
	SumWeight=0;
	SumDiff=0;
	Sum=0;
	Sum2=0;
        WorstCube=0;
        WorstMDCube=0;
	Real dWeight2_1=0;
	Real NumDiffGreaterThan_1=0;
	Real GreatestDifference=0;
        bool BoolPrint=false;
	for (m=0;m<CubeNum;m++)
	{
                if (m==282496) BoolPrint=true;
                //cout <<"m= "<<m<<endl;
		xgrid=int((Cubes1[m].x-MinX)/bin+0.5);
		ygrid=int((Cubes1[m].y-MinY)/bin+0.5);
		zgrid=int((Cubes1[m].z-MinZ)/bin+0.5);
		if (Exists[xgrid][ygrid][zgrid])
		{
                        if (BoolPrint) cout <<"At beginning"<<endl;
			CubeDensityXYZ=CubeDensity[xgrid][ygrid][zgrid];
                        if (BoolPrint) cout <<"n= "<<n<<" xgrid= "<<xgrid<<" ygrid= "<<ygrid<<" zgrid= "<<zgrid<<endl;
	                if (xgrid==xtest && ygrid==ytest && zgrid==ztest)
                        {
                                cout <<"TestCube= "<<TestCube<<endl;
                                PrintCubeInfo(Cubes2[TestCube]);
                                cout <<endl;
                                PrintCubeInfo(Cubes1[m]);
                        }
			/*
			if (Cubes1[m].AtomType!=Cubes2[n].AtomType || Cubes1[m].IntDist!=Cubes2[n].IntDist)
			{
				cout <<"Warning. AtomType1["<<m<<"]= "<<Cubes1[m].AtomType<<" AtomType2["<<n<<"]= "<<Cubes2[n].AtomType<<endl;
				cout <<"Warning. IntDist1["<<m<<"]= "<<Cubes1[m].IntDist<<" IntDist2["<<n<<"]= "<<Cubes2[n].IntDist<<endl;
				cout <<"x= "<<Cubes1[m].x<<" y= "<<Cubes1[m].y<<" z= "<<Cubes1[m].z<<endl;
				cout <<endl;
			}
			*/
			dWeight=Cubes1[m].density-CubeDensityXYZ;
                        if (m%100==0) 
                        {
                                //cout <<"m= "<<m<<" dWeight= "<<dWeight<<endl;
                        }
                        if (BoolPrint) cout <<"dWeight= "<<dWeight<<endl;
			//if (m%100==0)
			//{
			//	cout <<"dWeight= "<<dWeight<<" Cubes1["<<m<<"].density= "<<Cubes1[m].density<<" CubeDensityXYZ= "<<CubeDensityXYZ<<" Cubes2["<<n<<"].density= "<<Cubes2[n].density<<endl;
			//}
			dWeight2=dWeight*dWeight;
			SqrDiff+=dWeight2;
			if (GreatestDifference<abs(dWeight))
			{
				GreatestDifference=abs(dWeight);
				WorstCube=m;
				//WorstMDCube=n;
			}
			if (dWeight2>1.0)
			{
				//PrintCubeInfo(Cubes1[m]);
				//cout <<endl;
				//PrintCubeInfo(Cubes2[n]);
				//cout <<endl<<endl;
				dWeight2_1+=dWeight2;
				NumDiffGreaterThan_1+=1.0;
			}
			if (Cubes1[m].density!=0 && CubeDensityXYZ==0)
			{
				SolventMisidentifiedAsProtein+=1.0;
                                if (verbose)
                                {
                                        cout <<"Solvent misidentified as protein"<<endl;
                                        PrintCubeInfo(Cubes1[m]);
                                        cout <<endl;
                                }
                                DiffSolvent+=abs(dWeight);
				SqrDiffSolvent+=dWeight2;
			}
	                if (BoolPrint) cout <<"SqrDiff= "<<SqrDiff<<endl;
			if (Cubes1[m].density==0 && CubeDensityXYZ!=0)
			{
				ProteinMisidentifiedAsSolvent+=1.0;
				SqrDiffProtein+=dWeight2;
                                DiffProtein+=abs(dWeight);
                                cout <<"ProteinMisidentifiedAsSolvent"<<endl;
                                cout <<"m= "<<m<<endl;
                                PrintCubeInfo(Cubes1[m]);
                                cout <<endl;
			}
                        if (BoolPrint) cout <<"if Cubes1[m]"<<endl;
			if (Cubes1[m].density>0 || CubeDensityXYZ>0)
			{
				SumWeight+=abs(CubeDensityXYZ);
				SumDiff+=abs(dWeight);
				Sum+=abs(Cubes1[m].density+CubeDensityXYZ);
				Sum2+=abs(Cubes1[m].density+CubeDensityXYZ-0.33484*2.0);
				NonZeroGridPoints+=1.0;
				/*
                                SumWeightDist[Cubes1[m].AtomType]+=CubeDensityXYZ;
				SumDiffDist[Cubes1[m].AtomType]+=abs(dWeight);
				SumDist[Cubes1[m].AtomType]+=abs(Cubes1[m].density+CubeDensityXYZ);
				Sum2Dist[Cubes1[m].AtomType]+=abs(Cubes1[m].density+CubeDensityXYZ-0.33484*2.0);
				SumWeightType[Cubes1[m].IntDist]+=CubeDensityXYZ;
				SumDiffType[Cubes1[m].IntDist]+=abs(dWeight);
				SumType[Cubes1[m].IntDist]+=abs(Cubes1[m].density+CubeDensityXYZ);
				Sum2Type[Cubes1[m].IntDist]+=abs(Cubes1[m].density+CubeDensityXYZ-0.33484*2.0);
			        */
                        }
                        if (BoolPrint) cout <<"Near end"<<endl;
			Atom.x=Cubes1[m].x;
			Atom.y=Cubes1[m].y;
			Atom.z=Cubes1[m].z;
			Atom.weight=dWeight;
                        if (PrintDifferencePdb) SafePushBack(Atoms, Atom, "Atoms in ComparePDBs");
                        if (BoolPrint) cout <<"At end"<<endl;
                        BoolPrint=false;
		}
	}
        cout <<"About to delete memory"<<endl;
	delete [] SumWeightDist;
	delete [] SumDiffDist;
	delete [] SumDist;
	delete [] Sum2Dist;

	delete [] SumWeightType;
	delete [] SumDiffType;
	delete [] SumType;
	delete [] Sum2Type;
        RFactor2=SumDiff/Sum2;	
        cout <<"GreatestDifference= "<<GreatestDifference<<endl;
	cout <<"WorstCube= "<<WorstCube<<endl;
	cout <<"WorstMDCube= "<<WorstMDCube<<endl;
	PrintCubeInfo(Cubes1[WorstCube]);
	cout <<endl;
	PrintCubeInfo(Cubes2[WorstMDCube]);
	cout <<endl;
	cout <<"dWeight2_1= "<<dWeight2_1<<endl;
	cout <<"NumDiffGreaterThan_1= "<<NumDiffGreaterThan_1<<endl;
	cout <<"NonZeroGridPoints= "<<NonZeroGridPoints<<endl;
	cout <<"SumWeight= "<<SumWeight<<endl;
	cout <<"SqrDiff= "<<SqrDiff<<endl;
	cout <<"SumDiff= "<<SumDiff<<endl;
	cout <<"Sum= "<<Sum<<endl;
	cout <<"RFactor= "<<SumDiff/Sum<<endl;
	cout <<"RFactor2= "<<SumDiff/Sum2<<endl;
	cout <<"RMSD= "<<sqrt(SqrDiff*NonZeroGridPoints)/SumWeight<<endl;
	cout <<"Number of solvent cubes misidentified as protein= "<<SolventMisidentifiedAsProtein<<endl;
	cout <<"SqrDiff of those cubes= "<<SqrDiffSolvent<<endl;
        cout <<"Diff of those cubes= "<<DiffSolvent<<endl;
	cout <<"Number of protein cubes misidentified as solvent= "<<ProteinMisidentifiedAsSolvent<<endl;
	cout <<"SqrDiff of those cubes= "<<SqrDiffProtein<<endl;
        cout <<"Diff of thoese cubes= "<<DiffProtein<<endl;
	cout <<"SqrDiff of misidentified cubes= "<<SqrDiffSolvent+SqrDiffProtein<<endl;
	cout <<"Score= "<<SqrDiffSolvent+SqrDiffProtein<<endl;
	return SumDiff/Sum;
        //return SqrDiffSolvent+SqrDiffProtein;
}

Real ComparePDBs(Real contrast, vector<AtomStruct> &Atoms1, vector<AtomStruct> &Atoms2, vector<AtomStruct> &DiffAtoms, Real cubeSize, bool PrintDifferencePdb)
{
        vector<CubeStruct> Cubes1, Cubes2;

        ConvertAtomsToCubes(Atoms1, Cubes1);
        ConvertAtomsToCubes(Atoms2, Cubes2);

        return ComparePDBs(contrast, Cubes1, Cubes2, DiffAtoms, cubeSize, PrintDifferencePdb);

}

#endif
