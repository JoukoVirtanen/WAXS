#ifndef _ComparePDBs_included_
#define _ComparePDBs_included_
 
Real ComparePDBs(Real contrast, vector<CubeStruct> &Cubes1, vector<CubeStruct> &Cubes2, vector<AtomStruct> &Atoms, bool PrintDifferencePdb=false);
Real ComparePDBs(Real contrast, vector<AtomStruct> &Atoms1, vector<AtomStruct> &Atoms2, vector<AtomStruct> &DiffAtoms, bool PrintDifferencPdb=false);

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

Real IntersectionCubes(vector<CubeStruct> &Cubes1, vector<CubeStruct> &Cubes2, vector<AtomStruct> &IntersectionAtoms1, vector<AtomStruct> &IntersectionAtoms2, bool PrintDifferencePdb)
{
        bool verbose=false;
	string DistFile, TypeFile;
	int m, n, t;
	int CubeNum, CubeNum2;
        vector< vector< vector<int> > > conversion;
	vector< vector< vector<bool> > > Exists;
	int Xbins, Ybins, Zbins;
	int xgrid, ygrid, zgrid;
	int xtest, ytest, ztest;
        int MaxAtomType=0;
	int WorstCube, WorstMDCube;
	Real SumDiff, Sum2;
	Real MinX, MinY, MinZ;
	Real MaxX, MaxY, MaxZ;
	Real bin=0.5;
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
	AtomStruct Atom;
        IntersectionAtoms1.clear();
        IntersectionAtoms2.clear();
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
	Safe3DAlloc(conversion, Xbins, Ybins, Zbins, "conversion");
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
		conversion[xgrid][ygrid][zgrid]=n;
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

        bool BoolPrint=false;
	for (m=0;m<CubeNum;m++)
	{
                //cout <<"m= "<<m<<endl;
		xgrid=int((Cubes1[m].x-MinX)/bin+0.5);
		ygrid=int((Cubes1[m].y-MinY)/bin+0.5);
		zgrid=int((Cubes1[m].z-MinZ)/bin+0.5);
		if (Exists[xgrid][ygrid][zgrid])
		//if (true)
		{
                        if (BoolPrint) cout <<"At beginning"<<endl;
                        if (BoolPrint) cout <<"n= "<<n<<" xgrid= "<<xgrid<<" ygrid= "<<ygrid<<" zgrid= "<<zgrid<<endl;
	                if (xgrid==xtest && ygrid==ytest && zgrid==ztest)
                        {
                                cout <<"TestCube= "<<TestCube<<endl;
                                PrintCubeInfo(Cubes2[TestCube]);
                                cout <<endl;
                                PrintCubeInfo(Cubes1[m]);
                        }
	                if (BoolPrint) cout <<"SqrDiff= "<<SqrDiff<<endl;
                        if (BoolPrint) cout <<"if Cubes1[m]"<<endl;
			if (Cubes1[m].density!=0 && Cubes2[conversion[xgrid][ygrid][zgrid]].density!=0)
			{
			        //Atom.x=Cubes1[m].x;
			        //Atom.y=Cubes1[m].y;
			        //Atom.z=Cubes1[m].z;
			        //Atom.weight=Cubes1[m].density;
                                ConvertCubeToAtom(Cubes1[m], Atom);
                                SafePushBack(IntersectionAtoms1, Atom, "IntersectionAtoms1");
                                ConvertCubeToAtom(Cubes2[conversion[xgrid][ygrid][zgrid]], Atom);
                                SafePushBack(IntersectionAtoms2, Atom, "IntersectionAtoms2");
                        }
                        if (BoolPrint) cout <<"Near end"<<endl;
                        if (BoolPrint) cout <<"At end"<<endl;
                        BoolPrint=false;
		}
	}
	//cout <<"Score= "<<SqrDiffSolvent+SqrDiffProtein<<endl;
	return SumDiff/Sum;
        //return SqrDiffSolvent+SqrDiffProtein;
}

Real IntersectionCubes(vector<AtomStruct> &Atoms1, vector<AtomStruct> &Atoms2, vector<AtomStruct> &IntersectionAtoms1, vector<AtomStruct> &IntersectionAtoms2, bool PrintDifferencePdb)
{
        vector<CubeStruct> Cubes1, Cubes2;
        RemoveNonHetatm(Atoms1);
        RemoveNonHetatm(Atoms2);
        ConvertAtomsToCubes(Atoms1, Cubes1);
        ConvertAtomsToCubes(Atoms2, Cubes2);

        IntersectionCubes(Cubes1, Cubes2, IntersectionAtoms1, IntersectionAtoms2, PrintDifferencePdb);

}

#endif
