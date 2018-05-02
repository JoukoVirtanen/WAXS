#ifndef _AverageCubes_included_
#define _AverageCubes_included_

# include <vector>

# include "Structures.h"
# include "TypeDef.h"
# include "VectorManip.h"

void AverageCubes(vector<CubeStruct> &Cubes1, vector<CubeStruct> &Cubes2, vector<CubeStruct> &AverageCubes, Real NumStructures1, Real NumStructures2, Real bin)
{
	//Averages two density maps.
	int m, n, t;
	int CubeNum, CubeNum2;
	vector< vector< vector<bool> > > Exists1, Exists2;
	int Xbins, Ybins, Zbins;
	int xgrid, ygrid, zgrid;
	int xtest, ytest, ztest;
	vector< vector< vector<int> > > CubeConversion;
	Real SumDiff, Sum2;
	Real MinX, MinY, MinZ;
	Real MaxX, MaxY, MaxZ;
	vector< vector< vector<Real> > > CubeDensity1, CubeDensity2;
	Real CubeDensityXYZ;
	Real dWeight, dWeight2;
	Real MinCubeWeight, MinCubeWeight2;
	Real SqrDiffSolvent, SqrDiffProtein;
	Real RFactor2, MisidentifiedPoint;
	Real SolventMisidentifiedAsProtein;
	Real ProteinMisidentifiedAsSolvent;
	CubeStruct Cube;
	AtomStruct Atom;
	vector<AtomStruct> Atoms;

	Atom.AtomNumber=1;
	Atom.ResidueNum=1;
	Atom.weight=1.0;
	Atom.Occupancy=1.0;
	Atom.AtomName="O";
	Atom.ResidueName="HOH";
	Atom.ChainName="A";

	MinX=Cubes1[0].x;
	MinY=Cubes1[0].y;
	MinZ=Cubes1[0].z;
	MaxX=Cubes1[0].x;
	MaxY=Cubes1[0].y;
	MaxZ=Cubes1[0].z;
	CubeNum=Cubes1.size();
	CubeNum2=Cubes2.size();
	cout <<"CubeNum= "<<CubeNum<<" CubeNum2= "<<CubeNum2<<endl;
	for (n=0;n<CubeNum;n++)
	{
		if (Cubes1[n].x<MinX) MinX=Cubes1[n].x;
		if (Cubes1[n].y<MinY) MinY=Cubes1[n].y;
		if (Cubes1[n].z<MinZ) MinZ=Cubes1[n].z;
		if (Cubes1[n].x>MaxX) MaxX=Cubes1[n].x;
		if (Cubes1[n].y>MaxY) MaxY=Cubes1[n].y;
		if (Cubes1[n].z>MaxZ) MaxZ=Cubes1[n].z;
	}
	//cout <<"MinX= "<<MinX<<" MinY= "<<MinY<<" MinZ= "<<MinZ<<endl;
	for (n=0;n<CubeNum2;n++)
	{
		if (Cubes2[n].x<MinX) MinX=Cubes2[n].x;
		if (Cubes2[n].y<MinY) MinY=Cubes2[n].y;
		if (Cubes2[n].z<MinZ) MinZ=Cubes2[n].z;
		if (Cubes2[n].x>MaxX) MaxX=Cubes2[n].x;
		if (Cubes2[n].y>MaxY) MaxY=Cubes2[n].y;
		if (Cubes2[n].z>MaxZ) MaxZ=Cubes2[n].z;
	}
	cout <<"MinX= "<<MinX<<" MinY= "<<MinY<<" MinZ= "<<MinZ<<endl;
	Xbins=int(floor((MaxX-MinX)/bin))+10;
	Ybins=int(floor((MaxY-MinY)/bin))+10;
	Zbins=int(floor((MaxZ-MinZ)/bin))+10;
	//cout <<"Xbins= "<<Xbins<<" Ybins= "<<Ybins<<" Zbins= "<<Zbins<<endl;
	Safe3DAlloc(CubeDensity1, Xbins, Ybins, Zbins, "CubeDensity1");
	Safe3DAlloc(CubeDensity2, Xbins, Ybins, Zbins, "CubeDensity2");
	Safe3DAlloc(Exists1, Xbins, Ybins, Zbins, "Exists1");
	Safe3DAlloc(Exists2, Xbins, Ybins, Zbins, "Exists2");

	MinCubeWeight=Cubes1[0].density;
	MinCubeWeight2=Cubes2[0].density;
	//cout <<"Initialized Exist"<<endl;

	for (n=0;n<CubeNum;n++)
	{
		xgrid=int((Cubes1[n].x-MinX)/bin+0.5);
		ygrid=int((Cubes1[n].y-MinY)/bin+0.5);
		zgrid=int((Cubes1[n].z-MinZ)/bin+0.5);
		CubeDensity1[xgrid][ygrid][zgrid]=Cubes1[n].density;
		Exists1[xgrid][ygrid][zgrid]=true;
	}
	//cout <<"After exist1"<<endl;
	for (n=0;n<CubeNum2;n++)
	{
		xgrid=int((Cubes2[n].x-MinX)/bin+0.5);
		ygrid=int((Cubes2[n].y-MinY)/bin+0.5);
		zgrid=int((Cubes2[n].z-MinZ)/bin+0.5);
		CubeDensity2[xgrid][ygrid][zgrid]=Cubes2[n].density;
		Exists2[xgrid][ygrid][zgrid]=true;
	}
	SumDiff=0;
	Sum2=0;
	//cout <<"Assigned Exists"<<endl;
	n=0;
	for (m=0;m<CubeNum;m++)
	{
		xgrid=int((Cubes1[m].x-MinX)/bin+0.5);
		ygrid=int((Cubes1[m].y-MinY)/bin+0.5);
		zgrid=int((Cubes1[m].z-MinZ)/bin+0.5);
		Cube.x=Cubes1[m].x;
		Cube.y=Cubes1[m].y;
		Cube.z=Cubes1[m].z;
                Cube.IntDist=Cubes1[m].IntDist;
                Cube.AtomType=Cubes1[m].AtomType;
		if (Exists2[xgrid][ygrid][zgrid])
		{
			CubeDensityXYZ=CubeDensity2[xgrid][ygrid][zgrid];
			Cube.density=(Cubes1[m].density*NumStructures1+CubeDensityXYZ*NumStructures2)/(NumStructures1+NumStructures2);
                        //cout <<"Cube.density= "<<Cube.density<<" Cubes1["<<m<<"].density= "<<Cubes1[m].density<<" NumStructures1= "<<NumStructures1<<" CubeDensityXYZ= "<<CubeDensityXYZ<<" NumStructures2= "<<NumStructures2<<endl;
		}
		else Cube.density=Cubes1[m].density*NumStructures1/(NumStructures1+NumStructures2);
                SafePushBack(AverageCubes, Cube, "AverageCubes");
	}

	for (m=0;m<CubeNum2;m++)
	{
		xgrid=int((Cubes2[m].x-MinX)/bin+0.5);
		ygrid=int((Cubes2[m].y-MinY)/bin+0.5);
		zgrid=int((Cubes2[m].z-MinZ)/bin+0.5);
		if (!Exists1[xgrid][ygrid][zgrid])
		{
			Cube.x=Cubes2[m].x;
			Cube.y=Cubes2[m].y;
			Cube.z=Cubes2[m].z;
			Cube.density=Cubes2[m].density*NumStructures2/(NumStructures1+NumStructures2);
                        Cube.IntDist=Cubes2[m].IntDist;
                        Cube.AtomType=Cubes2[m].AtomType;
                        SafePushBack(AverageCubes, Cube, "AverageCubes");
		}
	}
}
#endif
