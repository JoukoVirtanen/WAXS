#ifndef _ComparePDBs_included_
#define _ComparePDBs_included_
 
# include "Structures.h"
# include "ConvertStructures.h"
# include "StringUtils.h"
# include "LinkedList.h"
# include "ReadPdb.h"
# include "TypeDef.h"
# include "VectorManip.h"
# include "WritePdb.h"

Real ComparePDBs(Real contrast, vector<CubeStruct> &Cubes1, vector<CubeStruct> &Cubes2, vector<AtomStruct> &Atoms)
{
	string DistFile, TypeFile;
	int m, n, t;
	int CubeNum, CubeNum2;
	vector< vector< vector<bool> > > Exists;
	int Xbins, Ybins, Zbins;
	int xgrid, ygrid, zgrid;
	int xtest, ytest, ztest;
	int WorstCube, WorstMDCube;
	vector< vector< vector<int> > > CubeConversion;
	Real SumDiff, Sum2;
	Real MinX, MinY, MinZ;
	Real MaxX, MaxY, MaxZ;
	vector< vector< vector<Real> > > CubeDensity;
	Real bin=0.5;
	Real CubeDensityXYZ;
	Real dWeight, dWeight2;
	Real MinCubeWeight, MinCubeWeight2;
	Real SqrDiffSolvent, SqrDiffProtein;
	Real SqrDiff, NonZeroGridPoints;
	Real SumWeight, Sum;
	Real RFactor2, MisidentifiedPoint;
	Real SolventMisidentifiedAsProtein;
	Real ProteinMisidentifiedAsSolvent;
	Real *SumWeightDist, *SumDiffDist;
	Real *SumDist, *Sum2Dist;
	Real *SumWeightType, *SumDiffType;
	Real *SumType, *Sum2Type;
	AtomStruct Atom;

	Atom.AtomNumber=1;
	Atom.ResidueNum=1;
	Atom.weight=1.0;
	Atom.Occupancy=1.0;
	Atom.AtomName="O";
	Atom.ResidueName="HOH";
	Atom.ChainName="A";

	SumWeightDist=new Real[1000];
	SumDiffDist=new Real[1000];
	SumDist=new Real[1000];
	Sum2Dist=new Real[1000];

	SumWeightType=new Real[1000];
	SumDiffType=new Real[1000];
	SumType=new Real[1000];
	Sum2Type=new Real[1000];
	cout <<"Allocated sum arrays"<<endl;
	for (n=0;n<1000;n++)
	{
		SumWeightDist[n]=0;
		SumDiffDist[n]=0;
		SumDist[n]=0;
		Sum2Dist[n]=0;

		SumWeightType[n]=0;
		SumDiffType[n]=0;
		SumType[n]=0;
		Sum2Type[n]=0;
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
        Create3DArray(CubeConversion, Xbins, Ybins, Zbins, "CubeConversion");
	Create3DArray(CubeDensity, Xbins, Ybins, Zbins, "CubeDensity");
	Create3DArray(Exists, Xbins, Ybins, Zbins, "Exists");

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

	for (n=0;n<CubeNum;n++)
	{
		Cubes1[n].density-=MinCubeWeight;
	}

	for (n=0;n<CubeNum2;n++)
	{
		Cubes2[n].density-=MinCubeWeight2;
	}
	cout <<"Initialized Exist"<<endl;
	xtest=int((-9.870-MinX)/bin+0.5);
	ytest=int((-8.001-MinY)/bin+0.5);
	ztest=int((-6.631-MinZ)/bin+0.5);
	for (n=0;n<CubeNum2;n++)
	{
		xgrid=int((Cubes2[n].x-MinX)/bin+0.5);
		ygrid=int((Cubes2[n].y-MinY)/bin+0.5);
		zgrid=int((Cubes2[n].z-MinZ)/bin+0.5);
		CubeDensity[xgrid][ygrid][zgrid]=Cubes2[n].density;
		//cout <<"CubeDensity["<<xgrid<<"]["<<ygrid<<"]["<<zgrid<<"]= "<<CubeDensity[xgrid][ygrid][zgrid]<<endl;
		Exists[xgrid][ygrid][zgrid]=true;
		CubeConversion[xgrid][ygrid][zgrid]=n;
	}
	cout <<"CubeDensity["<<xtest<<"]["<<ytest<<"]["<<ztest<<"]= "<<CubeDensity[xtest][ytest][ztest]<<endl;
	cout <<"xgrid_-9.870= "<<int((-9.870-MinX)/bin+0.5)<<endl;
	cout <<"xgrid_-8.001= "<<int((-8.001-MinY)/bin+0.5)<<endl;
	cout <<"xgrid_-6.631= "<<int((-6.631-MinZ)/bin+0.5)<<endl;
	//cout <<"xgrid_16.631= "<<int(floor((16.631-MinX+0.5)/bin));
	//cout <<"xgrid_16.630= "<<int(floor((16.630-MinX+0.5)/bin));
	SumDiff=0;
	Sum2=0;
	cout <<"Assigned Exists"<<endl;
	n=0;
	cout <<"Exists["<<xtest<<"]["<<ztest<<"]["<<ytest<<"]= "<<Exists[xtest][ytest][ztest]<<endl;

	SqrDiff=0;
	SolventMisidentifiedAsProtein=0;
	SqrDiffSolvent=0;
	ProteinMisidentifiedAsSolvent=0;
	NonZeroGridPoints=0;
	SqrDiffProtein=0;
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
                if (m==688381) BoolPrint=true;
                cout <<"m= "<<m<<endl;
		xgrid=int((Cubes1[m].x-MinX)/bin+0.5);
		ygrid=int((Cubes1[m].y-MinY)/bin+0.5);
		zgrid=int((Cubes1[m].z-MinZ)/bin+0.5);
		if (Exists[xgrid][ygrid][zgrid])
		{
                        if (BoolPrint) cout <<"At beginning"<<endl;
			CubeDensityXYZ=CubeDensity[xgrid][ygrid][zgrid];
			n=CubeConversion[xgrid][ygrid][zgrid];
                        if (BoolPrint) cout <<"n= "<<n<<" xgrid= "<<xgrid<<" ygrid= "<<ygrid<<" zgrid= "<<zgrid<<endl;
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
				WorstMDCube=n;
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
				SqrDiffSolvent+=dWeight2;
				n=CubeConversion[xgrid][ygrid][zgrid];
			}
	                if (BoolPrint) cout <<"SqrDiff= "<<SqrDiff<<endl;
			if (Cubes1[m].density==0 && CubeDensityXYZ!=0)
			{
				ProteinMisidentifiedAsSolvent+=1.0;
				SqrDiffProtein+=dWeight2;
				n=CubeConversion[xgrid][ygrid][zgrid];
                                cout <<"ProteinMisidentifiedAsSolvent"<<endl;
                                PrintCubeInfo(Cubes1[m]);
                                cout <<endl;
                                PrintCubeInfo(Cubes2[n]);
                                cout <<endl;
			}
                        if (BoolPrint) cout <<"if Cubes1[m]"<<endl;
			if (Cubes1[m].density>0 || CubeDensityXYZ>0)
			{
				SumWeight+=CubeDensityXYZ;
				SumDiff+=abs(dWeight);
				Sum+=abs(Cubes1[m].density+CubeDensityXYZ);
				Sum2+=abs(Cubes1[m].density+CubeDensityXYZ-0.33484*2.0);
				NonZeroGridPoints+=1.0;
				SumWeightDist[Cubes1[m].AtomType]+=CubeDensityXYZ;
				SumDiffDist[Cubes1[m].AtomType]+=abs(dWeight);
				SumDist[Cubes1[m].AtomType]+=abs(Cubes1[m].density+CubeDensityXYZ);
				Sum2Dist[Cubes1[m].AtomType]+=abs(Cubes1[m].density+CubeDensityXYZ-0.33484*2.0);
				SumWeightType[Cubes1[m].IntDist]+=CubeDensityXYZ;
				SumDiffType[Cubes1[m].IntDist]+=abs(dWeight);
				SumType[Cubes1[m].IntDist]+=abs(Cubes1[m].density+CubeDensityXYZ);
				Sum2Type[Cubes1[m].IntDist]+=abs(Cubes1[m].density+CubeDensityXYZ-0.33484*2.0);
			}
                        if (BoolPrint) cout <<"Near end"<<endl;
			Atom.x=Cubes1[m].x;
			Atom.y=Cubes1[m].y;
			Atom.z=Cubes1[m].z;
			Atom.weight=dWeight;
                        SafePushBack(Atoms, Atom, "Atoms");
                        if (BoolPrint) cout <<"At end"<<endl;
                        BoolPrint=false;
		}
	}
	delete [] SumWeightDist;
	delete [] SumDiffDist;
	delete [] SumDist;
	delete [] Sum2Dist;

	delete [] SumWeightType;
	delete [] SumDiffType;
	delete [] SumType;
	delete [] Sum2Type;
	
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
	cout <<"Number of protein cubes misidentified as solvent= "<<ProteinMisidentifiedAsSolvent<<endl;
	cout <<"SqrDiff of those cubes= "<<SqrDiffProtein<<endl;
	cout <<"SqrDiff of misidentified cubes= "<<SqrDiffSolvent+SqrDiffProtein<<endl;
	//cout <<"Score= "<<SqrDiffSolvent+SqrDiffProtein<<endl;
	return SumDiff/Sum;
}

#endif
