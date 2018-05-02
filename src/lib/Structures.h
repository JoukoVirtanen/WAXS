#ifndef _Structures_included_
#define _Structures_included_

# include <iostream>
# include <fstream>
# include <string>
# include <cstring>
# include <sstream>
# include <time.h>
# include <vector>
# include <iomanip>
# include <complex>
# include <numeric>

//using namespace std;

# include "Constants.h"
# include "TypeDef.h"

class AtomStruct
{
	public:
		bool HetAtom, solute;
		int AtomNumber, ResidueNum;
		int atomid, AtomType, AtomType2;
		int ParmType;
		int residueid;
                int xbin, ybin, zbin;
		//int atomtype;
		vector<int> connectivity;
                Real ResidencyTime;
		Real weight, Electrons, charge, mass;
		//Real partialcharge;
		Real vdw;
		Real epsilon, epsilon14, vdw14, HydrationRadius;
		Real x, y, z, Occupancy, BFactor;
		string AtomName, CharmmAtomName, duplicate, ResidueName, ChainName, SegID, ID;
		bool operator == (AtomStruct Atom);
		AtomStruct();
};

struct ProteinStruct
{
	string name;
        int FirstResidue, LastResidue, NumResidues;
	vector<int> ResidueType;
        Real XBoxLength, YBoxLength, ZBoxLength;
        vector<Real> phi, psi, omega, dihedral;
        Matrix DistMatrix;
	vector<AtomStruct> Atoms;
	ProteinStruct& operator = (const ProteinStruct &Protein);
};

ProteinStruct& ProteinStruct::operator = (const ProteinStruct &Protein)
{
	name=Protein.name;
        FirstResidue=Protein.FirstResidue;
        LastResidue=Protein.LastResidue;
        NumResidues=Protein.NumResidues;
        ResidueType=Protein.ResidueType;
        XBoxLength=Protein.XBoxLength;
        YBoxLength=Protein.YBoxLength;
        ZBoxLength=Protein.ZBoxLength;
        phi=Protein.phi;
        psi=Protein.psi;
        omega=Protein.omega;
        dihedral=Protein.dihedral;
        DistMatrix=Protein.DistMatrix;
        Atoms=Protein.Atoms;
        return *this;
}

bool AtomStruct::operator== (AtomStruct Atom)
{
	if (HetAtom!=Atom.HetAtom) return false;
	if (AtomNumber!=Atom.AtomNumber) return false;
	if (ResidueNum!=Atom.ResidueNum) return false;
	if (atomid!=Atom.atomid) return false;
	if (AtomType!=Atom.AtomType) return false;
	if (AtomType2!=Atom.AtomType2) return false;
	if (residueid!=Atom.residueid) return false;
	//if (atomtype!=Atom.atomtype) return false;
	if (weight!=Atom.weight) return false;
	if (Electrons!=Atom.Electrons) return false;
	if (charge!=Atom.charge) return false;
	if (mass!=Atom.mass) return false;
	//if (partialcharge!=Atom.partialcharge) return false;
	if (vdw!=Atom.vdw) return false;
	if (epsilon!=Atom.epsilon) return false;
	if (epsilon14!=Atom.epsilon14) return false;
	if (vdw14!=Atom.vdw14) return false;
	if (HydrationRadius!=Atom.HydrationRadius) return false;
	if (x!=Atom.x) return false;
	if (y!=Atom.y) return false;
	if (z!=Atom.z) return false;

	return true;
}

bool CompareAtoms(AtomStruct Atom1, AtomStruct Atom2)
{
	if (Atom1==Atom2) return true;
	else return false;
}

AtomStruct::AtomStruct(void)
{
	HetAtom=false;
	solute=true;
        AtomNumber=0;
	ResidueNum=0;
	atomid=0;
	AtomType=0;
	AtomType2=0;
	ParmType=0;
	residueid=0;
        ResidencyTime=0;
        xbin=-1;
        ybin=-1;
        zbin=-1;
	//atomtype=0;
	weight=1.0;
	Electrons=0;
	charge=0;
	mass=0;
	//partialcharge=0;
	vdw=0;
	epsilon=0;
	epsilon14=0;
	vdw14=0;
	HydrationRadius=0;
	x=0;
	y=0;
	z=0;
	Occupancy=0;
	BFactor=0;
	AtomName="";
	CharmmAtomName="";
	duplicate="";
	ResidueName="";
	ChainName="A";
	SegID="";
	ID="";
}

struct VectorStruct
{
	Real x, y, z;
	VectorStruct operator + (VectorStruct);
	VectorStruct operator - (VectorStruct);
	VectorStruct operator * (Real);
	Real operator * (VectorStruct);
	VectorStruct operator / (Real);
};

class CubeStruct
{
	public:
		//bool angular, phobic, concave, occupied;
		bool angular, phobic, concave;
		//Real x, y, z, density, entropy, MinDist, deg;
		Real x, y, z, density, MinDist, deg;
                #ifndef __DELPHI__
		int AllAtomType, AtomType, ElementType; 
                #endif
		int IntDist, LocatedIn, nearest, SecondNearest; 
		//vector<Real> IntOrientation, ResidencyTimes, ComponentDensity;
		//vector<Real> IntOrientation, ComponentDensity;
                //vector<Real> ComponentDensity; Commented out 6/21/2014
                //vector<Real> IntOrientation; Commented out 6/21/2014
                //Real NumAtoms[NumAtomTypes];
		CubeStruct& operator = (const CubeStruct &cube);
		CubeStruct();
};

CubeStruct& CubeStruct::operator = (const CubeStruct &cube)
{
	angular=cube.angular;
	phobic=cube.phobic;
	concave=cube.concave;
        //occupied=cube.occupied;
	x=cube.x;
	y=cube.y;
	z=cube.z;
	density=cube.density;
        //entropy=cube.entropy;
	MinDist=cube.MinDist;
        deg=cube.deg;
        #ifndef __DELPHI__
	AllAtomType=cube.AllAtomType;
	AtomType=cube.AtomType;
	ElementType=cube.ElementType;
        #endif 
	IntDist=cube.IntDist;
	LocatedIn=cube.LocatedIn;
	nearest=cube.nearest;
	SecondNearest=cube.SecondNearest;
        //IntOrientation=cube.IntOrientation;
        //ResidencyTimes=cube.ResidencyTimes;
	return *this;
}

CubeStruct::CubeStruct(void)
{
	angular=false;
	phobic=false;
	concave=false;
        //occupied=false;
	x=0;
	y=0;
	z=0;
	density=0;
	//entropy=0;
        MinDist=0;
        deg=0;
        #ifndef __DELPHI__
	AllAtomType=0;
        AtomType=0;
	ElementType=0;
        #endif
	IntDist=0;
	LocatedIn=InProtein;
	nearest=0;
	SecondNearest=0;
}

void InitializeCube(CubeStruct &Cube)
{
	Cube.angular=false;
	Cube.phobic=false;
	Cube.concave=false;
	Cube.x=0;
	Cube.y=0;
	Cube.z=0;
	Cube.density=0;
        //Cube.entropy=0;
	Cube.MinDist=1000.0;
        Cube.deg=0;
        #ifndef __DELPHI__
	Cube.AllAtomType=0;
	Cube.AtomType=0;
	Cube.ElementType=0;
        #endif
	Cube.IntDist=0;
	Cube.LocatedIn=InProtein;
	Cube.nearest=0;
	Cube.SecondNearest=0;
	//for (int i=0;i<NumAtomTypes;i++) Cube.NumAtoms[i]=0;
}

void PrintCubeInfo(CubeStruct &Cube)
{
	cout <<"x= "<<Cube.x<<" y= "<<Cube.y<<" z= "<<Cube.z<<endl;
	cout <<"density= "<<Cube.density<<" MinDist= "<<Cube.MinDist<<" Cube.IntDist= "<<Cube.IntDist<<endl;
	cout <<"deg= "<<Cube.deg<<endl;
        #ifndef __DELPHI__
        cout <<"AllAtomType= "<<Cube.AllAtomType<<" AtomType= "<<Cube.AtomType<<" ElementType= "<<Cube.ElementType<<endl;
        #endif
	cout <<"LocatedIn= "<<Cube.LocatedIn<<" nearest= "<<Cube.nearest<<" SecondNearest= "<<Cube.SecondNearest<<endl;
	cout <<"angular= "<<Cube.angular<<" phobic= "<<Cube.phobic<<" concave= "<<Cube.concave<<endl;
}

void CopyCube(CubeStruct Cube1, CubeStruct &Cube2)
{
	Cube2.angular=Cube1.angular;
	Cube2.phobic=Cube1.phobic;
	Cube2.concave=Cube1.concave;
	Cube2.x=Cube1.x;
	Cube2.y=Cube1.y;
	Cube2.z=Cube1.z;
	Cube2.density=Cube1.density;
	//Cube2.entropy=Cube1.entropy;
        Cube2.MinDist=Cube1.MinDist;
        Cube2.deg=Cube1.deg;
        #ifndef __DELPHI__
	Cube2.AllAtomType=Cube1.AllAtomType;
	Cube2.AtomType=Cube1.AtomType;
	Cube2.ElementType=Cube1.ElementType;
        #endif
	Cube2.IntDist=Cube1.IntDist;
	Cube2.nearest=Cube1.nearest;
	Cube2.SecondNearest=Cube1.SecondNearest;
        //Cube2.IntOrientation=Cube1.IntOrientation;
        //Cube2.ResidencyTimes=Cube1.ResidencyTimes;
}

typedef vector< vector< vector<CubeStruct> > > lattice;

struct ElementStruct
{
	Real HydrationRadius, SolventCorrectedElectrons;
        vector<Real> alpha, beta;
};

struct PrdfStruct
{
        Matrix density;
        Array3D density3d;
        bool AngularDependence, AngularDependence2, ConvexConcave;
        bool ConvexConcave2, PhobicPhilic, UseChargeRadii;
        bool SecondNearestNeighbor;
        Real AngleBin, AngularCutOff, CubeSize, ChargeRadiiScale;
        Real HydrationRadiusScale, RecBin;
        vector< vector<Real> > AllAtomEntropy, AtomTypeEntropy, ElementEntropy;
        vector< vector<Real> > AllAtomResidency, AtomTypeResidency, ElementResidency;
        vector< vector<Real> > AllAtomHydrogenBond, AtomTypeHydrogenBond, ElementHydrogenBond;
        vector< vector<Real> > NumAllAtomCubes, NumAtomTypeCubes, NumElementCubes;
        vector< vector<Real> > AllAtomMass, AtomTypeMass, ElementMass;
        vector< vector<Real> > AllAtomPrdf, AtomTypePrdf, ElementPrdf;
        vector< vector< vector<Real> > > NumAllAtomSolvent, NumAtomTypeSolvent, NumElementSolvent;
        vector< vector< vector<Real> > > NumAllAtomOriented, NumAtomTypeOriented, NumElementOriented;
        vector< vector< vector<Real> > > AllAtomSolvent, AtomTypeSolvent, ElementSolvent;
        PrdfStruct();
};

PrdfStruct::PrdfStruct(void)
{
        density.resize(0);
}

/*
struct PrStruct
{
	vector<Real> r, all;
};
*/
#endif
