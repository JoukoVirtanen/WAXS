#ifndef _AtomUtils_included_
#define _AtomUtils_included_

# include "AtomCode.h"
# include "AtomIDs.h"
# include "GetAtomType.h"
# include "GetAtomType2.h"
# include "MinMax.h"
# include "ReadPdb.h"
# include "ReadPrm.h"
# include "ReadPsf.h"
# include "ReadRtf.h"
# include "ResidueCode.h"
# include "Structures.h"
# include "VectorManip.h"

#ifndef _ReadPdb_included_
#error ReadPdb not included
#endif

const Real UNK_REAL=-666.0;

void CopyAtom(AtomStruct Atom, AtomStruct &Atom2)
{
	Atom2.AtomNumber=Atom.AtomNumber;
	Atom2.ResidueNum=Atom.ResidueNum;
	Atom2.atomid=Atom.atomid;
	Atom2.AtomType=Atom.AtomType;
	Atom2.ParmType=Atom.ParmType;
	Atom2.residueid=Atom.residueid;
	Atom2.weight=Atom.weight;
	Atom2.Electrons=Atom.Electrons;
	Atom2.charge=Atom.charge;
	Atom2.mass=Atom.mass;
	Atom2.x=Atom.x;
	Atom2.y=Atom.y;
	Atom2.z=Atom.z;
	Atom2.vdw=Atom.vdw;
	Atom2.epsilon=Atom.epsilon;
	Atom2.vdw14=Atom.vdw14;
	Atom2.epsilon14=Atom.epsilon14;
	Atom2.Occupancy=Atom.Occupancy;
	Atom2.BFactor=Atom.BFactor;
	Atom2.AtomName=Atom.AtomName;
	Atom2.CharmmAtomName=Atom.CharmmAtomName;
	Atom2.ResidueName=Atom.ResidueName;
	Atom2.ChainName=Atom.ChainName;
	Atom2.SegID=Atom.SegID;
	Atom2.ID=Atom.ID;
}

void InitializeAtom(AtomStruct &Atom)
{
	Atom.AtomNumber=0;
	Atom.ResidueNum=0;
	Atom.atomid=0;
	Atom.AtomType=0;
	Atom.AtomType2=0;
	Atom.ParmType=0;
	Atom.residueid=0;
	Atom.atomtype=0;
	Atom.weight=1.0;
	Atom.Electrons=0;
	Atom.charge=0;
	Atom.mass=0;
	Atom.x=0;
	Atom.y=0;
	Atom.z=0;
	Atom.vdw=0;
	Atom.epsilon=0;
	Atom.vdw14=0;
	Atom.epsilon14=0;
	Atom.Occupancy=0;
	Atom.BFactor=0;
	Atom.AtomName="";
	Atom.CharmmAtomName="";
	Atom.ResidueName="";
	Atom.ChainName="A";
	Atom.SegID="";
	Atom.ID="";
}

void MoveAtoms(vector<AtomStruct> &Atoms, Real dx, Real dy, Real dz)
{
	int i, natom=Atoms.size();
	for (i=0;i<natom;i++)
	{
		Atoms[i].x+=dx;
		Atoms[i].y+=dy;
		Atoms[i].z+=dz;
	}
}

void MergeProteins(vector<ProteinStruct> Proteins, vector<AtomStruct> &Atoms)
{
	int i, nprotein=Proteins.size();
	DeleteVector(Atoms);

	for (i=0;i<nprotein;i++) AppendVector(Proteins[i].Atoms, Atoms);
}

void CalcDistMatrix(vector< vector<Real> > &DistMatrix, vector<AtomStruct> &Atoms)
{
	int i, j, natom=Atoms.size();
	Real dist;
	Real dx, dy, dz;
	for (i=0;i<natom;i++)
	{
		for (j=i+1;j<natom;j++)
		{
			dx=Atoms[i].x-Atoms[j].x;
			dy=Atoms[i].y-Atoms[j].y;
			dz=Atoms[i].z-Atoms[j].z;
			dist=sqrt(dx*dx+dy*dy+dz*dz);
			DistMatrix[i][j]=dist;
			DistMatrix[j][i]=dist;
		}
	}
}

Real AtomDistance(AtomStruct Atom1, AtomStruct Atom2)
{
	Real dx, dy, dz;
	dx=Atom1.x-Atom2.x;
	dy=Atom1.y-Atom2.y;
	dz=Atom1.z-Atom2.z;
	return sqrt(dx*dx+dy*dy+dz*dz);
}

void RotateAroundX(AtomStruct &Atom, Real theta)
{
	Real tempY, tempZ;
	tempY=Atom.y*cos(theta)-Atom.z*sin(theta);
	tempZ=Atom.y*sin(theta)+Atom.z*cos(theta);
	Atom.y=tempY;
	Atom.z=tempZ;
}

void RotateAroundY(AtomStruct &Atom, Real theta)
{
	Real tempX, tempZ;
	tempX=Atom.x*cos(theta)-Atom.z*sin(theta);
	tempZ=Atom.x*sin(theta)+Atom.z*cos(theta);
	Atom.x=tempX;
	Atom.z=tempZ;
}

void RotateAroundZ(AtomStruct &Atom, Real theta)
{
	Real tempX, tempY;
	tempX=Atom.x*cos(theta)-Atom.y*sin(theta);
	tempY=Atom.x*sin(theta)+Atom.y*cos(theta);
	Atom.x=tempX;
	Atom.y=tempY;
}

/*
void CreateRotMatrix(vector<Real> Axis, Real angle, vector< vector<Real> > &RotMat)
{
	Real cosine=cos(angle);
	Real sine=sin(angle);
	Real one_minus_cosine=1-cosine;

	RotMat[X][X]=Sqr(Axis[X])+(1-Sqr(Axis[x]))*cosine;
	RotMat[X][Y]=Axis[X]*Axis[Y]*one_minus_cosine+Axis[Z]*sine;
	RotMat[X][Z]=Axis[X]*Axis[Z]*one_minus_cosine-Axis[Y]*sine;
}
*/
void RotateAtoms(vector<AtomStruct> &Atoms, Real theta, Real phi, Real psi)
{
	int i, natom=Atoms.size();
	for (i=0;i<natom;i++)
	{
		RotateAroundX(Atoms[i], theta);
		RotateAroundY(Atoms[i], phi);
		RotateAroundZ(Atoms[i], psi);
	}
} 

bool IsWater(AtomStruct Atom)
{
	if (Atom.ResidueName=="HOH") return true;
	if (Atom.ResidueName=="TIP") return true;
	if (Atom.ResidueName=="TIP3") return true;
	if (Atom.ResidueName=="TIP4") return true;
	if (Atom.ResidueName=="TIP5") return true;
	if (Atom.ResidueName=="WAT") return true;
	return false;
}

void RemoveWaters(vector<AtomStruct> &Atoms)
{
	vector<AtomStruct> TempAtoms;
	int i, natom=Atoms.size();
	for (i=0;i<natom;i++)
	{
		if (!IsWater(Atoms[i])) TempAtoms.push_back(Atoms[i]);
	}
	CopyVector(TempAtoms, Atoms);
}

void RemoveUnknownResidues(vector<AtomStruct> &Atoms)
{
	vector<AtomStruct> TempAtoms;
	int i, natom=Atoms.size();
	for (i=0;i<natom;i++)
	{
		if (Atoms[i].residueid!=UNK) TempAtoms.push_back(Atoms[i]);
	}
	CopyVector(TempAtoms, Atoms);
}

int FindNumProteinAtoms(vector<AtomStruct> Atoms)
{
	int n, TotalParticles=Atoms.size();
	for (n=0;n<TotalParticles;n++) 
	{
		if (IsWater(Atoms[n])) return n;
		if (Atoms[n].atomid==HydrationShell) return n;
		if (Atoms[n].atomid==ExcludedVolume) return n;
	}
	return TotalParticles;
}

void CalcCenter(vector<AtomStruct> &Atoms, Real &xave, Real &yave, Real &zave)
{
	Real xsum, ysum, zsum;
	Real NumAve=0;
	int i, natom=Atoms.size();
	xsum=ysum=zsum=0;
	for (i=0;i<natom;i++)
	{
		if (Atoms[i].AtomName=="CA")
		{
			xsum+=Atoms[i].x;
			ysum+=Atoms[i].y;
			zsum+=Atoms[i].z;
			NumAve+=1.0;
		}
	}
	xave=xsum/NumAve;
	yave=ysum/NumAve;
	zave=zsum/NumAve;
}

void center(vector<AtomStruct> &Atoms)
{
	Real xave, yave, zave;
	CalcCenter(Atoms, xave, yave, zave);
	MoveAtoms(Atoms, -xave, -yave, -zave);
}

void CenterAtomsInBox(vector<AtomStruct> &Atoms)
{
	Real XMax, YMax, ZMax;
	Real XMin, YMin, ZMin;
	MinMax(XMin, YMin, ZMin, XMax, YMax, ZMax, Atoms);
	cout <<"XMin= "<<XMin<<" YMin= "<<YMin<<" ZMin= "<<ZMin<<endl;
	cout <<"XMax= "<<XMax<<" YMax= "<<YMax<<" ZMax= "<<ZMax<<endl;
	//PrintPDB("/home2/jouko/Density/TempPdb/reference.pdb", Atoms);
	MoveAtoms(Atoms, -XMax+0.5*(XMax-XMin), -YMax+0.5*(YMax-YMin), -ZMax+0.5*(ZMax-ZMin));
	//PrintPDB("/home2/jouko/Density/TempPdb/reference.pdb", Atoms);
}

void ApplyPeriodicBoundaryConditions(vector<AtomStruct> &Atoms, Real XBoxLength, Real YBoxLength, Real ZBoxLength)
{
	int n, TotalParticles=Atoms.size();
	for (n=0;n<TotalParticles;n++)
	{
		if (Atoms[n].x>XBoxLength*0.5) Atoms[n].x-=XBoxLength;
		if (Atoms[n].x<-XBoxLength*0.5) Atoms[n].x+=XBoxLength;
		if (Atoms[n].y>YBoxLength*0.5) Atoms[n].y-=YBoxLength;
		if (Atoms[n].y<-YBoxLength*0.5) Atoms[n].y+=YBoxLength;
		if (Atoms[n].z>ZBoxLength*0.5) Atoms[n].z-=ZBoxLength;
		if (Atoms[n].z<-ZBoxLength*0.5) Atoms[n].z+=ZBoxLength;
	}
}

void SetCharmmAtomNames(vector<AtomStruct> CharmmAtoms, vector<AtomStruct> &Atoms)
{
	int i, j, NumCharmmAtoms=CharmmAtoms.size(), NumAtoms=Atoms.size();
	for (i=0;i<NumAtoms;i++)
	{
		for (j=0;j<NumCharmmAtoms;j++)
		{
			if (CharmmAtoms[j].ResidueName==Atoms[i].ResidueName)
			{
				if (CharmmAtoms[j].AtomName==Atoms[i].AtomName)
				{
					Atoms[i].CharmmAtomName=CharmmAtoms[j].CharmmAtomName;
					break;
				}
			}
		}
	}
}

void CharmmAtomNamesToVDW(vector<AtomStruct> CharmmAtoms, vector<AtomStruct> &Atoms)
{
	int i, j, NumCharmmAtoms=CharmmAtoms.size(), NumAtoms=Atoms.size();
	Real DefaultVDW=1.5;
	for (i=0;i<NumAtoms;i++)
	{
		Atoms[i].vdw=UNK_REAL;
		//if (Atoms[i].AtomName.substr(0,2)=="HT") Atoms[i].CharmmAtomName="HC";
		//if (Atoms[i].AtomName.substr(0,2)=="OT") Atoms[i].CharmmAtomName="OC";
		for (j=0;j<NumCharmmAtoms;j++)
		{
			if ((Atoms[i].AtomName=="HT1" || Atoms[i].AtomName=="HT2" || Atoms[i].AtomName=="HT3") && CharmmAtoms[j].AtomName=="HC")
			{
				Atoms[i].vdw=CharmmAtoms[j].vdw;
				Atoms[i].epsilon=CharmmAtoms[j].epsilon;
				Atoms[i].vdw14=CharmmAtoms[j].vdw14;
				Atoms[i].epsilon14=CharmmAtoms[j].epsilon14;
				Atoms[i].charge=CharmmAtoms[j].charge;
				break;
			}
			if ((Atoms[i].AtomName=="OT1" || Atoms[i].AtomName=="OT2") && CharmmAtoms[j].AtomName=="OC")
			{
				Atoms[i].vdw=CharmmAtoms[j].vdw;
				Atoms[i].epsilon=CharmmAtoms[j].epsilon;
				Atoms[i].vdw14=CharmmAtoms[j].vdw14;
				Atoms[i].epsilon14=CharmmAtoms[j].epsilon14;
				Atoms[i].charge=CharmmAtoms[j].charge;
				break;
			}
			if (CharmmAtoms[j].CharmmAtomName==Atoms[i].CharmmAtomName)
			{
				PrintAtomInfo(Atoms[i]);
				Atoms[i].vdw=CharmmAtoms[j].vdw;
				Atoms[i].epsilon=CharmmAtoms[j].epsilon;
				Atoms[i].vdw14=CharmmAtoms[j].vdw14;
				Atoms[i].epsilon14=CharmmAtoms[j].epsilon14;
				Atoms[i].charge=CharmmAtoms[j].charge;
				cout <<"vdw= "<<Atoms[i].vdw<<endl;
				cout <<endl;
				break;
			}
		}
		if (Atoms[i].vdw==UNK_REAL)
		{
			Atoms[i].vdw=DefaultVDW;
			cout <<"Warning.  Unable to set vdw for atom "<<i<<endl;
			cout <<"Setting vdw to "<<DefaultVDW<<endl;
		}
	}
}

void ScaleVDW(Real VDWScale, Real VDWOffSet, vector<AtomStruct> &Atoms)
{
	int i, natom=Atoms.size();
	for (i=0;i<natom;i++)
	{
		Atoms[i].vdw=Atoms[i].vdw*VDWScale+VDWOffSet;
	}
}

void SetVDW(string PrmFile, string RtfFile, vector<AtomStruct> &Atoms)
{
	vector<AtomStruct> CharmmAtoms, TempCharmmAtoms;
	ReadRtfFile(RtfFile, CharmmAtoms);
	ReadPrmFile(PrmFile, TempCharmmAtoms);
	SetCharmmAtomNames(CharmmAtoms, Atoms);
	CharmmAtomNamesToVDW(TempCharmmAtoms, Atoms);
}

void ParseAtoms(string PdbFile, vector<AtomStruct> &Atoms, bool BoolRemoveWaters, bool BoolRemoveUnknownResidues)
{
	string str;
	int position;
	position=PdbFile.rfind(".");
	if (position!=string::npos)
	{
		str=PdbFile.substr(position+1,4);
		if (str=="pdb") ReadPdb(PdbFile, Atoms);
		else ReadXYZ(PdbFile, Atoms);
	}
	else ReadPdb(PdbFile, Atoms);
	if (BoolRemoveWaters) RemoveWaters(Atoms);
	if (BoolRemoveUnknownResidues) RemoveUnknownResidues(Atoms);
	AssignAtomIDs(Atoms);
	GetAtomType(Atoms);
	GetAtomType2(Atoms);		
}
#endif
