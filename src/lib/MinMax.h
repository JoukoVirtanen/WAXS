#ifndef _MinMax_included_
#define _MinMax_included_

# include "Structures.h"
# include "TypeDef.h"

void MinMax(Real &XMin, Real &YMin, Real &ZMin, Real &XMax, Real &YMax, Real &ZMax, vector<AtomStruct> &Atoms);
void MinMaxProtein(Real &XMin, Real &YMin, Real &ZMin, Real &XMax, Real &YMax, Real &ZMax, vector<AtomStruct> Atoms);
template <class Vector_T>
void MinMax(Real &XMin, Real &YMin, Real &ZMin, Real &XMax, Real &YMax, Real &ZMax, Vector_T &pos);

# include "AtomUtils.h"

int findMax(vector<int> &v)
{
        int Size=v.size();
        int max=v[0];

        for (int i=1;i<Size;i++)
        {
                if (v[i]>max) max=v[i];
        }
        return max;
}

int findMin(vector<Real> &v, Real &min)
{
        int Size=v.size(), lowest;
        min=v[0];
        lowest=0;
        for (int i=1;i<Size;i++)
        {
                if (v[i]>min) 
                {
                        min=v[i];
                        lowest=i;
                }
        }
        return lowest;
}

void MinMax(vector<Real> &v, Real &min, Real &max)
{
        int Size=v.size();

        min=v[0];
        max=v[0];
        for (int i=0;i<Size;i++)
        {
                if (v[i]<min) min=v[i];
                else if (v[i]>max) max=v[i];
        }
}

void MinMax(Real &XMin, Real &YMin, Real &ZMin, Real &XMax, Real &YMax, Real &ZMax, vector<AtomStruct> &Atoms)
{
	int n, natom;
        cout <<"In MinMax"<<endl;
	natom=Atoms.size();
        cout <<"natom= "<<natom<<endl;
	if (natom==0)
        {
                cout <<"There are no atoms in the protein"<<endl;
                exit(EXIT_FAILURE);
        }
        XMin=Atoms[0].x;
	YMin=Atoms[0].y;
	ZMin=Atoms[0].z;
	XMax=Atoms[0].x;
	YMax=Atoms[0].y;
	ZMax=Atoms[0].z;
	for (n=0;n<natom;n++)
	{
                //cout <<"n= "<<n<<endl;
		if (Atoms[n].x<XMin) XMin=Atoms[n].x;
		if (Atoms[n].y<YMin) YMin=Atoms[n].y;
		if (Atoms[n].z<ZMin) ZMin=Atoms[n].z;
		if (Atoms[n].x>XMax) XMax=Atoms[n].x;
		if (Atoms[n].y>YMax) YMax=Atoms[n].y;
		if (Atoms[n].z>ZMax) ZMax=Atoms[n].z;
	}
}

void MinMaxProtein(Real &XMin, Real &YMin, Real &ZMin, Real &XMax, Real &YMax, Real &ZMax, vector<AtomStruct> Atoms)
{
	RemoveUnknownResidues(Atoms);
	RemoveWaters(Atoms);
	MinMax(XMin, YMin, ZMin, XMax, YMax, ZMax, Atoms);
	cout <<"XMin= "<<XMin<<" YMin= "<<YMin<<" ZMin= "<<ZMin<<endl;
	cout <<"XMax= "<<XMax<<" YMax= "<<YMax<<" ZMax= "<<ZMax<<endl;
}

template <class Vector_T>
void MinMax(Real &XMin, Real &YMin, Real &ZMin, Real &XMax, Real &YMax, Real &ZMax, Vector_T &pos)
{
	int n, natom;
	natom=pos.size();
	XMin=pos[0].x;
	YMin=pos[0].y;
	ZMin=pos[0].z;
	XMax=pos[0].x;
	YMax=pos[0].y;
	ZMax=pos[0].z;
	for (n=0;n<natom;n++)
	{
		if (pos[n].x<XMin) XMin=pos[n].x;
		if (pos[n].y<YMin) YMin=pos[n].y;
		if (pos[n].z<ZMin) ZMin=pos[n].z;
		if (pos[n].x>XMax) XMax=pos[n].x;
		if (pos[n].y>YMax) YMax=pos[n].y;
		if (pos[n].z>ZMax) ZMax=pos[n].z;
	}
}

void MinMax(Real &XMin, Real &YMin, Real &ZMin, Real &XMax, Real &YMax, Real &ZMax, lattice &cubes)
{
        int MaxXBin, MaxYBin, MaxZBin;
        Get3DVectorSize(cubes, MaxXBin, MaxYBin, MaxZBin, "cubes in MinMax");
	
        XMin=cubes[0][0][0].x;
	YMin=cubes[0][0][0].y;
	ZMin=cubes[0][0][0].z;
	XMax=cubes[0][0][0].x;
	YMax=cubes[0][0][0].y;
	ZMax=cubes[0][0][0].z;

        for (int i=0;i<MaxXBin;i++)
        {
                for (int j=0;j<MaxYBin;j++)
                {
                        for (int k=0;k<MaxZBin;k++)
                        {
		                if (cubes[i][j][k].x<XMin) XMin=cubes[i][j][k].x;
		                if (cubes[i][j][k].y<YMin) YMin=cubes[i][j][k].y;
		                if (cubes[i][j][k].z<ZMin) ZMin=cubes[i][j][k].z;
		                if (cubes[i][j][k].x>XMax) XMax=cubes[i][j][k].x;
		                if (cubes[i][j][k].y>YMax) YMax=cubes[i][j][k].y;
		                if (cubes[i][j][k].z>ZMax) ZMax=cubes[i][j][k].z;
                        }
                }
        }
}
#endif

