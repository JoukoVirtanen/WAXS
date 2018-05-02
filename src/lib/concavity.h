#ifndef _concavity_Included_
#define _concavity_Included_

# include "/home2/jouko/project/HeaderFiles/MathUtils.h"
# include "/home2/jouko/project/HeaderFiles/VectorManip.h"

struct PlaneStruct
{
	vector<Real> coef;
	string location;
};

void CalcPlane(vector<Real> &plane, AtomStruct &Atom1, AtomStruct &Atom2, AtomStruct &Atom3)
{
	vector< vector<Real> > coefficient;
	
	CreateMatrix(coefficient, 4, 3);

	coefficient[0][0]=Atom1.x; coefficient[1][0]=Atom1.y; coefficient[2][0]=1.0; coefficient[3][0]=Atom1.z;
	coefficient[0][1]=Atom2.x; coefficient[1][1]=Atom2.y; coefficient[2][1]=1.0; coefficient[3][1]=Atom2.z;
	coefficient[0][2]=Atom3.x; coefficient[1][2]=Atom3.y; coefficient[2][2]=1.0; coefficient[3][2]=Atom3.z;

	Equation(coefficient, plane);
}

string PlaneLocation(vector<Real> &plane, vector<AtomStruct> &Atoms)
{
	bool AllAbove=true, AllBelow=true;
	int natom=Atoms.size();
	Real z, dz=0.01;

	if (natom==0) return "middle";
        //cout <<endl;
	for (int i=0;i<natom;i++)
	{
		z=plane[0]*Atoms[i].x+plane[1]*Atoms[i].y+plane[2];
                //cout <<"z= "<<z<<" Atoms.z= "<<Atoms[i].z<<" AllBelow= "<<AllBelow<<" AllAbove= "<<AllAbove<<endl;
		if (z-dz>Atoms[i].z) AllBelow=false;
		if (z+dz<Atoms[i].z) AllAbove=false;
		if (!AllAbove && !AllBelow) return "middle";
	}
        //cout <<endl;
	if (AllAbove) return "AllAbove";
	if (AllBelow) return "AllBelow";
        cout <<"ERROR: Should not be at the end of PlaneLocation"<<endl;
        return "WTF";
}

void CheckPlane(vector<Real> &plane, AtomStruct &Atom1, AtomStruct &Atom2, AtomStruct &Atom3)
{
	Real z;

	z=plane[0]*Atom1.x+plane[1]*Atom1.y+plane[2];
	cout <<"z= "<<z<<" Atom.z= "<<Atom1.z<<endl;
	z=plane[0]*Atom2.x+plane[1]*Atom2.y+plane[2];
	cout <<"z= "<<z<<" Atom.z= "<<Atom2.z<<endl;
	z=plane[0]*Atom3.x+plane[1]*Atom3.y+plane[2];
	cout <<"z= "<<z<<" Atom.z= "<<Atom3.z<<endl;
	
}

void CalcPlanes(vector<AtomStruct> &Atoms, vector<PlaneStruct> &planes)
{
	int natom=Atoms.size();
	vector<Real> plane;
	PlaneStruct TempPlane;
        cout <<"In CalcPlanes. natom= "<<natom<<endl;
	plane.resize(3);
	for (int i=0;i<natom;i++)
	{
                cout <<"i= "<<i<<" natom= "<<natom<<endl;
		for (int j=i+1;j<natom;j++)
		{
			for (int k=j+1;k<natom;k++)
			{
				CalcPlane(plane, Atoms[i], Atoms[j], Atoms[k]);
				TempPlane.location=PlaneLocation(plane, Atoms);
				TempPlane.coef=plane;
				if (TempPlane.location=="AllBelow" || TempPlane.location=="AllAbove")
				{
					SafePushBack(planes, TempPlane, "planes");
                                        /*
                                        cout <<"TempPlane.coeff[0]= "<<TempPlane.coef[0]<<endl;
                                        cout <<"TempPlane.coeff[1]= "<<TempPlane.coef[1]<<endl;
                                        cout <<"TempPlane.coeff[2]= "<<TempPlane.coef[2]<<endl;
				        */
                                }
			}
		}
	}
}

#endif
