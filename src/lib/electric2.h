#ifndef _electric_included_
#define _electric_included_

# include "Structures.h"
# include "MathUtils.h"
# include "ReadPdb.h"
# include "ConvertStructures.h"
# include "TypeDef.h"
# include "VectorManip.h"
# include "Constants.h"

Real potential(Real x, Real y, Real z, AtomStruct &Atom)
{
        Real dx, dy, dz;
        Real r2, r, u;

        dx=x-Atom.x;
        dy=y-Atom.y;
        dz=z-Atom.z;

        r2=dx*dx+dy*dy+dz*dz;
        r=sqrt(r2);
        u=Atom.charge/r;
        return u;
}

Real potential(Real x, Real y, Real z, vector<AtomStruct> &Atoms)
{
        int natom=Atoms.size();
        Real u, ui;

        u=0;

        for (int i=0;i<natom;i++)
        {
                ui=potential(x, y, z, Atoms[i]);
                u+=ui;
        }
        return u;
}

Real columbPotential(Real q1, Real q2, Real r)
{
        return q1*q2/r;
}

Real columbPotential(Real q, Real r)
{
	return q/r;
}

Real columbPotential(AtomStruct &Atom1, AtomStruct &Atom2)
{
        Real dist=AtomDistance(Atom1, Atom2);
        return columbPotential(Atom1.charge, Atom2.charge, dist);
}

Real calcElectrostaticFreeEnergyOfAtom(AtomStruct &Atom, vector<CubeStruct> &Cubes, Real CubeSize)
{
        bool boolPrint=false;
        int ncube=Cubes.size();
        Real cubeCharge, CubeVolume, ePotential, freeEnergy=0;
        Real x, y, z;

        if (Atom.AtomName=="O1P" && Atom.ResidueName=="ADE")
        {
                boolPrint=false;
        }
        else boolPrint=false;

        

        CubeVolume=CubeSize*CubeSize*CubeSize;
        for (int j=0;j<ncube;j++)
        {
                x=Cubes[j].x;
                y=Cubes[j].y;
                z=Cubes[j].z;
                cubeCharge=Cubes[j].density*CubeVolume;  //density is charge density
                ePotential=potential(x, y, z, Atom);
                freeEnergy+=ePotential*cubeCharge;
                if (boolPrint)
                {
                cout <<endl;
                PrintCubeInfo(Cubes[j]);
                PrintAtomInfo(Atom);
                cout <<"Atom.charge= "<<Atom.charge<<endl;
                cout <<"cubeCharge= "<<cubeCharge<<" CubeSize= "<<CubeSize<<" ePotential= "<<ePotential<<endl;
                cout <<"freeEnergy= "<<freeEnergy<<endl;
                cout <<endl;
                }
        }
        //exit(EXIT_FAILURE);
        return freeEnergy;
}

Real calcElectrostaticFreeEnergyOfAtom(AtomStruct &Atom, lattice &Cubes, Real CubeSize)
{
        int MaxXBin, MaxYBin, MaxZBin;
        Real cubeCharge, CubeVolume, ePotential, freeEnergy;
        Real x, y, z;

        CubeVolume=CubeSize*CubeSize*CubeSize;
        Get3DVectorSize(Cubes, MaxXBin, MaxYBin, MaxZBin, "Cubes in calcElectrostaticFreeEnergyOfAtom");
        for (int xbin=0;xbin<MaxXBin;xbin++)
        {
                for (int ybin=0;ybin<MaxYBin;ybin++)
                {
                        for (int zbin=0;zbin<MaxZBin;zbin++)
                        {
                                x=Cubes[xbin][ybin][zbin].x;
                                y=Cubes[xbin][ybin][zbin].y;
                                z=Cubes[xbin][ybin][zbin].z;
                                cubeCharge=Cubes[xbin][ybin][zbin].density*CubeVolume;  //In this case density is charge density.
                                ePotential=potential(x, y, z, Atom);
                                freeEnergy+=ePotential*cubeCharge;
                        }
                }
        }
        return freeEnergy;
}

        template <class TCubes>
Real calcElectrostaticFreeEnergy(vector<AtomStruct> &Atoms, TCubes &Cubes, Real CubeSize)
{
        int natom=Atoms.size();
        Real FreeEnergy=0, sumFreeEnergy=0;
        vector<Real> SolventAtomCharge;
        cout <<"In calcElectrostaticFreeEnergy"<<endl;
        for (int i=0;i<natom;i++)
        {
                FreeEnergy=calcElectrostaticFreeEnergyOfAtom(Atoms[i], Cubes, CubeSize);
                //cout <<"FreeEnergy= "<<FreeEnergy<<endl;
                sumFreeEnergy+=FreeEnergy;
                //PrintAtomInfo(Atoms[i]);
        }
        sumFreeEnergy*=0.5;

        return sumFreeEnergy;

}

Real calcElectrostaticFreeEnergy(vector<AtomStruct> &Atoms, string CubeFile, Real CubeSize)
{
        vector<CubeStruct> Cubes;

        ReadPdb(CubeFile, Cubes);
        return calcElectrostaticFreeEnergy(Atoms, Cubes, CubeSize);
}

Real calcElectrostaticPotentialDueToAtoms(vector<AtomStruct> &Atoms, CubeStruct &Cube)
{
	return potential(Cube.x, Cube.y, Cube.z, Atoms);
}

Real calcElectrostaticPotentialDueToCubes(vector<CubeStruct> &Cubes, Real cubeSize, int index)
{
	int ncube=Cubes.size();
	Real u=0, r;
	Real dx, dy, dz;

	for (int i=0;i<ncube;i++)
	{
		if (i!=index) 
		{
			dx=Cubes[index].x-Cubes[i].x;
			dy=Cubes[index].y-Cubes[i].y;
			dz=Cubes[index].z-Cubes[i].z;
			r=sqrt(dx*dx+dy*dy+dz*dz);
			u+=Cubes[i].density/r;
		}
	}
	u*=cubeSize*cubeSize*cubeSize;
	return u;
}

Real calcSelfElectrostaticPotential(CubeStruct &Cube, Real cubeSize)
{
	Real r, q;

	r=calcPower(3.0/(4.0*pi), 1.0/3.0);
        r=1.290381*cubeSize;
	q=Cube.density*cubeSize*cubeSize*cubeSize;

	return q*3.0/(2.0*r);
}

Real calcElectrostaticPotentialOfCube(vector<AtomStruct> &Atoms, vector<CubeStruct> &Cubes, Real cubeSize, int index)
{
	Real u=0;

	u=calcElectrostaticPotentialDueToAtoms(Atoms, Cubes[index]);
	u+=calcElectrostaticPotentialDueToCubes(Cubes, cubeSize, index);
	u+=calcSelfElectrostaticPotential(Cubes[index], cubeSize);

        return u;
}

void calcElectrostaticPotentialEnergyMap(vector<AtomStruct> &Atoms, vector<CubeStruct> &Cubes, Real cubeSize)
{
	int ncube=Cubes.size();
	Real u;

	for (int i=0;i<ncube;i++)
	{
                if (Cubes[i].z==-3.8)
                {
		        u=calcElectrostaticPotentialOfCube(Atoms, Cubes, cubeSize, i);
		        Cubes[i].ComponentDensity[0]=u*SCALE_TO_KCAL;
                }
	}
}

#endif
