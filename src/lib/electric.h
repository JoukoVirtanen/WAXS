#ifndef _electric_included_
#define _electric_included_

# include "Structures.h"
# include "MathUtils.h"
# include "ReadPdb.h"
# include "ConvertStructures.h"
# include "TypeDef.h"
# include "VectorManip.h"

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

Real potential(Real r, Real q)
{
        return q/r;
}

Real potential(Real x1, Real y1, Real z1, Real x2, Real y2, Real z2, Real q)
{
        Real dx, dy, dz;
        Real r;

        dx=x2-x1;
        dy=y2-y1;
        dz=z2-z1;
        r=sqrt(dx*dx+dy*dy+dz*dz);
        return potential(r, q);
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

Real potential(Real x, Real y, Real z, CubeStruct Cube, Real CubeSize)
{
        Real q=Cube.density*CubeSize*CubeSize*CubeSize;

        return potential(x, y, z, Cube.x, Cube.y, Cube.z, q);
}

Real columbPotential(AtomStruct &Atom1, AtomStruct &Atom2)
{
        Real dist=AtomDistance(Atom1, Atom2);
        return columbPotential(Atom1.charge, Atom2.charge, dist);
}

Real calcElectrostaticPotential(Real x, Real y, Real z, vector<AtomStruct> &Atoms)
{
        int natom=Atoms.size();
        Real v=0;
        for (int i=0;i<natom;i++)
        {
                if (x!=Atoms[i].x || y!=Atoms[i].y || z!=Atoms[i].z)
                {
                        v+=potential(x, y, z, Atoms[i]);
                }
        }
        return v;
}

Real calcElectrostaticPotential(Real x, Real y, Real z, vector<CubeStruct> &Cubes, Real CubeSize)
{
        int ncube=Cubes.size();
        Real v=0;
        for (int i=0;i<ncube;i++)
        {
                if (x!=Cubes[i].x || y!=Cubes[i].y || z!=Cubes[i].z)
                {
                        if (Cubes[i].density!=0)
                        {
                                v+=potential(x, y, z, Cubes[i], CubeSize);
                        }
                }
        }
        return v;
}


Real calcElectrostaticPotential(Real x, Real y, Real z, vector<AtomStruct> &Atoms, vector<CubeStruct> &Cubes, Real CubeSize)
{
        Real vFromAtoms=calcElectrostaticPotential(x, y, z, Atoms);
        Real vFromCubes=calcElectrostaticPotential(x, y, z, Cubes, CubeSize);
        return vFromAtoms+vFromCubes;
}

Real calcElectrostaticPotential(CubeStruct &cube, vector<AtomStruct> &Atoms, vector<CubeStruct> &Cubes, Real CubeSize)
{
        Real r=CubeSize*CubeToSphere;
        Real v=calcElectrostaticPotential(cube.x, cube.y, cube.z, Atoms, Cubes, CubeSize);
        Real selfEnergy=cube.density*2.0*pi*r*r;  //Approximate cube as sphere
        return v+selfEnergy;
}

void calcElectrostaticPotentialEnergyMap(vector<AtomStruct> &Atoms, vector<CubeStruct> &Cubes, Real CubeSize)
{
        int ncube=Cubes.size();

        for (int i=0;i<ncube;i++)
        {
                if (Cubes[i].density!=0)
                {
                        //Cubes[i].ComponentDensity[0]=calcElectrostaticPotential(Cubes[i], Atoms, Cubes, CubeSize); //Hack to avoid allocating more memory.
                }
        }
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

Real calcElectrostaticFreeEnergy(vector<AtomStruct> &Atoms, CubeStruct &Cube, Real CubeSize)
{
        int natom=Atoms.size();
        Real cubeCharge, CubeVolume, ePotential, freeEnergy=0;
        Real x, y, z;

        CubeVolume=CubeSize*CubeSize*CubeSize;
        for (int j=0;j<natom;j++)
        {
                x=Cube.x;
                y=Cube.y;
                z=Cube.z;
                cubeCharge=Cube.density*CubeVolume;  //density is charge density
                ePotential=potential(x, y, z, Atoms[j]);
                freeEnergy+=ePotential*cubeCharge;
        }
        //exit(EXIT_FAILURE);
        return freeEnergy;
}

Real calcElectrostaticFreeEnergyOfAtom(AtomStruct &Atom, lattice &Cubes, Real CubeSize)
{
        int MaxXBin, MaxYBin, MaxZBin;
        Real cubeCharge, CubeVolume, ePotential, freeEnergy=0;
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
                sumFreeEnergy+=FreeEnergy;
        }
        sumFreeEnergy*=0.5;

        return sumFreeEnergy;

}
/**
Real calcElectrostaticFreeEnergy(vector<AtomStruct> &Atoms, string CubeFile, Real CubeSize)
{
        vector<CubeStruct> Cubes;

        ReadNonZeroPdb(CubeFile, Cubes);
        return calcElectrostaticFreeEnergy(Atoms, Cubes, CubeSize);
}
**/
Real calcElectrostaticFreeEnergy(vector<AtomStruct> &Atoms, string CubeFile, Real CubeSize)
{
        fstream file;
        string line;
        Real sumFreeEnergy=0;
        AtomStruct Atom;
        CubeStruct Cube;
        OpenFile(CubeFile, file, "CubeFile");

        getline(file, line);

        while (true)
        {
                if (file.eof()) break;
                if (line.substr(0,4)=="HETA") 
                {
                        PdbLineToAtom(line, Atom);
                        ConvertAtomToCube(Atom, Cube);
                        sumFreeEnergy+=calcElectrostaticFreeEnergy(Atoms, Cube, CubeSize);
                }
                getline(file, line);
        }
        sumFreeEnergy*=0.5;
        
        return sumFreeEnergy;
}
#endif
