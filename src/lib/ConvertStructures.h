#ifndef _ConvertStructures_included_
#define _ConvertStructures_included_

# include "AtomIDs.h"
# include "Cube.h"
# include "Constants.h"
# include "VectorManip.h"
# include "MemoryUsage.h"
# include "MinMax.h"
//# include "HydrationShell.h"

void ConvertCubeToAtom(CubeStruct &Cube, AtomStruct &Atom)
{
        Atom.AtomName="O";
        Atom.ResidueName="HOH";
        Atom.HetAtom=true;
        Atom.x=Cube.x;
        Atom.y=Cube.y;
        Atom.z=Cube.z;
        Atom.weight=Cube.density;
        Atom.BFactor=Cube.density;
        Atom.atomid=CUBE;
        Atom.AtomNumber=Cube.IntDist;
        Atom.ResidueNum=Cube.AtomType;
        if (Cube.concave) Atom.Occupancy=1.0;
        else Atom.Occupancy=0.0;
}

void ConvertCubesToAtoms(vector<CubeStruct> &Cubes, vector<AtomStruct> &Atoms)
{
        int i, ncubes=Cubes.size();
        AtomStruct Atom;
        Atom.AtomName="O";
        Atom.ResidueName="HOH";
        Atom.HetAtom=true;
        DeleteVector(Atoms);
        for (i=0;i<ncubes;i++)
        {
                Atom.x=Cubes[i].x;
                Atom.y=Cubes[i].y;
                Atom.z=Cubes[i].z;
                Atom.weight=Cubes[i].density;
                Atom.atomid=CUBE;
                Atom.AtomNumber=Cubes[i].IntDist;
                Atom.ResidueNum=Cubes[i].AtomType;
                SafePushBack(Atoms, Atom, "Atoms in ConvertCubesToAtoms");
        }
}


void LatticeDimensions(vector<CubeStruct> &cubes, Real &xmin, Real &ymin, Real &zmin, int &MaxXBin, int &MaxYBin, int &MaxZBin, Real CubeSize)
{
        Real xmax, ymax, zmax;
        cout <<"In LatticeDimensions"<<endl;
        MinMax(xmin, ymin, zmin, xmax, ymax, zmax, cubes);
        GetBins(xmax, ymax, zmax, xmin, ymin, zmin, MaxXBin, MaxYBin, MaxZBin, CubeSize);
        MaxXBin++;
        MaxYBin++;
        MaxZBin++;
}

void LatticeDimensions(Real xmin, Real ymin, Real zmin, Real xmax, Real ymax, Real zmax, int &MaxXBin, int &MaxYBin, int &MaxZBin, Real CubeSize)
{
        GetBins(xmax, ymax, zmax, xmin, ymin, zmin, MaxXBin, MaxYBin, MaxZBin, CubeSize);
        MaxXBin++;
        MaxYBin++;
        MaxZBin++;
}

void GetLatticeDimensions(lattice &Cubes, int &MaxXBin, int &MaxYBin, int &MaxZBin)
{
        MaxXBin=Cubes.size();
        if (MaxXBin!=0)
        {
                MaxYBin=Cubes[0].size();
                if (MaxYBin!=0)
                {
                        MaxZBin=Cubes[0][0].size();
                }
                else
                {
                        cout <<"ERROR: Lattice dimension=0"<<endl;
                        exit(EXIT_FAILURE);
                }
        }
        else
        {
                cout <<"ERROR: Lattice dimensions=0"<<endl;
                exit(EXIT_FAILURE);
        }
}

void SetCoordinates(lattice &CubeLattice, Real xmin, Real ymin, Real zmin, Real CubeSize)
{
        int MaxXBin, MaxYBin, MaxZBin;
        Real xgrid, ygrid, zgrid;
        MaxXBin=CubeLattice.size();
        MaxYBin=CubeLattice[0].size();
        MaxZBin=CubeLattice[0][0].size();
        xgrid=xmin-CubeSize;
        cout <<"In SetCoordinates"<<endl;;
        for (int j=0;j<MaxXBin;j++)
        {
                //cout <<"MaxXBin= "<<MaxXBin<<" MaxYBin= "<<MaxYBin<<" MaxZBin= "<<MaxZBin<<endl;
                xgrid+=CubeSize;
                ygrid=ymin-CubeSize;
                for (int k=0;k<MaxYBin;k++)
                {
                        //cout <<"MaxXBin= "<<MaxXBin<<" MaxYBin= "<<MaxYBin<<" MaxZBin= "<<MaxZBin<<endl;
                        ygrid+=CubeSize;
                        zgrid=zmin-CubeSize;
                        for (int l=0;l<MaxZBin;l++)
                        {
                                //cout <<"MaxXBin= "<<MaxXBin<<" MaxYBin= "<<MaxYBin<<" MaxZBin= "<<MaxZBin<<endl;
                                //cout <<"j= "<<j<<" k= "<<k<<" l= "<<l<<endl;
                                zgrid+=CubeSize;
                                CubeLattice[j][k][l].x=xgrid;
                                CubeLattice[j][k][l].y=ygrid;
                                CubeLattice[j][k][l].z=zgrid;
                                CubeLattice[j][k][l].LocatedIn=InSolution;
                        }
                }
        }
        cout <<"Leaving SetCoordinates"<<endl;
}

void ConvertCubeVectorToLattice(vector<CubeStruct> &cubes, lattice &CubeLattice, Real CubeSize)
{
        int NumCubes=cubes.size();
        int MaxXBin, MaxYBin, MaxZBin;
        int xbin, ybin, zbin;
        Real xmin, ymin, zmin;
        cout <<"In ConvertCubeVerctorToLattice"<<endl;
        LatticeDimensions(cubes, xmin, ymin, zmin, MaxXBin, MaxYBin, MaxZBin, CubeSize);
        InitializeLattice(CubeLattice, MaxXBin, MaxYBin, MaxZBin);
        cout <<"xmin= "<<xmin<<" ymin= "<<ymin<<" zmin= "<<zmin<<endl;
        SetCoordinates(CubeLattice, xmin, ymin, zmin, CubeSize);
        cout <<"After SetCoordinates.  About to enter Print3DVectorSize"<<endl;
        Print3DVectorSize(CubeLattice, "CubeLattice");
        cout <<"MaxXBin= "<<MaxXBin<<" MaxYBin= "<<MaxYBin<<" MaxZBin= "<<MaxZBin<<endl;
        //CubesToPdb(cubes, "/home2/jouko/project/WAXS/test/pdb/cubes.pdb");
        for (int j=0;j<MaxXBin;j++)
        {
                for (int k=0;k<MaxYBin;k++)
                {
                        for (int l=0;l<MaxZBin;l++)
                        {
                                CubeLattice[j][k][l].LocatedIn=InSolution;
                        }
                }
        }
        //cout <<"MaxXBin*MaxYBin*MaxZBin= "<<MaxXBin*MaxYBin*MaxZBin<<" NumCubes= "<<NumCubes<<endl;
        for (int n=NumCubes-1;n>=0;n--)
        {
                GetBins(cubes[n], xmin, ymin, zmin, xbin, ybin, zbin, CubeSize);
                //PrintCubeInfo(cubes[n]);
                //cout <<"xbin= "<<xbin<<" ybin= "<<ybin<<" zbin= "<<zbin<<endl;
                //cout <<endl;
                CubeLattice[xbin][ybin][zbin]=cubes[n];
                //cubes.pop_back();
        }
        for (int j=0;j<MaxXBin;j++)
        {
                for (int k=0;k<MaxYBin;k++)
                {
                        for (int l=0;l<MaxZBin;l++)
                        {
                                //cout <<"j= "<<j<<" MaxXBin= "<<MaxXBin<<endl;
                                //cout <<"k= "<<k<<" MaxYBin= "<<MaxYBin<<endl;
                                //cout <<"l= "<<l<<" MaxZBin= "<<MaxZBin<<endl;
                                if (CubeLattice[j][k][l].x==0 && CubeLattice[j][k][l].y==0 && CubeLattice[j][k][l].z==0)
                                {
                                        cout <<"j= "<<j<<" k= "<<k<<" l= "<<l<<endl;
                                        PrintCubeInfo(CubeLattice[j][k][l]);
                                }
                        }
                }
        }
        //DensityMapToPdb(CubeLattice, "/home2/jouko/project/WAXS/test/pdb/ConvertCubeVectorToLattice.pdb");
}


void ConvertAtomToCube(AtomStruct &Atoms, CubeStruct &Cubes)
{
        CubeStruct Cube;
        InitializeCube(Cube);
        Cubes.x=Atoms.x;
        Cubes.y=Atoms.y;
        Cubes.z=Atoms.z;
        Cubes.density=Atoms.BFactor;
        Cubes.AtomType=Atoms.ResidueNum;
        Cubes.IntDist=Atoms.AtomNumber;
}

void ConvertAtomsToCubes(vector<AtomStruct> &Atoms, vector<CubeStruct> &Cubes)
{
        int i, natom=Atoms.size();
        CubeStruct Cube;
        InitializeCube(Cube);
        SafeAlloc(Cubes, Cube, natom, "Cubes in Convert AtomsToCubes");
        for (i=0;i<natom;i++)
        {
                Cubes[i].x=Atoms[i].x;
                Cubes[i].y=Atoms[i].y;
                Cubes[i].z=Atoms[i].z;
                Cubes[i].density=Atoms[i].BFactor;
                Cubes[i].AtomType=Atoms[i].ResidueNum;
                Cubes[i].IntDist=Atoms[i].AtomNumber;
        }
}

void ReadPdb(string PdbFile, vector<CubeStruct> &Cubes)
{
        vector<AtomStruct> Atoms;
        ReadPdb(PdbFile, Atoms);
        RemoveNonHetatm(Atoms);
        ConvertAtomsToCubes(Atoms, Cubes);
}

void ReadNonZeroPdb(string PdbFile, vector<CubeStruct> &Cubes)
{
        vector<AtomStruct> Atoms;
        ReadNonZeroPdb(PdbFile, Atoms);
        RemoveNonHetatm(Atoms);
        ConvertAtomsToCubes(Atoms, Cubes);
}

void MinMaxPdb(string CubeFile, Real &xmin, Real &ymin, Real &zmin, Real &xmax, Real &ymax, Real &zmax)
{
        string line, ElementID, StrResidueID, resstrm, str1, str2, str3, StrAtomName;
        Real x, y, z;
        fstream protein;
        OpenFile(CubeFile, protein, "CubeFile");
        getline(protein, line);
        xmin=10e10;
        ymin=10e10;
        zmin=10e10;
        xmax=-10e10;
        ymax=-10e10;
        zmax=-10e10;
        while (true)
        {
                if (line.substr(0,3)=="END") break;
                if (protein.eof()) break;		
                if (line.substr(0,4)=="HETA")
                {
                        x=StrToFloat(line.substr(30,8));
                        y=StrToFloat(line.substr(38,8));
                        z=StrToFloat(line.substr(46,8));
                        if (x<xmin) xmin=x;
                        if (y<ymin) ymin=y;
                        if (z<zmin) zmin=z;
                        if (x>xmax) xmax=x;
                        if (y>ymax) ymax=y;
                        if (z>zmax) zmax=z;
                }

                getline(protein, line);
        }
        protein.close();
}

void ReadPdbCube(string PdbFile, lattice &Cubes, Real CubeSize)
{
        fstream file;
        string line;
        int xbin, ybin, zbin;
        int MaxXBin, MaxYBin, MaxZBin;
        Real density_temp;
        Real x, y, z;
        Real xmax, ymax, zmax;
        Real xmin, ymin, zmin;
        
        cout <<"In ReadPdbCube"<<endl;
        MinMaxPdb(PdbFile, xmin, ymin, zmin, xmax, ymax, zmax);
        cout <<"After MinMax"<<endl;
        GetBins(xmax, ymax, zmax, xmin, ymin, zmin, MaxXBin, MaxYBin, MaxZBin, CubeSize);
        cout <<"After GetBins"<<endl;
        InitializeGrid(MaxXBin+1, MaxYBin+1, MaxZBin+1, Cubes);
        cout <<"After InitializeGrid"<<endl;
        AssignXYZCoordinates(Cubes, xmin+CubeSize*0.5, ymin+CubeSize*0.5, zmin+CubeSize*0.5, CubeSize);
        cout <<"xmin= "<<xmin<<" ymin= "<<ymin<<" zmin= "<<zmin<<" xmax= "<<xmax<<" ymax= "<<ymax<<" zmax= "<<zmax<<endl;
        cout <<"MaxXBin= "<<MaxXBin<<" MaxYBin= "<<MaxYBin<<" MaxZBin= "<<MaxZBin<<endl;
        cout <<"CubeSize= "<<CubeSize<<endl;
        OpenFile(PdbFile, file, "PdbFile in ReadPdbCube");
        getline(file, line);
        while (true)
        {
                if (line.substr(0,3)=="END") break;
                if (file.eof()) break;		
                //if ((line.substr(0,4)=="HETA" || line.substr(0,4)=="ATOM") && line.substr(17,3)=="HOH")
                if ((line.substr(0,4)=="HETA") && line.substr(17,3)=="HOH")
                {
                        
                        x=StrToFloat(line.substr(30,8));
                        y=StrToFloat(line.substr(38,8));
                        z=StrToFloat(line.substr(46,8));
                        density_temp=StrToFloat(line.substr(61,10));
                        GetBins(x, y, z, xmin, ymin, zmin, xbin, ybin, zbin, CubeSize);
                        //cout <<"xbin= "<<xbin<<" ybin= "<<ybin<<" zbin= "<<zbin<<endl;
                        if (xbin<0 || xbin>=MaxXBin+1 || ybin<0 || ybin>=MaxYBin+1 || zbin<0 || zbin>=MaxZBin+1)
                        {
                                cout <<"xmin= "<<xmin<<" ymin= "<<ymin<<" zmin= "<<zmin<<" xmax= "<<xmax<<" ymax= "<<ymax<<" zmax= "<<zmax<<endl;
                                cout <<"MaxXBin= "<<MaxXBin<<" MaxYBin= "<<MaxYBin<<" MaxZBin= "<<MaxZBin<<endl;
                                cout <<"xbin= "<<xbin<<" ybin= "<<ybin<<" zbin= "<<zbin<<endl;
                                cout <<"ERROR: Grid point outside of box"<<endl;
                                exit(EXIT_FAILURE);
                        }
                        Cubes[xbin][ybin][zbin].IntDist=StrToInt(line.substr(7,5));
                        Cubes[xbin][ybin][zbin].x=x;
                        Cubes[xbin][ybin][zbin].y=y;
                        Cubes[xbin][ybin][zbin].z=z;
                        Cubes[xbin][ybin][zbin].density=density_temp;
                        if (density_temp==0) Cubes[xbin][ybin][zbin].LocatedIn=InProtein;
                        else Cubes[xbin][ybin][zbin].LocatedIn=InHydrationShell;
                        //if (density_temp==0) PrintPdbLine("/home2/jouko/Density/TempPdb/ZeroDensity.pdb", true, 0, "O", "HOH", 0, Cubes[xbin][ybin][zbin].x, Cubes[xbin][ybin][zbin].y, Cubes[xbin][ybin][zbin].z, 1.0, 1.0);
                }
                getline(file, line);
        }
        file.close();
        cout <<"Finished ReadingPdbCube"<<endl;
}
/*
void ReadPdb(string PdbFile, lattice &Cubes, Real CubeSize)
{
        vector<AtomStruct> Atoms;
        vector<CubeStruct> CubeVector;
        ReadPdb(PdbFile, Atoms);
        RemoveNonHetatm(Atoms);
        ConvertAtomsToCubes(Atoms, CubeVector);
        ConvertCubeVectorToLattice(CubeVector, Cubes, CubeSize);
}
*/
int CountNonZeroCubes(lattice &CubeLattice)
{
        int ncube=0;
        int MaxXBin, MaxYBin, MaxZBin;

        MaxXBin=CubeLattice.size();
        MaxYBin=CubeLattice[0].size();
        MaxZBin=CubeLattice[0][0].size();

        for (int j=0;j<MaxXBin;j++)
        {
                for (int k=0;k<MaxYBin;k++)
                {
                        for (int l=0;l<MaxZBin;l++)
                        {
                                if (CubeLattice[j][k][l].density!=0) ncube++;
                        }
                }
        }
        return ncube;
}

void ConvertCubeLatticeToCubeVector(lattice &CubeLattice, vector<CubeStruct> &cubes)
{
        int ncube=0, NonZeroCubes;
        int MaxXBin, MaxYBin, MaxZBin;
        CubeStruct cube;

        InitializeCube(cube);
        MaxXBin=CubeLattice.size();
        MaxYBin=CubeLattice[0].size();
        MaxZBin=CubeLattice[0][0].size();
        cubes.clear();
        NonZeroCubes=CountNonZeroCubes(CubeLattice);
        SafeAlloc(cubes, cube, MaxXBin*MaxYBin*MaxZBin, "cubes");
        cout <<"In ConvertCubeLatticeToCubeVector"<<endl;
        cout <<"MaxXBin= "<<MaxXBin<<" MaxYBin= "<<MaxYBin<<" MaxZBin= "<<MaxZBin<<endl;
        for (int j=0;j<MaxXBin;j++)
        {
                //cout <<"MaxXBin= "<<MaxXBin<<" j= "<<j<<endl;
                for (int k=0;k<MaxYBin;k++)
                {
                        for (int l=0;l<MaxZBin;l++)
                        {
                                if (CubeLattice[j][k][l].density!=InSolution)
                                {
                                        cubes[ncube]=CubeLattice[j][k][l];
                                        //cout <<"ncube= "<<ncube<<" NonZeroCubes= "<<NonZeroCubes<<endl;
                                        ncube++;
                                }
                        }
                }
        }
        cout <<"At end of ConvertCubeLatticeToCubeVector"<<endl;
}

#endif
