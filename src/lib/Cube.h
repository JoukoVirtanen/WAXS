#ifndef _Cubes_included_
#define _Cubes_included_

# include "/home2/jouko/project/HeaderFiles/VectorManip.h"
# include "/home2/jouko/project/HeaderFiles/MemoryUsage.h"

void AssignXYZCoordinates(lattice &Cubes, Real xmin, Real ymin, Real zmin, Real CubeSize)
{
        int MaxXBin, MaxYBin, MaxZBin;
        Real xgrid, ygrid, zgrid;

        Get3DVectorSize(Cubes, MaxXBin, MaxYBin, MaxZBin, "CubesAssingXYZCoordinates");
        cout <<"In AssignXYZCoordinates"<<endl;
        xgrid=xmin-CubeSize*0.5;      //Is this correct?
        cout <<"In AssignXYZCoordinates"<<endl;
        for (int i=0;i<MaxXBin;i++)
        {
                xgrid+=CubeSize;
                cout <<"xgrid= "<<xgrid<<endl;
                ygrid=ymin-CubeSize*0.5;
                for (int j=0;j<MaxYBin;j++)
                {
                        ygrid+=CubeSize;
                        zgrid=zmin-CubeSize*0.5;
                        for (int k=0;k<MaxZBin;k++)
                        {
                                zgrid+=CubeSize;
                                Cubes[i][j][k].x=xgrid;
                                Cubes[i][j][k].y=ygrid;
                                Cubes[i][j][k].z=zgrid;
                        }
                }
        }
}

void GetBins(Real xgrid, Real ygrid, Real zgrid, Real xmin, Real ymin, Real zmin, int &xbin, int &ybin, int &zbin, Real CubeSize)
{
        xbin=int((xgrid-xmin)/CubeSize+0.501);
        ybin=int((ygrid-ymin)/CubeSize+0.501);
        zbin=int((zgrid-zmin)/CubeSize+0.501);
}

void GetBins(CubeStruct &cube, Real xmin, Real ymin, Real zmin, int &xbin, int &ybin, int &zbin, Real CubeSize)
{
        GetBins(cube.x, cube.y, cube.z, xmin, ymin, zmin, xbin, ybin, zbin, CubeSize);
}

void InitializeLattice(lattice &CubeLattice, int MaxXBin, int MaxYBin, int MaxZBin)
{
        CubeStruct cube;
        vector<CubeStruct> CubeVector;
        vector< vector<CubeStruct> > CubeGrid;
        cout <<"In InitializeLattice"<<endl;
        InitializeCube(cube);
        printMemoryUsage();
        SafeAlloc(CubeVector, cube, MaxZBin, "CubeVector");
        cout <<"Allocated CubeVector"<<endl;
        printMemoryUsage();
        SafeAlloc(CubeGrid, CubeVector, MaxYBin, "CubeGrid");
        cout <<"Allocated CubeGrid"<<endl;
        printMemoryUsage();
        SafeAlloc(CubeLattice, CubeGrid, MaxXBin, "CubeLattice");
        cout <<"Initialized lattice"<<endl;
        printMemoryUsage();
}

void InitializeGrid(int MaxXBin, int MaxYBin, int MaxZBin, lattice &Cubes)
{
        CubeStruct Cube;
        vector<CubeStruct> CubeVector;
        vector< vector<CubeStruct> > CubeMatrix;
        InitializeCube(Cube);
        Cube.LocatedIn=InSolution;
        SafeAlloc(CubeVector, Cube, MaxZBin, "CubeVector");
        cout <<"After CubeVector"<<endl;
        SafeAlloc(CubeMatrix, CubeVector, MaxYBin, "CubeMatrix");
        cout <<"After CubeMatrix"<<endl;
        SafeAlloc(Cubes, CubeMatrix, MaxXBin, "Cubes");
        cout <<"After CubeMatrix"<<endl;
}

#endif
