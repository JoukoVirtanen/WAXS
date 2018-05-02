#ifndef _ReadGrasp_included_
#define _ReadGrasp_included_

# include "/home2/jouko/project/HeaderFiles/VectorManip.h"
# include "/home2/jouko/project/HeaderFiles/Cube.h"

const Real KT_TO_KCAL=0.550756879;

void ReadGrasp(string phiMap, Real xmid, Real ymid, Real zmid, Real cubeSize, lattice &Cubes)
{
        float Float;
        double Double;
        Real xmin, ymin, zmin;
        char *title;
        char charPhiMap[1000];
        int Int;
        int MaxXBin, MaxYBin, MaxZBin;
        int size, max=1000000;
        Real SCALE_TO_KCAL=332.0634;
        CubeStruct Cube;

        strcpy(charPhiMap, phiMap.c_str());
        ifstream file(charPhiMap, ios::in|ios::binary|ios::ate);

        MaxXBin=MaxYBin=MaxZBin=211; 
        xmin=xmid-Real(MaxXBin)*cubeSize*0.5;
        ymin=ymid-Real(MaxYBin)*cubeSize*0.5;
        zmin=zmid-Real(MaxZBin)*cubeSize*0.5;
        InitializeLattice(Cubes,  MaxXBin, MaxYBin, MaxZBin);
        cout <<"After InitializeLattice"<<endl;
        AssignXYZCoordinates(Cubes, xmin, ymin, zmin, cubeSize);
        cout <<"After AssignXYZCoordinates"<<endl;
        title=new char[1];
        if (file.is_open())
        {
                size=4;
                file.seekg(size, ios::beg);
                for (int i=0;i<max;i++)
                {
                        file.seekg(i);
                        file.read(reinterpret_cast<char *>(&Double), sizeof(Double));
                        file.seekg(i);
                        file.read(reinterpret_cast<char *>(&Float), sizeof(Float));
                        file.seekg(i);
                        file.read(reinterpret_cast<char *>(&Int), sizeof(Int));
                        file.seekg(i);
                        file.read(title, 1);
                        file.seekg(i);
                        //cout <<Double<<"\t"<<Float<<"\t"<<Int<<"\t"<<title[0]<<endl;
                }
                file.seekg(0);
                file.read(reinterpret_cast<char *>(&Int), sizeof(Int));
                //cout <<Int<<endl;
                for (int i=0;i<20;i++)
                {
                        file.read(title, 1);
                        //cout <<title[0]<<endl;
                }
                file.read(reinterpret_cast<char *>(&Int), sizeof(Int));
                //cout <<Int<<endl;
                file.read(reinterpret_cast<char *>(&Int), sizeof(Int));
                //cout <<Int<<endl;
                for (int i=0;i<70;i++)
                {
                        file.read(title, 1);
                        //cout <<title[0]<<endl;
                }
                file.read(reinterpret_cast<char *>(&Int), sizeof(Int));
                //cout <<Int<<endl;
                file.read(reinterpret_cast<char *>(&Double), sizeof(Double));
                //cout <<Double<<endl;
                for (int i=0;i<MaxXBin;i++)
                {
                        cout <<"i= "<<i<<endl;
                        for (int j=0;j<MaxYBin;j++)
                        {
                                for (int k=0;k<MaxZBin;k++)
                                {
                                        file.read(reinterpret_cast<char *>(&Float), sizeof(Float));
                                        //Cubes[k][j][i].density=Real(Float)*KT_TO_KCAL/SCALE_TO_KCAL;
                                        //Cubes[k][j][i].density=Real(Float)*KT_TO_KCAL;
                                        Cubes[k][j][i].density=Real(Float);
                                        //cout <<Cubes[i][j][k].density<<endl;
                                }
                        }
                }
                /*
                for (int i=0;i<max;i++)
                {
                        file.read(reinterpret_cast<char *>(&Float), sizeof(Float));
                        cout <<Float<<endl;
                }
                */
        }
        else
        {
                cout <<"ERROR: Unable to open phi file"<<endl;
                exit(EXIT_FAILURE);
        }
        cout <<"Read GraspFile"<<endl;
}

#endif
