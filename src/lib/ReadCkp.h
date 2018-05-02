#ifndef _ReadDCD_included_
#define _ReadDCD_included_
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

# include "StringUtils.h"
# include "Structures.h"

using namespace std;


void ReadCkp(string CkpFile, vector<AtomStruct> &Atoms, int nthstruct, int &nstructures, Real &XBoxLength, Real &YBoxLength, Real &ZBoxLength)
{
	bool verbose=false;
	ifstream::pos_type size;
	char *title, *memblock, CharCkpFile[1000];
	//int istrt;
	//int interval, zero, nfreat, natom, steps, ndegf, delta, cryst, vernum;
	int natom, ntitle1;
	int cont;
	int position;
        Real max=1000.0;
	vector<int> control;
	float x, y, z;
	int Int, MAX=30000;
        long int LongInt;
        long long int LongLongInt;
	//float Float;
	double Double;
	vector<double> xtal;
        AtomStruct TempAtom;
	size=4;
	cout <<"Reading structure number "<<nthstruct<<" of CkpFile "<<CkpFile<<endl;
	strcpy(CharCkpFile, CkpFile.c_str());
	ifstream file(CharCkpFile, ios::in|ios::binary|ios::ate);
	natom=Atoms.size();
        title=new char[1];
        if (file.is_open())
	{
		size=file.tellg();
                file.seekg(164);
                for (int i=0;i<MAX;i++)
                {
                        file.seekg(i);
                        file.read(reinterpret_cast<char *>(&Int), sizeof(int));
                        file.seekg(i);
                        file.read(reinterpret_cast<char *>(&x), sizeof(float));
                        file.seekg(i);
                        file.read(reinterpret_cast<char *>(&LongInt), sizeof(long int));
                        file.seekg(i);
                        file.read(reinterpret_cast<char *>(&LongLongInt), sizeof(long long int));
			file.seekg(i);
                        file.read(title, 1);
                        cout <<"i= "<<i<<"\t"<<Int<<"\t"<<x<<"\t"<<title[0]<<"\t"<<LongInt<<"\t"<<LongLongInt<<endl;
                }
                file.seekg(164);
                while (true)
                {
                        file.read(reinterpret_cast<char *>(&x), sizeof(float));
                        file.read(reinterpret_cast<char *>(&y), sizeof(float));
                        file.read(reinterpret_cast<char *>(&z), sizeof(float));
                        position=file.tellg();
                        for (int i=0;i<10;i++)
                        {
                                file.read(title, 1);
                                //cout <<title[0];
                        }
                        //cout <<endl;
                        file.seekg(position);
                        cout <<"x= "<<x<<" y= "<<y<<" z= "<<z<<endl;
                        //cout <<"size= "<<Atoms.size()<<endl;
                        cout <<"position= "<<position<<endl;
                        if (abs(x)<max && abs(y)<max && abs(z)<max || (x==0 && y==0 && z==0))
                        {
                                TempAtom.x=x;
                                TempAtom.y=y;
                                TempAtom.z=z;
                                SafePushBack(Atoms, TempAtom, "Atoms in ReadCkp");
                        }
                        else
                        {
                                //break;
                        }
                }
                /*
                for (int i=0;i<size;i++)
                {
                        file.seekg(i);
                        file.read(reinterpret_cast<char *>(&cont), sizeof(int));
                        file.seekg(i);
                        file.read(reinterpret_cast<char *>(&x), sizeof(float));
                        file.seekg(i);
                        file.read(reinterpret_cast<char *>(&Double), sizeof(double));
                        file.seekg(i);
                        cout <<cont<<"\t"<<x<<"\t"<<Double<<"\t"<<i<<"\t"<<size<<endl;
                }
                */
	}
	else cout <<"Unable to open DCD file "<<CkpFile;
}
#endif
