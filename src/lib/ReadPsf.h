#ifndef _ReadPSF_included_
#define _ReadPSF_included_
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

# include "Structures.h"
# include "StringUtils.h"
# include "IOUtils.h"

using namespace std;

void ReadPsf(string PsfFile, vector<AtomStruct> &Atoms)
{
	AtomStruct Atom;
	bool verbose=true;
	string line;
	size_t found;
	int i, natom, nbond;
	int BondsInLine, atom1, atom2;
	int AtomNumber=0;
	int SegID=9;
	int ResidueNum=14;
	int ResidueName=19;
	int AtomName=24;
	int CharmmAtomName=29;
	int Charge=35;
	int Mass=44;
	int Last=69;

	fstream file;
	OpenFile(PsfFile, file, "Psf file");
        Atoms.clear();	
	if (verbose) cout <<"Opened file"<<endl;
	if (verbose) cout <<"file.open"<<endl;
	getline(file, line);
	if (verbose) cout <<"line= "<<line<<endl;;	
	while (true)
	{
		//cout <<"In while true"<<endl;
		if (file.eof()) break;
		//cout <<"After if (file.eof())"<<endl;
		found=line.find("!NATOM");
		if (found!=string::npos)
		{
			natom=StrToInt(line.substr(0,found));
			getline(file, line);
			if (verbose) cout <<"natom= "<<natom<<endl;
			for (i=0;i<natom;i++)
			{
				if (file.eof()) break;
				Atom.AtomNumber=StrToInt(line.substr(AtomNumber, SegID-AtomNumber-1));
				Atom.SegID=line.substr(SegID, ResidueNum-SegID-1);
				Atom.ResidueNum=StrToInt(line.substr(ResidueNum, ResidueName-ResidueNum-1));
				Atom.ResidueName=line.substr(ResidueName, AtomName-ResidueName-1);
				Atom.AtomName=line.substr(AtomName, CharmmAtomName-AtomName-1);
				Atom.CharmmAtomName=line.substr(CharmmAtomName, Charge-CharmmAtomName-1);
				Atom.charge=StrToFloat(line.substr(Charge, Mass-Charge-1));
				Atom.mass=StrToFloat(line.substr(Mass, Last-Mass-1));
				Trim(Atom.ChainName);
				Trim(Atom.ResidueName);
				Trim(Atom.AtomName);
				Trim(Atom.CharmmAtomName);
				Atoms.push_back(Atom);
				getline(file, line);
			}
		}
		found=line.find("!NBOND:");
		if (found!=string::npos)
		{
			nbond=StrToInt(line.substr(0,found));
			getline(file, line);
			if (verbose) cout <<"nbond= "<<nbond<<endl;
			while (true)
			{
				if (file.eof()) break;
				if (line=="") break;
				if (line.length()<8) break;
				if (line.substr(10,7)=="!NTHETA:") break;
				BondsInLine=line.length()/16;
				for (i=0;i<BondsInLine;i++)
				{
					atom1=StrToInt(line.substr(i*16,8))-1;
					atom2=StrToInt(line.substr(i*16+8,8))-1;
					//cout <<"atom1= "<<atom1<<" atom2= "<<atom2<<endl;
					Atoms[atom1].connectivity.push_back(atom2);
					Atoms[atom2].connectivity.push_back(atom1);
                                }
				//cout <<endl;
				getline(file, line);
			}
		}
		getline(file, line);
	}
	file.close();
/*************
	for (i=0;i<1;i++)
	{
		cout <<Atoms[i].AtomNumber<<"\t"<<endl
		<<Atoms[i].ChainName<<"\t"<<endl
		<<Atoms[i].ResidueNum<<"\t"<<endl
		<<Atoms[i].ResidueName<<endl
		<<Atoms[i].AtomName<<endl
		<<Atoms[i].CharmmAtomName<<endl
		<<Atoms[i].charge<<endl
		<<Atoms[i].mass<<endl;
	}
*************/

}
#endif