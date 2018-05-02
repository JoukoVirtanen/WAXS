#ifndef _ReadPdb_included_
#define _ReadPdb_included_

# include "AssignAtomIDs.h"
# include "AtomCode.h"
# include "GetAtomType.h"
# include "GetAtomType2.h"
# include "ResidueCode.h"
# include "ResidueID.h"
# include "Structures.h"
# include "StringUtils.h"
# include "TypeDef.h"
# include "LinkedList.h"
# include "IOUtils.h"

void ReadPdb(string PdbFile, ProteinStruct &Protein)
{
	AtomStruct Atom;
	//bool verbose=false;
	string line, ElementID, StrResidueID, StrAtomName;
	int natom;
	fstream file;
	
	DeleteVector(Protein.Atoms);
	natom=0;
	OpenFile(PdbFile, file, "pdb file");
	cout <<"Reading pdb file "<<PdbFile<<endl;
	getline(file, line);

	while (true)
	{
		line=line+"                                         ";

		if (file.eof()) break;
		if (line.substr(0,4)=="CRYS")
		{
			Protein.XBoxLength=StrToFloat(GetWord(line,2));
			Protein.YBoxLength=StrToFloat(GetWord(line,3));
			Protein.ZBoxLength=StrToFloat(GetWord(line,4));
		}

                if (line.substr(0,3)=="END") break;

                //cout <<line.substr(0,4)<<endl;    //Atom or Hetatom


                if (line.substr(0,4)=="ATOM" || line.substr(0,4)=="HETA")
                {
                        Atom.AtomNumber=StrToInt(line.substr(7,5));
                        Atom.AtomName=line.substr(12,4);
                        Atom.ResidueName=line.substr(17,4);
                        Atom.ChainName=line.substr(21,1);
                        Atom.ResidueNum=StrToInt(line.substr(22,4));
                        Atom.x=StrToFloat(line.substr(31,7));
                        Atom.y=StrToFloat(line.substr(39,7));
                        Atom.z=StrToFloat(line.substr(47,7));
                        Atom.Occupancy=StrToFloat(line.substr(56,4));
                        Atom.BFactor=StrToFloat(line.substr(61,8));
                        Atom.SegID=line.substr(72,4);
                        Atom.ID=line.substr(73,3);
			Atom.weight=1.0;

			Trim(Atom.AtomName);
			Trim(Atom.ResidueName);
			Atom.residueid=GetResidueID(Atom.ResidueName);

			Protein.Atoms.push_back(Atom);
			natom++;
                }
		getline(file, line);
        }
	file.close();
	cout <<"natom= "<<natom<<endl;
}

int ReadPdb(string PdbFile, vector<AtomStruct> &Atoms)
{
	AtomStruct Atom;
	bool valid;
	//bool verbose=false;
	char CharPdbFile[100];
	string line, ElementID, StrResidueID, StrAtomName;
	int natom;
	Real XBoxLength, YBoxLength, ZBoxLength;
	DeleteVector(Atoms);
	natom=0;
	strcpy(CharPdbFile, PdbFile.c_str());
	valid=exist(CharPdbFile);
	if (!valid)
	{
		cout <<"Error could not open pdb file "<<CharPdbFile<<endl;
		exit(0);
	}

	fstream protein;
	protein.open(CharPdbFile, ios::in);

	getline(protein,line);

	while (true)
	{
		line=line+"                                         ";

		if (protein.eof()) break;
		if (line.substr(0,4)=="CRYS")
		{
			XBoxLength=StrToFloat(GetWord(line,2));
			YBoxLength=StrToFloat(GetWord(line,3));
			ZBoxLength=StrToFloat(GetWord(line,4));
		}

                if (line.substr(0,3)=="END") break;

                //cout <<line.substr(0,4)<<endl;    //Atom or Hetatom


                if (line.substr(0,4)=="ATOM" || line.substr(0,4)=="HETA")
                {
                        Atom.AtomNumber=StrToInt(line.substr(7,5));
                        Atom.AtomName=line.substr(12,4);
                        Atom.ResidueName=line.substr(17,4);
                        Atom.ChainName=line.substr(21,1);
                        Atom.ResidueNum=StrToInt(line.substr(22,4));
                        Atom.x=StrToFloat(line.substr(31,7));
                        Atom.y=StrToFloat(line.substr(39,7));
                        Atom.z=StrToFloat(line.substr(47,7));
                        Atom.Occupancy=StrToFloat(line.substr(56,4));
                        Atom.BFactor=StrToFloat(line.substr(61,8));
                        Atom.SegID=line.substr(72,4);
                        Atom.ID=line.substr(73,3);
			Atom.weight=1.0;

			Trim(Atom.AtomName);
			Trim(Atom.ResidueName);
			Atom.residueid=GetResidueID(Atom.ResidueName);
			Atoms.push_back(Atom);
			natom++;
                }
		getline(protein,line);
        }
	protein.close();
	cout <<"natom= "<<natom<<endl;
	return natom;
}

void ReadXYZ(string xyz, vector<AtomStruct> &Atoms)
{
	AtomStruct Atom;
	bool valid;
	char CharXYZ[100];
	string line;
	int n, natom;
	strcpy(CharXYZ, xyz.c_str());
	valid=exist(CharXYZ);
	if (!valid)
	{
		cout <<"Error unable to open xyz file "<<CharXYZ<<endl;
		exit(0);
	}

	int AtomNumber=6;
	int AtomName=12;
	int XPos=24;
	int YPos=36;
	int ZPos=48;
	int AtomType=54;
	int Connect1=60;
	int Connect2=66;
	int Connect3=72;
	int Connect4=78;
	
	fstream protein;
	protein.open(CharXYZ, ios::in);

	getline(protein, line);
	natom=StrToInt(GetWord2(line, 1));
	for (n=0;n<natom;n++)
	{
		if (protein.eof()) break;
		Atom.AtomNumber=StrToInt(line.substr(0, AtomNumber));
		Atom.AtomName=line.substr(AtomNumber, XPos-AtomName);
		Atom.x=StrToFloat(line.substr(AtomName, YPos-XPos));
		Atom.y=StrToFloat(line.substr(XPos, ZPos-XPos));
		Atom.z=StrToFloat(line.substr(ZPos, AtomType-ZPos));
		Atom.ParmType=StrToInt(line.substr(Connect1, Connect2-Connect1));
		Atom.connectivity.push_back(StrToInt(line.substr(Connect2, Connect3-Connect2)));
		Atom.connectivity.push_back(StrToInt(line.substr(Connect3, Connect4-Connect3)));
		Atom.connectivity.push_back(StrToInt(line.substr(Connect4, 100)));
		Atom.connectivity.push_back(StrToInt(line.substr(Connect4, 100)));
		Atoms.push_back(Atom);
		getline(protein, line);
	}
}

void ReadXYZLine(string line, AtomStruct Atom)
{
	int AtomNumber=6;
	int AtomName=12;
	int XPos=24;
	int YPos=36;
	int ZPos=48;
	int AtomType=54;
	int Connect1=60;
	int Connect2=66;
	int Connect3=72;
	int Connect4=78;

	Atom.AtomNumber=StrToInt(line.substr(0, AtomNumber));
	Atom.AtomName=line.substr(AtomNumber, XPos-AtomName);
	Atom.x=StrToFloat(line.substr(AtomName, YPos-XPos));
	Atom.y=StrToFloat(line.substr(XPos, ZPos-XPos));
	Atom.z=StrToFloat(line.substr(ZPos, AtomType-ZPos));
	Atom.ParmType=StrToInt(line.substr(Connect1, Connect2-Connect1));
	Atom.connectivity.push_back(StrToInt(line.substr(Connect2, Connect3-Connect2)));
	Atom.connectivity.push_back(StrToInt(line.substr(Connect3, Connect4-Connect3)));
	Atom.connectivity.push_back(StrToInt(line.substr(Connect4, 100)));
	Atom.connectivity.push_back(StrToInt(line.substr(Connect4, 100)));
}
//#else
//#error ReadPdb is already defined
#endif
