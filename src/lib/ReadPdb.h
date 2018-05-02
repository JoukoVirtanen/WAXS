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

int CountAtomsInPdb(string PdbFile)
{
        //bool verbose=false;
        string line, ElementID, StrResidueID, StrAtomName;
        int natom;
        natom=0;

        fstream protein;
        OpenFile(PdbFile, protein, "PdbFile");

        getline(protein,line);

        while (true)
        {

                if (protein.eof()) break;
                if (line.substr(0,3)=="END") break;


                if (line.substr(0,4)=="ATOM" || line.substr(0,4)=="HETA")
                {
                        natom++;
                }
                getline(protein,line);
        }
        protein.close();
        return natom;
}

int CountAtomsInPdb(vector<string> &lines)
{
        bool valid;
        //bool verbose=false;
        char CharPdbFile[1000];
        string line, ElementID, StrResidueID, StrAtomName;
        int natom, nlines;
        nlines=lines.size();
        natom=0;

        for (int i=0;i<nlines;i++)
        {
                if (lines[i].substr(0,4)=="ATOM" || lines[i].substr(0,4)=="HETA")
                {
                        natom++;
                }
        }
        return natom;
}

int CountNonZeroAtomsInPdb(string PdbFile)
{
        bool valid;
        //bool verbose=false;
        char CharPdbFile[1000];
        string line, ElementID, StrResidueID, StrAtomName;
        Real BFactor;
        int natom;
        natom=0;
        strcpy(CharPdbFile, PdbFile.c_str());
        valid=exist(CharPdbFile);
        if (!valid)
        {
                cout <<"Error could not open pdb file "<<CharPdbFile<<endl;
                exit(EXIT_FAILURE);
        }

        fstream protein;
        protein.open(CharPdbFile, ios::in);

        getline(protein,line);

        while (true)
        {

                if (protein.eof()) break;
                if (line.substr(0,3)=="END") break;


                BFactor=StrToFloat(line.substr(61,8));
                if ((line.substr(0,4)=="ATOM" || line.substr(0,4)=="HETA") && BFactor!=0)
                {
                        natom++;
                }
                getline(protein,line);
        }
        protein.close();
        return natom;
}

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
                        if (line.substr(0,4)=="ATOM") Atom.HetAtom=false;
                        else if (line.substr(0,4)=="HETA") Atom.HetAtom=true;
                        Atom.AtomNumber=StrToInt(line.substr(7,5));
                        Atom.AtomName=line.substr(12,4);
                        Atom.duplicate=line.substr(16,1);
                        Atom.ResidueName=line.substr(17,4);
                        Atom.ChainName=line.substr(21,1);
                        Atom.ResidueNum=StrToInt(line.substr(22,4));
                        Atom.x=StrToFloat(line.substr(30,8));
                        Atom.y=StrToFloat(line.substr(38,8));
                        Atom.z=StrToFloat(line.substr(46,8));
                        Atom.Occupancy=StrToFloat(line.substr(56,4));
                        Atom.BFactor=StrToFloat(line.substr(61,8));
                        Atom.SegID=line.substr(72,4);
                        Atom.ID=line.substr(73,3);
                        Atom.weight=1.0;

                        Trim(Atom.AtomName);
                        Trim(Atom.ResidueName);
                        DeleteStartingAt(Atom.ResidueName, " ");
                        Atom.residueid=GetResidueID(Atom.ResidueName);

                        SafePushBack(Protein.Atoms, Atom, "Protein.Atoms");
                        natom++;
                }
                getline(file, line);
        }
        file.close();
}

void PdbLineToAtom(string line, AtomStruct &Atom)
{
        if (line.substr(0,4)=="ATOM") Atom.HetAtom=false;
        else if (line.substr(0,4)=="HETA") Atom.HetAtom=true;
        Atom.AtomNumber=StrToInt(line.substr(7,5));
        Atom.AtomName=line.substr(12,4);
        Atom.duplicate=line.substr(16,1);
        Atom.ResidueName=line.substr(17,4);
        Atom.ChainName=line.substr(21,1);
        Atom.ResidueNum=StrToInt(line.substr(22,4));
        Atom.x=StrToFloat(line.substr(30,8));
        Atom.y=StrToFloat(line.substr(38,8));
        Atom.z=StrToFloat(line.substr(46,8));
        Atom.Occupancy=StrToFloat(line.substr(56,4));
        Atom.BFactor=StrToFloat(line.substr(61,8));
        Atom.SegID=line.substr(72,4);
        Atom.ID=line.substr(73,3);
        Atom.weight=1.0;
        //cout <<"Atoms[natom].AtomName= "<<Atoms[natom].AtomName<<"."<<endl;
        Trim(Atom.AtomName);
        //cout <<"Atoms[natom].AtomName= "<<Atoms[natom].AtomName<<"."<<endl;
        Trim(Atom.ResidueName);
        DeleteStartingAt(Atom.ResidueName, " ");
        Atom.residueid=GetResidueID(Atom.ResidueName);

}

void PdbLinesToAtoms(vector<string> &lines, vector<AtomStruct> &Atoms)
{
        int natom;
        AtomStruct Atom;
        natom=CountAtomsInPdb(lines);
        SafeAlloc(Atoms, Atom, natom, "in LinesToAtoms");
        for (int i=0;i<natom;i++)
        {
                if (lines[i].substr(0,4)=="ATOM" || lines[i].substr(0,4)=="HETA")
                {
                        PdbLineToAtom(lines[i], Atoms[i]);
                }
        }
}

int ReadPdb(string PdbFile, vector<AtomStruct> &Atoms)
{
        AtomStruct Atom;
        //bool verbose=false;
        string line, ElementID, StrResidueID, StrAtomName;
        int natom;
        Real XBoxLength, YBoxLength, ZBoxLength;
        DeleteVector(Atoms);
        natom=0;

        fstream protein;
        OpenFile(PdbFile, protein, "PdbFile");


        getline(protein,line);
        natom=CountAtomsInPdb(PdbFile);
        SafeAlloc(Atoms, Atom, natom, "Atoms in ReadPdb");
        natom=0;
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
                        if (line.substr(0,4)=="ATOM") Atoms[natom].HetAtom=false;
                        else if (line.substr(0,4)=="HETA") Atoms[natom].HetAtom=true;
                        Atoms[natom].AtomNumber=StrToInt(line.substr(7,5));
                        Atoms[natom].AtomName=line.substr(12,4);
                        Atoms[natom].duplicate=line.substr(16,1);
                        Atoms[natom].ResidueName=line.substr(17,4);
                        Atoms[natom].ChainName=line.substr(21,1);
                        Atoms[natom].ResidueNum=StrToInt(line.substr(22,4));
                        //Atoms[natom].x=StrToFloat(line.substr(31,7));
                        //Atoms[natom].y=StrToFloat(line.substr(39,7));
                        //Atoms[natom].z=StrToFloat(line.substr(47,7));
                        Atoms[natom].x=StrToFloat(line.substr(30,8));
                        Atoms[natom].y=StrToFloat(line.substr(38,8));
                        Atoms[natom].z=StrToFloat(line.substr(46,8));
                        Atoms[natom].Occupancy=StrToFloat(line.substr(56,4));
                        Atoms[natom].BFactor=StrToFloat(line.substr(61,8));
                        Atoms[natom].SegID=line.substr(72,4);
                        Atoms[natom].ID=line.substr(73,3);
                        Atoms[natom].weight=1.0;
                        //cout <<"Atoms[natom].AtomName= "<<Atoms[natom].AtomName<<"."<<endl;
                        Trim(Atoms[natom].AtomName);
                        //cout <<"Atoms[natom].AtomName= "<<Atoms[natom].AtomName<<"."<<endl;
                        Trim(Atoms[natom].ResidueName);
                        cout <<Atoms[natom].ResidueName<<endl;
                        DeleteStartingAt(Atoms[natom].ResidueName, " ");
                        cout <<Atoms[natom].ResidueName<<endl;
                        Atoms[natom].residueid=GetResidueID(Atoms[natom].ResidueName);
                        natom++;
                }
                getline(protein,line);
        }
        protein.close();
        return natom;
}

int ReadNonZeroPdb(string PdbFile, vector<AtomStruct> &Atoms)
{
        AtomStruct Atom;
        bool valid;
        //bool verbose=false;
        char CharPdbFile[1000];
        string line, ElementID, StrResidueID, StrAtomName;
        int natom;
        Real BFactor;
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
        natom=CountNonZeroAtomsInPdb(PdbFile);
        SafeAlloc(Atoms, Atom, natom, "Atoms in ReadPdb");
        natom=0;
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
                        BFactor=StrToFloat(line.substr(61,8));
                        if (BFactor!=0)
                        {
                                PdbLineToAtom(line, Atoms[natom]);
                        }
                        natom++;
                }
                getline(protein,line);
        }
        protein.close();
        return natom;
}

void ReadXYZ(string xyz, vector<AtomStruct> &Atoms)
{
        AtomStruct Atom;
        bool valid;
        char CharXYZ[1000];
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
                Atom.connectivity.clear();
                Atom.connectivity.push_back(StrToInt(line.substr(Connect2, Connect3-Connect2)));
                Atom.connectivity.push_back(StrToInt(line.substr(Connect3, Connect4-Connect3)));
                Atom.connectivity.push_back(StrToInt(line.substr(Connect4, 100)));
                Atom.connectivity.push_back(StrToInt(line.substr(Connect4, 100)));
                SafePushBack(Atoms, Atom, "Atoms");
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
