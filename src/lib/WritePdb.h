#ifndef _WritePdb_included_
#define _WritePdb_included_

# include "IOUtils.h"
# include "Structures.h"
# include "TypeDef.h"
# include "StringUtils.h"

void PrintPdbLine(string PdbDir, bool HetAtom, int AtomNumber, string AtomName, string ResidueName, int ResidueNum, Real x, Real y, Real z, Real Occupancy, Real BFactor)
{
	int chars, AtomChars;
	char CharPdbDir[1000];

	ResidueName=ResidueName+".";
	AtomName=AtomName+".";

	chars=ResidueName.find('.');
	AtomChars=AtomName.find('.');
        if (ResidueNum>999) ResidueNum=999;
	ResidueName=ResidueName.substr(0,chars);
	AtomName=AtomName.substr(0,AtomChars)+" ";
	strcpy(CharPdbDir, PdbDir.c_str());
	ofstream pdb(CharPdbDir, ios::app);
	pdb << setfill(' ')
	        << setw(6) << setiosflags(ios::left) << (HetAtom ? "HETATM" : "ATOM") << resetiosflags(ios::left)
		<< setw(5) << setiosflags(ios::right) << (AtomNumber>100000 ? hex : dec ) << AtomNumber << resetiosflags(ios::right)
		<< (AtomChars>3 ? " " : "  ")
		<< (AtomChars>3 ? setw(3) : setw(4))
		<< (AtomChars>3 ? setiosflags(ios::right) : setiosflags(ios::left))
		<< AtomName << resetiosflags(ios::left)
		<< setiosflags(ios::right);

	pdb <<setfill(' ')
		<< setiosflags(ios::right) << ResidueName << resetiosflags(ios::left);
	if (chars>3) pdb << setw(1) <<"A"<< resetiosflags(ios::right);
	else if (chars==3) pdb <<setw(2) <<" A"<< resetiosflags(ios::right);
	else if (chars==2) pdb <<setw(3) <<"  A"<< resetiosflags(ios::right);
	else if (chars==1) pdb <<setw(4) <<"   A"<< resetiosflags(ios::right);

	pdb << setw(4) << setiosflags(ios::right) << ResidueNum<<resetiosflags(ios::right);
	pdb	<< setw(4) << " ";
	pdb	<< setiosflags(ios::fixed);
	pdb	<< setw(8) << setprecision(3);
	pdb	<< x;
	pdb	<< setw(8) << setprecision(3);
	pdb	<< y;
	pdb	<< setw(8) << setprecision(3);
	pdb	<< z;
	pdb	<< setw(6) << setprecision(2);
	pdb	<< Occupancy;
	pdb	<< setw(9) << setprecision(5);
	pdb	<< BFactor;
	pdb	<< setw(14) << " ";
	pdb	<< resetiosflags(ios::fixed);
	pdb	<< resetiosflags(ios::right);
	pdb	<< setw(0) <<endl;
}

void PrintPdbLine(ofstream &pdb, bool HetAtom, int AtomNumber, string AtomName, string ResidueName, int ResidueNum, Real x, Real y, Real z, Real Occupancy, Real BFactor)
{
	int chars, AtomChars;
	//char CharPdbDir[1000];

	ResidueName=ResidueName+".";
	AtomName=AtomName+".";

	chars=ResidueName.find('.');
	AtomChars=AtomName.find('.');
        if (ResidueNum>999) ResidueNum=999;
	ResidueName=ResidueName.substr(0,chars);
	AtomName=AtomName.substr(0,AtomChars)+" ";
	//strcpy(CharPdbDir, PdbDir.c_str());
	//ofstream pdb(CharPdbDir, ios::app);
	pdb << setfill(' ')
	        << setw(6) << setiosflags(ios::left) << (HetAtom ? "HETATM" : "ATOM") << resetiosflags(ios::left)
		<< setw(5) << setiosflags(ios::right) << (AtomNumber>100000 ? hex : dec ) << AtomNumber << resetiosflags(ios::right)
		<< (AtomChars>3 ? " " : "  ")
		<< (AtomChars>3 ? setw(3) : setw(4))
		<< (AtomChars>3 ? setiosflags(ios::right) : setiosflags(ios::left))
		<< AtomName << resetiosflags(ios::left)
		<< setiosflags(ios::right);

	pdb <<setfill(' ')
		<< setiosflags(ios::right) << ResidueName << resetiosflags(ios::left);
        cout <<"chars= "<<chars<<endl;
	if (chars>3) pdb << setw(1) <<"A"<< resetiosflags(ios::right);
	else if (chars==3) pdb <<setw(2) <<" A"<< resetiosflags(ios::right);
	else if (chars==2) pdb <<setw(3) <<"  A"<< resetiosflags(ios::right);
	else if (chars==1) pdb <<setw(4) <<"   A"<< resetiosflags(ios::right);

	pdb << setw(4) << setiosflags(ios::right) << ResidueNum<<resetiosflags(ios::right);
	pdb	<< setw(4) << " ";
	pdb	<< setiosflags(ios::fixed);
	pdb	<< setw(8) << setprecision(3);
	pdb	<< x;
	pdb	<< setw(8) << setprecision(3);
	pdb	<< y;
	pdb	<< setw(8) << setprecision(3);
	pdb	<< z;
	pdb	<< setw(6) << setprecision(2);
	pdb	<< Occupancy;
	pdb	<< setw(9) << setprecision(5);
	pdb	<< BFactor;
	pdb	<< setw(14) << " ";
	pdb	<< resetiosflags(ios::fixed);
	pdb	<< resetiosflags(ios::right);
	pdb	<< setw(0) <<endl;
}

void PrintAtomToFile(string PdbFile, AtomStruct &Atom)
{
        PrintPdbLine(PdbFile, Atom.HetAtom, Atom.AtomNumber, Atom.AtomName, Atom.ResidueName, Atom.ResidueNum, Atom.x, Atom.y, Atom.z, Atom.Occupancy, Atom.weight);
}

void PrintAtomToFile(ofstream &file, AtomStruct &Atom)
{
        PrintPdbLine(file, Atom.HetAtom, Atom.AtomNumber, Atom.AtomName, Atom.ResidueName, Atom.ResidueNum, Atom.x, Atom.y, Atom.z, Atom.Occupancy, Atom.weight);
}

void AddAtomsToTrj(string PdbOutputFile, vector<AtomStruct> &Atoms)
{
	int m, n, natom;
	string PdbDir, StrFileIndex;
	char CharPdbDir[1000];
	int chars, AtomChars;
	//cout <<"In PrintPdb"<<endl;
	natom=Atoms.size();
	//cout <<"natom= "<<natom<<endl;
	m=1;
	cout <<"writing "<<PdbOutputFile<<endl;
        strcpy(CharPdbDir, PdbOutputFile.c_str());
        ofstream pdb(CharPdbDir, ios::app);

	cout <<"CharPdbDir="<<CharPdbDir<<endl;

	for (n=0;n<natom;n++)
	{
		if (Atoms[n].AtomNumber>9999) Atoms[n].AtomNumber=9999;
                //if (Atoms[n].ResidueNum>999) Atoms[n].ResidueNum=999;
		chars=Atoms[n].ResidueName.length();
		AtomChars=Atoms[n].AtomName.length();

		Atoms[n].ResidueName=Atoms[n].ResidueName.substr(0,chars);
		Atoms[n].AtomName=Atoms[n].AtomName.substr(0,AtomChars)+" ";
		pdb  << setfill(' ')
			<< setw(6) << setiosflags(ios::left) << (Atoms[n].HetAtom ? "HETATM" : "ATOM") << resetiosflags(ios::left)
			<< setw(5) << setiosflags(ios::right) << (Atoms[n].AtomNumber>=100000 ? hex : dec ) << Atoms[n].AtomNumber << resetiosflags(ios::right)
			<< (AtomChars>3 ? " " : "  ")
			<< (AtomChars>3 ?  setw(3) : setw(4))
			<< (AtomChars>3 ?  setiosflags(ios::right) : setiosflags(ios::left))
			<< Atoms[n].AtomName << resetiosflags(ios::left)
			//<< " "
			<< setiosflags(ios::right);
		//<< setw(4)

		pdb <<setfill(' ')
			<< setiosflags(ios::right)<< Atoms[n].ResidueName <<resetiosflags(ios::left);

                string chainNameBlanks;
                if (chars==3) chainNameBlanks=" "+Atoms[n].ChainName;
                if (chars==2) chainNameBlanks="  "+Atoms[n].ChainName;
                if (chars==1) chainNameBlanks="   "+Atoms[n].ChainName;
		if (chars>3) pdb << setw(1) <<Atoms[n].ChainName<< resetiosflags(ios::right);
	        else if (chars==3) pdb <<setw(2) <<chainNameBlanks<< resetiosflags(ios::right);
	        else if (chars==2) pdb <<setw(3) <<chainNameBlanks<< resetiosflags(ios::right);
	        else if (chars==1) pdb <<setw(4) <<chainNameBlanks<< resetiosflags(ios::right);

		pdb << setw(4) << setiosflags(ios::right)<<Atoms[n].ResidueNum<<resetiosflags(ios::right)
			<< setw(4) << " "
			<< setiosflags(ios::fixed)
			<< setw(8) << setprecision(3)
			<< Atoms[n].x
			<< setw(8) << setprecision(3)
			<< Atoms[n].y
			<< setw(8) << setprecision(3)
			<< Atoms[n].z
			<< setw(6) << setprecision(2)
			<< 1.0
			<< setw(9) << setprecision(5)
			<< Atoms[n].weight
			<< setw(14) << " "
			<< resetiosflags(ios::fixed)
			<< resetiosflags(ios::right)
			<< setw(0) <<"\n";                       

		if (m<9999) m++;
		Trim(Atoms[n].AtomName);
		Trim(Atoms[n].ResidueName);
	}

	pdb <<"ENDMDL"<<endl;


}

void PrintPdb(string PdbFile, vector<AtomStruct> &Atoms)
{
	AddIndexToFile(PdbFile);
        AddAtomsToTrj(PdbFile, Atoms);
}

void PrintAtomInfo(AtomStruct Atom)
{
	int m=1;
	int chars, AtomChars;



	chars=Atom.ResidueName.length();
	AtomChars=Atom.AtomName.length();

	Atom.ResidueName=Atom.ResidueName.substr(0,chars);
	Atom.AtomName=Atom.AtomName.substr(0,AtomChars)+" ";

	cout  << setfill(' ')
		<< setw(6) << setiosflags(ios::left) << (Atom.HetAtom ? "HETATM" : "ATOM") << resetiosflags(ios::left)
		<< setw(5) << setiosflags(ios::right) << (m>=100000 ? hex : dec ) << Atom.AtomNumber << resetiosflags(ios::right)
		<< (AtomChars>3 ? " " : "  ")
		<< (AtomChars>3 ?  setw(3) : setw(4))
		<< (AtomChars>3 ?  setiosflags(ios::right) : setiosflags(ios::left))
		<< Atom.AtomName << resetiosflags(ios::left)
		//<< " "
		<< setiosflags(ios::right);
	//<< setw(4)

	cout <<setfill(' ')
		<< setiosflags(ios::right)<< Atom.ResidueName <<resetiosflags(ios::left);

	if (chars>3) cout << setw(1) <<"A"<< resetiosflags(ios::right);
	else cout << setw(2) <<" A"<< resetiosflags(ios::right);

	cout << setw(4) << setiosflags(ios::right)<<Atom.ResidueNum<<resetiosflags(ios::right)
		<< setw(4) << " "
		<< setiosflags(ios::fixed)
		<< setw(8) << setprecision(3)
		<< Atom.x
		<< setw(8) << setprecision(3)
		<< Atom.y
		<< setw(8) << setprecision(3)
		<< Atom.z
		<< setw(6) << setprecision(2)
		<< 1.0
		<< setw(9) << setprecision(5)
		<< Atom.weight
		<< setw(14) << " "
		<< resetiosflags(ios::fixed)
		<< resetiosflags(ios::right)
		<< setw(0) <<endl;                       
}

#endif



