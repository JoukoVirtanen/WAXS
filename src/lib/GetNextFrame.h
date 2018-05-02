#ifndef _GetNextFrame_included_
#define _GetNextFrame_included_

# include "IOUtils.h"
# include "ReadDcd.h"
# include "ReadPdb.h"
# include "Structures.h"

//bool verbose=true;

bool GetDcdFileName(string &DcdFile, string DcdFilePaths)
{
	string line;
	fstream file;
	//if (verbose) 
	cout <<"In GetDcdFileName"<<endl;
	OpenFile(DcdFilePaths, file, "DcdFilePaths");
	getline(file, line);
	if (DcdFile=="") DcdFile=line;
	else
	{
		while (true)
		{
			if (file.eof()) break;
			if (line==DcdFile)
			{
				getline(file, line);
				if (file.eof()) return false;
				DcdFile=line;
				file.close();
				return true;
			}
			getline(file, line);
		}
	}
	file.close();
	return false;
}

bool GetNextDcdFrame(vector<AtomStruct> &Atoms, int &nthstruct, int &nInDcd, int &nthDcd, string &DCDFile, string DCDFilePaths, ParamStruct &params)
{
	bool IsStructure=true;
	//if (verbose) 
	cout <<"In GetNextDcdFrame"<<endl;
	cout <<"nthDcd= "<<nthDcd<<" nthstruct= "<<nthstruct<<endl;
	if (nthDcd==1 && nthstruct==1)
	{
		IsStructure=GetLineInFile(DCDFilePaths, DCDFile, nthDcd);
		ReadDcd(DCDFile, Atoms, nthstruct, nInDcd, params.XBoxLength, params.YBoxLength, params.ZBoxLength);
	}
	if (nthstruct>nInDcd) 
	{
		nthDcd++;
		//cout <<"nthstruct= "<<nthstruct<<" nInDcd= "<<nInDcd<<" nthDcd= "<<nthDcd<<endl;
		IsStructure=GetLineInFile(DCDFilePaths, DCDFile, nthDcd);
		//cout <<"IsStructure= "<<IsStructure<<" DCDFilePaths= "<<DCDFilePaths<<" DCDFile= "<<DCDFile<<endl;
		if (!IsStructure || DCDFile=="") 
                {
                        cout <<" In if (!IsStructure || DCDFile==)"<<endl;
                        return false;
                }
		nthstruct=1;
	}	
	ReadDcd(DCDFile, Atoms, nthstruct, nInDcd, params.XBoxLength, params.YBoxLength, params.ZBoxLength);
	cout <<"In GetNextDcdFrame"<<endl;
        cout <<"DCDFile= "<<DCDFile<<" IsStructrue= "<<IsStructure<<endl;
        cout <<"params.XBoxLength= "<<params.XBoxLength<<endl;
        cout <<"params.YBoxLength= "<<params.YBoxLength<<endl;
        cout <<"params.ZBoxLength= "<<params.ZBoxLength<<endl;
	PrintAtomInfo(Atoms[0]);
	PrintAtomInfo(Atoms[1]);
	PrintAtomInfo(Atoms[2]);
	cout <<"Leaving GetNextDcdFrame"<<endl;
	nthstruct++;
	return IsStructure;
}

bool GetNextPdbFrame(vector<AtomStruct> &Atoms, string PdbFilePaths, int &nthstruct, ParamStruct &params)
{
	bool IsStructure;
	ProteinStruct Protein;
	string PdbFilePath;
	//if (verbose) 
	IsStructure=GetLineInFile(PdbFilePaths, PdbFilePath, nthstruct);
	if (IsStructure) ReadPdb(PdbFilePath, Protein);
	CopyVector(Protein.Atoms, Atoms);
        if (params.XBoxLength==0)
        {
	        params.XBoxLength=Protein.XBoxLength;
	        params.YBoxLength=Protein.YBoxLength;
	        params.ZBoxLength=Protein.ZBoxLength;
	}
        nthstruct++;
	return IsStructure;
}

bool GetNextFrame(vector<AtomStruct> &Atoms, int &nthstruct, int &nInDcd, int &nthDcd, string &DCDFile, string DCDFilePaths, string PdbFilePaths, ParamStruct &params)
{
	bool IsStructure;
	//if (verbose) cout <<"in GetNextFrame. DCDFilePaths= "<<DCDFilePaths<<endl;
	cout <<"DCDFilePaths= "<<DCDFilePaths<<endl;
        if (DCDFilePaths!="") IsStructure=GetNextDcdFrame(Atoms, nthstruct, nInDcd, nthDcd, DCDFile, DCDFilePaths, params);
	else if (PdbFilePaths!="") IsStructure=GetNextPdbFrame(Atoms, PdbFilePaths, nthstruct, params);
	else
	{
		cout <<"Fatal error. No DCDFiles or PDBFiles given"<<endl;
		exit(0);
	}
        cout <<"params.XBoxLength= "<<params.XBoxLength<<endl;
        cout <<"params.YBoxLength= "<<params.YBoxLength<<endl;
        cout <<"params.ZBoxLength= "<<params.ZBoxLength<<endl;
	AtomAlias(Atoms);
        AssignAtomIDs(Atoms);
	//GetAtomType(Atoms, !params.UniformHydrationShell);
	GetAtomType(Atoms, true);
	GetAtomType2(Atoms);
	SetSolute(Atoms);
        MoveAtoms(Atoms, -params.XOrigin, -params.YOrigin, -params.ZOrigin);
	return IsStructure;
}
#endif
