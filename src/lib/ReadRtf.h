#ifndef _ReadRTF_included_
#define _ReadRTF_included_

# include "IOUtils.h"
# include "LinkedList.h"
# include "ReadPdb.h"
# include "ResidueCode.h"
# include "StringUtils.h"
# include "Structures.h"
# include "WritePdb.h"

struct PatchStruct
{
        string first, last, ResidueName;
};

void applyPatch(string patch, string ResidueName, vector<AtomStruct> &CharmmAtoms)
{
        int natom=CharmmAtoms.size();
        for (int i=0;i<natom;i++)
        {
                if (CharmmAtoms[i].ResidueName==patch)
                {
                        CharmmAtoms[i].ResidueName=ResidueName;
                }
        }
}

void applyPatches(vector<PatchStruct> &Patch, vector<AtomStruct> &CharmmAtoms)
{
        int npatch=Patch.size();
        for (int i=0;i<npatch;i++)
        {
                applyPatch(Patch[i].first, Patch[i].ResidueName, CharmmAtoms);
                applyPatch(Patch[i].last, Patch[i].ResidueName, CharmmAtoms);
        }
}

void ReadRtfFile(string RtfFile, vector<AtomStruct> &CharmmAtoms)
{
	bool valid, verbose=false;
	char CharRtfFile[1000];
	string line, first4, ResidueName;
        PatchStruct tempPatch;
        vector<PatchStruct> Patch;
	AtomStruct CharmmAtom;
	
	strcpy(CharRtfFile, RtfFile.c_str());
	valid=exist(CharRtfFile);
	if (!valid)
	{
		cout <<"Error unable to open Rtf file "<<RtfFile<<endl;
		exit(EXIT_FAILURE);
	}
	
	fstream File;
	File.open(CharRtfFile, ios::in);
	getline(File, line);
	while (true)
	{
		if (File.eof()) break;
                first4=line.substr(0,4);
		if (first4=="RESI" || first4=="PRES")
		{
                        ResidueName=GetWord2(line, 2);
			CharmmAtom.ResidueName=ResidueName;
			Trim(CharmmAtom.ResidueName);
		}
		if (first4=="ATOM")
		{
			CharmmAtom.AtomName=line.substr(5,5);
                        if (verbose) cout <<line.substr(5,4)<<"."<<endl;
			CharmmAtom.CharmmAtomName=line.substr(10,4);
			//CharmmAtom.charge=StrToFloat(line.substr(17,9));
                        CharmmAtom.charge=StrToFloat(GetWord2(line, 4));
                        Trim(CharmmAtom.AtomName);
			Trim(CharmmAtom.CharmmAtomName);
                        DeleteStartingAt(CharmmAtom.AtomName, "\t");
                        DeleteStartingAt(CharmmAtom.CharmmAtomName, "\t");
                        DeleteStartingAt(CharmmAtom.AtomName, " ");
                        DeleteStartingAt(CharmmAtom.CharmmAtomName, " ");
			CharmmAtoms.push_back(CharmmAtom);
		}
                if (first4=="PATC") //PATCH
                {
                        tempPatch.ResidueName=ResidueName;
                        tempPatch.first=GetWord2(line, 3);
                        tempPatch.last=GetWord2(line, 5);
                        Trim(tempPatch.first);
                        Trim(tempPatch.last);
                        DeleteChar(tempPatch.last, CARRIAGE_RETURN);
                        SafePushBack(Patch, tempPatch, "patch");
                }
		getline(File, line);
	}
	File.close();
        applyPatches(Patch, CharmmAtoms);
}
#endif
