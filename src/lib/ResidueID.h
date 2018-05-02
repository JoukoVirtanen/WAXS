#ifndef _ResidueID_included_
#define _ResidueID_included_

# include "ResidueCode.h"

//const string OneLetterCode[20]={"A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"};

string OneLetterCode(int ResidueID)
{
        if (ResidueID==GLY) return "G";
        else if (ResidueID==ALA) return "A";
        else if (ResidueID==VAL) return "V";
        else if (ResidueID==LEU) return "L";
        else if (ResidueID==ILE) return "I";
        else if (ResidueID==SER) return "S";
        else if (ResidueID==THR) return "T";
        else if (ResidueID==CYS) return "C";
        else if (ResidueID==MET) return "M";
        else if (ResidueID==ASP) return "D";
        else if (ResidueID==GLU) return "E";
        else if (ResidueID==ASN) return "N";
        else if (ResidueID==GLN) return "Q";
        else if (ResidueID==LYS) return "K";
        else if (ResidueID==ARG) return "R";
        else if (ResidueID==PHE) return "F";
        else if (ResidueID==TYR) return "Y";
        else if (ResidueID==PRO) return "P";
        else if (ResidueID==HIS) return "H";
        else if (ResidueID==TRP) return "W";
        else
        {
                cout <<"ERROR: Unknown ResidueID "<<ResidueID<<endl;
                exit(EXIT_FAILURE);
        }
}

/*
OneLetterCode[ALA]="A";
OneLetterCode[ARG]="R";
OneLetterCode[ASN]="N";
OneLetterCode[ASP]="D";
OneLetterCode[CYS]="C";
OneLetterCode[GLU]="E";
OneLetterCode[GLN]="Q";
OneLetterCode[GLY]="G";
OneLetterCode[HIS]="H";
OneLetterCode[ILE]="I";
OneLetterCode[LEU]="L";
OneLetterCode[LYS]="K";
OneLetterCode[MET]="M";
OneLetterCode[PHE]="F";
OneLetterCode[PRO]="P";
OneLetterCode[SER]="S";
OneLetterCode[THR]="T";
OneLetterCode[TRP]="W";
OneLetterCode[TYR]="Y";
OneLetterCode[VAL]="V";
*/
int GetResidueID(string StrResidueID)
{
	bool verbose=false;
	if (StrResidueID=="GLY") return GLY;
	else if (StrResidueID=="ALA") return ALA;
	else if (StrResidueID=="VAL") return VAL;
	else if (StrResidueID=="LEU") return LEU;
	else if (StrResidueID=="ILE") return ILE;
	else if (StrResidueID=="SER") return SER;
	else if (StrResidueID=="THR") return THR;
	else if (StrResidueID=="CYS") return CYS;
	else if (StrResidueID=="MET") return MET;
	else if (StrResidueID=="ASP") return ASP;
	else if (StrResidueID=="GLU") return GLU;
	else if (StrResidueID=="ASN") return ASN;
	else if (StrResidueID=="GLN") return GLN;
	else if (StrResidueID=="LYS") return LYS;
	else if (StrResidueID=="ARG") return ARG;
	else if (StrResidueID=="PHE") return PHE;
	else if (StrResidueID=="TYR") return TYR;
	else if (StrResidueID=="PRO") return PRO;
	else if (StrResidueID=="HIS") return HIS;
	else if (StrResidueID=="HSE") return HIS;
	else if (StrResidueID=="HSD") return HIS;
	else if (StrResidueID=="TRP") return TRP;
	else if (StrResidueID=="HEM") return HEM;
	else if (StrResidueID=="WAT") return WAT;
	else if (StrResidueID=="TIP3") return WAT;
	else if (StrResidueID=="TIP4") return WAT;
	else if (StrResidueID=="TIP5") return WAT;
	else
	{
		return UNK;
		if (verbose) cout <<"Warning unknown resiude "<<StrResidueID<<endl;
	}
}

void GetResidueIDs(vector<AtomStruct> &Atoms)
{
	int i, natom=Atoms.size();
	for (i=0;i<natom;i++) Atoms[i].residueid=GetResidueID(Atoms[i].ResidueName);
}

#endif
