#ifndef _GetAtomTypes_included
#define _GetAtomTypes_included

# include "AtomCode.h"
# include "LinkedList.h"
# include "StringUtils.h"
# include "Structures.h"

int GetSolventType(string resstr, string str)
{
	bool verbose=true, old_equivalence=true;
	string str1, str2, str2b, str3, str_trim;
	Trim(str);
	str1=str[0];
	str2=str.substr(0,2);
	str3=str.substr(0,3);
	//cout <<"str1= "<<str1<<" str2= "<<str2<<" str3= "<<str3<<endl;
	if (str2=="HT") return HT;
        else if (str=="H3T") return HT;
	else if (str2=="OT") return OT;

	else if (resstr=="ALA")
	{
		if (str=="N") return ALA_N;
		else if (str2=="HN") return ALA_HN;
		else if (str2=="CA") return ALA_CA;
		else if (str2=="HA") return ALA_HA;
		else if (str2=="CB") return ALA_CB;
		else if (str2=="HB") return ALA_HB;
		else if (str=="C") return ALA_C;
		else if (str=="O") return ALA_O;
		else
		{
			if (verbose) cout <<"Error unable to identify solvent type.  Residue= "<<resstr<<" str= "<<str<<endl;
			return UNK_ATOM_TYPE;
		}
	}

	else if (resstr=="ARG")
	{
		if (str=="N") return ARG_N;
		else if (str2=="HN") return ARG_HN;
		else if (str2=="CA") return ARG_CA;
		else if (str2=="HA") return ARG_HA;
		else if (str2=="CB") return ARG_CB;
		else if (str2=="HB") return ARG_HB;
		else if (str2=="CG") return ARG_CG;
		else if (str2=="HG") return ARG_HG;
		else if (str2=="CD") return ARG_CD;
		else if (str2=="HD") return ARG_HD;
		else if (str2=="NE") return ARG_NE;
		else if (str2=="HE") return ARG_HE;
		else if (str2=="CZ") return ARG_CZ;
		else if (str2=="NH") return ARG_NH;
		else if (str2=="HH") return ARG_HH;
		else if (str=="C") return ARG_C;
		else if (str=="O") return ARG_O;
		else
		{
			if (verbose) cout <<"Error unable to identify solvent type.  Residue= "<<resstr<<" str= "<<str<<endl;
			return UNK_ATOM_TYPE;
		}
	}

	else if (resstr=="ASN")
	{
		if (str=="N") return ASN_N;
		else if (str2=="HN") return ASN_HN;
		else if (str2=="CA") return ASN_CA;
		else if (str2=="HA") return ASN_HA;
		else if (str2=="CB") return ASN_CB;
		else if (str2=="HB") return ASN_HB;
		else if (str2=="CG") return ASN_CG;
		else if (str2=="OD") return ASN_OD;
		else if (str2=="ND") return ASN_ND;
		else if (str2=="HD") return ASN_HD;
		else if (str=="C") return ASN_C;
		else if (str=="O") return ASN_O;
		else
		{

			if (verbose) cout <<"Error unable to identify solvent type.  Residue= "<<resstr<<" str= "<<str<<endl;
			return UNK_ATOM_TYPE;
		}
	}

	else if (resstr=="ASP")
	{
		if (str=="N") return ASP_N;
		else if (str2=="HN") return ASP_HN;
		else if (str2=="CA") return ASP_CA;
		else if (str2=="HA") return ASP_HA;
		else if (str2=="CB") return ASP_CB;
		else if (str2=="HB") return ASP_HB;
		else if (str2=="CG") return ASP_CG;
		else if (str2=="OD") return ASP_OD;
		else if (str2=="ND") return ASP_ND;
		else if (str2=="HD") return ASP_HD;
		else if (str=="C") return ASP_C;
		else if (str=="O") return ASP_O;
		else
		{
			if (verbose) cout <<"Error unable to identify solvent type.  Residue= "<<resstr<<" str= "<<str<<endl;
			return UNK_ATOM_TYPE;
		}
	}

	else if (resstr=="CYS")
	{
		if (str=="N") return CYS_N;
		else if (str2=="HN") return CYS_HN;
		else if (str2=="CA") return CYS_CA;
		else if (str2=="HA") return CYS_HA;
		else if (str2=="CB") return CYS_CB;
		else if (str2=="HB") return CYS_HB;
		else if (str2=="SG") return CYS_SG;
		else if (str2=="HG") return CYS_HG;
		else if (str=="C") return CYS_C;
		else if (str=="O") return CYS_O;
		else
		{
			if (verbose) cout <<"Error unable to identify solvent type.  Residue= "<<resstr<<" str= "<<str<<endl;
			return UNK_ATOM_TYPE;
		}
	}

	else if (resstr=="GLU")
	{
		if (str=="N") return GLU_N;
		else if (str2=="HN") return GLU_HN;
		else if (str2=="CA") return GLU_CA;
		else if (str2=="HA") return GLU_HA;
		else if (str2=="CB") return GLU_CB;
		else if (str2=="HB") return GLU_HB;
		else if (str2=="CG") return GLU_CG;
		else if (str2=="HG") return GLU_HG;
		else if (str2=="CD") return GLU_CD;
		else if (str2=="OE") return GLU_OE;
		else if (str=="C") return GLU_C;
		else if (str=="O") return GLU_O;
		else
		{
			if (verbose) cout <<"Error unable to identify solvent type.  Residue= "<<resstr<<" str= "<<str<<endl;
			return UNK_ATOM_TYPE;
		}
	}

	else if (resstr=="GLN")
	{
		if (str=="N") return GLN_N;
		else if (str2=="HN") return GLN_HN;
		else if (str2=="CA") return GLN_CA;
		else if (str2=="HA") return GLN_HA;
		else if (str2=="CB") return GLN_CB;
		else if (str2=="HB") return GLN_HB;
		else if (str2=="CG") return GLN_CG;
		else if (str2=="HG") return GLN_HG;
		else if (str2=="CD") return GLN_CD;
		else if (str2=="OE") return GLN_OE;
		else if (str2=="NE") return GLN_NE;
		else if (str2=="HE") return GLN_HE;
		else if (str=="C") return GLN_C;
		else if (str=="O") return GLN_O;
		else
		{
			if (verbose) cout <<"Error unable to identify solvent type.  Residue= "<<resstr<<" str= "<<str<<endl;
			return UNK_ATOM_TYPE;
		}				
	}

	else if (resstr=="GLY")
	{
		if (str=="N") return GLY_N;
		else if (str2=="HN") return GLY_HN;
		else if (str2=="CA") return GLY_CA;
		else if (str2=="HA") return GLY_HA;
		else if (str=="C") return GLY_C;
		else if (str=="O") return GLY_O;
		else
		{
			if (verbose) cout <<"Error unable to identify solvent type.  Residue= "<<resstr<<" str= "<<str<<endl;
			return UNK_ATOM_TYPE;
		}
	}

	else if (resstr=="HIS" || resstr=="HSE" || resstr=="HSD")
	{
		if (str=="N") return HIS_N;
		else if (str2=="HN") return HIS_HN;
		else if (str2=="CA") return HIS_CA;
		else if (str2=="HA") return HIS_HA;
		else if (str2=="CB") return HIS_CB;
		else if (str2=="HB") return HIS_HB;
		else if (str2=="ND") return HIS_ND;
		else if (str3=="HD1") return HIS_HD1;
		else if (str2=="CG") return HIS_CG;
		else if (str2=="CE") return HIS_CE;
		else if (str2=="HE") return HIS_HE;
		else if (str2=="NE") return HIS_NE;
		else if (str2=="CD") return HIS_CD;
		else if (str3=="HD2") return HIS_HD2;
		else if (str=="C") return HIS_C;
		else if (str=="O") return HIS_O;
		else
		{
			if (verbose) cout <<"Error unable to identify solvent type.  Residue= "<<resstr<<" str= "<<str<<endl;
			return UNK_ATOM_TYPE;
		}					
	}

	else if (resstr=="ILE")
	{
		if (str=="N") return ILE_N;
		else if (str2=="HN") return ILE_HN;
		else if (str2=="CA") return ILE_CA;
		else if (str2=="HA") return ILE_HA;
		else if (str2=="CB") return ILE_CB;
		else if (str2=="HB") return ILE_HB;
		else if (str2=="CG") return ILE_CG;
		else if (str2=="HG") return ILE_HG;
		else if (str2=="CD") return ILE_CD;
		else if (str2=="HD") return ILE_HD;
		else if (str=="C") return ILE_C;
		else if (str=="O") return ILE_O;
		else
		{
			if (verbose) cout <<"Error unable to identify solvent type.  Residue= "<<resstr<<" str= "<<str<<endl;
			return UNK_ATOM_TYPE;
		}
	}

	else if (resstr=="LEU")
	{
		if (str=="N") return LEU_N;
		else if (str2=="HN") return LEU_HN;
		else if (str2=="CA") return LEU_CA;
		else if (str2=="HA") return LEU_HA;
		else if (str2=="CB") return LEU_CB;
		else if (str2=="HB") return LEU_HB;
		else if (str2=="CG") return LEU_CG;
		else if (str2=="HG") return LEU_HG;
		else if (str2=="CD" && old_equivalence) return LEU_CD;
		else if (str3=="CD1") return LEU_CD1;
		else if (str3=="CD2") return LEU_CD2;
		else if (str2=="HD" && old_equivalence) return LEU_HD;
		else if (str3=="HD1") return LEU_HD1;
		else if (str3=="HD2") return LEU_HD2;
		else if (str=="C") return LEU_C;
		else if (str=="O") return LEU_O;
		else
		{
			if (verbose) cout <<"Error unable to identify solvent type.  Residue= "<<resstr<<" str= "<<str<<endl;
			return UNK_ATOM_TYPE;
		}
	}

	else if (resstr=="LYS")
	{
		if (str=="N") return LYS_N;
		else if (str2=="HN") return LYS_HN;
		else if (str2=="CA") return LYS_CA;
		else if (str2=="HA") return LYS_HA;
		else if (str2=="CB") return LYS_CB;
		else if (str2=="HB") return LYS_HB;
		else if (str2=="CG") return LYS_CG;
		else if (str2=="HG") return LYS_HG;
		else if (str2=="CD") return LYS_CD;
		else if (str2=="HD") return LYS_HD;
		else if (str2=="CE") return LYS_CE;
		else if (str2=="HE") return LYS_HE;
		else if (str2=="NZ") return LYS_NZ;
		else if (str2=="HZ") return LYS_HZ;
		else if (str=="C") return LYS_C;
		else if (str=="O") return LYS_O;
		else
		{
			if (verbose) cout <<"Error unable to identify solvent type.  Residue= "<<resstr<<" str= "<<str<<endl;
			return UNK_ATOM_TYPE;
		}
	}

	else if (resstr=="MET")
	{
		if (str=="N") return MET_N;
		else if (str2=="HN") return MET_HN;
		else if (str2=="CA") return MET_CA;
		else if (str2=="HA") return MET_HA;
		else if (str2=="CB") return MET_CB;
		else if (str2=="HB") return MET_HB;
		else if (str2=="CG") return MET_CG;
		else if (str2=="HG") return MET_HG;
		else if (str2=="SD") return MET_SD;
		else if (str2=="CE") return MET_CE;
		else if (str2=="HE") return MET_HE;
		else if (str=="C") return MET_C;
		else if (str=="O") return MET_O;
		else
		{
			if (verbose) cout <<"Error unable to identify solvent type.  Residue= "<<resstr<<" str= "<<str<<endl;
			return UNK_ATOM_TYPE;
		}					
	}

	else if (resstr=="PHE")
	{
		if (str=="N") return PHE_N;
		else if (str2=="HN") return PHE_HN;
		else if (str2=="CA") return PHE_CA;
		else if (str2=="HA") return PHE_HA;
		else if (str2=="CB") return PHE_CB;
		else if (str2=="HB") return PHE_HB;
		else if (str2=="CG") return PHE_CG;
		else if (str2=="CD") return PHE_CD;
		else if (str2=="HD") return PHE_HD;
		else if (str2=="CE") return PHE_CE;
		else if (str2=="HE") return PHE_HE;
		else if (str2=="CZ") return PHE_CZ;
		else if (str2=="HZ") return PHE_HZ;
		else if (str=="C") return PHE_C;
		else if (str=="O") return PHE_O;
		else
		{
			if (verbose) cout <<"Error unable to identify solvent type.  Residue= "<<resstr<<" str= "<<str<<endl;
			return UNK_ATOM_TYPE;
		}					
	}

	else if (resstr=="PRO")
	{
		if (str=="N") return PRO_N;
		else if (str2=="CD") return PRO_CD;
		else if (str2=="HD") return PRO_HD;
		else if (str2=="CA") return PRO_CA;
		else if (str2=="HA") return PRO_HA;
		else if (str2=="CB") return PRO_CB;
		else if (str2=="HB") return PRO_HB;
		else if (str2=="CG") return PRO_CG;
		else if (str2=="HG") return PRO_HG;
		else if (str=="C") return PRO_C;
		else if (str=="O") return PRO_O;
		else
		{
			if (verbose) cout <<"Error unable to identify solvent type.  Residue= "<<resstr<<" str= "<<str<<endl;
			return UNK_ATOM_TYPE;
		}
	}

	else if (resstr=="SER")
	{
		if (str=="N") return SER_N;
		else if (str2=="HN") return SER_HN;
		else if (str2=="CA") return SER_CA;
		else if (str2=="HA") return SER_HA;
		else if (str2=="CB") return SER_CB;
		else if (str2=="HB") return SER_HB;
		else if (str2=="OG") return SER_OG;
		else if (str2=="HG") return SER_HG;
		else if (str=="C") return SER_C;
		else if (str=="O") return SER_O;
		else
		{
			if (verbose) cout <<"Error unable to identify solvent type.  Residue= "<<resstr<<" str= "<<str<<endl;
			return UNK_ATOM_TYPE;
		}
	}

	else if(resstr=="THR")
	{
		if (str=="N") return THR_N;
		else if (str2=="HN") return THR_HN;
		else if (str2=="CA") return THR_CA;
		else if (str2=="HA") return THR_HA;
		else if (str2=="CB") return THR_CB;
		else if (str2=="HB") return THR_HB;
		else if (str2=="OG") return THR_OG;
		else if (str3=="HG1") return THR_HG1;
		else if (str2=="CG") return THR_CG;
		else if (str3=="HG2") return THR_HG2;
		else if (str=="C") return THR_C;
		else if (str=="O") return THR_O;
		else
		{
			if (verbose) cout <<"Error unable to identify solvent type.  Residue= "<<resstr<<" str= "<<str<<endl;
			return UNK_ATOM_TYPE;
		}
	}

	else if (resstr=="TRP")
	{
		if (str=="N") return TRP_N;
		else if (str2=="HN") return TRP_HN;
		else if (str2=="CA") return TRP_CA;
		else if (str2=="HA") return TRP_HA;
		else if (str2=="CB") return TRP_CB;
		else if (str2=="HB") return TRP_HB;
		else if (str2=="CG") return TRP_CG;
		else if (str3=="CD1") return TRP_CD1;
		else if (str2=="HD") return TRP_HD;
		else if (str2=="NE") return TRP_NE;
		else if (str3=="HE1") return TRP_HE1;
		else if (str3=="CE2") return TRP_CE2;
		else if (str3=="CD2") return TRP_CD2;
		else if (str3=="CE3") return TRP_CE3;
		else if (str3=="HE3") return TRP_HE3;
		else if (str3=="CZ3") return TRP_CZ3;
		else if (str3=="HZ3") return TRP_HZ3;
		else if (str3=="CZ2") return TRP_CZ2;
		else if (str3=="HZ2") return TRP_HZ2;
		else if (str3=="CH2") return TRP_CH2;
		else if (str2=="HH") return TRP_HH;
		else if (str=="C") return TRP_C;
		else if (str=="O") return TRP_O;
		else
		{
			if (verbose) cout <<"Error unable to identify solvent type.  Residue= "<<resstr<<" str= "<<str<<endl;
			return UNK_ATOM_TYPE;
		}					
	}

	else if (resstr=="TYR")
	{
		if (str=="N") return TYR_N;
		else if (str2=="HN") return TYR_HN;
		else if (str2=="CA") return TYR_CA;
		else if (str2=="HA") return TYR_HA;
		else if (str2=="CB") return TYR_CB;
		else if (str2=="HB") return TYR_HB;
		else if (str2=="CG") return TYR_CG;
		else if (str2=="CD") return TYR_CD;
		else if (str2=="HD") return TYR_HD;
		else if (str2=="CE") return TYR_CE;
		else if (str2=="HE") return TYR_HE;
		else if (str2=="CZ") return TYR_CZ;
		else if (str2=="OH") return TYR_OH;
		else if (str2=="HH") return TYR_HH;
		else if (str=="C") return TYR_C;
		else if (str=="O") return TYR_O;
		else
		{
			if (verbose) cout <<"Error unable to identify solvent type.  Residue= "<<resstr<<" str= "<<str<<endl;
			return UNK_ATOM_TYPE;
		}					
	}

	else if (resstr=="VAL")
	{
		if (str=="N") return VAL_N;
		else if (str2=="HN") return VAL_HN;
		else if (str2=="CA") return VAL_CA;
		else if (str2=="HA") return VAL_HA;
		else if (str2=="CB") return VAL_CB;
		else if (str2=="HB") return VAL_HB;
		else if (str2=="CG") return VAL_CG;
		else if (str2=="HG") return VAL_HG;
		else if (str=="C") return VAL_C;
		else if (str=="O") return VAL_O;
		else
		{
			if (verbose) cout <<"Error unable to identify solvent type.  Residue= "<<resstr<<" str= "<<str<<endl;
			return UNK_ATOM_TYPE;
		}
	}

	else if (resstr=="HSP")
	{
		if (str=="N") return HSP_N;
		else if (str2=="HN") return HSP_HN;
		else if (str2=="CA") return HSP_CA;
		else if (str2=="HA") return HSP_HA;
		else if (str2=="CB") return HSP_CB;
		else if (str2=="HB") return HSP_HB;
		else if (str2=="CD") return HSP_CD;
		else if (str3=="HD2") return HSP_HD2;
		else if (str2=="CG") return HSP_CG;
		else if (str2=="NE") return HSP_NE;
		else if (str3=="HE2") return HSP_HE2;
		else if (str2=="ND") return HSP_ND;
		else if (str3=="HD1") return HSP_HD1;
		else if (str2=="CE") return HSP_CE;
		else if (str3=="HE1") return HSP_HE1;
		else if (str=="C") return HSP_C;
		else if (str=="O") return HSP_O;
		else
		{
			if (verbose) cout <<"Error unable to identify solvent type.  Residue= "<<resstr<<" str= "<<str<<endl;
			return UNK_ATOM_TYPE;
		}
	}
	else if (resstr=="HEM")
	{
		if (str2=="FE") return HEM_FE;
		else if (str2=="NA") return HEM_NA;
		else if (str2=="NB") return HEM_NB;
		else if (str2=="NC") return HEM_NC;
		else if (str2=="ND") return HEM_ND;
		else if (str2=="C1") return HEM_C1;
		else if (str2=="C2") return HEM_C2;
		else if (str2=="C3") return HEM_C3;
		else if (str2=="C4") return HEM_C4;
		else if (str2=="CH") return HEM_CH;
		else if (str2=="HA") return HEM_HA;
		else if (str2=="HB") return HEM_HB;
		else if (str2=="HC") return HEM_HC;
		else if (str2=="HD") return HEM_HD;
		else if (str2=="CM") return HEM_CM;
		else if (str2=="HM") return HEM_HM;
		else if (str2=="CA") return HEM_CA;
		else if (str2=="CB") return HEM_CB;
		else if (str2=="CG") return HEM_CG;
		else if (str2=="O1") return HEM_O1;
		else if (str2=="O2") return HEM_O2;
		else
		{

			if (verbose) cout <<"Error unable to identify solvent type.  Residue= "<<resstr<<" str= "<<str<<endl;
			return UNK_ATOM_TYPE;
		}
	}
	else
	{
		if (verbose) cout <<"Error unable to identify solvent type.  Residue= "<<resstr<<" str= "<<str<<endl;
		return UNK_ATOM_TYPE;
	}
}

void AtomTypeToResidueName(string ResidueName[])
{
	ResidueName[UNK_ATOM_TYPE]="UNK";
	ResidueName[HT]="UNK";
	ResidueName[OT]="UNK";
	ResidueName[ALA_N]="ALA";
	ResidueName[ALA_HN]="ALA";
	ResidueName[ALA_CA]="ALA";
	ResidueName[ALA_HA]="ALA";
	ResidueName[ALA_CB]="ALA";
	ResidueName[ALA_HB]="ALA";
	ResidueName[ALA_C]="ALA";
	ResidueName[ALA_O]="ALA";
	ResidueName[ARG_N]="ARG";
	ResidueName[ARG_HN]="ARG";
	ResidueName[ARG_CA]="ARG";
	ResidueName[ARG_HA]="ARG";
	ResidueName[ARG_CB]="ARG";
	ResidueName[ARG_HB]="ARG";
	ResidueName[ARG_CG]="ARG";
	ResidueName[ARG_HG]="ARG";
	ResidueName[ARG_CD]="ARG";
	ResidueName[ARG_HD]="ARG";
	ResidueName[ARG_NE]="ARG";
	ResidueName[ARG_HE]="ARG";
	ResidueName[ARG_CZ]="ARG";
	ResidueName[ARG_NH]="ARG";
	ResidueName[ARG_HH]="ARG";
	ResidueName[ARG_C]="ARG";
	ResidueName[ARG_O]="ARG";
	ResidueName[ASN_N]="ASN";
	ResidueName[ASN_HN]="ASN";
	ResidueName[ASN_CA]="ASN";
	ResidueName[ASN_HA]="ASN";
	ResidueName[ASN_CB]="ASN";
	ResidueName[ASN_HB]="ASN";
	ResidueName[ASN_CG]="ASN";
	ResidueName[ASN_OD]="ASN";
	ResidueName[ASN_ND]="ASN";
	ResidueName[ASN_HD]="ASN";
	ResidueName[ASN_C]="ASN";
	ResidueName[ASN_O]="ASN";
	ResidueName[ASP_N]="ASP";
	ResidueName[ASP_HN]="ASP";
	ResidueName[ASP_CA]="ASP";
	ResidueName[ASP_HA]="ASP";
	ResidueName[ASP_CB]="ASP";
	ResidueName[ASP_HB]="ASP";
	ResidueName[ASP_CG]="ASP";
	ResidueName[ASP_OD]="ASP";
	ResidueName[ASP_ND]="ASP";
	ResidueName[ASP_HD]="ASP";
	ResidueName[ASP_C]="ASP";
	ResidueName[ASP_O]="ASP";
	ResidueName[CYS_N]="CYS";
	ResidueName[CYS_HN]="CYS";
	ResidueName[CYS_CA]="CYS";
	ResidueName[CYS_HA]="CYS";
	ResidueName[CYS_CB]="CYS";
	ResidueName[CYS_HB]="CYS";
	ResidueName[CYS_SG]="CYS";
	ResidueName[CYS_HG]="CYS";
	ResidueName[CYS_C]="CYS";
	ResidueName[CYS_O]="CYS";
	ResidueName[GLU_N]="GLU";
	ResidueName[GLU_HN]="GLU";
	ResidueName[GLU_CA]="GLU";
	ResidueName[GLU_HA]="GLU";
	ResidueName[GLU_CB]="GLU";
	ResidueName[GLU_HB]="GLU";
	ResidueName[GLU_CG]="GLU";
	ResidueName[GLU_HG]="GLU";
	ResidueName[GLU_CD]="GLU";
	ResidueName[GLU_OE]="GLU";
	ResidueName[GLU_C]="GLU";
	ResidueName[GLU_O]="GLU";
	ResidueName[GLN_N]="GLN";
	ResidueName[GLN_HN]="GLN";
	ResidueName[GLN_CA]="GLN";
	ResidueName[GLN_HA]="GLN";
	ResidueName[GLN_CB]="GLN";
	ResidueName[GLN_HB]="GLN";
	ResidueName[GLN_CG]="GLN";
	ResidueName[GLN_HG]="GLN";
	ResidueName[GLN_CD]="GLN";
	ResidueName[GLN_OE]="GLN";
	ResidueName[GLN_NE]="GLN";
	ResidueName[GLN_HE]="GLN";
	ResidueName[GLN_C]="GLN";
	ResidueName[GLN_O]="GLN";
	ResidueName[GLY_N]="GLY";
	ResidueName[GLY_HN]="GLY";
	ResidueName[GLY_CA]="GLY";
	ResidueName[GLY_HA]="GLY";
	ResidueName[GLY_C]="GLY";
	ResidueName[GLY_O]="GLY";
	ResidueName[HIS_N]="HIS";
	ResidueName[HIS_HN]="HIS";
	ResidueName[HIS_CA]="HIS";
	ResidueName[HIS_HA]="HIS";
	ResidueName[HIS_CB]="HIS";
	ResidueName[HIS_HB]="HIS";
	ResidueName[HIS_ND]="HIS";
	ResidueName[HIS_HD1]="HIS";
	ResidueName[HIS_CG]="HIS";
	ResidueName[HIS_CE]="HIS";
	ResidueName[HIS_HE]="HIS";
	ResidueName[HIS_NE]="HIS";
	ResidueName[HIS_CD]="HIS";
	ResidueName[HIS_HD2]="HIS";
	ResidueName[HIS_C]="HIS";
	ResidueName[HIS_O]="HIS";
	ResidueName[ILE_N]="ILE";
	ResidueName[ILE_HN]="ILE";
	ResidueName[ILE_CA]="ILE";
	ResidueName[ILE_HA]="ILE";
	ResidueName[ILE_CB]="ILE";
	ResidueName[ILE_HB]="ILE";
	ResidueName[ILE_CG]="ILE";
	ResidueName[ILE_HG]="ILE";
	ResidueName[ILE_CD]="ILE";
	ResidueName[ILE_HD]="ILE";
	ResidueName[ILE_C]="ILE";
	ResidueName[ILE_O]="ILE";
	ResidueName[LEU_N]="LEU";
	ResidueName[LEU_HN]="LEU";
	ResidueName[LEU_CA]="LEU";
	ResidueName[LEU_HA]="LEU";
	ResidueName[LEU_CB]="LEU";
	ResidueName[LEU_HB]="LEU";
	ResidueName[LEU_CG]="LEU";
	ResidueName[LEU_HG]="LEU";
	ResidueName[LEU_CD]="LEU";
	ResidueName[LEU_CD1]="LEU";
	ResidueName[LEU_CD2]="LEU";
	ResidueName[LEU_HD]="LEU";
	ResidueName[LEU_HD1]="LEU";
	ResidueName[LEU_HD2]="LEU";
	ResidueName[LEU_C]="LEU";
	ResidueName[LEU_O]="LEU";
	ResidueName[LYS_N]="LYS";
	ResidueName[LYS_HN]="LYS";
	ResidueName[LYS_CA]="LYS";
	ResidueName[LYS_HA]="LYS";
	ResidueName[LYS_CB]="LYS";
	ResidueName[LYS_HB]="LYS";
	ResidueName[LYS_CG]="LYS";
	ResidueName[LYS_HG]="LYS";
	ResidueName[LYS_CD]="LYS";
	ResidueName[LYS_HD]="LYS";
	ResidueName[LYS_CE]="LYS";
	ResidueName[LYS_HE]="LYS";
	ResidueName[LYS_NZ]="LYS";
	ResidueName[LYS_HZ]="LYS";
	ResidueName[LYS_C]="LYS";
	ResidueName[LYS_O]="LYS";
	ResidueName[MET_N]="MET";
	ResidueName[MET_HN]="MET";
	ResidueName[MET_CA]="MET";
	ResidueName[MET_HA]="MET";
	ResidueName[MET_CB]="MET";
	ResidueName[MET_HB]="MET";
	ResidueName[MET_CG]="MET";
	ResidueName[MET_HG]="MET";
	ResidueName[MET_SD]="MET";
	ResidueName[MET_CE]="MET";
	ResidueName[MET_HE]="MET";
	ResidueName[MET_C]="MET";
	ResidueName[MET_O]="MET";
	ResidueName[PHE_N]="PHE";
	ResidueName[PHE_HN]="PHE";
	ResidueName[PHE_CA]="PHE";
	ResidueName[PHE_HA]="PHE";
	ResidueName[PHE_CB]="PHE";
	ResidueName[PHE_HB]="PHE";
	ResidueName[PHE_CG]="PHE";
	ResidueName[PHE_CD]="PHE";
	ResidueName[PHE_HD]="PHE";
	ResidueName[PHE_CE]="PHE";
	ResidueName[PHE_HE]="PHE";
	ResidueName[PHE_CZ]="PHE";
	ResidueName[PHE_HZ]="PHE";
	ResidueName[PHE_C]="PHE";
	ResidueName[PHE_O]="PHE";
	ResidueName[PRO_N]="PRO";
	ResidueName[PRO_CD]="PRO";
	ResidueName[PRO_HD]="PRO";
	ResidueName[PRO_CA]="PRO";
	ResidueName[PRO_HA]="PRO";
	ResidueName[PRO_CB]="PRO";
	ResidueName[PRO_HB]="PRO";
	ResidueName[PRO_CG]="PRO";
	ResidueName[PRO_HG]="PRO";
	ResidueName[PRO_C]="PRO";
	ResidueName[PRO_O]="PRO";
	ResidueName[SER_N]="SER";
	ResidueName[SER_HN]="SER";
	ResidueName[SER_CA]="SER";
	ResidueName[SER_HA]="SER";
	ResidueName[SER_CB]="SER";
	ResidueName[SER_HB]="SER";
	ResidueName[SER_OG]="SER";
	ResidueName[SER_HG]="SER";
	ResidueName[SER_C]="SER";
	ResidueName[SER_O]="SER";
	ResidueName[THR_N]="THR";
	ResidueName[THR_HN]="THR";
	ResidueName[THR_CA]="THR";
	ResidueName[THR_HA]="THR";
	ResidueName[THR_CB]="THR";
	ResidueName[THR_HB]="THR";
	ResidueName[THR_OG]="THR";
	ResidueName[THR_HG1]="THR";
	ResidueName[THR_CG]="THR";
	ResidueName[THR_HG2]="THR";
	ResidueName[THR_C]="THR";
	ResidueName[THR_O]="THR";
	ResidueName[TRP_N]="TRP";
	ResidueName[TRP_HN]="TRP";
	ResidueName[TRP_CA]="TRP";
	ResidueName[TRP_HA]="TRP";
	ResidueName[TRP_CB]="TRP";
	ResidueName[TRP_HB]="TRP";
	ResidueName[TRP_CG]="TRP";
	ResidueName[TRP_CD1]="TRP";
	ResidueName[TRP_HD]="TRP";
	ResidueName[TRP_NE]="TRP";
	ResidueName[TRP_HE1]="TRP";
	ResidueName[TRP_CE2]="TRP";
	ResidueName[TRP_CD2]="TRP";
	ResidueName[TRP_CE3]="TRP";
	ResidueName[TRP_HE3]="TRP";
	ResidueName[TRP_CZ3]="TRP";
	ResidueName[TRP_HZ3]="TRP";
	ResidueName[TRP_CZ2]="TRP";
	ResidueName[TRP_HZ2]="TRP";
	ResidueName[TRP_CH2]="TRP";
	ResidueName[TRP_HH]="TRP";
	ResidueName[TRP_C]="TRP";
	ResidueName[TRP_O]="TRP";
	ResidueName[TYR_N]="TYR";
	ResidueName[TYR_HN]="TYR";
	ResidueName[TYR_CA]="TYR";
	ResidueName[TYR_HA]="TYR";
	ResidueName[TYR_CB]="TYR";
	ResidueName[TYR_HB]="TYR";
	ResidueName[TYR_CG]="TYR";
	ResidueName[TYR_CD]="TYR";
	ResidueName[TYR_HD]="TYR";
	ResidueName[TYR_CE]="TYR";
	ResidueName[TYR_HE]="TYR";
	ResidueName[TYR_CZ]="TYR";
	ResidueName[TYR_OH]="TYR";
	ResidueName[TYR_HH]="TYR";
	ResidueName[TYR_C]="TYR";
	ResidueName[TYR_O]="TYR";
	ResidueName[VAL_N]="VAL";
	ResidueName[VAL_HN]="VAL";
	ResidueName[VAL_CA]="VAL";
	ResidueName[VAL_HA]="VAL";
	ResidueName[VAL_CB]="VAL";
	ResidueName[VAL_HB]="VAL";
	ResidueName[VAL_CG]="VAL";
	ResidueName[VAL_HG]="VAL";
	ResidueName[VAL_C]="VAL";
	ResidueName[VAL_O]="VAL";
	ResidueName[HSP_N]="HSP";
	ResidueName[HSP_HN]="HSP";
	ResidueName[HSP_CA]="HSP";
	ResidueName[HSP_HA]="HSP";
	ResidueName[HSP_CB]="HSP";
	ResidueName[HSP_HB]="HSP";
	ResidueName[HSP_CD]="HSP";
	ResidueName[HSP_HD2]="HSP";
	ResidueName[HSP_CG]="HSP";
	ResidueName[HSP_NE]="HSP";
	ResidueName[HSP_HE2]="HSP";
	ResidueName[HSP_ND]="HSP";
	ResidueName[HSP_HD1]="HSP";
	ResidueName[HSP_CE]="HSP";
	ResidueName[HSP_HE1]="HSP";
	ResidueName[HSP_C]="HSP";
	ResidueName[HSP_O]="HSP";
	ResidueName[HEM_FE]="HEM";
	ResidueName[HEM_NA]="HEM";
	ResidueName[HEM_NB]="HEM";
	ResidueName[HEM_NC]="HEM";
	ResidueName[HEM_ND]="HEM";
	ResidueName[HEM_C1]="HEM";
	ResidueName[HEM_C2]="HEM";
	ResidueName[HEM_C3]="HEM";
	ResidueName[HEM_C4]="HEM";
	ResidueName[HEM_CH]="HEM";
	ResidueName[HEM_HA]="HEM";
	ResidueName[HEM_HB]="HEM";
	ResidueName[HEM_HC]="HEM";
	ResidueName[HEM_HD]="HEM";
	ResidueName[HEM_CM]="HEM";
	ResidueName[HEM_HM]="HEM";
	ResidueName[HEM_CA]="HEM";
	ResidueName[HEM_CB]="HEM";
	ResidueName[HEM_CG]="HEM";
	ResidueName[HEM_O1]="HEM";
	ResidueName[HEM_O2]="HEM";
}

void AtomTypeToAtomName(string AtomName[])
{
	AtomName[UNK_ATOM_TYPE]="UNK";
	AtomName[HT]="UNK";
	AtomName[OT]="UNK";
	AtomName[ALA_N]="N";
	AtomName[ALA_HN]="HN";
	AtomName[ALA_CA]="CA";
	AtomName[ALA_HA]="HA";
	AtomName[ALA_CB]="CB";
	AtomName[ALA_HB]="HB";
	AtomName[ALA_C]="C";
	AtomName[ALA_O]="O";
	AtomName[ARG_N]="N";
	AtomName[ARG_HN]="HN";
	AtomName[ARG_CA]="CA";
	AtomName[ARG_HA]="HA";
	AtomName[ARG_CB]="CB";
	AtomName[ARG_HB]="HB";
	AtomName[ARG_CG]="CG";
	AtomName[ARG_HG]="HG";
	AtomName[ARG_CD]="CD";
	AtomName[ARG_HD]="HD";
	AtomName[ARG_NE]="NE";
	AtomName[ARG_HE]="HE";
	AtomName[ARG_CZ]="CZ";
	AtomName[ARG_NH]="NH";
	AtomName[ARG_HH]="HH";
	AtomName[ARG_C]="C";
	AtomName[ARG_O]="O";
	AtomName[ASN_N]="N";
	AtomName[ASN_HN]="HN";
	AtomName[ASN_CA]="CA";
	AtomName[ASN_HA]="HA";
	AtomName[ASN_CB]="CB";
	AtomName[ASN_HB]="HB";
	AtomName[ASN_CG]="CG";
	AtomName[ASN_OD]="OD";
	AtomName[ASN_ND]="ND";
	AtomName[ASN_HD]="HD";
	AtomName[ASN_C]="C";
	AtomName[ASN_O]="O";
	AtomName[ASP_N]="N";
	AtomName[ASP_HN]="HN";
	AtomName[ASP_CA]="CA";
	AtomName[ASP_HA]="HA";
	AtomName[ASP_CB]="CB";
	AtomName[ASP_HB]="HB";
	AtomName[ASP_CG]="CG";
	AtomName[ASP_OD]="OD";
	AtomName[ASP_ND]="ND";
	AtomName[ASP_HD]="HD";
	AtomName[ASP_C]="C";
	AtomName[ASP_O]="O";
	AtomName[CYS_N]="N";
	AtomName[CYS_HN]="HN";
	AtomName[CYS_CA]="CA";
	AtomName[CYS_HA]="HA";
	AtomName[CYS_CB]="CB";
	AtomName[CYS_HB]="HB";
	AtomName[CYS_SG]="SG";
	AtomName[CYS_HG]="HG";
	AtomName[CYS_C]="C";
	AtomName[CYS_O]="O";
	AtomName[GLU_N]="N";
	AtomName[GLU_HN]="HN";
	AtomName[GLU_CA]="CA";
	AtomName[GLU_HA]="HA";
	AtomName[GLU_CB]="CB";
	AtomName[GLU_HB]="HB";
	AtomName[GLU_CG]="CG";
	AtomName[GLU_HG]="HG";
	AtomName[GLU_CD]="CD";
	AtomName[GLU_OE]="OE";
	AtomName[GLU_C]="C";
	AtomName[GLU_O]="O";
	AtomName[GLN_N]="N";
	AtomName[GLN_HN]="HN";
	AtomName[GLN_CA]="CA";
	AtomName[GLN_HA]="HA";
	AtomName[GLN_CB]="CB";
	AtomName[GLN_HB]="HB";
	AtomName[GLN_CG]="CG";
	AtomName[GLN_HG]="HG";
	AtomName[GLN_CD]="CD";
	AtomName[GLN_OE]="OE";
	AtomName[GLN_NE]="NE";
	AtomName[GLN_HE]="HE";
	AtomName[GLN_C]="C";
	AtomName[GLN_O]="O";
	AtomName[GLY_N]="N";
	AtomName[GLY_HN]="HN";
	AtomName[GLY_CA]="CA";
	AtomName[GLY_HA]="HA";
	AtomName[GLY_C]="C";
	AtomName[GLY_O]="O";
	AtomName[HIS_N]="N";
	AtomName[HIS_HN]="HN";
	AtomName[HIS_CA]="CA";
	AtomName[HIS_HA]="HA";
	AtomName[HIS_CB]="CB";
	AtomName[HIS_HB]="HB";
	AtomName[HIS_ND]="ND";
	AtomName[HIS_HD1]="HD1";
	AtomName[HIS_CG]="CG";
	AtomName[HIS_CE]="CE";
	AtomName[HIS_HE]="HE";
	AtomName[HIS_NE]="NE";
	AtomName[HIS_CD]="CD";
	AtomName[HIS_HD2]="HD2";
	AtomName[HIS_C]="C";
	AtomName[HIS_O]="O";
	AtomName[ILE_N]="N";
	AtomName[ILE_HN]="HN";
	AtomName[ILE_CA]="CA";
	AtomName[ILE_HA]="HA";
	AtomName[ILE_CB]="CB";
	AtomName[ILE_HB]="HB";
	AtomName[ILE_CG]="CG";
	AtomName[ILE_HG]="HG";
	AtomName[ILE_CD]="CD";
	AtomName[ILE_HD]="HD";
	AtomName[ILE_C]="C";
	AtomName[ILE_O]="O";
	AtomName[LEU_N]="N";
	AtomName[LEU_HN]="HN";
	AtomName[LEU_CA]="CA";
	AtomName[LEU_HA]="HA";
	AtomName[LEU_CB]="CB";
	AtomName[LEU_HB]="HB";
	AtomName[LEU_CG]="CG";
	AtomName[LEU_HG]="HG";
	AtomName[LEU_CD]="CD";
	AtomName[LEU_CD1]="CD1";
	AtomName[LEU_CD2]="CD2";
	AtomName[LEU_HD]="HD";
	AtomName[LEU_HD1]="HD1";
	AtomName[LEU_HD2]="HD2";
	AtomName[LEU_C]="C";
	AtomName[LEU_O]="O";
	AtomName[LYS_N]="N";
	AtomName[LYS_HN]="HN";
	AtomName[LYS_CA]="CA";
	AtomName[LYS_HA]="HA";
	AtomName[LYS_CB]="CB";
	AtomName[LYS_HB]="HB";
	AtomName[LYS_CG]="CG";
	AtomName[LYS_HG]="HG";
	AtomName[LYS_CD]="CD";
	AtomName[LYS_HD]="HD";
	AtomName[LYS_CE]="CE";
	AtomName[LYS_HE]="HE";
	AtomName[LYS_NZ]="NZ";
	AtomName[LYS_HZ]="HZ";
	AtomName[LYS_C]="C";
	AtomName[LYS_O]="O";
	AtomName[MET_N]="N";
	AtomName[MET_HN]="HN";
	AtomName[MET_CA]="CA";
	AtomName[MET_HA]="HA";
	AtomName[MET_CB]="CB";
	AtomName[MET_HB]="HB";
	AtomName[MET_CG]="CG";
	AtomName[MET_HG]="HG";
	AtomName[MET_SD]="SD";
	AtomName[MET_CE]="CE";
	AtomName[MET_HE]="HE";
	AtomName[MET_C]="C";
	AtomName[MET_O]="O";
	AtomName[PHE_N]="N";
	AtomName[PHE_HN]="HN";
	AtomName[PHE_CA]="CA";
	AtomName[PHE_HA]="HA";
	AtomName[PHE_CB]="CB";
	AtomName[PHE_HB]="HB";
	AtomName[PHE_CG]="CG";
	AtomName[PHE_CD]="CD";
	AtomName[PHE_HD]="HD";
	AtomName[PHE_CE]="CE";
	AtomName[PHE_HE]="HE";
	AtomName[PHE_CZ]="CZ";
	AtomName[PHE_HZ]="HZ";
	AtomName[PHE_C]="C";
	AtomName[PHE_O]="O";
	AtomName[PRO_N]="N";
	AtomName[PRO_CD]="CD";
	AtomName[PRO_HD]="HD";
	AtomName[PRO_CA]="CA";
	AtomName[PRO_HA]="HA";
	AtomName[PRO_CB]="CB";
	AtomName[PRO_HB]="HB";
	AtomName[PRO_CG]="CG";
	AtomName[PRO_HG]="HG";
	AtomName[PRO_C]="C";
	AtomName[PRO_O]="O";
	AtomName[SER_N]="N";
	AtomName[SER_HN]="HN";
	AtomName[SER_CA]="CA";
	AtomName[SER_HA]="HA";
	AtomName[SER_CB]="CB";
	AtomName[SER_HB]="HB";
	AtomName[SER_OG]="OG";
	AtomName[SER_HG]="HG";
	AtomName[SER_C]="C";
	AtomName[SER_O]="O";
	AtomName[THR_N]="N";
	AtomName[THR_HN]="HN";
	AtomName[THR_CA]="CA";
	AtomName[THR_HA]="HA";
	AtomName[THR_CB]="CB";
	AtomName[THR_HB]="HB";
	AtomName[THR_OG]="OG";
	AtomName[THR_HG1]="HG1";
	AtomName[THR_CG]="CG";
	AtomName[THR_HG2]="HG2";
	AtomName[THR_C]="C";
	AtomName[THR_O]="O";
	AtomName[TRP_N]="N";
	AtomName[TRP_HN]="HN";
	AtomName[TRP_CA]="CA";
	AtomName[TRP_HA]="HA";
	AtomName[TRP_CB]="CB";
	AtomName[TRP_HB]="HB";
	AtomName[TRP_CG]="CG";
	AtomName[TRP_CD1]="CD1";
	AtomName[TRP_HD]="HD";
	AtomName[TRP_NE]="NE";
	AtomName[TRP_HE1]="HE1";
	AtomName[TRP_CE2]="CE2";
	AtomName[TRP_CD2]="CD2";
	AtomName[TRP_CE3]="CE3";
	AtomName[TRP_HE3]="HE3";
	AtomName[TRP_CZ3]="CZ3";
	AtomName[TRP_HZ3]="HZ3";
	AtomName[TRP_CZ2]="CZ2";
	AtomName[TRP_HZ2]="HZ2";
	AtomName[TRP_CH2]="CH2";
	AtomName[TRP_HH]="HH";
	AtomName[TRP_C]="C";
	AtomName[TRP_O]="O";
	AtomName[TYR_N]="N";
	AtomName[TYR_HN]="HN";
	AtomName[TYR_CA]="CA";
	AtomName[TYR_HA]="HA";
	AtomName[TYR_CB]="CB";
	AtomName[TYR_HB]="HB";
	AtomName[TYR_CG]="CG";
	AtomName[TYR_CD]="CD";
	AtomName[TYR_HD]="HD";
	AtomName[TYR_CE]="CE";
	AtomName[TYR_HE]="HE";
	AtomName[TYR_CZ]="CZ";
	AtomName[TYR_OH]="OH";
	AtomName[TYR_HH]="HH";
	AtomName[TYR_C]="C";
	AtomName[TYR_O]="O";
	AtomName[VAL_N]="N";
	AtomName[VAL_HN]="HN";
	AtomName[VAL_CA]="CA";
	AtomName[VAL_HA]="HA";
	AtomName[VAL_CB]="CB";
	AtomName[VAL_HB]="HB";
	AtomName[VAL_CG]="CG";
	AtomName[VAL_HG]="HG";
	AtomName[VAL_C]="C";
	AtomName[VAL_O]="O";
	AtomName[HSP_N]="N";
	AtomName[HSP_HN]="HN";
	AtomName[HSP_CA]="CA";
	AtomName[HSP_HA]="HA";
	AtomName[HSP_CB]="CB";
	AtomName[HSP_HB]="HB";
	AtomName[HSP_CD]="CD";
	AtomName[HSP_HD2]="HD2";
	AtomName[HSP_CG]="CG";
	AtomName[HSP_NE]="NE";
	AtomName[HSP_HE2]="HE2";
	AtomName[HSP_ND]="ND";
	AtomName[HSP_HD1]="HD1";
	AtomName[HSP_CE]="CE";
	AtomName[HSP_HE1]="HE1";
	AtomName[HSP_C]="C";
	AtomName[HSP_O]="O";
	AtomName[HEM_FE]="FE";
	AtomName[HEM_NA]="NA";
	AtomName[HEM_NB]="NB";
	AtomName[HEM_NC]="NC";
	AtomName[HEM_ND]="ND";
	AtomName[HEM_C1]="C1";
	AtomName[HEM_C2]="C2";
	AtomName[HEM_C3]="C3";
	AtomName[HEM_C4]="C4";
	AtomName[HEM_CH]="CH";
	AtomName[HEM_HA]="HA";
	AtomName[HEM_HB]="HB";
	AtomName[HEM_HC]="HC";
	AtomName[HEM_HD]="HD";
	AtomName[HEM_CM]="CM";
	AtomName[HEM_HM]="HM";
	AtomName[HEM_CA]="CA";
	AtomName[HEM_CB]="CB";
	AtomName[HEM_CG]="CG";
	AtomName[HEM_O1]="O1";
	AtomName[HEM_O2]="O2";
}

void GetAtomType(vector<AtomStruct> &Atoms)
{
	int i, natom=Atoms.size();
	for (i=0;i<natom;i++)
	{
		if (Atoms[i].ResidueName!="TIP3")
		{
			Atoms[i].AtomType=GetSolventType(Atoms[i].ResidueName, Atoms[i].AtomName);
		}
	}
}

#endif
