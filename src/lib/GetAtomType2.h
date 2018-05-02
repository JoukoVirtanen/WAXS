#ifndef _GetAtomTypes2_included
#define _GetAtomTypes2_included

# include "AtomCode.h"
# include "StringUtils.h"
# include "LinkedList.h"

int GetSolventType2(string resstr, string str)
{
        bool verbose=true;
        string str1, str2, str2b, str3, str_trim;
        Trim(str);
        str1=str[0];
        str2=str.substr(0,2);
        str3=str.substr(0,3);
        if (str1=="H")
        {

                if (str2=="HA") return HC;
                else if (str2=="HB") return HC;
                else if (str2=="HN") return HN;
                else if (str2=="HT") return HN;
                else if (str3=="H5T") return HO;
                else if (str3=="H3T") return HO;
                else if (resstr=="ARG")
                {
                        if (str2=="HG") return HC;
                        else if (str2=="HD") return HC;
                        else if (str2=="HE") return HN;
                        else if (str2=="HH") return HN;
                        else
                        {
                                cout <<"Error unable to identify solvent type"<<endl;
                                cout <<resstr<<"	"<<str<<endl;
                                return HC;
                        }
                }

                else if (resstr=="ASN")
                {
                        if (str2=="HD") return HN;
                        else
                        {
                                cout <<"Error unable to identify solvent type"<<endl;
                                cout <<resstr<<"	"<<str<<endl;
                                return HC;
                        }
                }

                else if (resstr=="CYS")
                {
                        if (str2=="HG") return HS;
                        else 
                        {
                                cout <<"Error unable to identify solvent type"<<endl;
                                cout <<resstr<<"	"<<str<<endl;
                                return HC;
                        }
                }
                else if (resstr=="GLN")
                {
                        if (str2=="HG") return HC;
                        else if (str2=="HE") return HN;
                        else
                        {
                                cout <<"Error unable to identify solvent type"<<endl;
                                cout <<resstr<<"	"<<str<<endl;
                                return 0;
                        }
                }

                else if (resstr=="GLU")
                {
                        if (str2=="HG") return HC;
                        else
                        {
                                cout <<"Error unable to identify solvent type"<<endl;
                                cout <<resstr<<"	"<<str<<endl;
                                return HC;
                        }
                }

                else if (resstr=="HSP" || resstr=="HSE" || resstr=="HSD" || resstr=="HIE")
                {
                        if (str3=="HD2") return HC;
                        else if (str3=="HE2") return HN;
                        else if (str3=="HD1") return HN;
                        else if (str3=="HE1") return HC;
                        else
                        {
                                cout <<"Error unable to identify solvent type"<<endl;
                                cout <<resstr<<"	"<<str<<endl;
                                return HC;
                        }
                }

                else if (resstr=="ILE")
                {
                        if (str2=="NH") return HN;
                        else if (str2=="HG") return HC;
                        else if (str2=="HD") return HC;
                        else
                        {
                                cout <<"Error unable to identify solvent type"<<endl;
                                cout <<resstr<<"	"<<str<<endl;
                                return HC;
                        }
                }

                else if (resstr=="LEU")
                {
                        if (str2!="NH") return HC;
                        else
                        {
                                cout <<"Error unable to identify solvent type"<<endl;
                                cout <<resstr<<"	"<<str<<endl;
                                return HC;
                        }
                }

                else if (resstr=="LYS")
                {
                        if (str2=="HG") return HC;
                        else if (str2=="HD") return HC;
                        else if (str2=="HE") return HC;
                        else if (str2=="HZ") return HN;
                        else
                        {
                                cout <<"Error unable to identify solvent type"<<endl;
                                cout <<resstr<<"	"<<str<<endl;
                                return HC;
                        }
                }

                else if (resstr=="MET")
                {
                        if (str2=="HT") return HN;
                        else if (str2=="HA") return HC;
                        else if (str2=="HB") return HC;
                        else if (str2=="HE") return HC;
                        else if (str2=="HG") return HC;
                        //Needs to be fixed.  Deleted CE.
                        else
                        {
                                cout <<"Error unable to identify solvent type"<<endl;
                                cout <<resstr<<"	"<<str<<endl;
                                return HC;
                        }
                }

                else if (resstr=="PHE")
                {
                        if (str2=="HN") return HN;
                        else if (str2=="HD") return HC;
                        else if (str2=="HE") return HC;
                        else if (str2=="HZ") return HC;
                        else
                        {
                                cout <<"Error unable to identify solvent type"<<endl;
                                cout <<resstr<<"	"<<str<<endl;
                                return HC;
                        }
                }

                else if (resstr=="PRO")
                {
                        if (str2!="HN") return HC;
                        else
                        {
                                cout <<"Error unable to identify solvent type"<<endl;
                                cout <<resstr<<"	"<<str<<endl;
                                return HC;
                        }
                }

                else if (resstr=="THR")
                {
                        if (str3=="HG1") return HO;
                        else if (str3=="HG2") return HC;
                        else
                        {
                                cout <<"Error unable to identify solvent type"<<endl;
                                cout <<resstr<<"	"<<str<<endl;
                                return HC;
                        }
                }

                else if (resstr=="VAL")
                {
                        if (str2!="HN") return HC;
                        else
                        {
                                cout <<"Error unable to identify solvent type"<<endl;
                                cout <<resstr<<"	"<<str<<endl;
                                return HC;
                        }
                }

                else if (resstr=="SER")
                {
                        if (str2=="HG") return HO;
                        else
                        {
                                cout <<"Error unable to identify solvent type"<<endl;
                                cout <<resstr<<"	"<<str<<endl;
                                return HC;
                        }
                }

                else if (resstr=="TRP")
                {
                        if (str2=="HA") return HC;
                        else if (str2=="HD") return HC;
                        else if (str3=="HE1") return HN;
                        else if (str3=="HE3") return HC;
                        else if (str2=="HZ") return HC;
                        else if (str2=="HH") return HC;
                        else
                        {
                                cout <<"Error unable to identify solvent type"<<endl;
                                cout <<resstr<<"	"<<str<<endl;
                                return 0;
                        }
                }

                else if (resstr=="TYR")
                {
                        if (str2=="HD") return HC;
                        else if (str2=="HE") return HC;
                        else if (str2=="HH") return HO;
                        else
                        {
                                cout <<"Error unable to identify solvent type"<<endl;
                                cout <<resstr<<"	"<<str<<endl;
                                return HC;
                        }
                }
                else if (resstr=="HEM")
                {
                        if (str2=="HC") return HC;
                        else if (str2=="HD") return HC;
                        else if (str2=="HM") return HC;
                        else
                        {
                                cout <<"Error unable to identify solvent type"<<endl;
                                cout <<resstr<<"	"<<str<<endl;
                                return HC;
                        }
                }
                else if (resstr=="GUA")
                {
                        //Check to make sure that this is correct
                        if (str=="H5'") return HC;
                        else if (str=="H5''") return HC;
                        else if (str=="H3'") return HC;
                        else if (str=="H4'") return HC;
                        else if (str=="H1'") return HC;
                        else if (str=="H8") return HC;
                        else if (str=="H1") return HN;
                        else if (str=="H21") return HN;
                        else if (str=="H22") return HN;
                        else if (str=="H3") return HC;
                        else if (str=="H2''") return HC;
                        else if (str=="H2'") return HO;
                        else
                        {

                                if (verbose) cout <<"Error unable to identify solvent type.  Residue= "<<resstr<<" str= "<<str<<endl;
                                return UNK_ATOM_TYPE;
                        }
                }
                else if (resstr=="URA")
                {
                        if (str=="H5'") return HC;
                        else if (str=="H5''") return HC;
                        else if (str=="H4'") return HC;
                        else if (str=="H1'") return HC;
                        else if (str=="H6") return HC;
                        else if (str=="H5") return HC;
                        else if (str=="H3") return HN;
                        else if (str=="H3'") return HC;
                        else if (str=="H2''") return HC;
                        else if (str=="H2'") return HO;
                        else
                        {

                                if (verbose) cout <<"Error unable to identify solvent type.  Residue= "<<resstr<<" str= "<<str<<endl;
                                return UNK_ATOM_TYPE;
                        }
                }
                else if (resstr=="ADE")
                {
                        if (str=="H5'") return HC;
                        else if (str=="H5''") return HC;
                        else if (str=="H4'") return HC;
                        else if (str=="H1'") return HC;
                        else if (str=="H8") return HC;
                        else if (str=="H61") return HN;
                        else if (str=="H62") return HN;
                        else if (str=="H2") return HC;
                        else if (str=="H3'") return HC;
                        else if (str=="H2''") return HC;
                        else if (str=="H2'") return HO;
                        else
                        {

                                if (verbose) cout <<"Error unable to identify solvent type.  Residue= "<<resstr<<" str= "<<str<<endl;
                                return UNK_ATOM_TYPE;
                        }
                }
                else if (resstr=="CYT")
                {
                        if (str=="H5'") return HC;
                        else if (str=="H5''") return HC;
                        else if (str=="H4'") return HC;
                        else if (str=="H1'") return HC;
                        else if (str=="H6") return HC;
                        else if (str=="H5") return HC;
                        else if (str=="H41") return HN;
                        else if (str=="H42") return HN;
                        else if (str=="H3'") return HC;
                        else if (str=="H2''") return HC;
                        else if (str=="H2'") return HO;
                        else
                        {

                                if (verbose) cout <<"Error unable to identify solvent type.  Residue= "<<resstr<<" str= "<<str<<endl;
                                return UNK_ATOM_TYPE;
                        }
                }
                else
                {
                        if (verbose) cout <<"Error unable to identify atomtype2."
                                <<resstr<<"	"<<str<<endl;
                        return UNK_ATOM_TYPE;
                }
        }

        else if (str1=="C") return C;
        else if (str1=="N") return N;
        else if (str1=="O") return O;
        else if (str1=="S") return S;
        else if (str1=="P") return P;
        else
        {
                if (verbose) cout <<"Error unable to identify solvent type"<<endl;
                if (verbose) cout <<resstr<<"	"<<str<<endl;
                return UNK_ATOM_TYPE;
        }

}

void ElementTypeToElementName(string ElementName[])
{
        ElementName[HC]="HC";
        ElementName[HN]="HN";
        ElementName[HO]="HO";
        ElementName[HS]="HS";
        ElementName[C]="C";
        ElementName[N]="N";
        ElementName[O]="O";
        ElementName[S]="S";
}

void GetAtomType2(vector<AtomStruct> &Atoms)
{
        int i, natom=Atoms.size();
        for (i=0;i<natom;i++)
        {
                if (Atoms[i].ResidueName!="TIP3")
                {
                        Atoms[i].AtomType2=GetSolventType2(Atoms[i].ResidueName, Atoms[i].AtomName);
                }
        }
}

#endif
