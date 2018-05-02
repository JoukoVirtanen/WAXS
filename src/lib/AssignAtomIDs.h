#ifndef _AssignAtomIDs_included_
#define _AssignAtomIDs_included_

# include <vector>

# include "Structures.h"
# include "AtomIDs.h"
# include "StringUtils.h"
# include "WritePdb.h"

void AssignAtomIDs(vector<AtomStruct> &Atoms)
{
	string ElementID, AtomName;
	int n, natom;

	natom=Atoms.size();

	for (n=0;n<natom;n++)
	{
		Trim(Atoms[n].AtomName);
		ElementID=Atoms[n].AtomName.substr(0,1);
		if (ElementID=="H") Atoms[n].atomid=HYDROGEN;
		else if (ElementID=="C" && Atoms[n].AtomName!="CLA") Atoms[n].atomid=CARBON;
		else if (ElementID=="N") Atoms[n].atomid=NITROGEN;
		else if (ElementID=="O") Atoms[n].atomid=OXYGEN;
		else if (ElementID=="S" && Atoms[n].ResidueName!="SOD") Atoms[n].atomid=SULFUR;
		else if (ElementID=="F") Atoms[n].atomid=IRON;
		else if (ElementID=="P") Atoms[n].atomid=PHOSPHORUS;
		else if (Atoms[n].ResidueName=="SOD" || Atoms[n].ResidueName=="NA0") Atoms[n].atomid=NA_Plus;
		else if (Atoms[n].AtomName=="IP") Atoms[n].atomid=NA_Plus;
		else if (Atoms[n].AtomName=="CLA") Atoms[n].atomid=CL;
                else if (Atoms[n].AtomName=="IM") Atoms[n].atomid=CL;
                else if (Atoms[n].AtomName=="MG") Atoms[n].atomid=MG;
                else if (Atoms[n].AtomName=="M") Atoms[n].atomid=MG;
                else
		{
			cout <<"Warning unrecognized atom "<<ElementID<<" will be treated as a hydrogen."<<endl;
			PrintAtomInfo(Atoms[n]);
                        Atoms[n].atomid=HYDROGEN;
		}
	}
	return;
}

string AtomIDToStr(int AtomID)
{
        if (AtomID==HYDROGEN) return "H";
        else if (AtomID==CARBON) return "C";
        else if (AtomID==NITROGEN) return "N";
        else if (AtomID==OXYGEN) return "O";
        else if (AtomID==SULFUR) return "S";
        else if (AtomID==IRON) return "FE";
        else if (AtomID==ExcludedVolume) return "ExcludedVolume";
        else if (AtomID==HydrationShell) return "HydrationShell";
        else if (AtomID==PHOSPHORUS) return "P";
        else if (AtomID==CL) return "CL";
        else if (AtomID==NA_Plus) return "NA";
        else if (AtomID==WaterSphere) return "WaterSphere";
        else if (AtomID==CUBE) return "CUBE";
        else return "UNK";
}
#endif

