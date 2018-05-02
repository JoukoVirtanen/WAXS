#ifndef _ReadPRM_included_
#define _ReadPRM_included_

# include "StringUtils.h"
# include "LinkedList.h"
# include "IOUtils.h"
# include "ResidueCode.h"
# include "ReadPdb.h"
# include "WritePdb.h"
# include "Structures.h"

void ReadPrmFile(string PrmFile, vector<AtomStruct> &TempCharmmAtoms)
{
	bool verbose=false;
	string line;
	vector<string> str;
	AtomStruct CharmmAtom;
	cout <<"In ReadPrmFile"<<endl;
	fstream File;
	OpenFile(PrmFile, File, "PrmFile");
	//cout <<"Opened file "<<PrmFile<<endl;
	getline(File, line);
	while (true)
	{
		//cout <<line<<endl;
		if (File.eof()) break;
		if (line.substr(0,9)=="NONBONDED")
		{
			while (true)
			{
				if (File.eof()) break;
				if (line.substr(0,5)=="HBOND") break;
				Trim(line);
				DeleteStartingAt(line, "!");
				if (line.substr(0,1)!="!" && line!="" && line.substr(0,5)!="cutnb" && line.substr(0,9)!="NONBONDED")
				{
					DeleteStartingAt(line, "!");
					Tokenize2(line, " ", str);
					if (str.size()>3)
					{
						//cout <<"After Tokenize2"<<endl;
						//CharmmAtom.CharmmAtomName=line.substr(0,4);
						CharmmAtom.CharmmAtomName=str[0];
                                                CharmmAtom.epsilon=-StrToFloat(str[2]);
						CharmmAtom.vdw=StrToFloat(str[3]);
                                                //cout <<"epsilon= "<<CharmmAtom.epsilon<<endl;
						//cout <<"vdw= "<<CharmmAtom.vdw<<endl;
						if (str.size()>6)
						{
							CharmmAtom.epsilon14=-StrToFloat(str[5]);
							//cout <<"epsilon14= "<<CharmmAtom.epsilon14<<endl;
							CharmmAtom.vdw14=StrToFloat(str[6]);
							//cout <<"vdw14= "<<CharmmAtom.vdw14<<endl;
						}
						Trim(CharmmAtom.CharmmAtomName);
						TempCharmmAtoms.push_back(CharmmAtom);
					}
				}
				getline(File, line);
			}
		}
		getline(File, line);
	}
	File.close();
}
#endif
