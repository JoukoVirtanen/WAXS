#ifndef _ReadGofRFile_included_
#define _ReadGofRFile_included_

# include "IOUtils.h"
# include "VectorManip.h"
# include "TypeDef.h"
# include "IOUtils.h"

const Real UNK_DENSITY=-666.0;

void ReadGofRFile(string GofRFile, Real contrast, Matrix &DensityMatrix, Real RecBin)
{
	bool valid;
	char CharLine[1000], chargofrfile[1000];
	char *Space;
        int previousDistance=0;
	int NumSolventTypes, MaxRecord;
	int distance, n;
	string line, str;
	Real dist;
        vector<Real> tempDensity;

	cout <<"In ReadGofRFile"<<endl;

	NumSolventTypes=2000;
	MaxRecord=1000;
	strcpy(chargofrfile, GofRFile.c_str());
	valid=exist(chargofrfile);
	if (!valid)
	{
		cout <<"Error. GofRFile "<<GofRFile<<" cannot be opened."<<endl;
		exit(EXIT_FAILURE);
	}

	fstream gofr;
	gofr.open(chargofrfile,ios::in);
	cout <<"chargofrfile= "<<chargofrfile<<endl;
	getline(gofr, line);
	while (true)
	{
		if (line.substr(0,3)=="END") break;
		if (gofr.eof()) break;
		if (line.substr(0,15)=="NumSolventTypes")
		{
			NumSolventTypes=StrToInt(line.substr(16,10));
			cout <<"NumSolventTypes= "<<NumSolventTypes<<endl;
		}

		if (line.substr(0,9)=="MaxRecord")
		{
			MaxRecord=StrToInt(line.substr(10,10));
			cout <<"NumSolventTypes= "<<NumSolventTypes<<endl;
			cout <<"MaxRecord= "<<MaxRecord<<endl;
			CreateMatrix(DensityMatrix, NumSolventTypes, 0);
			cout <<"After CreateMatrix"<<endl;
		}
		//cout <<"After Initialized DensityMatrix"<<endl;
		if (line.substr(0,8)=="AtomType")
		{
			n=StrToInt(line.substr(9,15));
			str=line.substr(0,10);
			distance=0;
			//cout <<line<<endl;
			//cout <<line.substr(9,15)<<endl;
			//cout <<"n= "<<n<<endl;
                        previousDistance=0;
			while (true)
			{
				getline(gofr,line);
				if (line.substr(0,1)==" ") cout <<"Break found"<<endl;
				if (line.substr(0,8)=="AtomType")
				{
					n++;
					break;
				}
				strcpy(CharLine, line.c_str());
				Space=strchr(CharLine,' ');
				if (Space==NULL) break;
				if (line.substr(0,1)==" " || line.substr(0,3)=="END") break;
				line=" " + line + " ";
				dist=StrToFloat(GetWord(line, 1));
				distance=int(dist/RecBin+0.5);
				str=GetWord(line, 2);
                                tempDensity.resize(distance+1);
                                for (int i=previousDistance+1;i<distance;i++) tempDensity[i]=UNK_DENSITY;
				tempDensity[distance]=StrToFloat(str)-contrast;
                                previousDistance=distance;
				//cout <<dist<<" "<<str<<endl;
				//cout <<"DensityMatrix["<<n<<"]["<<distance<<"]= "<<DensityMatrix[n][distance]<<endl;
			}
                        DensityMatrix[n]=tempDensity;
                        tempDensity.clear();
		}
		getline(gofr, line);
	}
	for (n=0;n<NumSolventTypes;n++)
	{
		for (distance=0;distance<MaxRecord;distance++)
		{
			//cout <<"DensityMatrix["<<n<<"]["<<distance<<"]= "<<DensityMatrix[n][distance]<<endl;
		}
	}
	gofr.close();	
        cout <<"Leaving ReadGofRFile"<<endl;	
}
/*
void ReadGofRFile(string GofRFile, Real contrast, PrdfStruct &prdf)
{
        fstream gofr;
	//char *Space;
        int previousDistance=0;
	int NumSolventTypes, MaxRecord;
	int distance, n;
        vector<string> str;
	string line;
	Real dist;
        vector<Real> tempDensity;

	cout <<"In ReadGofRFile"<<endl;

	NumSolventTypes=2000;
	MaxRecord=1000;
	
        OpenFile(GofRFile, gofr, "prdfFile");
        
	cout <<"GofRFile= "<<GofRFile<<endl;
	getline(gofr, line);
	while (true)
	{
		if (line.substr(0,3)=="END") break;
		if (gofr.eof()) break;
		if (line.substr(0,15)=="NumSolventTypes")
		{
			NumSolventTypes=StrToInt(line.substr(16,10));
			cout <<"NumSolventTypes= "<<NumSolventTypes<<endl;
		}

		if (line.substr(0,9)=="MaxRecord")
		{
			MaxRecord=StrToInt(line.substr(10,10));
			cout <<"NumSolventTypes= "<<NumSolventTypes<<endl;
			cout <<"MaxRecord= "<<MaxRecord<<endl;
			Safe2DAlloc(prdf.density, NumSolventTypes, 0, "density");
			//CreateMatrix(prdf.density, NumSolventTypes, 0);
                        cout <<"After CreateMatrix"<<endl;
                        prdf.AngularDependence=getValueOfBool("AngularDependence", line);
                        cout <<"After AngularDependence"<<endl;
                        prdf.AngularDependence2=getValueOfBool("AngularDependence2", line);
                        prdf.ConvexConcave=getValueOfBool("ConvexConcave", line);
                        prdf.ConvexConcave2=getValueOfBool("ConvexConcave2", line);
                        prdf.PhobicPhilic=getValueOfBool("PhobicPhilic", line);
                        prdf.SecondNearestNeighbor=getValueOfBool("SecondNearestNeighbor", line);
                        prdf.UseChargeRadii=getValueOfBool("UseChargeRadii", line);
                        prdf.AngleBin=getValueOfReal("AngleBin", line);
                        prdf.AngularCutOff=getValueOfReal("AngularCutOff", line);
                        prdf.ChargeRadiiScale=getValueOfReal("ChargeRadiiScale", line);
                        prdf.CubeSize=getValueOfReal("CubeSize", line);
                        prdf.HydrationRadiusScale=getValueOfReal("HydrationRadiusScale", line);
                        prdf.RecBin=getValueOfReal("RecBin", line);
			cout <<"After RecBin"<<endl;
		}
		//cout <<"After Initialized DensityMatrix"<<endl;
		if (line.substr(0,8)=="AtomType")
		{
			n=StrToInt(line.substr(9,15));
			//str=line.substr(0,10);
			distance=0;
			//cout <<line<<endl;
			//cout <<line.substr(9,15)<<endl;
			//cout <<"n= "<<n<<endl;
                        previousDistance=0;
			while (true)
			{
				getline(gofr,line);
				if (line.substr(0,1)==" ") cout <<"Break found"<<endl;
				if (line.substr(0,8)=="AtomType")
				{
					n++;
					break;
				}
				//strcpy(CharLine, line.c_str());
				//Space=strchr(CharLine,' ');
				//if (Space==NULL) break;
				if (line.substr(0,1)==" " || line.substr(0,3)=="END") break;
                                //line+=" ";
                                cout <<line<<endl;
                                Tokenize2(line, " ", str);
				cout <<"str.size()= "<<str.size()<<endl;
                                if (str.size()>2)
                                {
                                        dist=StrToFloat(str[0]);
				        distance=int(dist/prdf.RecBin+0.5);
                                        tempDensity.resize(distance+1);
                                        for (int i=previousDistance+1;i<distance;i++) tempDensity[i]=UNK_DENSITY;
				        tempDensity[distance]=StrToFloat(str[1])-contrast;
                                        previousDistance=distance;
				        //cout <<dist<<" "<<str<<endl;
				        //cout <<"DensityMatrix["<<n<<"]["<<distance<<"]= "<<DensityMatrix[n][distance]<<endl;
                                }
                        }
                        prdf.density[n]=tempDensity;
                        tempDensity.clear();
		}
		getline(gofr, line);
	}
	for (n=0;n<NumSolventTypes;n++)
	{
		for (distance=0;distance<MaxRecord;distance++)
		{
			//cout <<"DensityMatrix["<<n<<"]["<<distance<<"]= "<<DensityMatrix[n][distance]<<endl;
		}
	}
	gofr.close();	
        cout <<"Leaving ReadGofRFile"<<endl;	
}
*/

void ReadPRDFFile(string prdfFile, Real contrast, PrdfStruct &prdf)
{
        fstream file;
        vector<string> str;
        string line;
        int previousDistance=0, dist;
        int atomType, MaxRecord, NumSolventTypes;
        vector<Real> tempDensity;
        Real distance;

        OpenFile(prdfFile, file, prdfFile);

        getline(file, line);
        if (line.substr(0, 15)=="NumSolventTypes")
        {
                NumSolventTypes=StrToInt(line.substr(16));
        }
        getline(file, line);
        Tokenize2(line, " ", str);
        MaxRecord=StrToInt(str[1]);

        prdf.AngularDependence=getValueOfBool("AngularDependence=", line);
        cout <<"After AngularDependence"<<endl;
        cout <<"line= "<<line<<endl;
        prdf.AngularDependence2=getValueOfBool("AngularDependence2=", line);
        prdf.ConvexConcave=getValueOfBool("ConvexConcave=", line);
        prdf.ConvexConcave2=getValueOfBool("ConvexConcave2=", line);
        prdf.PhobicPhilic=getValueOfBool("PhobicPhilic=", line);
        prdf.SecondNearestNeighbor=getValueOfBool("SecondNearestNeighbor=", line);
        prdf.UseChargeRadii=getValueOfBool("UseChargeRadii=", line);
        prdf.AngleBin=getValueOfReal("AngleBin=", line);
        prdf.AngularCutOff=getValueOfReal("AngularCutOff=", line);
        prdf.ChargeRadiiScale=getValueOfReal("ChargeRadiiScale=", line);
        prdf.CubeSize=getValueOfReal("CubeSize=", line);
        prdf.HydrationRadiusScale=getValueOfReal("HydrationRadiusScale=", line);
        prdf.RecBin=getValueOfReal("RecBin=", line);
        cout <<"Got prdf parameters"<<endl;
        cout <<"prdf.RecBin= "<<prdf.RecBin<<endl;
        while (true)
        {
                if (file.eof()) break;
                Tokenize2(line, " ", str);
                if (str.size()>0)
                {
                        if (str[0]=="AtomType") atomType=StrToInt(str[1]);
                        else if (str.size()>1)
                        {
                                distance=StrToReal(str[0]);
                                dist=int(distance/prdf.RecBin+0.5);
                                tempDensity.resize(dist+1);
                                for (int i=previousDistance+1;i<dist;i++) 
                                {
                                        tempDensity[i]=UNK_DENSITY;
                                }
				tempDensity[dist]=StrToReal(str[1])-contrast;
                                previousDistance=dist;
                        }
                }
                else
                {
                        previousDistance=0;
                        SafePushBack(prdf.density, tempDensity, "prdf.density");
                        tempDensity.clear();
                }
                getline(file, line);
        }
        cout <<"Leaving ReadPRDFFile"<<endl;
}

void ReadGofRFile(string GofRFile, Array3D &DensityMatrix, Real RecBin)
{
	char CharLine[1000];
	char *Space;
	int NumSolventTypes, MaxRecord;
	int type, distance, Size;
	string line;
        vector<string> str;
        vector<Real> v;
	Real dist;
        fstream gofr;

	cout <<"In ReadGofRFile"<<endl;

	NumSolventTypes=2000;
	MaxRecord=1000;
	type=0;

        OpenFile(GofRFile, gofr, "3D gofr file"); 
        
	getline(gofr, line);
	while (true)
	{
		if (line.substr(0,3)=="END") break;
		if (gofr.eof()) break;
		if (line.substr(0,15)=="NumSolventTypes")
		{
			NumSolventTypes=StrToInt(line.substr(16,10));
			cout <<"NumSolventTypes= "<<NumSolventTypes<<endl;
		}

		if (line.substr(0,9)=="MaxRecord")
		{
			MaxRecord=StrToInt(line.substr(10,10));
			cout <<"NumSolventTypes= "<<NumSolventTypes<<endl;
			cout <<"MaxRecord= "<<MaxRecord<<endl;
			Safe3DAlloc(DensityMatrix, NumSolventTypes, MaxRecord, 0, "3D density matrix");
			cout <<"After CreateMatrix"<<endl;
		}
		//cout <<"After Initialized DensityMatrix"<<endl;
		if (line.substr(0,8)=="AtomType")
		{
			type=StrToInt(line.substr(9,15));
			//str=line.substr(0,10);
			distance=0;
			//cout <<line<<endl;
			//cout <<line.substr(9,15)<<endl;
			//cout <<"n= "<<n<<endl;
			while (true)
			{
				getline(gofr,line);
				if (line.substr(0,1)==" ") cout <<"Break found"<<endl;
				if (line.substr(0,8)=="AtomType")
				{
					type++;
					break;
				}
				strcpy(CharLine, line.c_str());
				Space=strchr(CharLine,' ');
				if (Space==NULL) break;
				if (line.substr(0,1)==" " || line.substr(0,3)=="END") break;
				line=" " + line + " ";
                                Tokenize2(line, " ", str);
				dist=StrToFloat(str[0]);
				distance=int(dist/RecBin+0.5);
                                Size=str.size();
                                v.clear();
				for (int n=1;n<Size-2;n++) SafePushBack(v, Real(StrToFloat(str[n])/StrToFloat(str[Size-2])), "v in ReadGofRFile");
                                //PrintVector(v);
                                //PrintVector(str);
                                //cout <<line<<endl;
				DensityMatrix[type][distance]=v;
				//cout <<line<<endl;
				//cout <<dist<<" "<<str<<endl;
				//cout <<"DensityMatrix["<<type<<"]["<<distance<<"]= "<<DensityMatrix[type][distance]<<endl;
			}
		}
		getline(gofr, line);
	}
	gofr.close();
        cout <<"Leaving ReadGofRFile"<<endl;	
}

void ReadGofRFile3D(string GofRFile, PrdfStruct &prdf)
{
	char CharLine[1000];
	char *Space;
	int NumSolventTypes, MaxRecord;
	int type, distance, Size;
	string line;
        vector<string> str;
        vector<Real> v;
	Real dist;
        fstream gofr;

	cout <<"In ReadGofRFile"<<endl;

	NumSolventTypes=2000;
	MaxRecord=1000;
	type=0;

        OpenFile(GofRFile, gofr, "3D gofr file"); 
        
	getline(gofr, line);
	while (true)
	{
		if (line.substr(0,3)=="END") break;
		if (gofr.eof()) break;
		if (line.substr(0,15)=="NumSolventTypes")
		{
			NumSolventTypes=StrToInt(line.substr(16,10));
			cout <<"NumSolventTypes= "<<NumSolventTypes<<endl;
		}

		if (line.substr(0,9)=="MaxRecord")
		{
			MaxRecord=StrToInt(line.substr(10,10));
			cout <<"NumSolventTypes= "<<NumSolventTypes<<endl;
			cout <<"MaxRecord= "<<MaxRecord<<endl;
			Safe3DAlloc(prdf.density3d, NumSolventTypes, MaxRecord, 0, "3D density matrix");
                        prdf.AngularDependence=getValueOfBool("AngularDependence", line);
                        prdf.AngularDependence2=getValueOfBool("AngularDependence2", line);
                        prdf.ConvexConcave=getValueOfBool("ConvexConcave", line);
                        prdf.ConvexConcave2=getValueOfBool("ConvexConcave2", line);
                        prdf.PhobicPhilic=getValueOfBool("PhobicPhilic", line);
                        prdf.SecondNearestNeighbor=getValueOfBool("SecondNearestNeighbor", line);
                        prdf.UseChargeRadii=getValueOfBool("UseChargeRadii", line);
                        prdf.AngleBin=getValueOfReal("AngleBin", line);
                        prdf.AngularCutOff=getValueOfReal("AngularCutOff", line);
                        prdf.ChargeRadiiScale=getValueOfReal("ChargeRadiiScale", line);
                        prdf.CubeSize=getValueOfReal("CubeSize", line);
                        prdf.HydrationRadiusScale=getValueOfReal("HydrationRadiusScale", line);
                        prdf.RecBin=getValueOfReal("RecBin", line);
			cout <<"After CreateMatrix"<<endl;
		}
		//cout <<"After Initialized DensityMatrix"<<endl;
		if (line.substr(0,8)=="AtomType")
		{
			type=StrToInt(line.substr(9,15));
			//str=line.substr(0,10);
			distance=0;
			//cout <<line<<endl;
			//cout <<line.substr(9,15)<<endl;
			//cout <<"n= "<<n<<endl;
			while (true)
			{
				getline(gofr,line);
				if (line.substr(0,1)==" ") cout <<"Break found"<<endl;
				if (line.substr(0,8)=="AtomType")
				{
					type++;
					break;
				}
				strcpy(CharLine, line.c_str());
				Space=strchr(CharLine,' ');
				if (Space==NULL) break;
				if (line.substr(0,1)==" " || line.substr(0,3)=="END") break;
				line=" " + line + " ";
                                Tokenize2(line, " ", str);
				dist=StrToFloat(str[0]);
				distance=int(dist/prdf.RecBin+0.5);
                                Size=str.size();
                                v.clear();
				for (int n=1;n<Size-2;n++) SafePushBack(v, Real(StrToFloat(str[n])/StrToFloat(str[Size-2])), "v in ReadGofRFile");
                                //PrintVector(v);
                                //PrintVector(str);
                                //cout <<line<<endl;
				prdf.density3d[type][distance]=v;
				//cout <<line<<endl;
				//cout <<dist<<" "<<str<<endl;
				//cout <<"DensityMatrix["<<type<<"]["<<distance<<"]= "<<DensityMatrix[type][distance]<<endl;
			}
		}
		getline(gofr, line);
	}
	gofr.close();
        cout <<"Leaving ReadGofRFile"<<endl;	
}

#endif
