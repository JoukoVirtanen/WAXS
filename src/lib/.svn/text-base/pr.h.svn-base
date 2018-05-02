#ifndef _pr_included_
#define _pr_included_

# include "Structures.h"
# include "TypeDef.h"

Real FindGreatestDistance(vector<AtomStruct> &Atoms)
{
	Real dx, dy, dz, dist;
	Real xmax, ymax, zmax;
	Real xmin, ymin, zmin;
	MinMax(xmin, ymin, zmin, xmax, ymax, zmax, Atoms);
	dx=xmax-xmin;
	dy=ymax-ymin;
	dz=zmax-zmin;
	dist=sqrt(dx*dx+dy*dy+dz*dz);
	return dist;
}

void InitializePr(PrStruct &pr, int GreatestDistance, Real bin)
{
	Real dist=0;
	CreateVector(pr.r, GreatestDistance);
	CreateVector(pr.all, GreatestDistance);
	for (int i=0;i<GreatestDistance;i++)
	{
		pr.r[i]=dist;
		dist+=bin;	
	}
}

void PrintPr(PrStruct pr, string PrOutFile)
{
	char CharPrOutFile[1000];
	
	strcpy(CharPrOutFile, PrOutFile.c_str());
	ofstream out;
	out.open(CharPrOutFile, ios::app);
	for (int i=0;i<pr.r.size();i++)
	{
		out <<pr.r[i]<<"\t"<<pr.all[i]<<endl;
	}
	out.close();
}

void NormalizePr(PrStruct &pr, Real bin)
{
	Real total=0;
	for (int i=0;i<pr.r.size();i++) total+=pr.all[i];
	total*=bin;
	for (int i=0;i<pr.r.size();i++) pr.all[i]/=total;
}

void FindPr(PrStruct &pr, vector<AtomStruct> Atoms, ElementStruct ElementData[], Real bin)
{
	int GreatestDistance, IntDist, natom=Atoms.size();
	Real dist, InvBin, mass1, mass2, MassProduct;
	InvBin=1.0/bin;

	dist=FindGreatestDistance(Atoms);
	GreatestDistance=int(dist*InvBin+0.5);
	InitializePr(pr, GreatestDistance, bin);
	cout <<"In FindPr. natom= "<<natom;
	cout <<"GreatestDistance= "<<GreatestDistance<<endl;
	for (int i=0;i<natom;i++)
	{
		mass1=Atoms[i].weight*ElementData[Atoms[i].atomid].SolventCorrectedElectrons;
		for (int j=i+1;j<natom;j++)
		{
			dist=AtomDistance(Atoms[i], Atoms[j]);
			IntDist=int(dist*InvBin+0.5);
			mass2=Atoms[j].weight*ElementData[Atoms[j].atomid].SolventCorrectedElectrons;
			pr.all[IntDist]+=mass1*mass2;
			//cout <<"pr.all["<<IntDist<<"]= "<<pr.all[IntDist]<<endl;
		}
	}
}
#endif
