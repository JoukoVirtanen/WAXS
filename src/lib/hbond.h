#ifndef _hbond_included_
#define _hbond_included_

# include "TypeDef.h"
# include "AtomUtils.h"

void CalcThetaPsi(AtomStruct &o1, AtomStruct &h1, AtomStruct &h2, AtomStruct &o2, Real &theta, Real &psi)
{
        
        theta=CalcAngle(o1, h1, o2)*360.0/(2.0*pi);
        //cout <<"theta= "<<theta<<endl;
        psi=CalcDihedralAngle(o1, h1, o2, h2)*360.0/(2.0*pi);
        //cout <<"psi= "<<psi<<endl;
        if (psi>180.0) psi-=180.0;
}

bool CheckHBondAngles(AtomStruct &o1, AtomStruct &h1, AtomStruct &h2, AtomStruct &o2)
{
        Real thetaCutoff=30.0;
        Real theta;
        
        //PrintAtomInfo(o1);
        //PrintAtomInfo(h1);
        //PrintAtomInfo(o2);
        theta=CalcAngle(o1, h1, o2)*360.0/(2.0*pi);
        //cout <<"theta= "<<theta<<endl;
        if (theta>thetaCutoff) return false;
        //psi=CalcDihedralAngle(o1, h1, o2, h2)*360.0/(2.0*pi);
        //cout <<"psi= "<<psi<<endl;
        //if (psi>180.0) psi-=180.0;
        //if (psi>psiMin || psi<psiMax) return false;
        return true;
}

bool WaterWaterHBond(AtomStruct &h1a, AtomStruct &o1, AtomStruct &h1b, AtomStruct &h2a, AtomStruct &o2, AtomStruct &h2b)
{
        //J. Mol. Biol. (2003) 326, 1239-1259
        Real distCutoff=3.5;
        Real dist, distCutoff2;

        distCutoff2=distCutoff*distCutoff;
        //Real chiMin=
        //Real chiMax=         //Not currently using chi
        //cout <<endl;
        //PrintAtomInfo(h1a);
        //PrintAtomInfo(o1);
        //PrintAtomInfo(h1b);
        //PrintAtomInfo(h2a);
        //PrintAtomInfo(o2);
        //PrintAtomInfo(h2b);
        //cout <<endl;

        dist=AtomSquareDistance(o1, h2a);
        //cout <<"dist= "<<dist<<endl;
        if (dist<distCutoff2) return CheckHBondAngles(o1, h2a, h1a, o2);
        dist=AtomSquareDistance(o1, h2b);
        //cout <<"dist= "<<dist<<endl;
        //cout <<endl;
        if (dist<distCutoff2) return CheckHBondAngles(o1, h2b, h1a, o2);

        return false;
}

void WatersToAtoms(vector<AtomStruct> &WaterMolecule1, vector<AtomStruct> &WaterMolecule2, AtomStruct &h1a, AtomStruct &o1, AtomStruct &h1b, AtomStruct &h2a, AtomStruct &o2, AtomStruct &h2b)
{
        bool BoolPrint=false;
        int natom1=WaterMolecule1.size();
        int natom2=WaterMolecule2.size();
        if (natom1!=3 || natom2!=3)
        {
                BoolPrint=true;
                cout <<"WaterMolecule1.size()= "<<natom1<<endl;
                cout <<"WaterMolecule2.size()= "<<natom2<<endl;
        }
        for (int i=0;i<3;i++)
        {
                if (BoolPrint)
                {
                        PrintAtomInfo(WaterMolecule1[i]);
                        PrintAtomInfo(WaterMolecule2[i]);
                }
                if (WaterMolecule1[i].AtomName=="H1") h1a=WaterMolecule1[i];
                if (WaterMolecule1[i].AtomName=="OH2") o1=WaterMolecule1[i];
                if (WaterMolecule1[i].AtomName=="H2") h1b=WaterMolecule1[i];
                if (WaterMolecule2[i].AtomName=="H1") h2a=WaterMolecule2[i];
                if (WaterMolecule2[i].AtomName=="OH2") o2=WaterMolecule2[i];
                if (WaterMolecule2[i].AtomName=="H2") h2b=WaterMolecule2[i];
        }
}

void WatersToAtoms(vector<AtomStruct> &WaterMolecule, AtomStruct &ha, AtomStruct &o, AtomStruct &hb)
{
        for (int i=0;i<3;i++)
        {
                if (WaterMolecule[i].AtomName=="H1") ha=WaterMolecule[i];
                if (WaterMolecule[i].AtomName=="OH2") o=WaterMolecule[i];
                if (WaterMolecule[i].AtomName=="H2") hb=WaterMolecule[i];
        }
}

bool WaterWaterHBond(vector<AtomStruct> &WaterMolecule1, vector<AtomStruct> &WaterMolecule2)
{
        AtomStruct h1a, o1, h1b;
        AtomStruct h2a, o2, h2b;

        WatersToAtoms(WaterMolecule1, WaterMolecule2, h1a, o1, h1b, h2a, o2, h2b);
        //cout <<endl;
        //PrintAtomInfo(h1a);
        //PrintAtomInfo(o1);
        //PrintAtomInfo(h1b);
        //PrintAtomInfo(h2a);
        //PrintAtomInfo(o2);
        //PrintAtomInfo(h2b);
        //cout <<endl;
        return WaterWaterHBond(h1a, o1, h1b, h2a, o2, h2b);
}

int CalcHBonds(vector<AtomStruct> &Atoms, int index)
{
        int natom=Atoms.size();
        int numHBonds=0;
        vector<AtomStruct> WaterMolecule1, WaterMolecule2;

        GetWaterMolecule(Atoms, Atoms[index].ResidueNum, index, WaterMolecule1);

        for (int i=0;i<natom;i++)
        {
                if (i!=index && isWaterOxygen(Atoms[i]))
                {
                        GetWaterMolecule(Atoms, Atoms[i].ResidueNum, i, WaterMolecule2);
                        if (WaterWaterHBond(WaterMolecule1, WaterMolecule2)) numHBonds++;
                        if (WaterWaterHBond(WaterMolecule2, WaterMolecule1)) numHBonds++;
                }
        }
        return numHBonds;
}

#endif
