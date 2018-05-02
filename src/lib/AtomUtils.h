#ifndef _AtomUtils_included_
#define _AtomUtils_included_

# include "AtomCode.h"
# include "AtomIDs.h"
# include "Constants.h"
# include "GetAtomType.h"
# include "GetAtomType2.h"
# include "MathUtils.h"
# include "ReadPdb.h"
# include "ReadPrm.h"
# include "ReadPsf.h"
# include "ReadRtf.h"
# include "ResidueCode.h"
# include "Structures.h"
# include "VectorManip.h"

bool IsProtein(AtomStruct &Atom);
void RemoveUnknownResidues(vector<AtomStruct> &Atoms);
void RemoveWaters(vector<AtomStruct> &Atoms);
void GetFirstLastResidue(vector<AtomStruct> &Atoms, int &FirstResidue, int &LastResidue);
void MoveAtoms(AtomStruct &Atoms, Real dx, Real dy, Real dz, bool boolPrint=true);
bool isNA(AtomStruct &Atom);
# include "MinMax.h"


void CopyAtom(const AtomStruct &Atom, AtomStruct &Atom2)
{
        Atom2.HetAtom=Atom.HetAtom;
        Atom2.solute=Atom.solute;
        Atom2.AtomNumber=Atom.AtomNumber;
        Atom2.ResidueNum=Atom.ResidueNum;
        Atom2.atomid=Atom.atomid;
        Atom2.AtomType=Atom.AtomType;
        Atom2.ParmType=Atom.ParmType;
        Atom2.residueid=Atom.residueid;
        Atom2.weight=Atom.weight;
        Atom2.Electrons=Atom.Electrons;
        Atom2.charge=Atom.charge;
        Atom2.mass=Atom.mass;
        Atom2.x=Atom.x;
        Atom2.y=Atom.y;
        Atom2.z=Atom.z;
        Atom2.vdw=Atom.vdw;
        Atom2.epsilon=Atom.epsilon;
        Atom2.vdw14=Atom.vdw14;
        Atom2.epsilon14=Atom.epsilon14;
        Atom2.Occupancy=Atom.Occupancy;
        Atom2.BFactor=Atom.BFactor;
        Atom2.AtomName=Atom.AtomName;
        Atom2.CharmmAtomName=Atom.CharmmAtomName;
        Atom2.ResidueName=Atom.ResidueName;
        Atom2.ChainName=Atom.ChainName;
        Atom2.SegID=Atom.SegID;
        Atom2.ID=Atom.ID;
}

void InitializeAtom(AtomStruct &Atom)
{
        Atom.HetAtom=false;
        Atom.solute=true;
        Atom.AtomNumber=0;
        Atom.ResidueNum=0;
        Atom.atomid=0;
        Atom.AtomType=0;
        Atom.AtomType2=0;
        Atom.ParmType=0;
        Atom.residueid=0;
        //Atom.atomtype=0;
        Atom.weight=1.0;
        Atom.Electrons=0;
        Atom.charge=0;
        Atom.mass=0;
        Atom.x=0;
        Atom.y=0;
        Atom.z=0;
        Atom.vdw=0;
        Atom.epsilon=0;
        Atom.vdw14=0;
        Atom.epsilon14=0;
        Atom.Occupancy=0;
        Atom.BFactor=0;
        Atom.AtomName="";
        Atom.CharmmAtomName="";
        Atom.ResidueName="";
        Atom.ChainName="A";
        Atom.SegID="";
        Atom.ID="";
}

void MoveAtoms(vector<AtomStruct> &Atoms, Real dx, Real dy, Real dz, bool boolPrint=true)
{
        int i, natom=Atoms.size();
        if (boolPrint) cout <<"MoveX= "<<dx<<" MoveY= "<<dy<<" MoveZ= "<<dz<<endl;
        for (i=0;i<natom;i++)
        {
                Atoms[i].x+=dx;
                Atoms[i].y+=dy;
                Atoms[i].z+=dz;
        }
}

void MoveAtom(AtomStruct &Atom, Real dx, Real dy, Real dz, bool boolPrint=true)
{
        if (boolPrint) cout <<"MoveX= "<<dx<<" MoveY= "<<dy<<" MoveZ= "<<dz<<endl;
        Atom.x+=dx;
        Atom.y+=dy;
        Atom.z+=dz;
}

void MergeProteins(const vector<ProteinStruct> &Proteins, vector<AtomStruct> &Atoms)
{
        int i, nprotein=Proteins.size();
        DeleteVector(Atoms);

        for (i=0;i<nprotein;i++) AppendVector(Proteins[i].Atoms, Atoms);
}

void CalcDistMatrix(vector< vector<Real> > &DistMatrix, vector<AtomStruct> &Atoms)
{
        int i, j, natom=Atoms.size();
        Real dist;
        Real dx, dy, dz;
        for (i=0;i<natom;i++)
        {
                for (j=i+1;j<natom;j++)
                {
                        dx=Atoms[i].x-Atoms[j].x;
                        dy=Atoms[i].y-Atoms[j].y;
                        dz=Atoms[i].z-Atoms[j].z;
                        dist=sqrt(dx*dx+dy*dy+dz*dz);
                        DistMatrix[i][j]=dist;
                        DistMatrix[j][i]=dist;
                }
        }
}

Real AtomSquareDistance(AtomStruct &Atom1, AtomStruct &Atom2)
{
        Real dx, dy, dz;
        dx=Atom1.x-Atom2.x;
        dy=Atom1.y-Atom2.y;
        dz=Atom1.z-Atom2.z;
        return dx*dx+dy*dy+dz*dz;
}

Real AtomDistance(AtomStruct &Atom1, AtomStruct &Atom2)
{
        Real dx, dy, dz;
        dx=Atom1.x-Atom2.x;
        dy=Atom1.y-Atom2.y;
        dz=Atom1.z-Atom2.z;
        return sqrt(dx*dx+dy*dy+dz*dz);
}

Real AtomDistance(Real x, Real y, Real z, AtomStruct &Atom)
{
        Real dx, dy, dz;
        dx=x-Atom.x;
        dy=y-Atom.y;
        dz=z-Atom.z;
        return sqrt(dx*dx+dy*dy+dz*dz);
}

Real calcDistanceMatrix(vector<AtomStruct> &Atoms, Matrix &DistanceMatrix)
{
        int natom=Atoms.size();
        int Size1, Size2;
        Real dist;

        Get2DVectorSize(DistanceMatrix, Size1, Size2);
        if (Size1!=natom || Size2!=natom)
        {
                Safe2DAlloc(DistanceMatrix, natom, natom, "DistanceMatrix");
        }
        for (int i=0;i<natom;i++)
        {
                for (int j=i+1;j<natom;j++)
                {
                        dist=AtomDistance(Atoms[i], Atoms[j]);
                        DistanceMatrix[i][j]=dist;
                        DistanceMatrix[j][i]=dist;
                }
        }
}

void updateDistanceMatrix(vector<AtomStruct> &Atoms, Matrix &DistanceMatrix, int start1, int end1, int start2, int end2)
{
        Real dist;
        for (int i=start1;i<=end1;i++)
        {
                for (int j=start2;j<=end2;j++)
                {
                        dist=AtomDistance(Atoms[i], Atoms[j]);
                        DistanceMatrix[i][j]=dist;
                        DistanceMatrix[j][i]=dist;
                }
        }
}

void updateDistanceMatrix(ProteinStruct &Protein, vector<int> &start, vector<int> &end)
{
        int nstart=start.size();
        int nend=end.size();

        if (nstart!=nend)
        {
                cout <<"ERROR: nstart!=nend"<<endl;
                exit(EXIT_FAILURE);
        }
        for (int i=0;i<nstart;i++)
        {
                for (int j=i+1;j<nstart;j++)
                {
                        updateDistanceMatrix(Protein.Atoms, Protein.DistMatrix, start[i], end[i], start[j], end[j]);
                }
        }
}

void updateAtomDistances(vector<AtomStruct> &Atoms, Matrix &DistMatrix, int atomNum)
{
        int natom=Atoms.size();
        Real dist;
        for (int i=0;i<natom;i++)
        {
                if (i!=atomNum)
                {
                        dist=AtomDistance(Atoms[i], Atoms[atomNum]);
                        DistMatrix[i][atomNum]=dist;
                        DistMatrix[atomNum][i]=dist;
                }
        }
}

void updateLocalDistanceMatrix(ProteinStruct &Protein, int residue1, int residue2)
{
        int natom=Protein.Atoms.size();

        for (int i=0;i<natom;i++)
        {
                if (Protein.Atoms[i].ResidueNum>=residue1 && Protein.Atoms[i].ResidueNum<=residue2)
                {
                        updateAtomDistances(Protein.Atoms, Protein.DistMatrix, i);
                }
        }
}

void RotateAroundX(AtomStruct &Atom, Real theta)
{
        Real tempY, tempZ;
        tempY=Atom.y*cos(theta)-Atom.z*sin(theta);
        tempZ=Atom.y*sin(theta)+Atom.z*cos(theta);
        Atom.y=tempY;
        Atom.z=tempZ;
}

void RotateAroundY(AtomStruct &Atom, Real theta)
{
        Real tempX, tempZ;
        tempX=Atom.x*cos(theta)-Atom.z*sin(theta);
        tempZ=Atom.x*sin(theta)+Atom.z*cos(theta);
        Atom.x=tempX;
        Atom.z=tempZ;
}

void RotateAroundZ(AtomStruct &Atom, Real theta)
{
        Real tempX, tempY;
        tempX=Atom.x*cos(theta)-Atom.y*sin(theta);
        tempY=Atom.x*sin(theta)+Atom.y*cos(theta);
        Atom.x=tempX;
        Atom.y=tempY;
}

void CalcRotationMatrix(VectorStruct &v, Real angle, Real RotationMatrix[3][3])
{
        RotationMatrix[0][0]=cos(angle)+v.x*v.x*(1.0-cos(angle));
        RotationMatrix[0][1]=v.x*v.y*(1.0-cos(angle))-v.z*sin(angle);
        RotationMatrix[0][2]=v.x*v.z*(1.0-cos(angle))+v.y*sin(angle);

        RotationMatrix[1][0]=v.x*v.y*(1.0-cos(angle))+v.z*sin(angle);
        RotationMatrix[1][1]=cos(angle)+v.y*v.y*(1.0-cos(angle));
        RotationMatrix[1][2]=v.y*v.z*(1.0-cos(angle))-v.x*sin(angle);

        RotationMatrix[2][0]=v.x*v.z*(1.0-cos(angle))-v.y*sin(angle);
        RotationMatrix[2][1]=v.y*v.z*(1.0-cos(angle))+v.x*sin(angle);
        RotationMatrix[2][2]=cos(angle)+v.z*v.z*(1.0-cos(angle));
}

void ApplyRotationMatrix(Real &x, Real &y, Real &z, Real RotationMatrix[3][3])
{
        Real tempX, tempY, tempZ;

        tempX=RotationMatrix[0][0]*x+RotationMatrix[0][1]*y+RotationMatrix[0][2]*z;
        tempY=RotationMatrix[1][0]*x+RotationMatrix[1][1]*y+RotationMatrix[1][2]*z;
        tempZ=RotationMatrix[2][0]*x+RotationMatrix[2][1]*y+RotationMatrix[2][2]*z;

        x=tempX;
        y=tempY;
        z=tempZ;
}

void ApplyRotationMatrixToVector(VectorStruct &v, Real RotationMatrix[3][3])
{
        ApplyRotationMatrix(v.x, v.y, v.z, RotationMatrix);
}

VectorStruct AtomVector(AtomStruct &Atom1, AtomStruct &Atom2)
{
        VectorStruct v1;
        //PrintAtomInfo(Atom1);
        //PrintAtomInfo(Atom2);
        v1.x=Atom1.x-Atom2.x;
        v1.y=Atom1.y-Atom2.y;
        v1.z=Atom1.z-Atom2.z;

        return v1;
}

void ApplyRotationMatrixToAtoms(vector<AtomStruct> &Atoms, int Atom, Real RotationMatrix[3][3])
{
        int natom=Atoms.size();
        VectorStruct v;

        for (int i=Atom+1;i<natom;i++)
        {
                v=AtomVector(Atoms[i], Atoms[Atom]);
                ApplyRotationMatrixToVector(v, RotationMatrix);
                Atoms[i].x=Atoms[Atom].x+v.x;
                Atoms[i].y=Atoms[Atom].y+v.y;
                Atoms[i].z=Atoms[Atom].z+v.z;
        }
}

void ApplyRotationMatrixToAtoms(vector<AtomStruct> &Atoms, Real x, Real y, Real z, Real RotationMatrix[3][3])
{
        int natom=Atoms.size();
        VectorStruct v;

        for (int i=0;i<natom;i++)
        {
                v.x=Atoms[i].x-x;
                v.y=Atoms[i].y-y;
                v.z=Atoms[i].z-z;
                ApplyRotationMatrixToVector(v, RotationMatrix);
                Atoms[i].x=x+v.x;
                Atoms[i].y=y+v.y;
                Atoms[i].z=z+v.z;
        }
}

void RotateAtomsAroundVector(vector<AtomStruct> &Atoms, VectorStruct &v, Real x, Real y, Real z, Real angle)
{
        Real RotationMatrix[3][3];

        CalcRotationMatrix(v, angle, RotationMatrix);
        ApplyRotationMatrixToAtoms(Atoms, x, y, z, RotationMatrix);
}

void RotateAtomsAroundVector(vector<AtomStruct> &Atoms, VectorStruct &v, int Atom, Real angle)
{
        Real RotationMatrix[3][3];

        CalcRotationMatrix(v, angle, RotationMatrix);
        ApplyRotationMatrixToAtoms(Atoms, Atom, RotationMatrix);
}

VectorStruct RotateAroundVector(VectorStruct v, Real angle, VectorStruct &rotationV)
{
        Real RotationMatrix[3][3];
        VectorStruct vout;

        vout.x=v.x;
        vout.y=v.y;
        vout.z=v.z;
        cout <<"vout= "<<VectorToStr(vout)<<endl;
        CalcRotationMatrix(rotationV, angle, RotationMatrix);
        ApplyRotationMatrixToVector(vout, RotationMatrix);
        return vout;
}

Real CalcAngle(AtomStruct &Atom1, AtomStruct &Atom2, AtomStruct &Atom3)
{
        VectorStruct v1, v2;
        Real DotProduct;

        v1=AtomVector(Atom2, Atom1);
        v2=AtomVector(Atom3, Atom2);
        //cout <<"v1.x= "<<v1.x<<" v1.y= "<<v1.y<<" v1.z= "<<v1.z<<endl;
        //cout <<"v2.x= "<<v2.x<<" v2.y= "<<v2.y<<" v2.z= "<<v2.z<<endl;
        DotProduct=v1*v2;
        //cout <<"DotProduct= "<<DotProduct<<endl;
        return acos(DotProduct/(VectorLength(v1)*VectorLength(v2)));
}

VectorStruct CalcCrossProduct(AtomStruct &Atom1, AtomStruct &Atom2, AtomStruct &Atom3)
{
        VectorStruct v1, v2, product;

        v1.x=Atom2.x-Atom1.x;
        v1.y=Atom2.y-Atom1.y;
        v1.z=Atom2.z-Atom1.z;

        v2.x=Atom3.x-Atom2.x;
        v2.y=Atom3.y-Atom2.y;
        v2.z=Atom3.z-Atom2.z;

        //PrintAtomInfo(Atom1);
        //PrintAtomInfo(Atom2);
        //PrintAtomInfo(Atom3);

        product=CalcCrossProduct(v1, v2);

        return product;
}

Real CalcSign(VectorStruct &v1, VectorStruct &plane1, VectorStruct &plane2)
{
        VectorStruct CrossProduct;
        Real DotProduct;
        CrossProduct=CalcCrossProduct(plane1, plane2);
        DotProduct=v1*CrossProduct;
        if (DotProduct>0) return 1.0;
        else return -1.0;
}

Real CalcDihedralAngle(VectorStruct &v1, VectorStruct &v2, VectorStruct &v3)
{
        VectorStruct plane1, plane2;
        Real DotProduct, costheta, sign;
        plane1=CalcCrossProduct(v2, v1);
        plane2=CalcCrossProduct(v3, v2);
        DotProduct=plane1*plane2;
        sign=CalcSign(v2, plane1, plane2);

        costheta=DotProduct/(VectorLength(plane1)*VectorLength(plane2));
        if (sign>0) return acos(costheta);
        if (sign<=0) return -acos(costheta);
        cout <<"ERROR: Should not be here at the end of CalcDihedralAngle"<<endl;
        return 0;
}

Real CalcDihedralAngle(AtomStruct &Atom1, AtomStruct &Atom2, AtomStruct &Atom3, AtomStruct &Atom4)
{
        VectorStruct v1, v2, v3;

        v1=AtomVector(Atom2, Atom1);
        v2=AtomVector(Atom3, Atom2);
        v3=AtomVector(Atom4, Atom3);

        return CalcDihedralAngle(v1, v2, v3);
}

Real CalcDihedralAngle(int Atom1, int Atom2, int Atom3, int Atom4, vector<AtomStruct> &Atoms)
{
        return CalcDihedralAngle(Atoms[Atom1], Atoms[Atom2], Atoms[Atom3], Atoms[Atom4]);
}

Real CalcPhi(vector<AtomStruct> &Atoms, int ResidueNum)
{
        int natom=Atoms.size();
        int c, n, ca, c2;
        if (ResidueNum<2)
        {
                cout <<"ERROR: No phi angle for residue"<<endl;
                exit(EXIT_FAILURE);
        }
        for (int i=0;i<natom;i++)
        {
                if (Atoms[i].ResidueNum==ResidueNum-1 && Atoms[i].AtomName=="C") c=i;
                else if (Atoms[i].ResidueNum==ResidueNum && Atoms[i].AtomName=="N") n=i;
                else if (Atoms[i].ResidueNum==ResidueNum && Atoms[i].AtomName=="CA") ca=i;
                else if (Atoms[i].ResidueNum==ResidueNum && Atoms[i].AtomName=="C") c2=i;
        }
        return CalcDihedralAngle(c, n, ca, c2, Atoms);
}

Real CalcPsi(vector<AtomStruct> &Atoms, int ResidueNum)
{
        int natom=Atoms.size();
        int n, ca, c, n2;
        for (int i=0;i<natom;i++)
        {
                if (Atoms[i].ResidueNum==ResidueNum && Atoms[i].AtomName=="N") n=i;
                else if (Atoms[i].ResidueNum==ResidueNum && Atoms[i].AtomName=="CA") ca=i;
                else if (Atoms[i].ResidueNum==ResidueNum && Atoms[i].AtomName=="C") c=i;
                else if (Atoms[i].ResidueNum==ResidueNum+1 && Atoms[i].AtomName=="N") n2=i;
        }
        return CalcDihedralAngle(n, ca, c, n2, Atoms);
}

Real CalcOmega(vector<AtomStruct> &Atoms, int ResidueNum)
{
        int natom=Atoms.size();
        int ca, c, n, ca2;
        for (int i=0;i<natom;i++)
        {
                if (Atoms[i].ResidueNum==ResidueNum && Atoms[i].AtomName=="CA") ca=i;
                else if (Atoms[i].ResidueNum==ResidueNum && Atoms[i].AtomName=="C") c=i;
                else if (Atoms[i].ResidueNum==ResidueNum+1 && Atoms[i].AtomName=="N") n=i;
                else if (Atoms[i].ResidueNum==ResidueNum+1 && Atoms[i].AtomName=="CA") ca2=i;
        }
        return CalcDihedralAngle(ca, c, n, ca2, Atoms);
}

Real CalcDihedralAngle(vector<AtomStruct> &Atoms, int AngleType, int ResidueNum)
{
        if (AngleType==PHI) return CalcPhi(Atoms, ResidueNum);
        else if (AngleType==PSI) return CalcPsi(Atoms, ResidueNum);
        else if (AngleType==OMEGA) return CalcOmega(Atoms, ResidueNum);
        else
        {
                cout <<"ERROR: Unknown AngleType "<<AngleType<<endl;
                cout <<"Acceptable values are "<<PHI<<", "<<PSI<<", and "<<OMEGA<<endl;
                exit(EXIT_FAILURE);
        }
}

void RotateAtomsAroundVector(vector<AtomStruct> &Atoms, int Atom1, int Atom2, Real angle)
{
        VectorStruct v;

        v.x=Atoms[Atom2].x-Atoms[Atom1].x;
        v.y=Atoms[Atom2].y-Atoms[Atom1].y;
        v.z=Atoms[Atom2].z-Atoms[Atom1].z;
        VectorNormalize(v);

        RotateAtomsAroundVector(Atoms, v, Atom2, angle);
}

void SetDihedralAngle(vector<AtomStruct> &Atoms, int Atom1, int Atom2, int Atom3, int Atom4, Real DihedralAngle)
{
        Real OldDihedralAngle=CalcDihedralAngle(Atom1, Atom2, Atom3, Atom4, Atoms);
        Real deltaAngle=DihedralAngle-OldDihedralAngle;

        RotateAtomsAroundVector(Atoms, Atom2, Atom3, deltaAngle);
}

void rotateDihedralAngleBy(vector<AtomStruct> &Atoms, int atom1, int atom2, Real deltaAngle)
{
        RotateAtomsAroundVector(Atoms, atom1, atom2, deltaAngle);
}

void GetPhiDihedralAtoms(vector<AtomStruct> &Atoms, int &Atom1, int &Atom2, int &Atom3, int &Atom4, int ResidueNum)
{
        int natom=Atoms.size();
        if (ResidueNum<2)
        {
                cout <<"ERROR: No phi angle for residue"<<endl;
                exit(EXIT_FAILURE);
        }

        for (int i=0;i<natom;i++)
        {
                if (Atoms[i].ResidueNum==ResidueNum-1 && Atoms[i].AtomName=="C") Atom1=i;
                else if (Atoms[i].ResidueNum==ResidueNum && Atoms[i].AtomName=="N") Atom2=i;
                else if (Atoms[i].ResidueNum==ResidueNum && Atoms[i].AtomName=="CA") Atom3=i;
                else if (Atoms[i].ResidueNum==ResidueNum && Atoms[i].AtomName=="C") Atom4=i;
        }
}

void GetPsiDihedralAtoms(vector<AtomStruct> &Atoms, int &Atom1, int &Atom2, int &Atom3, int &Atom4, int ResidueNum)
{
        int natom=Atoms.size();

        for (int i=0;i<natom;i++)
        {
                if (Atoms[i].ResidueNum==ResidueNum && Atoms[i].AtomName=="N") Atom1=i;
                else if (Atoms[i].ResidueNum==ResidueNum && Atoms[i].AtomName=="CA") Atom2=i;
                else if (Atoms[i].ResidueNum==ResidueNum && Atoms[i].AtomName=="C") Atom3=i;
                else if (Atoms[i].ResidueNum==ResidueNum+1 && Atoms[i].AtomName=="N") Atom4=i;
        }
}

void GetOmegaDihedralAtoms(vector<AtomStruct> &Atoms, int &Atom1, int &Atom2, int &Atom3, int &Atom4, int ResidueNum)
{
        int natom=Atoms.size();

        for (int i=0;i<natom;i++)
        {
                if (Atoms[i].ResidueNum==ResidueNum && Atoms[i].AtomName=="CA") Atom1=i;
                else if (Atoms[i].ResidueNum==ResidueNum && Atoms[i].AtomName=="C") Atom2=i;
                else if (Atoms[i].ResidueNum==ResidueNum+1 && Atoms[i].AtomName=="N") Atom3=i;
                else if (Atoms[i].ResidueNum==ResidueNum+1 && Atoms[i].AtomName=="CA") Atom4=i;
        }
}

void GetDihedralAtoms(vector<AtomStruct> &Atoms, int &Atom1, int &Atom2, int &Atom3, int &Atom4, int AngleType, int ResidueNum)
{
        if (AngleType==PHI) GetPhiDihedralAtoms(Atoms, Atom1, Atom2, Atom3, Atom4, ResidueNum);
        else if (AngleType==PSI) GetPsiDihedralAtoms(Atoms, Atom1, Atom2, Atom3, Atom4, ResidueNum);
        else if (AngleType==OMEGA) GetOmegaDihedralAtoms(Atoms, Atom1, Atom2, Atom3, Atom4, ResidueNum);
        else
        {
                cout <<"ERROR: Unknown dihedral angle type"<<endl;
                exit(EXIT_FAILURE);
        }
}

void SetDihedralAngle(vector<AtomStruct> &Atoms, int AngleType, int ResidueNum, Real DihedralAngle)
{
        int Atom1, Atom2, Atom3, Atom4;
        GetDihedralAtoms(Atoms, Atom1, Atom2, Atom3, Atom4, AngleType, ResidueNum);
        SetDihedralAngle(Atoms, Atom1, Atom2, Atom3, Atom4, DihedralAngle);
}

void DihedralIndexToResidueAndAngleType(int dihedralIndex, int &residue, int &angleType)
{
        residue=int((dihedralIndex+1)/3)+1;
        if (dihedralIndex%3==0) angleType=PSI;
        if (dihedralIndex%3==1) angleType=OMEGA;
        if (dihedralIndex%3==2) angleType=PHI;
}

void SetDihedralAngle(vector<AtomStruct> &Atoms, int dihedralIndex, Real DihedralAngle)
{
        int angleType, residue;
        DihedralIndexToResidueAndAngleType(dihedralIndex, residue, angleType);
        SetDihedralAngle(Atoms, angleType, residue, DihedralAngle);
}

void rotateDihedralAngleBy(vector<AtomStruct> &Atoms, int dihedralIndex, Real deltaAngle)
{
        int atom1, atom2, atom3, atom4;
        int angleType, residue;
        DihedralIndexToResidueAndAngleType(dihedralIndex, residue, angleType);
        GetDihedralAtoms(Atoms, atom1, atom2, atom3, atom4, angleType, residue);
        rotateDihedralAngleBy(Atoms, atom2, atom3, deltaAngle); 
}

Real CalcDihedralAngle(vector<AtomStruct> &Atoms, int dihedralIndex)
{
        int residue, angleType;
        DihedralIndexToResidueAndAngleType(dihedralIndex, residue, angleType);
        return CalcDihedralAngle(Atoms, angleType, residue); 
}

void SetDihedralAngles(vector<AtomStruct> &Atoms, vector<Real> &dihedrals)
{
        int FirstResidue, LastResidue;
        int numDihedral=0;
        GetFirstLastResidue(Atoms, FirstResidue, LastResidue);
        for (int i=FirstResidue;i<=LastResidue;i++)
        {
                if (i!=FirstResidue)
                {
                        SetDihedralAngle(Atoms, PHI, i, dihedrals[numDihedral]);
                        numDihedral++;
                }
                if (i!=LastResidue)
                {
                        SetDihedralAngle(Atoms, PSI, i, dihedrals[numDihedral]);
                        numDihedral++;
                        SetDihedralAngle(Atoms, OMEGA, i, dihedrals[numDihedral]);
                        numDihedral++;
                }
        }
}
/*
   void CreateRotMatrix(vector<Real> Axis, Real angle, vector< vector<Real> > &RotMat)
   {
   Real cosine=cos(angle);
   Real sine=sin(angle);
   Real one_minus_cosine=1-cosine;

   RotMat[X][X]=Sqr(Axis[X])+(1-Sqr(Axis[x]))*cosine;
   RotMat[X][Y]=Axis[X]*Axis[Y]*one_minus_cosine+Axis[Z]*sine;
   RotMat[X][Z]=Axis[X]*Axis[Z]*one_minus_cosine-Axis[Y]*sine;
   }
   */
void RotateAtoms(vector<AtomStruct> &Atoms, Real theta, Real phi, Real psi)
{
        int i, natom=Atoms.size();
        for (i=0;i<natom;i++)
        {
                RotateAroundX(Atoms[i], theta);
                RotateAroundY(Atoms[i], phi);
                RotateAroundZ(Atoms[i], psi);
        }
} 

void RotateAtom(AtomStruct &Atom, Real theta, Real phi, Real psi)
{
        RotateAroundX(Atom, theta);
        RotateAroundY(Atom, phi);
        RotateAroundZ(Atom, psi);
} 
/*
   void allignBondWithZAxis(vector<AtomStruct> &Atoms, BondStruct &bond)
   {
   angleWithRespectoToXAxis=calcAngleWithRespectToXAxis(Atoms[bond.atom1], Atoms[bond.atom2]);
   RotateAroundZ(Atoms, angleWithRespectToXAxis);
   anngleWithRespectToZAxis=calcAngleWithRespectToZAxis(Atoms[bond.atom1], Atoms[bond.atom2]);
   RotateAroundX(Atoms, angleWithRespectToZAxis);
   }

   Real calcProjection(vector<AtomStruct> &Atoms, BondStruct &bond1, BondStruct &bond2)
   {
   VectorStruct v1, v2;

   v1=AtomVector(Atoms[bond1.atom1], Atoms[bond1.atom2]);
   v2=AtomVector(Atoms[bond2.atom1], Atoms[bond2.atom2]);
   return calcVectorAngle(v1, v2);
   }

   Real calcProjection(vector<AtomStruct> &Atoms, BondStruct &mainBond, BondStruct &bond1, BondStruct &bond2)
   {
   Real projection;
   allignBondWithZAxis(Atoms, mainBond);
   projection=calcProjection(Atoms, bond1, bond2);
   return projection;
   }

   void MakeProjectionsOfBondsPerpendicular(vector<AtomStruct> &Atoms, BondStruct &mainBond, BondStruct &bond1, BondStruct &bond2)
   {
   projection=calcProjection(Atoms, mainBond, bond1, bond2);

   }
   */

void doubleCrank(vector<AtomStruct> &Atoms, int Residue, Real dAngle)
{
        int FirstResidue, LastResidue;
        Real phi, psi;
        GetFirstLastResidue(Atoms, FirstResidue, LastResidue);
        if (Residue!=FirstResidue)
        {
                phi=CalcPhi(Atoms, Residue);
                SetDihedralAngle(Atoms, PHI, Residue, phi+dAngle);
                psi=CalcPsi(Atoms, Residue-1);
                SetDihedralAngle(Atoms, PSI, Residue-1, psi-dAngle);
        }
}

void pivot(vector<AtomStruct> &Atoms, int Residue, int AngleType, Real dAngle)
{
        Real angle;
        angle=CalcDihedralAngle(Atoms, AngleType, Residue);
        SetDihedralAngle(Atoms, AngleType, Residue, angle+dAngle);
}

int getAtomIndex(vector<AtomStruct> &Atoms, int residue, string AtomName)
{
        int natom=Atoms.size();

        for (int i=0;i<natom;i++)
        {
                if (Atoms[i].ResidueNum==residue && Atoms[i].AtomName==AtomName)
                {
                        return i;
                }
        }
        return 0;
}

int GetNumResidues(vector<AtomStruct> &Atoms)
{
        int natom=Atoms.size();
        return Atoms[natom-1].ResidueNum-Atoms[0].ResidueNum+1;
}

void GetResidues(vector<AtomStruct> &Atoms, vector<int> &ResidueTypes)
{
        int natom=Atoms.size();
        int ResidueType;
        for (int i=0;i<natom;i++)
        {
                if (Atoms[i].AtomName=="CA")
                {
                        ResidueType=GetResidueID(Atoms[i].ResidueName);                   
                        SafePushBack(ResidueTypes, ResidueType, "ResidueTypes");
                }
        }
}

void GetWaterMolecule(vector<AtomStruct> &Atoms, int ResidueNum, int index, vector<AtomStruct> &WaterMolecule)
{
        bool BoolPrint=false;
        int natom=Atoms.size();
        int min, max;
        //if (ResidueNum==9250) BoolPrint=true;
        WaterMolecule.clear();
        min=index-3;
        max=index+3;
        if (min<0) min=0;
        if (max>=natom) max=natom;
        if (BoolPrint) cout <<endl;
        //cout <<"In GetWaterMolecule"<<endl;
        for (int i=min;i<max;i++)
        {
                //cout <<"i= "<<i<<" ResidueNum= "<<ResidueNum<<endl;
                if (BoolPrint) PrintAtomInfo(Atoms[i]);
                if (Atoms[i].ResidueNum==ResidueNum)
                {
                        SafePushBack(WaterMolecule, Atoms[i], "WaterMolecule in GetWaterMolecule");
                        if (BoolPrint) cout <<"Added atom to water molecule"<<endl;
                }
        }
        if (BoolPrint) cout <<endl;
}

bool IsWater(AtomStruct &Atom)
{
        if (Atom.ResidueName=="HOH") return true;
        if (Atom.ResidueName=="TIP") return true;
        if (Atom.ResidueName=="TIP3") return true;
        if (Atom.ResidueName=="TIP4") return true;
        if (Atom.ResidueName=="TIP5") return true;
        if (Atom.ResidueName=="WAT") return true;
        return false;
}

bool isWaterOxygen(AtomStruct &Atom)
{
        if (Atom.atomid==OXYGEN && IsWater(Atom)) return true;
        else return false;
}

bool isIon(AtomStruct &Atom)
{
        if (Atom.ResidueName=="SOD") return true;
        else if (Atom.ResidueName=="CLA") return true;
        else if (Atom.ResidueName=="IP") return true;
        else if (Atom.ResidueName=="IM") return true;
        else if (Atom.ResidueName=="MG") return true;
        return false;
}

bool IsProtein(AtomStruct &Atom)
{
        if (Atom.ResidueName=="GLY") return true;
        else if (Atom.ResidueName=="ALA") return true;
        else if (Atom.ResidueName=="VAL") return true;
        else if (Atom.ResidueName=="LEU") return true;
        else if (Atom.ResidueName=="ILE") return true;
        else if (Atom.ResidueName=="SER") return true;
        else if (Atom.ResidueName=="THR") return true;
        else if (Atom.ResidueName=="CYS") return true;
        else if (Atom.ResidueName=="MET") return true;
        else if (Atom.ResidueName=="ASP") return true;
        else if (Atom.ResidueName=="GLU") return true;
        else if (Atom.ResidueName=="ASN") return true;
        else if (Atom.ResidueName=="GLN") return true;
        else if (Atom.ResidueName=="LYS") return true;
        else if (Atom.ResidueName=="ARG") return true;
        else if (Atom.ResidueName=="PHE") return true;
        else if (Atom.ResidueName=="TYR") return true;
        else if (Atom.ResidueName=="PRO") return true;
        else if (Atom.ResidueName=="HIS") return true;
        else if (Atom.ResidueName=="HSE") return true;
        else if (Atom.ResidueName=="HIE") return true;
        else if (Atom.ResidueName=="HSD") return true;
        else if (Atom.ResidueName=="TRP") return true;
        else if (Atom.ResidueName=="HEM") return true;
        else return false;
}

bool IsDuplicate(AtomStruct &Atom, vector<AtomStruct> &Atoms, int index)
{
        int natom=Atoms.size();
        int NumAppearances=0;
        for (int i=index+1;i<natom;i++)
        {
                if (Atoms[i].ResidueNum==Atom.ResidueNum && Atoms[i].AtomName==Atom.AtomName)
                {
                        if (Atoms[i].ChainName==Atom.ChainName && Atoms[i].ID==Atom.ID && Atoms[i].SegID==Atom.SegID)
                        {
                                if (Atoms[i].ResidueName==Atom.ResidueName)
                                {
                                        NumAppearances++;
                                        //cout <<endl;
                                        //PrintAtomInfo(Atom);
                                        //PrintAtomInfo(Atoms[i]);
                                        //cout <<endl;
                                        return true;
                                }
                        }
                }
        }
        return false;
}

int countHetatm(vector<AtomStruct> &Atoms)
{
        int natom=Atoms.size();
        int nhet=0;
        
        for (int i=0;i<natom;i++)
        {
                if (Atoms[i].HetAtom) nhet++;
        }
        return nhet;
}

void RemoveNonHetatm(vector<AtomStruct> &Atoms)
{
        int natom=Atoms.size(), nhet, count=0;
        AtomStruct tempAtom;
        vector<AtomStruct> TempAtoms;
        nhet=countHetatm(Atoms);
        SafeAlloc(TempAtoms, tempAtom, nhet, "TempAtoms");
        for (int i=0;i<natom;i++)
        {
                if (i%100000==0)
                {
                        PrintAtomInfo(Atoms[i]);
                        cout <<"Atoms.HetAtom= "<<Atoms[i].HetAtom<<endl;
                }
                if (Atoms[i].HetAtom) 
                {
                        TempAtoms[count]=Atoms[i];
                        count++;
                }
        }
        Atoms=TempAtoms;
}

void RemoveNonCA(vector<AtomStruct> &Atoms)
{
        int natom=Atoms.size();
        vector<AtomStruct> TempAtoms;

        for (int i=0;i<natom;i++)
        {
                if (Atoms[i].AtomName=="CA")
                {
                        SafePushBack(TempAtoms, Atoms[i], "in RemoveNonCa");
                }
        }
        Atoms=TempAtoms;
}

void RemoveDuplicateAtoms(vector<AtomStruct> &Atoms)
{
        int natom=Atoms.size();
        vector<AtomStruct> TempAtoms;
        for (int i=0;i<natom;i++)
        {
                if (!IsDuplicate(Atoms[i], Atoms, i)) 
                {
                        SafePushBack(TempAtoms, Atoms[i], "TempAtoms in RemoveDuplicateAtoms");
                }
        }
        Atoms=TempAtoms;
}

void RemoveProteinAtoms(vector<AtomStruct> &Atoms)
{
        vector<AtomStruct> TempAtoms;
        int i, natom=Atoms.size();
        for (i=0;i<natom;i++)
        {
                if (IsWater(Atoms[i])) SafePushBack(TempAtoms, Atoms[i], "TempAtoms in RemoveProteinAtoms");
        }
        Atoms=TempAtoms;
}

void RemoveNonProteinAtoms(vector<AtomStruct> &Atoms)
{
        vector<AtomStruct> TempAtoms;
        int natom=Atoms.size();
        for (int i=0;i<natom;i++)
        {
                if (IsProtein(Atoms[i])) SafePushBack(TempAtoms, Atoms[i], "TempAtoms in RemoveNonProteinAtoms");
        }
        Atoms=TempAtoms;
}

void RemoveNonRnaAtoms(vector<AtomStruct> &Atoms)
{
        vector<AtomStruct> TempAtoms;
        int natom=Atoms.size();
        for (int i=0;i<natom;i++)
        {
                if (isNA(Atoms[i])) SafePushBack(TempAtoms, Atoms[i], "TempAtoms in RemoveNonProteinAtoms");
        }
        Atoms=TempAtoms;
}

void RemoveWaters(vector<AtomStruct> &Atoms)
{
        vector<AtomStruct> TempAtoms;
        int i, natom=Atoms.size();
        for (i=0;i<natom;i++)
        {
                if (!IsWater(Atoms[i])) SafePushBack(TempAtoms, Atoms[i], "TempAtoms in RemoveWaters");
        }
        Atoms=TempAtoms;
}

void RemoveUnknownResidues(vector<AtomStruct> &Atoms)
{
        vector<AtomStruct> TempAtoms;
        int i, natom=Atoms.size();
        for (i=0;i<natom;i++)
        {
                if (IsProtein(Atoms[i]) || IsWater(Atoms[i]) || isNA(Atoms[i])) SafePushBack(TempAtoms, Atoms[i], "TempAtoms");
        }
        Atoms=TempAtoms;
}

void RemoveHydrogens(vector<AtomStruct> &Atoms1, vector<AtomStruct> &Atoms2)
{
        int natom=Atoms1.size();

        Atoms2.clear();
        for (int i=0;i<natom;i++)
        {
                if (Atoms1[i].AtomName.substr(0,1)!="H")
                {
                        SafePushBack(Atoms2, Atoms1[i], "RemoveHydrogens");
                }
        }
}
/*
   int FindNumProteinAtoms(vector<AtomStruct> &Atoms)
   {
   int n, TotalParticles=Atoms.size();
   for (n=0;n<TotalParticles;n++) 
   {
//PrintAtomInfo(Atoms[n]);
//cout <<"atomid= "<<Atoms[n].atomid<<endl;
if (IsWater(Atoms[n])) return n;
if (!IsProtein(Atoms[n])) return n;
if (Atoms[n].atomid==HydrationShell) return n;
if (Atoms[n].atomid==ExcludedVolume) return n;
}
return TotalParticles;
}
*/
int FindNumProteinAtoms(vector<AtomStruct> &Atoms)
{
        int TotalParticles=Atoms.size();
        int numProteinAtoms=0;
        for (int n=0;n<TotalParticles;n++) 
        {
                if (Atoms[n].solute) numProteinAtoms++;
        }
        return numProteinAtoms;
}

int countIons(vector<AtomStruct> &Atoms)
{
        int natom=Atoms.size();
        int nion=0;
        for (int i=0;i<natom;i++)
        {
                if (isIon(Atoms[i])) nion++;
        }
        return nion;
}

void SetSolute(vector<AtomStruct> &Atoms)
{
        int TotalParticles=Atoms.size();
        for (int n=0;n<TotalParticles;n++) 
        {
                Atoms[n].solute=true;
                if (IsWater(Atoms[n])) Atoms[n].solute=false;
                if (!IsProtein(Atoms[n]) && !isNA(Atoms[n])) Atoms[n].solute=false;
                if (Atoms[n].atomid==HydrationShell) Atoms[n].solute=false;
                if (Atoms[n].atomid==ExcludedVolume) Atoms[n].solute=false;
        }
}

void SetAllToSolute(vector<AtomStruct> &Atoms)
{
        int natom=Atoms.size();

        for (int i=0;i<natom;i++)
        {
                Atoms[i].solute=true;
        }
}

bool CalcCenterCa(vector<AtomStruct> &Atoms, Real &xave, Real &yave, Real &zave)
{
        Real xsum, ysum, zsum;
        Real NumAve=0;
        int i, natom=Atoms.size();
        xsum=ysum=zsum=0;
        //cout <<"natom= "<<natom;
        for (i=0;i<natom;i++)
        {
                //cout <<Atoms[i].AtomName<<endl;
                if (Atoms[i].AtomName=="CA")
                {
                        xsum+=Atoms[i].x;
                        ysum+=Atoms[i].y;
                        zsum+=Atoms[i].z;
                        NumAve+=1.0;
                }
        }
        if (NumAve==0)
        {
                cout <<"ERROR: In CalcCenter.  No CA atoms."<<endl;
                xave=0;
                yave=0;
                zave=0;
                return false;
        }
        xave=xsum/NumAve;
        yave=ysum/NumAve;
        zave=zsum/NumAve;

        return true;
}

bool CalcCenter(vector<AtomStruct> &Atoms, Real &xave, Real &yave, Real &zave)
{
        Real xsum, ysum, zsum;
        Real NumAve=0;
        int i, natom=Atoms.size();
        xsum=ysum=zsum=0;
        //cout <<"natom= "<<natom<<endl;
        //cout <<"In CalcCenter"<<endl;
        PrintAtomInfo(Atoms[0]);
        for (i=0;i<natom;i++)
        {
                //PrintAtomInfo(Atoms[i]);
                if (Atoms[i].AtomName=="CA" || Atoms[i].AtomName=="O1P")
                {
                        //cout <<Atoms[i].AtomName<<"."<<endl;
                        xsum+=Atoms[i].x;
                        ysum+=Atoms[i].y;
                        zsum+=Atoms[i].z;
                        NumAve+=1.0;
                        //cout <<"NumAve= "<<NumAve<<endl;
                }
        }
        if (NumAve==0)
        {
                cout <<"ERROR: In CalcCenter.  No CA or O1P atoms."<<endl;
                xave=0;
                yave=0;
                zave=0;
                return false;
        }
        xave=xsum/NumAve;
        yave=ysum/NumAve;
        zave=zsum/NumAve;

        return true;
}

bool CalcCenterSelection(vector<AtomStruct> &Atoms, Real &xave, Real &yave, Real &zave, string selection)
{
        Real xsum, ysum, zsum;
        Real NumAve=0;
        int i, natom=Atoms.size();
        xsum=ysum=zsum=0;
        //cout <<"natom= "<<natom;
        for (i=0;i<natom;i++)
        {
                //cout <<Atoms[i].AtomName<<endl;
                if (Atoms[i].AtomName==selection)
                {
                        xsum+=Atoms[i].x;
                        ysum+=Atoms[i].y;
                        zsum+=Atoms[i].z;
                        NumAve+=1.0;
                }
        }
        if (NumAve==0)
        {
                cout <<"ERROR: In CalcCenter.  No "<<selection<<" atoms."<<endl;
                xave=0;
                yave=0;
                zave=0;
                return false;
        }
        xave=xsum/NumAve;
        yave=ysum/NumAve;
        zave=zsum/NumAve;

        return true;
}

void CalcCenterAll(vector<AtomStruct> &Atoms, Real &xave, Real &yave, Real &zave)
{
        Real xsum, ysum, zsum;
        Real TotalMass=0;
        int natom=Atoms.size();
        xsum=ysum=zsum=0;
        for (int i=0;i<natom;i++)
        {
                xsum+=Atoms[i].x*Atoms[i].Electrons;
                ysum+=Atoms[i].y*Atoms[i].Electrons;
                zsum+=Atoms[i].z*Atoms[i].Electrons;
                TotalMass+=Atoms[i].Electrons;
        }
        if (TotalMass==0)
        {
                cout <<"ERROR: In CalcCenterAll.  TotalMass=0."<<endl;
                exit(EXIT_FAILURE);
        }
        xave=xsum/TotalMass;
        yave=ysum/TotalMass;
        zave=zsum/TotalMass;

}

void center(vector<AtomStruct> &Atoms)
{
        Real xave, yave, zave;
        CalcCenter(Atoms, xave, yave, zave);
        MoveAtoms(Atoms, -xave, -yave, -zave);
}

void selectAtomsByName(vector<AtomStruct> &Atoms, vector<AtomStruct> &SelectedAtoms, string SelectionName)
{
        int natom=Atoms.size();
        SelectedAtoms.clear();
        for (int i=0;i<natom;i++)
        {
                if (Atoms[i].AtomName==SelectionName)
                {
                        SafePushBack(SelectedAtoms, Atoms[i], "SelectedAtoms");
                }
        }
}

void selectAtomsByName(vector<AtomStruct> &Atoms, string SelectionName)
{
        vector<AtomStruct> tempAtoms;

        selectAtomsByName(Atoms, tempAtoms, SelectionName);
        Atoms=tempAtoms;
}

void removeSelection(vector<AtomStruct> &Atoms, string selection)
{
        int natom=Atoms.size();
        vector<AtomStruct> tempAtoms;

        for (int i=0;i<natom;i++)
        {
                Trim(Atoms[i].AtomName);
                //PrintAtomInfo(Atoms[i]);
                if (Atoms[i].AtomName!=selection)
                {
                        //cout <<"Added Atom"<<endl;
                        SafePushBack(tempAtoms, Atoms[i], "tempAtoms in removeSelection");
                }
        }
        Atoms=tempAtoms;
}

void removeResidue(vector<AtomStruct> &Atoms, string selection)
{
        int natom=Atoms.size();
        vector<AtomStruct> tempAtoms;

        for (int i=0;i<natom;i++)
        {
                Trim(Atoms[i].ResidueName);
                if (Atoms[i].AtomName!=selection)
                {
                        SafePushBack(tempAtoms, Atoms[i], "tempAtoms in removeResidue");
                }
        }
        Atoms=tempAtoms;
}

void CenterAtomsInBox(vector<AtomStruct> &Atoms, string SelectionName)
{
        Real XMax, YMax, ZMax;
        Real XMin, YMin, ZMin;
        Real XMove, YMove, ZMove;
        vector<AtomStruct> tempAtoms;
        selectAtomsByName(Atoms, tempAtoms, SelectionName);
        MinMaxProtein(XMin, YMin, ZMin, XMax, YMax, ZMax, tempAtoms);
        cout <<"XMin= "<<XMin<<" YMin= "<<YMin<<" ZMin= "<<ZMin<<endl;
        cout <<"XMax= "<<XMax<<" YMax= "<<YMax<<" ZMax= "<<ZMax<<endl;
        XMove=-(XMax+XMin)*0.5;
        YMove=-(YMax+YMin)*0.5;
        ZMove=-(ZMax+ZMin)*0.5;
        //PrintPdb("/home2/jouko/Density/TempPdb/reference.pdb", Atoms);
        MoveAtoms(Atoms, XMove, YMove, ZMove);
        cout <<"XMove= "<<XMove<<" YMove= "<<YMove<<" ZMove= "<<ZMove<<endl;
        //PrintPdb("/home2/jouko/Density/TempPdb/referenceb.pdb", Atoms);
}

void CenterAtomsInBox(vector<AtomStruct> &Atoms)
{
        Real XMax, YMax, ZMax;
        Real XMin, YMin, ZMin;
        Real XMove, YMove, ZMove;
        MinMaxProtein(XMin, YMin, ZMin, XMax, YMax, ZMax, Atoms);
        cout <<"XMin= "<<XMin<<" YMin= "<<YMin<<" ZMin= "<<ZMin<<endl;
        cout <<"XMax= "<<XMax<<" YMax= "<<YMax<<" ZMax= "<<ZMax<<endl;
        XMove=-(XMax+XMin)*0.5;
        YMove=-(YMax+YMin)*0.5;
        ZMove=-(ZMax+ZMin)*0.5;
        //PrintPdb("/home2/jouko/Density/TempPdb/reference.pdb", Atoms);
        MoveAtoms(Atoms, XMove, YMove, ZMove);
        cout <<"XMove= "<<XMove<<" YMove= "<<YMove<<" ZMove= "<<ZMove<<endl;
        //PrintPdb("/home2/jouko/Density/TempPdb/referenceb.pdb", Atoms);
}

void AllignCenters(vector<AtomStruct> &Atoms1, vector<AtomStruct> &Atoms2)
{
        Real xave1, yave1, zave1;
        Real xave2, yave2, zave2;
        Real XMove, YMove, ZMove;
        CalcCenter(Atoms1, xave1, yave1, zave1);
        CalcCenter(Atoms2, xave2, yave2, zave2);
        PrintAtomInfo(Atoms1[0]);
        PrintAtomInfo(Atoms2[0]);
        XMove=xave1-xave2;
        YMove=yave1-yave2;
        ZMove=zave1-zave2;
        cout <<"xave1= "<<xave1<<" yave1= "<<yave1<<" zave1= "<<zave1<<endl;
        cout <<"xave2= "<<xave2<<" yave2= "<<yave2<<" zave2= "<<zave2<<endl;
        cout <<"XMove= "<<XMove<<" YMove= "<<YMove<<" ZMove= "<<ZMove<<endl;
        MoveAtoms(Atoms2, XMove, YMove, ZMove);
}

Real CalcRg(vector<AtomStruct> &Atoms)
{
        int natom=Atoms.size();
        Real dist2, distsqr, TotalMass, rg;
        Real dx, dy, dz;
        Real xave, yave, zave;

        TotalMass=0;
        CalcCenterAll(Atoms, xave, yave, zave);
        cout <<"xave= "<<xave<<" yave= "<<yave<<" zave= "<<zave<<endl;
        for (int i=0;i<natom;i++)
        {
                dx=xave-Atoms[i].x;
                dy=yave-Atoms[i].y;
                dz=zave-Atoms[i].z;
                dist2=dx*dx+dy*dy+dz*dz;
                distsqr+=dist2*Atoms[i].Electrons;
                TotalMass+=Atoms[i].Electrons;
        }

        if (distsqr/TotalMass>0) rg=sqrt(distsqr/TotalMass);
        else rg=0;

        return rg;
}

void SetNumElectrons(Real NumElectrons[])
{
        for (int i=0;i<NumAtomTypes;i++) NumElectrons[i]=0;
        NumElectrons[HYDROGEN]=1.0;
        NumElectrons[CARBON]=6.0;
        NumElectrons[NITROGEN]=7.0;
        NumElectrons[OXYGEN]=8.0;
        NumElectrons[SULFUR]=16.0;
        NumElectrons[IRON]=26.0;
        NumElectrons[ExcludedVolume]=0.0;
        NumElectrons[HydrationShell]=0.0;
        NumElectrons[PHOSPHORUS]=15.0;
        NumElectrons[CL]=18.0;
        NumElectrons[NA_Plus]=10.0;
}

void SetAtomElectrons(vector<AtomStruct> &Atoms, Real NumElectrons[])
{
        int natom=Atoms.size();

        for (int i=0;i<natom;i++)
        {
                Atoms[i].Electrons=NumElectrons[Atoms[i].atomid];
        }
}

void ApplyPeriodicBoundaryConditions(vector<AtomStruct> &Atoms, Real xmin, Real ymin, Real zmin, Real xmax, Real ymax, Real zmax)
{
        int TotalParticles=Atoms.size();
        Real XBoxLength, YBoxLength, ZBoxLength;

        XBoxLength=xmax-xmin;
        YBoxLength=ymax-ymin;
        ZBoxLength=zmax-zmin;
        if (XBoxLength<=0 || YBoxLength<=0 || ZBoxLength<=0)
        {
                cout <<"ERROR: One of the box dimensions is 0 or less."<<endl;
                cout <<"XBoxLength= "<<XBoxLength<<" YBoxLength= "<<YBoxLength<<" ZBoxLength= "<<ZBoxLength<<endl;
                exit(EXIT_FAILURE);
        }
        for (int n=0;n<TotalParticles;n++)
        {
                if (Atoms[n].x>xmax) Atoms[n].x-=XBoxLength;
                if (Atoms[n].x<xmin) Atoms[n].x+=XBoxLength;
                if (Atoms[n].y>ymax) Atoms[n].y-=YBoxLength;
                if (Atoms[n].y<ymin) Atoms[n].y+=YBoxLength;
                if (Atoms[n].z>zmax) Atoms[n].z-=ZBoxLength;
                if (Atoms[n].z<zmin) Atoms[n].z+=ZBoxLength;
        }
}

void ApplyPeriodicBoundaryConditions(vector<AtomStruct> &Atoms, Real XBoxLength, Real YBoxLength, Real ZBoxLength)
{
        Real xmin, ymin, zmin;
        Real xmax, ymax, zmax;

        xmin=-0.5*XBoxLength;
        ymin=-0.5*YBoxLength;
        zmin=-0.5*ZBoxLength;
        xmax=0.5*XBoxLength;
        ymax=0.5*YBoxLength;
        zmax=0.5*ZBoxLength;

        ApplyPeriodicBoundaryConditions(Atoms, xmin, ymin, zmin, xmax, ymax, zmax);
}

void GetAtomsWithinRadius(vector<AtomStruct> &Atoms, Real xgrid, Real ygrid, Real zgrid, Real radius)
{
        int natom=Atoms.size();
        Real dist, r2;
        Real dx, dy, dz;
        vector<AtomStruct> TempAtoms;

        r2=radius*radius;
        cout <<"xgrid= "<<xgrid<<" ygrid= "<<ygrid<<" zgrid= "<<zgrid<<" r2= "<<r2<<endl;
        for (int i=0;i<natom;i++)
        {
                dx=xgrid-Atoms[i].x;
                dy=ygrid-Atoms[i].y;
                dz=zgrid-Atoms[i].z;
                dist=dx*dx+dy*dy+dz*dz;
                if (dist<r2) SafePushBack(TempAtoms, Atoms[i], "in GetAtomsWithinRadius");
        }
        Atoms=TempAtoms;
}

void Image(vector<AtomStruct> &Atoms, Real &XBoxLength, Real &YBoxLength, Real &ZBoxLength, int NumPerEdge)
{
        //Takes a system and replicates it.
        AtomStruct Atom;
        int TotalParticles2;
        int TotalParticles=Atoms.size();
        int a, b, c;
        Real Xlength, Ylength, Zlength;
        vector<AtomStruct> TempAtoms;

        TotalParticles2=NumPerEdge*NumPerEdge*NumPerEdge*TotalParticles;
        Xlength=Real(NumPerEdge)*XBoxLength;
        Ylength=Real(NumPerEdge)*YBoxLength;
        Zlength=Real(NumPerEdge)*ZBoxLength;

        for (a=1;a<=NumPerEdge;a++)
        {
                for (b=1;b<=NumPerEdge;b++)
                {
                        for (c=1;c<=NumPerEdge;c++)
                        {
                                for (int n=0;n<TotalParticles;n++)
                                {
                                        CopyAtom(Atoms[n], Atom);
                                        Atom.x=Atoms[n].x+Real(a)*XBoxLength;
                                        Atom.y=Atoms[n].y+Real(b)*YBoxLength;
                                        Atom.z=Atoms[n].z+Real(c)*ZBoxLength;
                                        SafePushBack(TempAtoms, Atom, "TempAtoms in Image");
                                }
                        }
                }
        }
        Atoms=TempAtoms;

        cout <<"XBoxLength="<<XBoxLength<<" YBoxLength="<<YBoxLength<<" ZBoxLength="<<ZBoxLength<<endl;
        MoveAtoms(Atoms, -Real(NumPerEdge-1)*XBoxLength*0.5, -Real(NumPerEdge-1)*YBoxLength*0.5, -Real(NumPerEdge-1)*ZBoxLength*0.5);
        ApplyPeriodicBoundaryConditions(Atoms, -0.5*Xlength, -0.5*Ylength, -0.5*Zlength, Xlength, Ylength, Zlength);
        XBoxLength=Xlength;
        YBoxLength=Ylength;
        ZBoxLength=Zlength;
        cout <<"XBoxLength="<<XBoxLength<<" YBoxLength="<<YBoxLength<<" ZBoxLength="<<ZBoxLength<<endl;
        cout <<"TotalParticles= "<<Atoms.size()<<endl;
}

void Image(vector<AtomStruct> &Atoms, int axis, Real &Length, int NumPerEdge)
{
        //Takes a system and replicates it.
        AtomStruct Atom;
        int TotalParticles2;
        int TotalParticles=Atoms.size();
        int natom;
        Real NewLength, MoveBy;
        vector<AtomStruct> TempAtoms;

        TotalParticles2=NumPerEdge*TotalParticles;
        NewLength=Real(NumPerEdge)*Length;
        SafeAlloc(TempAtoms, Atom, TotalParticles2, "in Image");
        natom=0;
        for (int a=0;a<NumPerEdge;a++)
        {
                for (int n=0;n<TotalParticles;n++)
                {
                        CopyAtom(Atoms[n], Atom);
                        if (axis==X) Atom.x=Atoms[n].x+Real(a)*Length;
                        if (axis==Y) Atom.y=Atoms[n].y+Real(a)*Length;
                        if (axis==Z) Atom.z=Atoms[n].z+Real(a)*Length;
                        TempAtoms[natom]=Atom;
                        natom++;
                }
        }
        Atoms=TempAtoms;
        //PrintPDB("/home2/jouko/Density/Pdb/BeforeMove.pdb", Atoms);
        MoveBy=-Real(NumPerEdge-1)*Length*0.5;
        if (axis==X) MoveAtoms(Atoms, MoveBy, 0, 0);
        if (axis==Y) MoveAtoms(Atoms, 0, MoveBy, 0);
        if (axis==Z) MoveAtoms(Atoms, 0, 0, MoveBy);
        //PrintPDB("/home2/jouko/Density/Pdb/AfterMove.pdb", Atoms);
        Length=NewLength;
        cout <<"TotalParticles= "<<Atoms.size()<<endl;
}

void CharmmAlias(AtomStruct &Atom)
{
        /*
        if (Atom.ResidueName=="HIS") Atom.ResidueName="HSE";
        if (Atom.ResidueName=="HIE") 
        {
                Atom.ResidueName="HSE";
        }
        */
        if (Atom.AtomName=="HW1") Atom.AtomName="HT1";
        if (Atom.AtomName=="HW2") Atom.AtomName="HT2";
        if (Atom.AtomName=="HW3") Atom.AtomName="HT3";
        if (Atom.AtomName=="OW1") Atom.AtomName="OT1";
        if (Atom.AtomName=="OW2") Atom.AtomName="OT2";
        if (Atom.AtomName=="OXT") Atom.AtomName="OT";
        if (Atom.ResidueName=="ILE" && Atom.AtomName=="CD1")
        {
                Atom.AtomName="CD";
        }
}

void AtomAlias(vector<AtomStruct> &Atoms)
{
        int natom=Atoms.size();

        for (int i=0;i<natom;i++)
        {
                CharmmAlias(Atoms[i]);
        }
}

void SetCharmmAtomNames(vector<AtomStruct> &CharmmAtoms, vector<AtomStruct> &Atoms)
{
        bool BoolPrint=false;
        int NumCharmmAtoms=CharmmAtoms.size(), NumAtoms=Atoms.size();
        cout <<"In SetCharmmAtomNames"<<endl;
        cout <<"NumAtoms= "<<NumAtoms<<endl;
        for (int i=0;i<NumAtoms;i++)
        {
                //if (Atoms[i].AtomName=="H5\'\'") BoolPrint=true;
                //else BoolPrint=false;
                //PrintAtomInfo(Atoms[i]);
                CharmmAlias(Atoms[i]);
                //if (Atoms[i].AtomName=="H5T") BoolPrint=true;
                //else BoolPrint=false;
                //PrintAtomInfo(Atoms[i]);
                for (int j=0;j<NumCharmmAtoms;j++)
                {
                        if (BoolPrint)
                        {
                                cout <<endl;
                                cout <<"CharmmAtomInfo"<<endl;
                                PrintAtomInfo(CharmmAtoms[j]);
                                PrintAtomInfo(Atoms[i]);
                                cout <<"CharmmAtomName= "<<Atoms[i].CharmmAtomName<<endl;
                                cout <<endl;
                        }
                        if (CharmmAtoms[j].ResidueName==Atoms[i].ResidueName)
                        {
                                if (WildCardCompare(Atoms[i].AtomName, CharmmAtoms[j].AtomName))
                                {
                                        Atoms[i].CharmmAtomName=CharmmAtoms[j].CharmmAtomName;
                                        if (Atoms[i].charge==0) Atoms[i].charge=CharmmAtoms[j].charge;
                                        if (BoolPrint) PrintAtomInfo(Atoms[i]);
                                        if (BoolPrint) cout <<"CharmmAtomName= "<<Atoms[i].CharmmAtomName<<endl;
                                        break;
                                }
                        }
                        if (BoolPrint)
                        {
                                cout <<"About to leave SetCharmmAtomNames"<<endl;
                                PrintAtomInfo(Atoms[i]);
                                cout <<"CharmmAtomName= "<<Atoms[i].CharmmAtomName<<endl;
                        }
                }
        }
        cout <<"Leaving SetCharmmAtomNames"<<endl;
}

void CharmmAtomNameToVDW(vector<AtomStruct> &CharmmAtoms, AtomStruct &Atom)
{
        bool BoolPrint=false;
        int NumCharmmAtoms;
        Real DefaultVDW=1.5;
        Atom.vdw=UNK_REAL;
        //if (Atom.AtomName.substr(0,2)=="HT") Atom.CharmmAtomName="HC";
        //if (Atom.AtomName.substr(0,2)=="OT") Atom.CharmmAtomName="OC";
        NumCharmmAtoms=CharmmAtoms.size();
        //cout <<"In CharmmAtomNameToVDW"<<endl;
        //if (Atom.ResidueName=="TRP" && Atom.AtomName=="CD1") BoolPrint=true;
        //else BoolPrint=false;
        
        for (int j=0;j<NumCharmmAtoms;j++)
        {
                if (BoolPrint) 
                {
                        PrintAtomInfo(Atom);
                        cout <<"Atom.CharmmAtomName= "<<Atom.CharmmAtomName<<endl;
                        PrintAtomInfo(CharmmAtoms[j]);
                        cout <<"CharmmAtoms["<<j<<"].CharmmAtomName= "<<CharmmAtoms[j].CharmmAtomName<<endl;
                        cout <<endl;
                }

                //cout <<"Atom.CharmmAtomName= "<<Atom.CharmmAtomName<<" CharmmAtoms["<<j<<"].CharmmAtomName== "<<CharmmAtoms[j].CharmmAtomName<<endl;
                if (WildCardCompare(CharmmAtoms[j].CharmmAtomName, Atom.CharmmAtomName))
                {
                        Atom.vdw=CharmmAtoms[j].vdw;
                        Atom.epsilon=CharmmAtoms[j].epsilon;
                        Atom.vdw14=CharmmAtoms[j].vdw14;
                        Atom.epsilon14=CharmmAtoms[j].epsilon14;
                        //Atom.charge=CharmmAtoms[j].charge;
                        break;
                }
        }
        if (Atom.vdw==UNK_REAL)
        {
                Atom.vdw=DefaultVDW;
                cout <<"Warning.  Unable to set vdw for atom "<<endl;
                PrintAtomInfo(Atom);
                cout <<"Atom.CharmmAtomName= "<<Atom.CharmmAtomName<<endl;
                cout <<"Setting vdw to "<<DefaultVDW<<endl;
        }
}

void CharmmAtomNamesToVDW(vector<AtomStruct> &CharmmAtoms, vector<AtomStruct> &Atoms)
{
        int NumAtoms=Atoms.size();
        for (int i=0;i<NumAtoms;i++)
        {
                CharmmAtomNameToVDW(CharmmAtoms, Atoms[i]);
                //cout <<"vdw= "<<Atoms[i].vdw<<endl;
                //PrintAtomInfo(Atoms[i]);
        }
}

void ScaleVDW(Real VDWScale, Real VDWOffSet, vector<AtomStruct> &Atoms)
{
        int i, natom=Atoms.size();
        for (i=0;i<natom;i++)
        {
                Atoms[i].vdw=Atoms[i].vdw*VDWScale+VDWOffSet;
        }
}

void AddChargeRadii(vector<AtomStruct> &Atoms, Real ChargeRadiiScale)
{
        int natom=Atoms.size();
        //cout <<"In AddChargeRadii"<<endl;
        for (int i=0;i<natom;i++)
        {
                Atoms[i].HydrationRadius+=abs(Atoms[i].charge)*ChargeRadiiScale;
                //cout <<"Atoms["<<i<<"].charge= "<<Atoms[i].charge<<" ChargeRadiiScale= "<<ChargeRadiiScale<<endl;
                //cout <<"Atoms["<<i<<"].HydrationRadius= "<<Atoms[i].HydrationRadius<<endl;
        }
}

void SetVDW(string PrmFile, string RtfFile, vector<AtomStruct> &Atoms)
{
        vector<AtomStruct> CharmmAtoms, TempCharmmAtoms;
        ReadRtfFile(RtfFile, CharmmAtoms);
        ReadPrmFile(PrmFile, TempCharmmAtoms);
        SetCharmmAtomNames(CharmmAtoms, Atoms);
        CharmmAtomNamesToVDW(TempCharmmAtoms, Atoms);
}

void GetCharges(vector<AtomStruct> &PsfAtoms, vector<AtomStruct> &Atoms)
{
        int natom=Atoms.size();
        int npsf=PsfAtoms.size();
        if (natom>npsf)
        {
                cout <<"ERROR: Atom number mismatch."<<endl;
                cout <<"natom= "<<natom<<" PsfAtoms.size()= "
                        <<PsfAtoms.size()<<endl;
                exit(EXIT_FAILURE);
        }

        for (int i=0;i<natom;i++)
        {
                Atoms[i].charge=PsfAtoms[i].charge;
        }
}

void GetCharges(string PsfFile, vector<AtomStruct> &Atoms)
{
        vector<AtomStruct> PsfAtoms;
        ReadPsf(PsfFile, PsfAtoms);
        GetCharges(PsfAtoms, Atoms);
}

void GetCharmmAtomNames(vector<AtomStruct> &PsfAtoms, vector<AtomStruct> &Atoms)
{
        int natom=Atoms.size();
        for (int i=0;i<natom;i++)
        {
                Atoms[i].CharmmAtomName=PsfAtoms[i].CharmmAtomName;
        }
}

void GetCharmmAtomNames(string PsfFile, vector<AtomStruct> &Atoms)
{
        vector<AtomStruct> PsfAtoms;
        ReadPsf(PsfFile, PsfAtoms);
        GetCharmmAtomNames(PsfAtoms, Atoms);
}

void GetConnectivity(vector<AtomStruct> &PsfAtoms, vector<AtomStruct> &Atoms)
{
        int natom=Atoms.size();
        //cout <<"natom= "<<natom<<endl;
        for (int i=0;i<natom;i++)
        {
                Atoms[i].connectivity=PsfAtoms[i].connectivity;
        }
        cout <<"Leaving GetConnectivity"<<endl;
}

void GetConnectivity(string PsfFile, vector<AtomStruct> &Atoms)
{
        vector<AtomStruct> PsfAtoms;
        cout <<"About to enter ReadPsf"<<endl;
        ReadPsf(PsfFile, PsfAtoms);
        cout <<"Read PsfFile"<<endl;
        GetConnectivity(PsfAtoms, Atoms);
}

void GetPsfInfo(string PsfFile, vector<AtomStruct> &Atoms)
{
        vector<AtomStruct> PsfAtoms;
        ReadPsf(PsfFile, PsfAtoms);
        GetConnectivity(PsfAtoms, Atoms);
        GetCharmmAtomNames(PsfAtoms, Atoms);
        GetCharges(PsfAtoms, Atoms);
}

Real calcTotalCharge(vector<AtomStruct> &Atoms)
{
        int natom=Atoms.size();
        Real totalCharge=0;
        for (int i=0;i<natom;i++)
        {
                totalCharge+=Atoms[i].charge;
        }
        return totalCharge;
}

void GetFirstLastResidue(vector<AtomStruct> &Atoms, int &FirstResidue, int &LastResidue)
{
        int natom=Atoms.size();
        FirstResidue=Atoms[0].ResidueNum;
        LastResidue=Atoms[natom-1].ResidueNum;
}

bool isNA(AtomStruct &Atom)
{
        if (Atom.ResidueName=="URA") return true;
        else if (Atom.ResidueName=="ADE") return true;
        else if (Atom.ResidueName=="CYT") return true;
        else if (Atom.ResidueName=="GUA") return true;
        else if (Atom.ResidueName=="THY") return true;
        else if (Atom.ResidueName=="U") return true;
        else if (Atom.ResidueName=="A") return true;
        else if (Atom.ResidueName=="C") return true;
        else if (Atom.ResidueName=="G") return true;
        else if (Atom.ResidueName=="T") return true;

        return false;
}

bool hasNA(vector<AtomStruct> &Atoms)
{
        int natom=Atoms.size();

        for (int i=0;i<natom;i++)
        {
                if (Atoms[i].AtomName=="O1P") return true;
        }
        return false;
}

bool isCompleteResidue(vector<AtomStruct> &Atoms, int ResidueNum)
{
        bool hasCA=false, hasN=false, hasC=false;
        int natom=Atoms.size();

        for (int i=0;i<natom;i++)
        {
                if (Atoms[i].ResidueNum==ResidueNum)
                {
                        if (Atoms[i].AtomName=="CA") hasCA=true;
                        if (Atoms[i].AtomName=="N") hasN=true;
                        if (Atoms[i].AtomName=="C") hasC=true;
                }
        }
        cout <<"hasCA= "<<hasCA<<" hasN= "<<hasN<<" hasC= "<<hasC<<endl;
        return (hasCA && hasN && hasC);
}

bool hasChainBreaks(vector<AtomStruct> &Atoms)
{
        int natom=Atoms.size();
        int FirstResidue, LastResidue;

        GetFirstLastResidue(Atoms, FirstResidue, LastResidue);
        for (int i=1;i<natom;i++)
        {
                if (Atoms[i].ResidueNum-Atoms[i-1].ResidueNum>1)
                {
                        return true;
                }
        }
        for (int i=FirstResidue;i<=LastResidue;i++)
        {
                if (!isCompleteResidue(Atoms, i)) return true;
        }
        return false;
}

bool isValidProtein(vector<AtomStruct> &Atoms)
{
        if (Atoms.size()==0) return false;
        else if (hasChainBreaks(Atoms)) return false;

        return true;
}

void RenumberResidues(vector<AtomStruct> &Atoms)
{
        int natom=Atoms.size();
        int FirstResidue=Atoms[0].ResidueNum;

        for (int i=0;i<natom;i++)
        {
                Atoms[i].ResidueNum-=(FirstResidue-1);
        }

}

void ParseAtoms(string PdbFile, vector<AtomStruct> &Atoms, bool BoolRemoveWaters, bool BoolRemoveUnknownResidues)
{
        bool BoolRenumberResidues=true, BoolRemoveDuplicates=true;
        string str;
        size_t position;
        position=PdbFile.rfind(".");
        cout <<"In ParseAtoms"<<endl;
        cout <<"PdbFile= "<<PdbFile<<endl;
        cout <<"BoolRemoveWaters= "<<BoolRemoveWaters<<endl;
        cout <<"BoolRemoveUnknownResidues= "<<BoolRemoveUnknownResidues<<endl;
        if (position!=string::npos)
        {
                str=PdbFile.substr(position+1,4);
                if (str=="pdb") ReadPdb(PdbFile, Atoms);
                else ReadXYZ(PdbFile, Atoms);
        }
        else ReadPdb(PdbFile, Atoms);
        cout <<"Atoms.size()= "<<Atoms.size()<<endl;
        if (BoolRemoveWaters) RemoveWaters(Atoms);
        cout <<"Atoms.size()= "<<Atoms.size()<<endl;
        if (BoolRemoveUnknownResidues) RemoveUnknownResidues(Atoms);
        if (BoolRemoveDuplicates) RemoveDuplicateAtoms(Atoms);
        if (BoolRenumberResidues) RenumberResidues(Atoms);
        cout <<"Atoms.size()= "<<Atoms.size()<<endl;
        AssignAtomIDs(Atoms);
        GetAtomType(Atoms);
        GetAtomType2(Atoms);		
}
#endif
