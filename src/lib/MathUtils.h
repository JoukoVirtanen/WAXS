#ifndef _MathUtils_included_
#define _MathUtils_included_

template <typename T>
void MakeDiagonalNonzero(vector< vector<T> > &coefficient, int n);
template <typename T>
void Equation(vector< vector<T> > &coefficient, vector<T> &variable);

# include "TypeDef.h"
# include "VectorManip.h"
# include "Structures.h"
# include "Constants.h"
# include "minimize.h"
# include "IOUtils.h"

struct xyStruct
{
        vector<Real> x, y;
};
/*
template <class Matrix_T>
void PrintMatrix(Matrix_T Matrix)
{
	int i, j, Size1, Size2;
	Size1=Matrix.size();
	Size2=Matrix[0].size();
	for (i=0;i<Size1;i++)
	{
		for (j=0;j<Size2;j++) cout <<Matrix[i][j]<<"\t";
		cout <<endl;
	}
	cout <<endl;
}
*/

void PrintMatrix(vector< vector<Real> > &Mat)
{
	int i, j, Size1, Size2;
	Size1=Mat.size();
	Size2=Mat[0].size();
	for (i=0;i<Size1;i++)
	{
                for (j=0;j<Size2;j++)
                {
                        cout <<setiosflags(ios::left)<<setw(10)<<setprecision(5)<<Mat[i][j]<<"\t";
		}
                cout <<endl;
	}
	cout <<endl;
}

string VectorToStr(VectorStruct &v)
{
        return "("+RealToStr(v.x)+", "+RealToStr(v.y)+", "+RealToStr(v.z)+")";
}

template <class Matrix_T>
void CreateMatrix(Matrix_T &Matrix, int Size1, int Size2)
{
	int i;
	vector<Real> v;
	for (i=0; i<Size2; i++) v.push_back(0);
	for (i=0; i<Size1; i++) Matrix.push_back(v);
}

template <class Matrix_T>
void DeleteMatrix(Matrix_T &Matrix)
{
	int i, Size=Matrix.size();
	for (i=0; i<Size; i++) Matrix.pop_back();
}

bool isEven(int x)
{
        if (x%2==0) return true;
        else return false;
}

bool isOdd(int x)
{
        if (x%2==0) return false;
        else return true;
}

bool isPositiveDefinite(Matrix &mat)
{
        //This does not actually determine if the matrix is positive definite.
        int Size=mat.size();
        Real sum=0;
        for (int i=0;i<Size;i++)
        {
                if (mat[i][i]<0) return false;
        }
        return true;
}

bool QuadraticEquation(Real a, Real b, Real c, Real &low, Real &high)
{
	Real sq;
	sq=b*b-4.0*a*c;
	if (sq<0) return false;
	low=(-b-sqrt(sq))/(2.0*a);
	high=(-b+sqrt(sq))/(2.0*a);
	return true;
}

Real solveTrig(Real a, Real b, Real c)
{
        //Solves a*cos(theta)+b*cos(theta)=c
        //a*cos(theta)+b*cos(theta)=d*sin(theta+phi)
        //d*sin(theta+phi)=d*cos(theta)*sin(phi)+d*sin(theta)*cos(phi)
        //a=d*sin(phi)  b=d*cos(phi)    a/b=tan(phi)
        Real d, phi, theta;
        phi=atan(a/b);
        d=sqrt(a*a+b*b);
        theta=asin(c/d)-phi;   //asin(c/d)+/-d 
        return theta;
}

VectorStruct VectorStruct::operator+ (VectorStruct v)
{
	VectorStruct temp;
	temp.x=x+v.x;
	temp.y=y+v.y;
	temp.z=z+v.z;
	return (temp);
}

VectorStruct VectorStruct::operator- (VectorStruct v)
{
	VectorStruct temp;
	temp.x=x-v.x;
	temp.y=y-v.y;
	temp.z=z-v.z;
	return (temp);
}

VectorStruct VectorStruct::operator* (Real scalar)
{
	VectorStruct temp;
	temp.x=x*scalar;
	temp.y=y*scalar;
	temp.z=z*scalar;
	return (temp);
}

Real VectorStruct::operator* (VectorStruct v)
{
	Real scalar;
	scalar=v.x*x+v.y*y+v.z*z;
	return scalar;
}

VectorStruct VectorStruct::operator/ (Real scalar)
{
	VectorStruct temp;
	temp.x=x/scalar;
	temp.y=y/scalar;
	temp.z=z/scalar;
	return (temp);
}

Real VectorLength(VectorStruct &v)
{
	return sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
}

void VectorNormalize(VectorStruct &v)
{
        Real length=VectorLength(v);
        v=v/length;
}

VectorStruct CalcCrossProduct(VectorStruct &v1, VectorStruct &v2)
{
	VectorStruct product;

	product.x=v1.y*v2.z-v1.z*v2.y;
	product.y=v2.x*v1.z-v2.z*v1.x;
	product.z=v1.x*v2.y-v1.y*v2.x;

	return product;
}

Real calcVectorAngle(VectorStruct &v1, VectorStruct &v2)
{
        Real dotProduct=v1*v2;
        return acos(dotProduct/(VectorLength(v1)*VectorLength(v2)));
}

void CalcEulerAngles(Real x, Real y, Real z, Real &theta, Real &phi)
{
        //Gives the Euler angles to rotate a vector the xyz coordinates given from the z-axis.
        theta=acos(-y/sqrt(1.0-z*z));
        phi=acos(z);
}

template<class pos>
void RotateAroundX(pos &coord, Real theta)
{
        Real tempY, tempZ;
        tempY=coord.y*cos(theta)-coord.z*sin(theta);
        tempZ=coord.y*sin(theta)+coord.z*cos(theta);
        coord.y=tempY;
        coord.z=tempZ;
}

template<class pos>
void RotateAroundY(pos &coord, Real theta)
{
        Real tempX, tempZ;
        tempX=coord.x*cos(theta)-coord.z*sin(theta);
        tempZ=coord.x*sin(theta)+coord.z*cos(theta);
        coord.x=tempX;
        coord.z=tempZ;
}

template<class pos>
void RotateAroundZ(pos &coord, Real theta)
{
        Real tempX, tempY;
        tempX=coord.x*cos(theta)-coord.y*sin(theta);
        tempY=coord.x*sin(theta)+coord.y*cos(theta);
        coord.x=tempX;
        coord.y=tempY;
}

Real coth(Real x)
{
	Real exp2;
	exp2=exp(2.0*x);
	return (exp2+1.0)/(exp2-1.0);
}

Real csch(Real x)
{
	//cout <<"exp("<<x<<")= "<<exp(x)<<endl;
	return 2.0/(exp(x)-exp(-x));
}

/*
void Equation(long double **coefficient, long double variable[], int nequations)
{
	//Finds the solution to a system of equations using gaussian elimination.
        bool verbose;
	int i, m, n;
        long double coefficientnn;
        long double scale;

	verbose=false;

	if (verbose)
	{
		cout <<"nequations= "<<nequations<<endl;
		for (n=0;n<nequations;n++)
		{
			for (m=0;m<=nequations;m++)
			{
				cout <<coefficient[m][n]<<"\t";
			}
			cout <<endl;
		}
		cout <<endl;
	}

	for (n=0;n<nequations-1;n++)
        {
                coefficientnn=coefficient[n][n];
                if (coefficientnn!=0)
                {
                        for (m=n+1;m<nequations;m++)
                        {
                                scale=coefficient[n][m]/coefficientnn;
                                for (i=0;i<nequations+1;i++)
                                {
                                        coefficient[i][m]=coefficient[i][m]-coefficient[i][n]*scale;
                                }
                        }
                }
        }
	
	if (verbose)
	{
        	for (n=0;n<nequations;n++)
        	{
                	for (m=0;m<=nequations;m++)
                	{
                        	cout <<coefficient[m][n]<<"\t";
                	}
                	cout <<endl;
        	}
		cout <<endl;
	}

        for (n=nequations-1;n>0;n=n-1)
        {
                coefficientnn=coefficient[n][n];
                if (coefficientnn!=0)
                {
			for (m=n-1;m>-1;m=m-1)
                        {
                                scale=coefficient[n][m]/coefficientnn;
				//cout <<"coefficient["<<n<<"]["<<m<<"]= "<<coefficient[n][m]<<" coefficientnn= "<<coefficientnn<<" scale= "<<scale<<endl;
                                for (i=n;i<nequations+1;i++)
                                {
                                        coefficient[i][m]=coefficient[i][m]-coefficient[i][n]*scale;
                                }

				if (verbose)
				{
        				for (int n2=0;n2<nequations;n2++)
        				{
                				for (int m2=0;m2<=nequations;m2++)
                				{
                        				cout <<coefficient[m2][n2]<<"\t";
                				}
                				cout <<endl;
        				}
					cout <<endl;
				}
                        }
                }
        }
	if (verbose)
	{
        	for (n=0;n<nequations;n++)
        	{
                	for (m=0;m<=nequations;m++)
                	{
                        	cout <<coefficient[m][n]<<"\t";
                	}
                	cout <<endl;
        	}
		cout <<endl;
	}

        for (n=0;n<nequations;n++)
        {
                variable[n]=coefficient[nequations][n]/coefficient[n][n];
		if (verbose) cout <<"Variable "<<n+1<<" is "<<variable[n]<<"\n";
        }
}

void Equation2(long double coefficient[][500], long double variable[], int nequations)
{
	//Finds the solution to a system of equations using gaussian elimination.
        int i, m, n;
        long double coefficientnn;
        long double scale;

	//cout <<"nequations= "<<nequations<<endl;
	//for (n=0;n<nequations;n++)
	//{
	//	for (m=0;m<=nequations;m++)
	//	{
	//		cout <<coefficient[m][n]<<"\t";
	//	}
	//	cout <<endl;
	//}
	//cout <<endl;

	for (n=0;n<nequations-1;n++)
        {
                coefficientnn=coefficient[n][n];
                if (coefficientnn!=0)
                {
                        for (m=n+1;m<nequations;m++)
                        {
                                scale=coefficient[n][m]/coefficientnn;
                                for (i=0;i<nequations+1;i++)
                                {
                                        coefficient[i][m]=coefficient[i][m]-coefficient[i][n]*scale;
                                }
                        }
                }
        }
        for (n=0;n<nequations;n++)
        {
                for (m=0;m<=nequations;m++)
                {
                        cout <<coefficient[m][n]<<"\t";
                }
                cout <<endl;
        }
	cout <<endl;
        for (n=nequations;n>0;n=n-1)
        {
                coefficientnn=coefficient[n][n];
                if (coefficientnn!=0)
                {
			for (m=n-1;m>-1;m=m-1)
                        {
                                scale=coefficient[n][m]/coefficientnn;
                                for (i=n;i<nequations+1;i++)
                                {
                                        coefficient[i][m]=coefficient[i][m]-coefficient[i][n]*scale;
                                }
                        }
                }
        }
        for (n=0;n<nequations;n++)
        {
                for (m=0;m<=nequations;m++)
                {
                        cout <<coefficient[m][n]<<"\t";
                }
                cout <<endl;
        }
	cout <<endl;
        for (n=0;n<nequations;n++)
        {
                variable[n]=coefficient[nequations][n]/coefficient[n][n];
                cout <<"Variable "<<n+1<<" is "<<variable[n]<<"\n";
        }
}
*/

template <typename T>
void SwapEquations(vector< vector<T> > &coefficient, int i, int j)
{
	int Size=coefficient.size();

	for (int k=0;k<Size;k++)
	{
		swap(coefficient[k][i], coefficient[k][j]);
	}
}

template <typename T>
void MakeDiagonalNonzero(vector< vector<T> > &coefficient, int n)
{
	int nequations=coefficient[0].size();

	for (int i=n+1;i<nequations;i++)
	{
		if (coefficient[i][n]!=0)
		{
			SwapEquations(coefficient, i, n);
			break;
		}
	}
}

template <typename T>
void Equation(vector< vector<T> > &coefficient, vector<T> &variable)
{
	//Finds the solution to a system of equations using gaussian elimination.
        bool verbose;
	int i, m, n;
        T coefficientnn;
        T scale;
        int nequations=variable.size();

	verbose=true;

	if (verbose)
	{
		cout <<"nequations= "<<nequations<<endl;
		for (n=0;n<nequations;n++)
		{
			for (m=0;m<=nequations;m++)
			{
				cout <<coefficient[m][n]<<"\t";
			}
			cout <<endl;
		}
		cout <<endl;
	}

	for (n=0;n<nequations-1;n++)
        {
		if (coefficient[n][n]==0) MakeDiagonalNonzero(coefficient, n);
                coefficientnn=coefficient[n][n];
                if (coefficientnn!=0)
                {
                        for (m=n+1;m<nequations;m++)
                        {
                                scale=coefficient[n][m]/coefficientnn;
                                for (i=0;i<nequations+1;i++)
                                {
                                        coefficient[i][m]=coefficient[i][m]-coefficient[i][n]*scale;
                                }
                        }
                }
        }
	
	if (verbose)
	{
        	for (n=0;n<nequations;n++)
        	{
                	for (m=0;m<=nequations;m++)
                	{
                        	cout <<coefficient[m][n]<<"\t";
                	}
                	cout <<endl;
        	}
		cout <<endl;
	}

        for (n=nequations-1;n>0;n=n-1)
        {
                coefficientnn=coefficient[n][n];
                if (coefficientnn!=0)
                {
			for (m=n-1;m>-1;m=m-1)
                        {
                                scale=coefficient[n][m]/coefficientnn;
				//cout <<"coefficient["<<n<<"]["<<m<<"]= "<<coefficient[n][m]<<" coefficientnn= "<<coefficientnn<<" scale= "<<scale<<endl;
                                for (i=n;i<nequations+1;i++)
                                {
                                        coefficient[i][m]=coefficient[i][m]-coefficient[i][n]*scale;
                                }

				if (verbose)
				{
        				for (int n2=0;n2<nequations;n2++)
        				{
                				for (int m2=0;m2<=nequations;m2++)
                				{
                        				cout <<coefficient[m2][n2]<<"\t";
                				}
                				cout <<endl;
        				}
					cout <<endl;
				}
                        }
                }
        }
	if (verbose)
	{
        	for (n=0;n<nequations;n++)
        	{
                	for (m=0;m<=nequations;m++)
                	{
                        	cout <<coefficient[m][n]<<"\t";
                	}
                	cout <<endl;
        	}
		cout <<endl;
	}

        for (n=0;n<nequations;n++)
        {
                variable[n]=coefficient[nequations][n]/coefficient[n][n];
		if (verbose) cout <<"Variable "<<n+1<<" is "<<variable[n]<<"\n";
        }
}

void GaussJordanElimination(long double **coefficient, long double **InverseMatrix, int nequations)
{
	//Finds the inverse matrix using Gauss-Jordan Elimination
        bool verbose;
	int i, m, n;
        long double coefficientnn;
        long double scale;

	verbose=false;

	if (verbose)
	{
		cout <<"nequations= "<<nequations<<endl;
		for (n=0;n<nequations;n++)
		{
			for (m=0;m<=nequations;m++)
			{
				cout <<coefficient[m][n]<<"\t";
			}
			cout <<endl;
		}
		cout <<endl;
	}

	for (m=0;m<nequations;m++)
	{
		for (n=0;n<nequations;n++)
		{
			if (m==n) InverseMatrix[m][n]=1.0;
			else InverseMatrix[m][n]=0.0;
		}
	}

	for (n=0;n<nequations-1;n++)
        {
                coefficientnn=coefficient[n][n];
                if (coefficientnn!=0)
                {
                        for (m=n+1;m<nequations;m++)
                        {
                                scale=coefficient[n][m]/coefficientnn;
                                for (i=0;i<nequations;i++)
                                {
                                        coefficient[i][m]=coefficient[i][m]-coefficient[i][n]*scale;
					InverseMatrix[i][m]=InverseMatrix[i][m]-InverseMatrix[i][n]*scale;
                                }
                        }
                }
        }
	
	if (verbose)
	{
        	for (n=0;n<nequations;n++)
        	{
                	for (m=0;m<=nequations;m++)
                	{
                        	cout <<coefficient[m][n]<<"\t";
                	}
                	cout <<endl;
        	}
		cout <<endl;
	}

        for (n=nequations-1;n>0;n=n-1)
        {
                coefficientnn=coefficient[n][n];
                if (coefficientnn!=0)
                {
			for (m=n-1;m>-1;m=m-1)
                        {
                                scale=coefficient[n][m]/coefficientnn;
				//cout <<"coefficient["<<n<<"]["<<m<<"]= "<<coefficient[n][m]<<" coefficientnn= "<<coefficientnn<<" scale= "<<scale<<endl;
                                for (i=n;i<nequations;i++)
                                {
                                        coefficient[i][m]=coefficient[i][m]-coefficient[i][n]*scale;
					InverseMatrix[i][m]=InverseMatrix[i][m]-InverseMatrix[i][n]*scale;
                                }

				if (verbose)
				{
        				for (int n2=0;n2<nequations;n2++)
        				{
                				for (int m2=0;m2<=nequations;m2++)
                				{
                        				cout <<coefficient[m2][n2]<<"\t";
                				}
                				cout <<endl;
        				}
					cout <<endl;
				}
                        }
                }
        }
	if (verbose)
	{
        	for (n=0;n<nequations;n++)
        	{
                	for (m=0;m<=nequations;m++)
                	{
                        	cout <<coefficient[m][n]<<"\t";
                	}
                	cout <<endl;
        	}
		cout <<endl;
	}
}

template <class Matrix_T>
void GaussJordanElimination(Matrix_T coefficient, Matrix_T &InverseMatrix)
{
	//Finds the inverse matrix using Gauss-Jordan Elimination
        bool verbose;
	int i, m, n;
	int nequations=coefficient.size();
        long double coefficientnn;
        long double scale;

	verbose=false;

	if (verbose)
	{
		cout <<"nequations= "<<nequations<<endl;
		PrintMatrix(coefficient);
	}

	for (m=0;m<nequations;m++)
	{
		for (n=0;n<nequations;n++)
		{
			if (m==n) InverseMatrix[m][n]=1.0;
			else InverseMatrix[m][n]=0.0;
		}
	}
	
	for (n=0;n<nequations-1;n++)
        {
                coefficientnn=coefficient[n][n];
                if (coefficientnn!=0)
                {
                        for (m=n+1;m<nequations;m++)
                        {
                                scale=coefficient[m][n]/coefficientnn;
                                for (i=0;i<nequations;i++)
                                {
                                        coefficient[m][i]-=coefficient[n][i]*scale;
					InverseMatrix[m][i]-=InverseMatrix[n][i]*scale;
                                }
                        }
                }
        }
	
	if (verbose) PrintMatrix(coefficient);

        for (n=nequations-1;n>0;n=n-1)
        {
                coefficientnn=coefficient[n][n];
                if (coefficientnn!=0)
                {
			for (m=n-1;m>-1;m=m-1)
                        {
                                scale=coefficient[m][n]/coefficientnn;
				//cout <<"coefficient["<<n<<"]["<<m<<"]= "<<coefficient[n][m]<<" coefficientnn= "<<coefficientnn<<" scale= "<<scale<<endl;
                                for (i=0;i<nequations;i++)
                                {
                                        coefficient[m][i]-=coefficient[n][i]*scale;
					InverseMatrix[m][i]-=InverseMatrix[n][i]*scale;
                                }

				if (verbose) PrintMatrix(coefficient);
                        }
                }
        }
	if (verbose) PrintMatrix(coefficient);
	for (n=0;n<nequations;n++)
	{
		for (m=0;m<nequations;m++) InverseMatrix[n][m]/=coefficient[n][n];
	}
}

long double GetMatrixElement(long double **Matrix1, long double **Matrix2, int m, int n, int Size)
{
	int i;
	Real sum=0;
	for (i=0;i<Size;i++) sum+=Matrix1[m][i]*Matrix2[i][n];
	return sum;
}

void MatrixMultiply(long double **Matrix1, long double **Matrix2, long double **MatrixOut, int Size1, int Size2)
{
	int m, n;

	for (m=0;m<Size1;m++)
	{
		for (n=0;n<Size2;n++)
		{
			MatrixOut[m][n]=GetMatrixElement(Matrix1, Matrix2, m, n, Size1);
		}
	}
}

template<class T>
T GetMatrixElement(vector< vector<T> > Matrix1, vector< vector<T> > Matrix2, int m, int n)
{
	int i, Size=Matrix2.size();
	Real sum=0;
	for (i=0;i<Size;i++) sum+=Matrix1[m][i]*Matrix2[i][n];
	
	return sum;
}

template<class Matrix_T>
void MatrixMultiply(Matrix_T Matrix1, Matrix_T Matrix2, Matrix_T &MatrixOut)
{
	int m, n, Size1, Size2;
	Size1=Matrix1.size();
	Size2=Matrix1[0].size();
	for (m=0;m<Size1;m++)
	{
		for (n=0;n<Size2;n++)
		{
			MatrixOut[m][n]=GetMatrixElement(Matrix1, Matrix2, m, n);
		}
	}
}

template<class Matrix_T>
void Transpose(Matrix_T Matrix, Matrix_T &MatrixTranspose)
{
	int i, j;
	int Size1, Size2;

	Size1=Matrix.size();
	Size2=Matrix[0].size();

	if (MatrixTranspose.size()!=Size2 || MatrixTranspose[0].size()!=Size1)
	{
		DeleteMatrix(MatrixTranspose);
		CreateMatrix(MatrixTranspose, Size2, Size1);
	}

	for (i=0;i<Size1;i++)
	{
		for (j=0;j<Size2;j++) MatrixTranspose[j][i]=Matrix[i][j];
	}
}

template<class Matrix_T>
void Transpose(Matrix_T &Matrix)
{
	int i, j;
	int Size1, Size2;
	Matrix_T MatrixTranspose;

	Size1=Matrix.size();
	Size2=Matrix[0].size();

	CreateMatrix(MatrixTranspose, Size2, Size1);

	for (i=0;i<Size1;i++)
	{
		for (j=0;j<Size2;j++) MatrixTranspose[j][i]=Matrix[i][j];
	}
	CopyVector(MatrixTranspose, Matrix);
}

void swap(Real &x, Real &y)
{
	Real temp;

	temp=x;
	x=y;
	y=temp;
}

template <class Matrix_T>
void Transpose(Matrix_T &Matrix, int Size1, int Size2)
{
	//Only works for square matrix.

	for (int i=0;i<Size1;i++)
	{
		for (int j=i;j<Size2;j++) 
		{
			swap(Matrix[i][j], Matrix[j][i]);
		}
	}
}

template <class Matrix_T>
void Transpose2(Matrix_T &Matrix, int Size1, int Size2)
{
	//Only works for square matrix.

	for (int i=0;i<Size1;i++)
	{
		for (int j=i;j<Size2;j++) 
		{
			for (int k=0;k<Size1;k++)
			{
				swap(Matrix[i][k][j], Matrix[j][k][i]);
			}
		}
	}
}

template<class T>
T calcSum(vector<T> &x)
{
        int Size=x.size();
        Real sum=0;

        for (int i=0;i<Size;i++) sum+=x[i];

        return sum;
}

template<class T>
T calcSumSqr(vector<T> &x)
{
        int Size=x.size();
        Real sum=0;

        for (int i=0;i<Size;i++) sum+=x[i]*x[i];

        return sum;
}

template<class T>
T calcAverage(vector<T> &x)
{
	int i, Size=x.size();
	Real sum=0;
        if (Size==0) return 0;
	for (i=0;i<Size;i++) sum+=x[i];

	return sum/Real(Size);
}

template<class T>
T calcAverage(vector<T> &x, int first, int last)
{
	int i, Size=x.size();
	Real sum=0;
        if (Size==0) return 0;
	for (i=first;i<=last;i++) sum+=x[i];

	return sum/Real(last-first+1);
}

template<class T>
T calcVariance(vector<T> &x)
{
        int Size=x.size();
        Real sum, sum_sqr;
        sum=0;
        sum_sqr=0;
        for (int i=0;i<Size;i++)
        {
                sum+=x[i];
                sum_sqr+=x[i]*x[i];
        }
        cout <<"Size= "<<Size<<endl;
        cout <<"sum_sqr= "<<sum_sqr<<" sum= "<<sum<<endl;
        return (sum_sqr/Real(Size)-sum*sum/(Real(Size)*Real(Size)));
}

template<class T>
T StandardDeviation(vector<T> &x)
{
        return sqrt(calcVariance(x));
}

template<class T>
T calcError(vector<T> &x)
{
        return StandardDeviation(x)/sqrt(Real(x.size()));
}

Real variance(vector< vector<Real> > &x)
{
	int i, j, Size1, Size2;
	Real sum, sum_sqr, NumElements;
	Size1=x.size();
	Size2=x[0].size();
	sum=0;
	sum_sqr=0;
	for (i=0;i<Size1;i++)
	{
		for (j=0;j<Size2;j++)
		{
			sum+=x[i][j];
			sum_sqr+=x[i][j]*x[i][j];
		}
	}

	NumElements=Real(Size1*Size2);
	return (sum_sqr/NumElements-sum*sum/(NumElements*NumElements));
}

template <class T>
Real covariance(vector<T> &x, vector<T> &y)
{
	int n, Size=x.size();
	Real sum=0;
	Real AverageX, AverageY;

	AverageX=average(x);
	AverageY=average(y);

	for (n=0;n<Size;n++) sum+=(x[n]-AverageX)*(y[n]-AverageY);

	return sum;
}

template<class Matrix_T>
void CovarianceMatrix(Matrix_T Matrix, Matrix_T &CovarianceMatrix)
{
	int i, j;
	int Size1, Size2;

	Size1=Matrix.size();
	Size2=Matrix[0].size();

	if (CovarianceMatrix.size()!=Size1 || CovarianceMatrix[0].size()!=Size1)
	{
		DeleteMatrix(CovarianceMatrix);
		CreateMatrix(CovarianceMatrix, Size1, Size1);
	}

	for (i=0;i<Size1;i++)
	{
		for (j=i;j<Size1;j++)
		{
			CovarianceMatrix[i][j]=covariance(Matrix[j], Matrix[i]);
			CovarianceMatrix[j][i]=CovarianceMatrix[i][j];
		}
	}
}

void normalize(vector<Real> &v)
{
        int Size=v.size();
        Real sum=calcSum(v);
        for (int i=0;i<Size;i++) v[i]/=sum;
}

void setAverageTo(vector<Real> &v, Real x)
{
        int Size=v.size();
        Real sum=calcSum(v);
        for (int i=0;i<Size;i++) v[i]=v[i]*x/sum;
}

Real getMin(Real x1, Real x2)
{
        if (x1<=x2) return x1;
        else return x2;
}

Real getMax(Real x1, Real x2)
{
        if (x1<=x2) return x2;
        else return x1;
}

int getMax(int x1, int x2)
{
        if (x1<=x2) return x2;
        else return x1;
}

Real Monotonicity(vector<Real> x, vector<Real> y)
{
	//Calculates the probability that if y decreases x decreases.
	int m, n, Size;
	Real monatonicity;
	Size=x.size();
	monatonicity=0;
	for (m=0;m<Size-1;m++)
	{
		for (n=m+1;n<Size;n++)
		{
			if (y[m]<y[n] && x[m]<y[n]) monatonicity+=1.0;
		}
	}
	monatonicity=monatonicity*2.0/Real((Size-1)*Size);
	return monatonicity;
}
/*
template <class Vector_T>
bool IsSorted(Vector_T &v)
{
	int n;
	int Size=v.size();
	if (Size>1)
	{
		for (n=1;n<Size;n++) 
		{
			if (v[n]<v[n-1]) return false;
		}
	}

	return true;

}

template <class Vector_T>
void Sort(Vector_T &v)
{
	bool Sorted1, Sorted2;
	int n;
	int Size, Size1, Size2;
	Real first, last;
	vector<Real> v1, v2;

	Size=v.size();
	first=v[0];
	for (n=0;n<Size;n++)
	{
		if (v[n]<first)
		{
			first=v[n];
			break;
		}
		else if (v[n]>first) break;
	}

	for (n=0;n<Size;n++)
	{
		if (v[n]<=first) v1.push_back(v[n]);
		else v2.push_back(v[n]);
	}

	Size1=v1.size();
	Size2=v2.size();
	//for (n=0;n<Size1;n++) cout <<"v1["<<n<<"]= "<<v1[n]<<endl;
	//cout <<endl;
	//for (n=0;n<Size2;n++) cout <<"v2["<<n<<"]= "<<v2[n]<<endl;
	//cout <<endl;
	Sorted1=IsSorted(v1);
	Sorted2=IsSorted(v2);
	if (Size1>1 && !Sorted1) Sort(v1);
	if (Size2>1 && !Sorted2) Sort(v2);

	for (n=0;n<Size;n++) v.pop_back();

	for (n=0;n<Size1;n++)
	{
		v.push_back(v1[n]);
	}

	for (n=0;n<Size2;n++)
	{
		v.push_back(v2[n]);
	}

}
*/

template <class Vector_T>
bool IsSorted(Vector_T &v)
{
	int n;
	int Size=v.size();
	if (Size>1)
	{
		for (n=1;n<Size;n++) 
		{
			if (v[n]<v[n-1]) return false;
		}
	}

	return true;

}

template <class T>
void Sort(vector<T> &v)
{
	bool Sorted1, Sorted2;
	int n;
	int Size, Size1, Size2;
	T first, last;
	vector<T> v1, v2;

	Size=v.size();
	first=v[0];
	for (n=0;n<Size;n++)
	{
		if (v[n]<first)
		{
			first=v[n];
			break;
		}
		else if (v[n]>first) break;
	}

	for (n=0;n<Size;n++)
	{
		if (v[n]<=first) SafePushBack(v1, v[n], "v1");
		else SafePushBack(v2, v[n], "v2");
	}

	Size1=v1.size();
	Size2=v2.size();
	//for (n=0;n<Size1;n++) cout <<"v1["<<n<<"]= "<<v1[n]<<endl;
	//cout <<endl;
	//for (n=0;n<Size2;n++) cout <<"v2["<<n<<"]= "<<v2[n]<<endl;
	//cout <<endl;
	Sorted1=IsSorted(v1);
	Sorted2=IsSorted(v2);
	if (Size1>1 && !Sorted1) Sort(v1);
	if (Size2>1 && !Sorted2) Sort(v2);

        v.clear();

	for (n=0;n<Size1;n++)
	{
                SafePushBack(v, v1[n], "v");
	}

	for (n=0;n<Size2;n++)
	{
                SafePushBack(v, v2[n], "v");
	}
}

Real calcGaussian(Real x, Real stdDev)
{
        return exp(-x*x/stdDev)*sqrt(abs(stdDev)/pi);
}

Real calcGaussian(Real x, Real mean, Real stdDev)
{
        return calcGaussian(x-mean, stdDev);
}

int RandInt(int Max)
{
	int Rand;
	Rand=int(Real(rand())*Real(Max)/Real(RAND_MAX)+0.5);
	return Rand;
}

int RandInt(int Min, int Max)
{
        return RandInt(Max-Min)+Min;
}

Real RandDouble(Real Max)
{
	Real Rand;
	Rand=2.0*Max*Real(rand())/Real(RAND_MAX)-Max;
	return Rand;
}

Real RandDouble(Real Min, Real Max)
{
        return (Max-Min)*Real(rand())/Real(RAND_MAX)+Min;
}

Real RandGaussian(Real stdDev)
{
        if (isnan(stdDev)) return RandDouble(-1.0, 1.0);
        Real max=4.0*stdDev;
        Real gaussian, randNum, randNum2;
        randNum=RandDouble(-max, max);
        gaussian=calcGaussian(randNum, stdDev);
        randNum2=RandDouble(0, 1.0);
        if (gaussian>randNum2) return randNum;
        else return RandGaussian(stdDev);
}

Real RandGaussian(Real mean, Real stdDev)
{
        return RandGaussian(stdDev)+mean;
}

template <typename T>
void LinearRegression(vector<T> x, vector<T> y, T &slope, T &yintercept)
{
	int m, n, Size;
	vector<long double> variable;
	vector< vector<long double> > coefficient;
	//cout <<"In FitLine"<<endl;
	CreateVector(variable, 2);
	CreateMatrix(coefficient, 3, 2);
	//cout <<"Initialized coefficeint"<<endl;
	Size=x.size();
	coefficient[0][0]=(long double)(Size);
	for (n=0;n<Size;n++) coefficient[1][0]+=x[n];
	for (n=0;n<Size;n++) coefficient[2][0]+=y[n];
	coefficient[0][1]=coefficient[1][0];
	for (n=0;n<Size;n++) coefficient[1][1]+=x[n]*x[n];
	for (n=0;n<Size;n++) coefficient[2][1]+=x[n]*y[n];
	//cout <<"Assigned coefficients"<<endl;
	Equation(coefficient, variable);
	//cout <<"Solved equation"<<endl;
	yintercept=variable[0];
	slope=variable[1];
}

Real calcAbsSum(vector<Real> &v)
{
        int Size=v.size();
        Real sum=0;

        for (int i=0;i<Size;i++)
        {
                sum+=abs(v[i]);
        }
        return sum;
}

void abs(vector<Real> &v)
{
        int Size=v.size();
        for (int i=0;i<Size;i++)
        {
                v[i]=abs(v[i]);
        }
}

Real calcAverage(vector<Real> &x, vector<Real> &y)
{
        int Size=x.size();
        Real sum=0;

        for (int i=0;i<Size;i++)
        {
                sum+=x[i]*abs(y[i]);
        }
        return sum/calcAbsSum(y);
}

void calcStdDevAndAverage(vector<Real> &x, vector<Real> &y, Real &ave, Real &stdDev)
{
        int Size=x.size();
        Real sum=0, sumSqr=0, sumY;
        Real aveSqr, sqrAve;

        for (int i=0;i<Size;i++)
        {
                sumSqr+=x[i]*x[i]*abs(y[i]);
                sum=x[i]*abs(y[i]);
        }
        sumY=calcAbsSum(y);
        ave=sum/sumY;
        aveSqr=sumSqr/sumY;
        sqrAve=ave*ave;
        stdDev=sqrt(aveSqr-sqrAve);
}

Real calcScaleFactor(vector<Real> &y1, vector<Real> &y2)
{
        //Minimize sum((scale*y1-y2)^2)
        int Size=y1.size();
        Real sumCrossTerm=0, sumSqr=0;

        for (int i=0;i<Size;i++)
        {
                sumCrossTerm+=y1[i]*y2[i];
                sumSqr+=y1[i]*y1[i];
                cout <<"y1= "<<y1[i]<<" y2= "<<y2[i]<<endl;
                cout <<"sumCrossTerm= "<<sumCrossTerm<<" sumSqr= "<<sumSqr<<endl;
        }
        return sumCrossTerm/sumSqr;
}

void calcGaussian(vector<Real> &x, vector<Real> &gaussian, Real mean, Real width)
{
        int Size=x.size();
        cout <<"in calcGaussian vector"<<endl;
        if (gaussian.size()==0) SafeAlloc(gaussian, Size, "gaussian");
        for (int i=0;i<Size;i++)
        {
                gaussian[i]=calcGaussian(x[i], mean, width);
                //cout <<"gaussian["<<i<<"]= "<<gaussian[i]<<" x= "<<x[i]
                //<<" mean= "<<mean<<" width= "<<width<<endl;
        }
}

Real calcScaleForGaussian(vector<Real> &x, vector<Real> &y, Real mean, Real width)
{
        vector<Real> gaussian;
        calcGaussian(x, gaussian, mean, width); 
        return calcScaleFactor(gaussian, y);
}

Real calcSqrDiff(vector<Real> &y1, vector<Real> &y2)
{
        int Size=y1.size();
        Real sum=0;
        for (int i=0;i<Size;i++)
        {
                sum+=(y2[i]-y1[i])*(y2[i]-y1[i]);
                cout <<"sum= "<<sum<<" y1["<<i<<"]= "<<y1[i]<<" y2= "<<y2[i]<<endl;
        }
        return sum;
}

void vectorMultiply(Real scale, vector<Real> &v)
{
        int Size=v.size();
        for (int i=0;i<Size;i++)
        {
                v[i]*=scale;
        }
}

Real gaussianFitScore(vector<Real> &x, vector<Real> &y, Real mean, Real width)
{
        Real scale;
        vector<Real> gaussian;
        cout <<"in gaussianFitScore"<<endl;
        calcGaussian(x, gaussian, mean, width);
        scale=calcScaleFactor(gaussian, y);
        cout <<"scale= "<<scale<<endl;
        vectorMultiply(scale, gaussian);
        return calcSqrDiff(gaussian, y);
}

Real gaussianFitScore(vector<Real> coefficient, xyStruct xy)
{
        Real fitScore=gaussianFitScore(xy.x, xy.y, coefficient[0], coefficient[1]);
        cout <<"fitScore= "<<fitScore<<endl;
        if (isnan(fitScore)) exit(EXIT_FAILURE);
        return fitScore;
}

void fitGaussian(vector<Real> &x, vector<Real> &y, Real &scale, Real &mean, Real &width)
{
        //y=scale*exp(-(x-mean)^2/width^0.5)
        Real dScoreCriteria=0.0001;
        vector<Real> coefficients;
        xyStruct xy;
        Real (*func_to_minimize)(vector<Real>, xyStruct)=gaussianFitScore;
        coefficients.resize(2);
        xy.x=x;
        xy.y=y;
        calcStdDevAndAverage(x, y, mean, width);
        coefficients[0]=mean;
        coefficients[1]=width;
        cout <<" mean= "<<mean<<" width= "<<width<<endl;
        nonlinearConjugateGradient(coefficients, xy, func_to_minimize, dScoreCriteria);
        mean=coefficients[0];
        width=coefficients[1];
}

void fitGaussian(vector<Real> &x, vector<Real> &y, Real &mean, Real &width)
{
        Real scale;
        fitGaussian(x, y, scale, mean, width);
}

vector<Real> calcLog(vector<Real> &v)
{
        int Size=v.size();
        vector<Real> logv;
        SafeAlloc(logv, Size, "logv");
        for (int i=0;i<Size;i++)
        {
                logv[i]=log(v[i]);
        }
        return logv;
}

Real calcPower(Real x, int degree)
{
        Real power=1;
        for (int i=0;i<degree;i++)
        {
                power*=x;
        }
        return power;
}

vector<Real> calcPower(vector<Real> &x, int degree)
{
        int Size=x.size();
        vector<Real> power;
        SafeAlloc(power, Size, "power");
        for (int i=0;i<Size;i++)
        {
                power[i]=calcPower(x[i], degree);
        }
        return power;
}

vector<Real> fitCurves(Matrix &curves, vector<Real> &y)
{
        int Size1, Size2, npoint=y.size();
        vector<Real> variable;
        Matrix coefficients;

        Get2DVectorSize(curves, Size1, Size2);
        SafeAlloc(variable, Size2, "variables in fitCurves");
        Safe2DAlloc(coefficients, Size2+1, Size2, "in fitCurves");
        for (int i=0;i<npoint;i++)
        {
                for (int j=0;j<Size2;j++)
                {
                        for (int k=j;k<Size2;k++)
                        {
                                coefficients[k][j]+=curves[i][j]*curves[i][k];
                                coefficients[j][k]=coefficients[k][j];
                        }
                        coefficients[Size2][j]+=y[i]*curves[i][j];
                }
        }
        Equation(coefficients, variable);
        return variable;
}

vector<Real> fitPolynomial(vector<Real> &x, vector<Real> &y, int degree)
{
        vector<Real> powers, coefficients;
        Matrix polynomials;
        for (int i=0;i<degree+1;i++)
        {
                powers=calcPower(x, i);
                SafePushBack(polynomials, powers, "polynomials");
        }
        Transpose(polynomials);
        coefficients=fitCurves(polynomials, y);
        return coefficients;
}

void fitGaussian2(vector<Real> &x, vector<Real> &y, Real &scale, Real &mean, Real &width)
{
        int degree=2;
        vector<Real> logy, coefficients;

        logy=calcLog(y);
        coefficients=fitPolynomial(x, logy, degree);
        width=-coefficients[2];
        mean=-coefficients[1]*0.5/coefficients[2];
        scale=exp(coefficients[0]-mean*mean*width);
        width=1.0/width;
}

void fitGaussian2(vector<Real> &x, vector<Real> &y, Real &mean, Real &width)
{
        Real scale=0;
        fitGaussian2(x, y, scale, mean, width);
}


void FFT(Real InputReal[], Real InputImag[], Real OutputReal[], Real OutputImag[], int max)
{
	int n, m;
	int halfmax;
	Real DoubleHalfMax, CosReal, CosImag, SinReal, SinImag;
	Real *EvenInputReal, *EvenInputImag, *OddInputReal, *OddInputImag;
	Real *EvenOutputReal, *EvenOutputImag, *OddOutputReal, *OddOutputImag;

	DoubleHalfMax=Real(max)*0.5;
	halfmax=int(floor(DoubleHalfMax));

	EvenInputReal=new Real[halfmax];
	EvenInputImag=new Real[halfmax];
	OddInputReal=new Real[halfmax];
	OddInputImag=new Real[halfmax];

	EvenOutputReal=new Real[halfmax];
	EvenOutputImag=new Real[halfmax];
	OddOutputReal=new Real[halfmax];
	OddOutputImag=new Real[halfmax];

	for (int i=0;i<halfmax;i++)
	{
		EvenInputReal[i]=0;
		EvenInputImag[i]=0;
		OddInputReal[i]=0;
		OddInputImag[i]=0;
		EvenOutputReal[i]=0;
		EvenOutputImag[i]=0;
		OddOutputReal[i]=0;
		OddOutputImag[i]=0;
	}

	n=0;
	m=0;

	if (max>1)
	{
		while (n<max)
		{
			EvenInputReal[m]=InputReal[n];
			EvenInputImag[m]=InputImag[n];
			n++;
			OddInputReal[m]=InputReal[n];
			OddInputImag[m]=InputImag[n];
			n++;
			m++;
		}

		FFT(EvenInputReal, EvenInputImag, EvenOutputReal, EvenOutputImag, halfmax);
		FFT(OddInputReal, OddInputImag, OddOutputReal, OddOutputImag, halfmax);

		for (int n=0;n<halfmax;n++)
                {
                        CosReal=cos(2.0*pi*Real(n)/Real(max))*OddOutputReal[n];
                        SinImag=sin(2.0*pi*Real(n)/Real(max))*OddOutputImag[n];
                        CosImag=cos(2.0*pi*Real(n)/Real(max))*OddOutputImag[n];
                        SinReal=sin(2.0*pi*Real(n)/Real(max))*OddOutputReal[n];
			OutputReal[n]=EvenOutputReal[n]+CosReal+SinImag;
			OutputImag[n]=EvenOutputImag[n]+CosImag-SinReal;
		        OutputReal[n+halfmax]=EvenOutputReal[n]-CosReal-SinImag;
                        OutputImag[n+halfmax]=EvenOutputImag[n]-CosImag+SinReal;
                }
	}
	else
	{
		OutputReal[0]=InputReal[0];
		OutputImag[0]=InputImag[0];
	}

	delete [] EvenInputReal;
	delete [] EvenInputImag;
	delete [] OddInputReal;
	delete [] OddInputImag;

	delete [] EvenOutputReal;
	delete [] EvenOutputImag;
	delete [] OddOutputReal;
	delete [] OddOutputImag;
}

void PrintArray(Real **output, int MAX)
{
	for (int i=0;i<MAX;i++)
	{
		for (int j=0;j<MAX;j++)
		{
			cout <<setw(15)<<output[i][j];
		}
		cout <<endl;
	}
}

void PrintTensor(Real ***output, int MAX)
{
	for (int i=0;i<MAX;i++)
	{
		PrintArray(output[i], MAX);
		cout <<endl;
	}
}

void FFT2D(Real **InputReal, Real **InputImag, Real **OutputReal, Real **OutputImag, int max)
{
	//Only works for square matrix.
	Real **OutputRealTemp, **OutputImagTemp;

	OutputRealTemp=new Real *[max];
	OutputImagTemp=new Real *[max];
	for (int i=0;i<max;i++)
	{
		*(OutputRealTemp+i)=new Real[max];
		*(OutputImagTemp+i)=new Real[max];
	}

	for (int i=0;i<max;i++)
	{
		FFT(InputReal[i], InputImag[i], OutputRealTemp[i], OutputImagTemp[i], max);
	}
	Transpose(OutputRealTemp, max, max);
	Transpose(OutputImagTemp, max, max);	
	for (int i=0;i<max;i++)
	{
		FFT(OutputRealTemp[i], OutputImagTemp[i], OutputReal[i], OutputImag[i], max);
	}
	Transpose(OutputReal, max, max);
	Transpose(OutputImag, max, max);
	for (int i=0;i<max;i++)
	{
		delete [] OutputRealTemp[i];
		delete [] OutputImagTemp[i];
	}
	delete [] OutputRealTemp;
	delete [] OutputImagTemp;
}

void FFT3D(Real ***InputReal, Real ***InputImag, Real ***OutputReal, Real ***OutputImag, int max)
{
	//Only works for cubic input.
	Real ***OutputRealTemp, ***OutputImagTemp;

	OutputRealTemp=new Real **[max];
	OutputImagTemp=new Real **[max];
	for (int i=0;i<max;i++)
	{
		*(OutputRealTemp+i)=new Real *[max];
		*(OutputImagTemp+i)=new Real *[max];
		for (int j=0;j<max;j++)
		{
			*(*(OutputRealTemp+i)+j)=new Real[max];
			*(*(OutputImagTemp+i)+j)=new Real[max];
		}
	}

	for (int i=0;i<max;i++)
	{
		FFT2D(InputReal[i], InputImag[i], OutputRealTemp[i], OutputImagTemp[i], max);
	}
	//PrintTensor(OutputRealTemp, max);
	Transpose(OutputRealTemp, max, max);
	//PrintTensor(OutputRealTemp, max);
	Transpose(OutputImagTemp, max, max);	
	for (int i=0;i<max;i++)
	{
		Transpose(OutputRealTemp[i], max, max);
		Transpose(OutputImagTemp[i], max, max);
	}
	for (int i=0;i<max;i++)
	{
		for (int j=0;j<max;j++)
		{
			FFT(OutputRealTemp[i][j], OutputImagTemp[i][j], OutputReal[i][j], OutputImag[i][j], max);
		}
	}

	for (int i=0;i<max;i++)
	{
		Transpose(OutputRealTemp[i], max, max);
		Transpose(OutputImagTemp[i], max, max);
	}
	Transpose(OutputReal, max, max);
	Transpose(OutputImag, max, max);
	Transpose2(OutputReal, max, max);
	Transpose2(OutputImag, max, max);
	//Transpose(OutputReal, max, max);
	//Transpose(OutputImag, max, max);
	for (int i=0;i<max;i++)
	{
		for (int j=0;j<max;j++)
		{
			delete [] OutputRealTemp[i][j];
			delete [] OutputImagTemp[i][j];
		}
		delete [] OutputRealTemp[i];
		delete [] OutputImagTemp[i];
	}
	delete [] OutputRealTemp;
	delete [] OutputImagTemp;
}

void CalcFTElement(Real ***InputReal, Real ***InputImag, Real &OutputReal, Real &OutputImag, int i, int j, int k, int max)
{
	Real theta;
	for (int a=0;a<max;a++)
	{
		for (int b=0;b<max;b++)
		{
			for (int c=0;c<max;c++)
			{
				theta=2.0*pi*Real(i*a+j*b+k*c)/Real(max);
				OutputReal+=InputReal[a][b][c]*cos(theta)-InputImag[a][b][c]*sin(theta);
				OutputImag+=InputImag[a][b][c]*cos(theta)+InputReal[a][b][c]*sin(theta);
			}
		}
	}
}

void FT3D(Real ***InputReal, Real ***InputImag, Real ***OutputReal, Real ***OutputImag, int max)
{
	for (int i=0;i<max;i++)
	{
		for (int j=0;j<max;j++)
		{
			for (int k=0;k<max;k++)
			{
				OutputReal[i][j][k]=0;
				OutputImag[i][j][k]=0;
				CalcFTElement(InputReal, InputImag, OutputReal[i][j][k], OutputImag[i][j][k], i, j, k, max);
			}
		}
	}
}

void FitLine(vector<Real> &x, vector<Real> &y, Real &a, Real &b)
{
        //There may be something wrong with function.
        int Size;
        vector<Real> variable;
        vector< vector<Real> > coefficient;

        SafeAlloc(variable, 2, "variable in FitLine");
        Safe2DAlloc(coefficient, 3, 2, "coefficient in FitLine");
        Size=x.size();
        coefficient[0][0]=Real(Size);
        for (int i=0;i<Size;i++)
        {
                coefficient[1][0]+=x[i];
                coefficient[2][0]+=y[i];
                coefficient[1][1]+=x[i]*x[i];
                coefficient[2][1]+=x[i]*y[i];
        }
        coefficient[0][1]=coefficient[1][0];
        Equation(coefficient, variable);

        a=variable[0];
        b=variable[1];
}

//template <typename T>
void MultiVariableRegression(vector<Real> &z, vector< vector<Real> > &data, vector<Real> &variable)
{
        int nvariable=variable.size();
        int npoint=z.size();
        vector< vector<Real> > coefficient;

        Safe2DAlloc(coefficient, nvariable+1, nvariable, "in MultiVariableRegression");
        //Print2DVectorSize(coefficient, "coefficeint");
        //Print2DVectorSize(data, "data");
        for (int i=0;i<nvariable-1;i++)
        {
                for (int j=0;j<nvariable-1;j++)
                {
                        for (int k=0;k<npoint;k++)
                        {
                                coefficient[i][j]+=data[i][k]*data[j][k];
                        }
                }
        }
        //PrintMatrix(coefficient);
        //cout <<endl;
        for (int i=0;i<nvariable-1;i++)
        {
                for (int k=0;k<npoint;k++)
                {
                        coefficient[i][nvariable-1]+=data[i][k];
                }
        }
        //PrintMatrix(coefficient);
        //cout <<endl;
        for (int i=0;i<nvariable;i++)
        {
                for (int j=0;j<nvariable;j++)
                {
                        coefficient[j][i]=coefficient[i][j];
                }
        }
        coefficient[nvariable-1][nvariable-1]=Real(npoint);
        //PrintMatrix(coefficient);
        //cout <<endl;
        
        for (int i=0;i<nvariable-1;i++)
        {
                for (int j=0;j<npoint;j++)
                {
                        coefficient[nvariable][i]+=z[j]*data[i][j];
                }
        }
        for (int j=0;j<npoint;j++)
        {
                coefficient[nvariable][nvariable-1]+=z[j];
        }
        //PrintMatrix(coefficient);
        //cout <<endl;
        Equation(coefficient, variable);
}

Real Sxy(vector<Real> &x, vector<Real> &y, Real a, Real b)
{
        int Size;
        Real diff, sum, sxy;

        Size=x.size();

        sum=0;

        for (int i=0;i<Size;i++)
        {
                diff=(y[i]-a-b*x[i]);
                sum+=diff*diff;
        }
        sxy=sum/Real(Size-1);
        return sxy;
}

Real Sy(vector<Real> &x, vector<Real> &y)
{
        int Size;
        Real AverageY, sum, diff, sy;

        Size=x.size();
        AverageY=calcAverage(y);

        sum=0;
        for (int i=0;i<Size;i++)
        {
                diff=y[i]-AverageY;
                sum+=diff*diff;
        }

        sy=sum/Real(Size-1);
        return sy;
}

Real correlation(vector<Real> &x, vector<Real> &y)
{
        Real a, b;
        Real sxy, sy, r2;

        a=0;
        b=0;

        //FitLine(x, y, a, b);
        LinearRegression(x, y, a, b);
        sxy=Sxy(x, y, a, b);
        sy=Sy(x, y);
        r2=1.0-sxy/sy;
        if (b<0) r2=-r2;        
        return r2;
}

//template <typename T>
Real Sxy(vector<Real> &z, vector< vector<Real> > &data, vector<Real> &variable)
{
        int Size;
        int nvariable=variable.size();
        Real diff, sum, sxy;

        Size=z.size();

        sum=0;

        for (int i=0;i<Size;i++)
        {
                diff=z[i]-variable[nvariable-1];
                for (int j=0;j<nvariable-1;j++)
                {
                        diff-=variable[j]*data[j][i];
                }
                sum+=diff*diff;
        }
        sxy=sum/Real(Size-1);
        return sxy;
}

//template <typename T>
Real Sy(vector<Real> &y)
{
        int Size;
        Real AverageY, sum, diff, sy;

        Size=y.size();
        AverageY=calcAverage(y);

        sum=0;
        for (int i=0;i<Size;i++)
        {
                diff=y[i]-AverageY;
                sum+=diff*diff;
        }

        sy=sum/Real(Size-1);
        return sy;
}

Real CalcCorrelationFromFit(vector<Real> &z, vector< vector<Real> > &data, vector<Real> &variable)
{
        Real sxy, sy, r2;
        sxy=Sxy(z, data, variable);
        sy=Sy(z);
        r2=1.0-sxy/sy;

        return r2;
}

Real linearRegression(vector<Real> &x, vector<Real> &y, Real &slope, Real &yintercept, Real &slopeError, Real &yError)
{
        Real r, sx, sxx, sy, syy, se2, sa2, sb2;
        Real Size=Real(x.size());
        LinearRegression(x, y, slope, yintercept);
        sxx=calcSumSqr(x);
        sx=calcSum(x);
        syy=calcSumSqr(y);
        sy=calcSum(y);
        se2=(Size*syy-sy*sy-slope*slope*(Size*sxx-sx*sx))/(Size*(Size-2));
        sb2=Size*se2/(Size*sxx-sx*sx);
        sa2=sb2*sxx/Size;
        cout <<"sa2= "<<sa2<<endl;
        cout <<"sb2= "<<sb2<<endl;
        yError=sqrt(sa2);
        slopeError=sqrt(sb2);
        return correlation(x, y);
}

//template <typename T>
Real CalcCorrelation(vector<Real> &z, vector< vector<Real> > &data, vector<Real> &variable)
{
        Real r2;

        MultiVariableRegression(z, data, variable);
        r2=CalcCorrelationFromFit(z, data, variable);

        return r2;
}

Real Round(Real x, Real nearest)
{
        return Real(nearest*int(x/nearest+0.5));
}

#endif
