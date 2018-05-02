#ifndef _VectorManip_included_
#define _VectorManip_included_

# include <iostream>
# include <fstream>
# include <string>
# include <cstring>
# include <cmath>
# include <sstream>
# include <time.h>
# include <iomanip>
# include <vector>

        template <class T>
void Print3DVectorSize(vector< vector< vector<T> > > &ThreeDArray, string name)
{
        int Size1=0, Size2=0, Size3=0;
        Size1=ThreeDArray.size();
        if (Size1!=0) Size2=ThreeDArray[0].size();
        if (Size2!=0) Size3=ThreeDArray[0][0].size();
        cout <<name<<"["<<Size1<<"]["<<Size2<<"]["<<Size3<<"]"<<endl;
}

        template <class T>
void Print2DVectorSize(vector< vector<T> >  &TwoDArray, string name)
{
        int Size1=0, Size2=0;
        Size1=TwoDArray.size();
        cout <<"In Print2DVectorSize"<<endl;
        if (Size1!=0) Size2=TwoDArray[0].size();
        cout <<name<<"["<<Size1<<"]["<<Size2<<"]"<<endl;
}

        template <class T>
void Get2DVectorSize(vector< vector<T> > &TwoDArray, int &Size1, int &Size2, string name)
{
        Size1=TwoDArray.size();
        if (Size1!=0) Size2=TwoDArray[0].size();
        else
        {
                cout <<"Vector "<<name<<" is empty"<<endl;
                exit(EXIT_FAILURE);
        }
}

        template <class T>
void Get2DVectorSize(vector< vector<T> > &TwoDArray, int &Size1, int &Size2)
{
        Size1=TwoDArray.size();
        if (Size1!=0) Size2=TwoDArray[0].size();
        else
        {
                Size2=0;
        }
}
        template <class T>
void Get3DVectorSize(vector< vector< vector<T> > > &ThreeDArray, int &Size1, int &Size2, int &Size3, string name)
{
        Size1=0; Size2=0; Size3=0;
        Size1=ThreeDArray.size();
        if (Size1!=0) Size2=ThreeDArray[0].size();
        if (Size2!=0) Size3=ThreeDArray[0][0].size();
        else
        {
                cout <<"Vector "<<name<<" is empty"<<endl;
                exit(EXIT_FAILURE);
        }
}

        template <class T>
int countNumZeroElements(vector< vector<T> > &TwoDArray)
{
        int Size1, Size2;
        int numZero=0;
        Get2DVectorSize(TwoDArray, Size1, Size2);

        for (int i=0;i<Size1;i++)
        {
                for (int j=0;j<Size2;j++)
                {
                        if (TwoDArray[i][j]==0) numZero++;
                }
        }
        return numZero;
}

        template <class T>
int countNumPositiveElements(vector< vector<T> > &TwoDArray)
{
        int Size1, Size2;
        int numPositive=0;
        Get2DVectorSize(TwoDArray, Size1, Size2);

        for (int i=0;i<Size1;i++)
        {
                for (int j=0;j<Size2;j++)
                {
                        if (TwoDArray[i][j]>0) numPositive++;
                }
        }
        return numPositive;
}

        template <class T>
int countNumNegativeElements(vector< vector<T> > &TwoDArray)
{
        int Size1, Size2;
        int numNegative=0;
        Get2DVectorSize(TwoDArray, Size1, Size2);

        for (int i=0;i<Size1;i++)
        {
                for (int j=0;j<Size2;j++)
                {
                        if (TwoDArray[i][j]<0) numNegative++;
                }
        }
        return numNegative;
}

        template <class T>
int countNumNanElements(vector< vector<T> > &TwoDArray)
{
        int Size1, Size2;
        int numNan=0;
        Get2DVectorSize(TwoDArray, Size1, Size2);

        for (int i=0;i<Size1;i++)
        {
                for (int j=0;j<Size2;j++)
                {
                        if (isnan(TwoDArray[i][j])) numNan++;
                }
        }
        return numNan;
}

        template <class T>
int countNumInfElements(vector< vector<T> > &TwoDArray)
{
        int Size1, Size2;
        int numInf=0;
        Get2DVectorSize(TwoDArray, Size1, Size2);

        for (int i=0;i<Size1;i++)
        {
                for (int j=0;j<Size2;j++)
                {
                        if (isinf(TwoDArray[i][j])) numInf++;
                }
        }
        return numInf;
}

        template <class T>
void check2DVector(vector< vector<T> > &TwoDArray, string name)
{
        int numZeroElements, numPositiveElements, numNegativeElements;
        int numNanElements, numInfElements;
        Print2DVectorSize(TwoDArray, name);
        numZeroElements=countNumZeroElements(TwoDArray);
        numPositiveElements=countNumPositiveElements(TwoDArray);
        numNegativeElements=countNumNegativeElements(TwoDArray);
        numNanElements=countNumNanElements(TwoDArray);
        numInfElements=countNumInfElements(TwoDArray);
        cout <<"numZeroElements= "<<numZeroElements
        <<" numPositiveElements= "<<numPositiveElements
        <<" numNegativeElements= "<<numNegativeElements
        <<" numNanElements= "<<numNanElements
        <<" numInfElements= "<<numInfElements<<endl;
}


        template <class Vector_T>
void CreateVector(Vector_T &v, int Size)
{
        int i;
        for (i=0;i<Size;i++) v.push_back(0);
}

        template <class T>
bool SafeAlloc(vector<T> &v, T fill, int Size)
{
        bool MemoryAllocated=true;
        v.clear();
        try
        {
                v.resize(Size, fill);
        }
        catch(const exception& e)
        {
                cout <<e.what()<<endl;
                cout <<"In SafeAlloc"<<endl;
                MemoryAllocated=false;
        }
        return MemoryAllocated;
}

        template <class T>
bool SafeAlloc(vector<T> &v, int Size)
{
        return SafeAlloc(v, T(0), Size);
}

        template <class T>
void SafeAlloc(vector<T> &v, T fill, int Size, string name)
{
        if (!SafeAlloc(v, fill, Size)) 
        {
                cout <<"ERROR: Unable to allocate memory for vector "<<name<<endl;
                exit(EXIT_FAILURE);
        }
}

        template <class T>
void SafeAlloc(vector<T> &v, int Size, string name)
{
        SafeAlloc(v, T(0), Size, name);
}

        template <class T>
bool Safe2DAlloc(vector< vector<T> > &m, T fill, int Size1, int Size2)
{
        vector<T> v;
        cout <<"m.size()= "<<m.size()<<endl;
        m.clear();
        bool MemoryAllocated=true;
        try
        {
                cout <<"Size2= "<<Size2<<" Size1= "<<Size1<<endl;
                v.resize(Size2, fill);
                m.resize(Size1, v);
        }
        catch(bad_alloc&)
        {
                MemoryAllocated=false;
        }
        return MemoryAllocated;
}

        template <class T>
void Safe2DAlloc(vector< vector<T> > &m, T fill, int Size1, int Size2, string name)
{
        if (!Safe2DAlloc(m, fill, Size1, Size2)) 
        {
                cout <<"ERROR: Unable to allocate memory for vector "<<name<<endl;
                exit(EXIT_FAILURE);
        }
}

        template <class T>
void Safe2DAlloc(vector< vector<T> > &m, int Size1, int Size2, string name)
{
        Safe2DAlloc(m, T(0), Size1, Size2, name);
}

        template <class T>
bool SafePushBack(vector<T> &v, T element)
{
        bool MemoryAllocated=true;
        try
        {
                v.push_back(T(element));
        }
        catch(bad_alloc&)
        {
                MemoryAllocated=false;
        }
        return MemoryAllocated;
}

        template <class T>
void SafePushBack(vector<T> &v, T element, string name)
{
        if (!SafePushBack(v, T(element)))
        {
                cout <<"ERROR: Unable to push back to vector "<<name<<endl;
                exit(EXIT_FAILURE);
        }
}

        template <class Vector_T>
void DeleteVector(Vector_T &v)
{
        int i, Size=v.size();
        for (i=0;i<Size;i++) v.pop_back();
}

        template <class T>
void ZeroVector(vector<T> &v)
{
        int i, Size=v.size();
        for (i=0;i<Size;i++) v[i]=0;
}

        template <class T>
bool SafeArrayAlloc(T *&array, T fill, int Size)
{
        try 
        {
                array=new T[Size];
        }
        catch(const exception &e) 
        {
                cout <<e.what()<<endl;
                return false;
        }
        for (int i=0;i<Size;i++) array[i]=fill;
        return true;
}

        template <class T>
bool SafeArrayAlloc(T *&array, int Size)
{
        return SafeArrayAlloc(array, T(0), Size);
}

        template <class T>
void SafeArrayAlloc(T *&array, int Size, string name)
{
        if (!SafeArrayAlloc(array, Size))
        {
                cout <<"ERROR: Unable to allocate memory for array "<<name<<endl;
                exit(EXIT_FAILURE);
        }
}

        template <class T>
void SafeArrayAlloc(T *&array, T fill, int Size, string name)
{
        if (!SafeArrayAlloc(array, fill, Size))
        {
                cout <<"ERROR: Unable to allocate memory for array "<<name<<endl;
                exit(EXIT_FAILURE);
        }
}

        template <class T>
bool Safe2DArrayAlloc(T **&array, int Size1, int Size2)
{
        try
        {
                array=new T *[Size1];
                for (int i=0;i<Size1;i++)
                {
                        *(array+i)=new T[Size2];
                }
        }
        catch(bad_alloc&)
        {
                return false;
        }
        return true;
}

        template <class T>
void Safe2DArrayAlloc(T **&array, int Size1, int Size2, string name)
{
        if (!Safe2DArrayAlloc(array, Size1, Size2))
        {
                cout <<"ERROR: Unable to allocate memory for array "<<name<<endl;
                exit(EXIT_FAILURE);
        }
}

        template <class T>
bool Safe3DArrayAlloc(T ***&array, int Size1, int Size2, int Size3)
{
        try
        {
                array=new T **[Size1];
                for (int i=0;i<Size1;i++)
                {
                        *(array+i)=new T *[Size2];
                        for (int j=0;j<Size2;j++)
                        {
                                *(*(array+i)+j)=new T[Size3];
                        }
                }
        }
        catch(bad_alloc&)
        {
                return false;
        }
        return true;
}

        template <class T>
void Safe3DArrayAlloc(T ***&array, int Size1, int Size2, int Size3, string name)
{
        if (!Safe3DArrayAlloc(array, Size1, Size2, Size3))
        {
                cout <<"ERROR: Unable to allocate memory for array "<<name<<endl;
                exit(EXIT_FAILURE);
        }
}

        template <class T>
void Delete3DArray(T ***&array, int Size1, int Size2, int Size3)
{
        for (int i=0;i<Size1;i++)
        {
                for (int j=0;j<Size2;j++)
                {
                        delete [] array[i][j];
                }
                delete [] array[i];
        }
        delete [] array;
}


/*
   template <class T>
   bool Create3DArray(vector< vector< vector<T> > > &array, int size1, int size2, int size3)
   {
   bool MemoryAllocated=true;
   int i;
   vector<T> v;
   vector< vector<T> > matrix;

   try
   {
   array.clear();
   SafeAlloc(v, T(0), size3);
   SafeAlloc(matrix, v, size2);
   MemoryAllocated=SafeAlloc(array, matrix, size1);
   }
   catch (bad_alloc&)
   {
   return false;
   }
   cout <<"MemoryAllocated= "<<MemoryAllocated<<endl;
   return MemoryAllocated;
   }
   */

        template <class T>
void EraseVector(vector<T> &v)
{
        cout <<"In EraseVector"<<endl;
        vector<T>().swap(v);
        cout <<"Leaving EraseVector"<<endl;
}

        template <class T>
void SafeResize(vector<T> &v, string name)
{
        vector<T> temp;
        SafeAlloc(temp, v[0], v.size(), name);
        temp=v;
        temp.swap(v);
}

        template <class T>
bool Safe3DAlloc(vector< vector< vector<T> > > &array, T &fill, int size1, int size2, int size3)
{
        bool MemoryAllocated=true;
        vector<T> v, temp;
        vector< vector<T> > matrix;
        MemoryAllocated=SafeAlloc(temp, size1*size2*size3);
        if (MemoryAllocated) cout <<"Allocated test 1D vector"<<endl;
        else cout <<"Unable to allocate test 1D vector"<<endl;
        temp.clear();
        vector<T>().swap(temp);
        if (MemoryAllocated)
        {
                cout <<"Attempting to allocate 3D array"<<endl;
                array.clear();
                try
                {
                        v.resize(size3, fill);
                        matrix.resize(size2, v);
                        array.resize(size1, matrix);
                }
                catch (const exception& e)
                {
                        cout <<e.what()<<endl;
                        MemoryAllocated=false;
                }
        }
        return MemoryAllocated;
}

        template <class T>
bool Safe3DAlloc(vector< vector< vector<T> > > &array, int size1, int size2, int size3)
{
        bool MemoryAllocated=true;
        vector<T> v, temp;
        vector< vector<T> > matrix;
        MemoryAllocated=SafeAlloc(temp, size1*size2*size3);
        if (MemoryAllocated) cout <<"Allocated test 1D vector"<<endl;
        else cout <<"Unable to allocate test 1D vector"<<endl;
        temp.clear();
        vector<T>().swap(temp);
        if (MemoryAllocated)
        {
                cout <<"Attempting to allocate 3D array"<<endl;
                array.clear();
                try
                {
                        v.resize(size3);
                        matrix.resize(size2, v);
                        array.resize(size1, matrix);
                }
                catch (const exception& e)
                {
                        cout <<e.what()<<endl;
                        MemoryAllocated=false;
                }
        }
        return MemoryAllocated;
}

        template <class T>
void Safe3DAlloc(vector< vector< vector<T> > > &array, int size1, int size2, int size3, string name)
{
        T fill=0;
        if (!Safe3DAlloc(array, fill, size1, size2, size3))
        {
                cout <<"ERROR: Unable to allocate memory for 3D array "<<name<<endl;
                exit(EXIT_FAILURE);
        }
}

        template <class T>
void Safe3DAlloc(vector< vector< vector<T> > > &array, T &fill, int size1, int size2, int size3, string name)
{
        if (!Safe3DAlloc(array, fill, size1, size2, size3))
        {
                cout <<"ERROR: Unable to allocate memory for 3D array "<<name<<endl;
                exit(EXIT_FAILURE);
        }
}

        template <class Vector_T>
void UnsafeCreate3DArray(vector< vector< Vector_T > > &array, int size1, int size2, int size3)
{
        int i;
        Vector_T v;
        vector< Vector_T > matrix;

        array.clear();
        for (int i=0;i<size3;i++) v.push_back(0);
        for (int i=0;i<size2;i++) matrix.push_back(v);
        for (int i=0;i<size1;i++) array.push_back(matrix);
}

/*
   template <class Vector_T>
   void Create3DArray(vector< vector<Vector_T> > &array, int size1, int size2, int size3)
   {
   bool MemoryAllocated=true;
   int i;
   Vector_T v;
   vector<Vector_T> matrix;

   DeleteVector(array);
   v.resize(size3, 0);
   matrix.resize(size2, v);
   array.resize(size1, matrix);
   }
   */
        template <class Vector_T>
void CopyVector(const Vector_T &v1, Vector_T &v2)
{
        int i, Size=v1.size();
        DeleteVector(v2);
        for (i=0;i<Size;i++) v2.push_back(v1[i]);
}

        template <class Vector_T>
void AppendVector(const Vector_T &v1, Vector_T &v2)
{
        int i, Size=v1.size();
        cout <<"Size= "<<Size<<endl;
        for (i=0;i<Size;i++) v2.push_back(v1[i]);
        cout <<"Size2= "<<v2.size()<<endl;
}

        template <class T>
bool SafeAppendVector(vector<T> &v1, vector<T> &v2)
{
        int Size=v1.size();
        try
        {
                for (int i=0;i<Size;i++) v2.push_back(v1[i]);
        }
        catch(bad_alloc&)
        {
                return false;
        }
        return true;
}

        template <class T>
void SafeAppendVector(vector<T> &v1, vector<T> &v2, string name)
{
        if (!SafeAppendVector(v1, v2))
        {
                cout <<"ERROR: Unable to append vector "<<name<<endl;
                exit(EXIT_FAILURE);
        }
}

        template <class T>
bool SafeCatVector(vector<T> &v1, vector<T> &v2)
{
        int Size1=v1.size();
        int Size2=v2.size();
        try
        {
                v2.resize(Size1+Size2);
                for (int i=0;i<Size1;i++) 
                {
                        v2[Size2+i]=v1[i];
                }
        }
        catch(bad_alloc&)
        {
                return false;
        }
        return true;
}

        template <class T>
void SafeCatVector(vector<T> &v1, vector<T> &v2, string name)
{
        if (!SafeCatVector(v1, v2))
        {
                cout <<"ERROR: Unable to concatenate vectors "<<name<<endl;
                exit(EXIT_FAILURE);
        }
}

        template <class T>
bool SafeCatVector2(vector<T> &v1, vector<T> &v2)
{
        int Size1=v1.size();
        int Size2=v2.size();
        try
        {
                for (int i=Size1-1;i<=0;i--) 
                {
                        v2.resize(v2.size()+1);
                        v2[v2.size()-2]=v1[i];
                        v1.resize(v1.size()-1);
                        if (i%1000==0)
                        {
                                cout <<"v1.capacity()= "<<v1.capacity()<<" v2.capacity()= "<<v2.capacity()<<endl;
                        }
                        //cout <<"v1.size()= "<<v1.size()<<" v2.size()= "<<v2.size()<<endl;
                }
        }
        catch(bad_alloc&)
        {
                return false;
        }
        return true;
}

        template <class T>
void SafeCatVector2(vector<T> &v1, vector<T> &v2, string name)
{
        if (!SafeCatVector2(v1, v2))
        {
                cout <<"ERROR: Unable to concatenate vectors "<<name<<endl;
                exit(EXIT_FAILURE);
        }
}

        template <class T>
bool SafeCatVector3(vector<T> &v1, vector<T> &v2)
{
        int Size1=v1.size();
        int Size2=v2.size();
        cout <<"Size1= "<<Size1<<endl;
        cout <<"Size2= "<<Size2<<endl;
        try
        {
                v1.resize(Size1+Size2);
                cout <<"v1.size()= "<<v1.size()<<endl;
                for (int i=Size1-1;i>=0;i--) v1[i+Size2]=v1[i];
                for (int i=0;i<Size2;i++) v1[i]=v2[i];
        }
        catch(bad_alloc&)
        {
                return false;
        }
        return true;
}

        template <class T>
void SafeCatVector3(vector<T> &v1, vector<T> &v2, string name)
{
        if (!SafeCatVector3(v1, v2))
        {
                cout <<"ERROR: Unable to concatenate vectors "<<name<<endl;
                exit(EXIT_FAILURE);
        }
}

        template <class Vector_T>
void CreateMatrix(vector<Vector_T> &m, int size1, int size2)
{
        Vector_T v;
        DeleteVector(m);
        v.resize(size2);
        m.resize(size1, v);
}

template <class T>
void PrintVector(vector<T> &v, string name)
{
        int Size=v.size();
        for (int i=0;i<Size;i++) cout <<name<<"["<<i<<"]= "<<v[i]<<endl;
}

        template <class Vector_T>
void PrintVector(Vector_T &v)
{
        int i, Size=v.size();
        for (i=0;i<Size;i++) cout <<v[i]<<endl;
}

        template <class T>
void PrintVectorHorizontal(vector<T> &v)
{
        int i, Size=v.size();
        for (i=0;i<Size;i++) cout <<v[i]<<"\t";
        cout <<endl;
}
#endif
