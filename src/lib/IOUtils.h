#ifndef _IOUtils_included_
#define _IOUtils_included_

# include "StringUtils.h"
# include "LinkedList.h"
# include "Constants.h"
# include "VectorManip.h"
# include "TypeDef.h"

bool exist(char charfile[1000])
{
        bool n;

        fstream file;
        file.open(charfile,ios::in);

        if (!file) n=false;
        else n=true;

        return n;
}

bool exist(string File)
{
        bool n;
        char charfile[1000];

        strcpy(charfile, File.c_str());
        fstream file;
        file.open(charfile,ios::in);

        if (!file) n=false;
        else n=true;

        return n;
}

void FileNotValid(string file, string FileType)
{
        bool valid=false;

        valid=exist(file);

        while (!valid)
        {
                cin.ignore();
                cout << "The "<<FileType<<" "<<file<<" cannot be opended."<<endl;
                cout << "Enter file path again: ";
                getline(cin, file);
                valid=exist(file);
        }
}

void OpenFile(string File, fstream &file, string FileType)
{
        //char CharFile[1000];
        char * CharFile = new char [File.length()+1];
        //Trim(File);
        cout <<"File= "<<File<<endl;
        strcpy(CharFile, File.c_str());
        if (!exist(CharFile))
        {
                cout <<"Error. "<<FileType<<" "<<File<<" cannot be opened."<<endl;
                exit(EXIT_FAILURE);
        }
        file.open(CharFile, ios::in);
}

void OpenFileForWriting(string File, ofstream &file)
{
        char CharFile[1000];
        Trim(File);
        strcpy(CharFile, File.c_str());
        file.open(CharFile, ios::app);
}

void ReadLines(string File, vector<string> &lines, string FileType)
{
        fstream file;
        string line;

        OpenFile(File, file, FileType);
        getline(file, line);
        while (true)
        {
                if (file.eof()) break;
                SafePushBack(lines, line, "lines");
                getline(file, line);
        }
        file.close();
}

string GetWord(string line, int WordNum)
{
        int n, position, WordCount;
        string str;

        WordCount=count(line.begin(),line.end(),' ');

        if (WordNum<WordCount)
        {
                for (n=0;n<WordNum;n++)
                {

                        while (true)
                        {
                                str=line[0];

                                if (str==" ") line=line.substr(1,100);
                                else break;

                                WordCount=count(line.begin(),line.end(),' ');
                                if (WordCount==0) break;
                        }

                        if (WordCount>0)
                        {
                                position=line.find(" ");
                                str=line.substr(0,position);
                                line=line.substr(position, 100);
                        }
                        else str=line.substr(0,100);
                }
        }

        return str;
}

bool GetLineInFile(string PdbFilePaths, string &line, int nthstruct)
{
        fstream file;
        OpenFile(PdbFilePaths, file, "StructureFilePaths");
        cout <<"PdbFilePaths= "<<PdbFilePaths<<endl;
        getline(file, line);
        for (int i=0;i<nthstruct-1;i++)
        {
                getline(file, line);
                if (file.eof()) 
                {
                        file.close();
                        return false;
                }
                cout <<"line= "<<line<<endl;
        }
        file.close();
        return true;
}

bool GetLinesInFile(string FilePath, vector<string> &lines, int start, int end)
{
        string line;
        fstream file;
        cout <<"In GetLinesInFile start= "<<start<<" end= "<<end<<endl;
        OpenFile(FilePath, file, "StructureFilePaths");
        //cout <<"PdbFilePaths= "<<PdbFilePaths<<endl;
        lines.clear();
        getline(file, line);
        for (int i=1;i<=end;i++)
        {
                if (file.eof()) 
                {
                        file.close();
                        return false;
                }
                //if (i%100==0) cout <<"i= "<<i<<" start= "<<start<<" end= "<<end<<endl;
                if (i>=start)
                {
                        SafePushBack(lines, line, "in GetLinesInFile");
                }
                getline(file, line);
                //cout <<"line= "<<line<<endl;
        }
        file.close();
        cout <<"Leaving GetLinesInFile"<<endl;
        return true;
}

string IntToStr(int x)
{
        stringstream out;
        out << x;
        return out.str();
}

string RealToStr(Real x)
{
        stringstream out;
        out << x;
        return out.str();
}

bool isNumber(string str)
{
        int length=str.length();
        if (length==0) return false;
        if (length==1 && isdigit(str[0])) return true;
        if (length>1 && (isdigit(str[0]) || isdigit(str[1]))) return true;
        return false;
}

void AddIndexToFile(string &str)
{
        char CharFile[1000];
        string prefix, suffix;
        string StrFileIndex;
        bool valid;
        int FileIndex, position;

        strcpy(CharFile, str.c_str());
        valid=exist(CharFile);
        FileIndex=2;

        position=str.rfind(".");
        prefix=str.substr(0, position);
        suffix=str.substr(position+1, 4);

        while (valid)
        {
                stringstream out;
                out << FileIndex;
                StrFileIndex=out.str();
                str=prefix+StrFileIndex+"."+suffix;
                strcpy(CharFile, str.c_str());
                valid=exist(CharFile);
                FileIndex++;
        }
}

string GetWord2(string line, int WordNum)
{
        int n, position, tab, space, WordCount;
        string str;
        //cout <<"In GetWord2"<<endl;
        line="\t"+line+"\t";
        line=" "+line+" ";

        WordCount=count(line.begin(),line.end(),'\t');
        WordCount+=count(line.begin(),line.end(),' ');
        //cout <<"WordCount= "<<WordCount<<endl;
        if (WordNum<WordCount)
        {
                for (n=0;n<WordNum;n++)
                {

                        while (true)
                        {
                                str=line[0];

                                if (str=="\t" || str==" ") line=line.substr(1,1000);
                                else break;

                                WordCount=count(line.begin(),line.end(),'\t');
                                WordCount+=count(line.begin(),line.end(),' ');
                                //cout <<"WordCount2= "<<WordCount<<endl;
                                if (WordCount==0) break;
                        }

                        if (WordCount>0)
                        {
                                tab=line.find("\t");
                                space=line.find(" ");
                                //cout <<"tab= "<<tab<<" space= "<<space<<endl;
                                if (space<tab && space>-1) position=space;
                                else position=tab;
                                //cout <<"position= "<<position<<endl;
                                str=line.substr(0,position);
                                //cout <<"str= "<<str<<endl;
                                line=line.substr(position, 1000);
                                //cout <<"line= "<<line<<endl;
                        }
                        else str=line.substr(0,1000);
                }
        }

        return str;
}

string RemoveExtension(string file)
{
        size_t position;

        position=file.rfind(".");
        if (position!=string::npos)
        {
                return file.substr(0,position);
        }
        return file;
}

string GetBase(string file)
{
        size_t position;

        position=file.rfind("/");
        if (position!=string::npos)
        {
                return file.substr(position+1, file.length()-position-1);
        }
        return file;
}

bool StrToBool(const string &str)
{
        return (str == BOOL_TRUE_STR) || (str == BOOL_T_STR) ||
                (str == BOOL_YES_STR) || (str == BOOL_Y_STR) ||
                (str == BOOL_true_STR) || (str == BOOL_t_STR) ||
                (str == BOOL_yes_STR) || (str == BOOL_y_STR) ||
                (str == BOOL_1_STR);
}

string BoolToStr(bool b)
{
        if (b) return "yes";
        else return "no";
}

Real StrToReal(string str)
{
        return Real(StrToFloat(str));
}

void Tokenize(string line, string str, vector<string> &v)
{
        size_t found;
        DeleteVector(v);
        while (true)
        {
                found=line.find(str);
                if (found!=string::npos)
                {
                        if (found>0) SafePushBack(v, line.substr(0,found), "v");
                        line=line.substr(found+1,line.length()-found-1);
                }
                else break;
        }
        SafePushBack(v, line, "v");
}

void DeleteFromLeft(string &line, string str)
{
        while (true)
        {
                if (line.substr(0,str.length())==str) line=line.substr(str.length(), line.length()-str.length());
                else break;
        }
}

void Tokenize2(string line, string str, vector<string> &v)
{
        size_t found, start, end;
        string line2;
        DeleteVector(v);
        if (line!="")
        {
                while (true)
                {
                        DeleteFromLeft(line, str);
                        found=line.find(str);
                        if (found!=string::npos)
                        {
                                if (found>0) SafePushBack(v, line.substr(0,found), "v");
                                start=found+1;
                                end=line.length()-found-1;
                                if (start<line.length() && end>0) 
                                {
                                        line2=line.substr(start);
                                        line=line2;
                                }
                                else line="";
                        }
                        else break;
                }
                SafePushBack(v, line, "v");
        }
}

void Delete(string &line, string str)
{
        int i, Size;
        vector<string> v;
        Tokenize(line, str, v);
        Size=v.size();
        line="";
        for (i=0;i<Size;i++) line+=v[i];
}

void DeleteStartingAt(string &line, string str)
{
        size_t found;
        //cout <<line<<endl;
        if (line!="")
        {
                found=line.find(str);
                if (found!=string::npos) 
                {
                        if (found==0) line="";
                        else line=line.substr(0,found);
                }
        }
}

void strip(string &line)
{
        if (line!="")
        {
                while (true)
                {
                        if (line.substr(line.length()-1, 1)==" ")
                        {
                                line=line.substr(0, line.length()-1);
                        }
                        else break;
                }
        }
}

void DeleteChar(string &str, const char &chr)
{
        string::size_type p = str.find(chr, 0);
        while (p != string::npos) { str.erase(p); p = str.find(chr, 0); }
}

void ProcessLinesForReadingParameters(vector<string> &lines)
{
        int NumLines=lines.size();
        vector<string> temp, str;
        size_t space;

        for (int i=0;i<NumLines;i++)
        {
                Delete(lines[i], "\n");
                DeleteChar(lines[i], CARRIAGE_RETURN);
                DeleteStartingAt(lines[i], "#");
                strip(lines[i]);
                space=lines[i].find(" ");
                if (space!=string::npos) 
                {
                        Tokenize2(lines[i], " ", str);
                        if (str.size()>1) SafePushBack(temp, lines[i], "temp");		
                }
        }
        lines=temp;
}

void readMatrix(string InputFile, Matrix &matrix)
{
        int ncolumn, nrow;
        vector<string> lines, str;
        ReadLines(InputFile, lines, "MatrixFile");
        nrow=lines.size();
        Tokenize2(lines[0], "\t", str);
        ncolumn=str.size();
        Safe2DAlloc(matrix, nrow, ncolumn, "matrix");
        for (int i=0;i<nrow;i++)
        {
                Tokenize2(lines[i], "\t", str);
                for (int j=0;j<ncolumn;j++)
                {
                        matrix[i][j]=StrToFloat(str[j]);
                }
        }
}

bool WildCardCompare2(string str1, string str2)
{
        bool BoolPrint=false;
        int Length;

        //if (str1=="CW") BoolPrint=true;
        //else BoolPrint=false;

        if (BoolPrint) cout <<endl;
        if (BoolPrint) cout <<"str1= "<<str1<<" str2= "<<str2<<endl;
        DeleteStartingAt(str2, "*");
        Length=str2.length();
        if (BoolPrint) cout <<"str2= "<<str2<<" str1.substr= "<<str1.substr(0, Length)<<endl;
        if (str1.substr(0, Length)==str2) return true;
        else return false;
}

bool WildCardCompare(string str1, string str2)
{
        bool BoolPrint=false;
        string::size_type found1, found2;
        //if (str1=="CW") BoolPrint=true;
        //else BoolPrint=false;
        if (BoolPrint) cout <<"In WildCardCompare"<<endl;
        found1=str1.find("*");
        found2=str2.find("*");
        if (BoolPrint) cout <<"found1= "<<found1<<" found2= "<<found2<<endl;
        if (BoolPrint) cout <<"str1= "<<str1<<" str2= "<<str2<<endl;

        if (found1==string::npos && found2==string::npos)
        {
                return str1==str2;
        }
        if (found1!=string::npos && found2==string::npos)
        {
                return WildCardCompare2(str2, str1);
        }
        if (found1==string::npos && found2!=string::npos)
        {
                return WildCardCompare2(str1, str2);
        }
        if (found1!=string::npos && found2!=string::npos)  //This might need to be corrected
        {
                return str1==str2;
        }
        return false;
}

Real getValueOfReal(string variableName, string str)
{
        string phrase, strVal;
        size_t position, length;
        str+=" ";
        length=str.length();
        position=str.find(variableName);
        if (position!=string::npos)
        {
                phrase=str.substr(position+variableName.length(), length-position-variableName.length()-1);
                strVal=GetWord2(phrase, 1);
                return StrToFloat(strVal);
        }
        return UNK_REAL;
}

bool getValueOfBool(string variableName, string str)
{
        string phrase, strVal;
        size_t position, length;
        length=str.length();
        position=str.find(variableName);
        if (position!=string::npos)
        {
                phrase=str.substr(position+variableName.length(), length-position-variableName.length()-1);
                strVal=GetWord2(phrase, 1);
                return StrToBool(strVal);
        }
        return false;
}

int getValueOfInt(string variableName, string str)
{
        string phrase, strVal;
        size_t position, length;
        length=str.length();
        position=str.find(variableName);
        if (position!=string::npos)
        {
                phrase=str.substr(position+variableName.length(), length-position-variableName.length()-1);
                strVal=GetWord2(phrase, 1);
                return StrToInt(strVal);
        }
        return false;
}

#endif

