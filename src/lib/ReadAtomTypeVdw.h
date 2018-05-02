#ifndef _ReadAtomTypeVdw_included_
#define _ReadAtomTypeVdw_included_

void ReadAtomTypeVdw(vector<Real> &AtomTypeVdw, string vdwFile)
{
        vector<string> lines;
        int Size;

        ReadLines(vdwFile, lines, "vdwFile");

        Size=lines.size();

        for (int i=0;i<Size;i++)
        {
                AtomTypeVdw[i]=StrToFloat(GetWord2(lines[i], 1));
        }
}
#endif
