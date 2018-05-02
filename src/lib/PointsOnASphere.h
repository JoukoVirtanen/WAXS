#ifndef _PointsOnASphere_included_
#define _PointsOnASphere_included_

# include "TypeDef.h"
# include "MathUtils.h"
# include "VectorManip.h"

void MakePointsOnAHemiSphere(Real ThetaInc, Real NumSlices, vector<VectorStruct> &v)
{
        Real cosphi=1.0/NumSlices-1.0;
        Real sinphi;
        Real theta=pi/ThetaInc;
        VectorStruct temp;
        while (true)
        {
                sinphi=sqrt(1.0-cosphi*cosphi);
                temp.x=sinphi*cos(theta);
                temp.y=sinphi*sin(theta);
                temp.z=cosphi;
                SafePushBack(v, temp, "in MakePointsOnAHemiSphere");
                theta+=2.0*pi/ThetaInc;
                if (theta>2.0*pi)
                {
                        theta=pi/ThetaInc;
                        cosphi+=2.0/NumSlices;
                        if (cosphi>-0.9/NumSlices) break;
                }
        }
}

void MakePointsOnASphere(Real ThetaInc, Real NumSlices, vector<VectorStruct> &v)
{
        Real cosphi=1.0/NumSlices-1.0;
        Real sinphi;
        Real theta=pi/ThetaInc;
        VectorStruct temp;
        while (true)
        {
                sinphi=sqrt(1.0-cosphi*cosphi);
                temp.x=sinphi*cos(theta);
                temp.y=sinphi*sin(theta);
                temp.z=cosphi;
                SafePushBack(v, temp, "in MakePointsOnAHemiSphere");
                theta+=2.0*pi/ThetaInc;
                if (theta>2.0*pi)
                {
                        theta=pi/ThetaInc;
                        cosphi+=2.0/NumSlices;
                        if (cosphi>1.0-0.9/NumSlices) break;
                }
        }
}

#endif
