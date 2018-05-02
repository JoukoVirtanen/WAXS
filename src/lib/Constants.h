#ifndef _Constants_included_
#define _Constants_included_

# include "TypeDef.h"

const Real pi=3.1415926535389793;
const long double LongPi=3.14159265358979323846264338327950288419716939937510;
const Real pi32=5.568327997;
const Real packing=1.105339; //(3*2^0.5/pi)^(1/3)
const Real GaussianToHardSphere=1.099543; //(4/(3*pi^0.5))^(1/3)
const Real HardSphereToGaussianSphere=1.0/GaussianToHardSphere;
const Real CubeToSphere=0.62035; //(3/(4pi))^(1/3)
const Real RAD_TO_DEGREE=360.0/(2.0*pi);
const Real DEGREE_TO_RAD=2.0*pi/360.0;
const Real Kb=1.3806504e-23;
const Real NA=6.02214179e23;
const Real ATOMS_PER_A3_TO_MOLAR=1660.540187;
const Real SCALE_TO_KCAL=332.0634;
const Real COULOMB=6.24150965e18;
const int InProtein=0;
const int InHydrationShell=1;
const int InSolution=2;

const int HYDROPHOBIC=0;
const int HYDROPHILIC=1;

const int ATOM_TYPES=0;
const int ELEMENTS=1;
const int NORMAL_ATOM_TYPES=2;
const int NORMAL_ELEMENTS=3;

const int X=0;
const int Y=1;
const int Z=2;

const int PHI=0;
const int PSI=1;
const int OMEGA=2;

const int NUM_RES_TYPES=20;

const Real UNK_REAL=-666666.0;
const int UNK_INT=-666666;
#endif
