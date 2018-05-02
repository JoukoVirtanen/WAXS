#ifndef _LJPotential_included_
#define _LJPotential_included_

# include "TypeDef.h"

Real LJPotential(Real epsilon, Real VDWRadius, Real r)
{
	Real vr, vr2, vr6;
	vr=VDWRadius/r;
	vr2=vr*vr;
	vr6=vr2*vr2*vr2;
	return 4.0*epsilon*(vr6*vr6-2.0*vr6);
}

Real LJPotential(AtomStruct &Atom1, AtomStruct &Atom2)
{
        Real epsilon12=sqrt(Atom1.epsilon*Atom2.epsilon);
        Real VDWRadius12=0.5*(Atom1.vdw+Atom2.vdw);
        Real r=AtomDistance(Atom1, Atom2);
        return LJPotential(epsilon12, VDWRadius12, r);
}

Real InvertLJPotential(Real epsilon, Real VDWRadius, Real energy)
{
	Real r, r6, sq;
	sq=4.0+energy/epsilon;
	if (sq<0)
	{
		cout <<"Error there is no distance that corresponds to that LJ potential"<<endl;
		return 0;
	}
	r6=2.0/(2.0+sqrt(sq));
	r=VDWRadius*exp(log(r6)/6.0);
	return r;
}

#endif
