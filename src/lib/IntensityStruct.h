#ifndef _included_IntensityStruct
#define _included_IntensityStruct

struct IntensityStruct
{
	vector<Real> calc, CalcError, error, expi, s;
	vector<Real> vacuum, NoHydration;
	vector< vector< vector<Real> > > point, CrossTerm;
	vector< vector<Real> > f;
	vector< vector<Real> > sinc;
};

#endif
