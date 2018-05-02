#ifndef _overlap_included_
#define _overlap_included_

# include "TypeDef.h"

Real ConeVolume(Real r2, Real h)
{
	return pi*r2*h/3.0;
}

Real SphericalArc(Real r, Real ra, Real rb)
{
	return 2.0*pi*ra*ra*ra*(1.0-(ra*ra+r*r-rb*rb/(2.0*r*ra)))/3.0;
}

Real section(Real r, Real ra, Real rb)
{
	Real arc, cone, x=(ra*ra+r*r-rb*rb)/(2.0*r);
	arc=SphericalArc(r, ra, rb);
	cone=ConeVolume(ra*ra-x*x, x);	
	return arc-cone;
}

Real SphereVolume(Real r)
{
	return 4.0*pi*r*r*r/3.0;
}

Real overlap(Real r, Real ra, Real rb)
{
	Real volume=0;
	if (ra+rb>r)
	{
		if (ra>r+rb) return SphereVolume(ra);
		else if (rb>r+ra) return SphereVolume(rb);
		else
		{
			volume=section(r, ra, rb)+section(r, rb, ra);
			return volume;
		}
	}
	else return 0;
}

Real OverlapVolumeElement(Real r, Real ra, Real rb, Real dr)
{
	Real volume=0;
	volume=overlap(r, ra+dr, rb+dr);
	volume+=overlap(r, ra, rb);
	volume-=overlap(r, ra+dr, rb);
	volume-=overlap(r, ra, rb+dr);
	return volume;	
}

Real FreedOverlap(Real r, Real u1, Real u2, Real v1, Real v2)
{
	//return r*r*r*pi*((u2*u2*u2-u1*u1*u1)*(v2-v1)+(u2-u1)*(v2*v2*v2-v1*v1*v1))/12.0;
	return r*r*r*pi*((u2*u2*u2-u1*u1*u1)/3.0+v1*v1*(u1-u2))/2.0;
}

Real FreedOverlap(Real r, Real ra, Real rb, Real dr)
{
	Real integral1, integral2;
	Real u1, u2, u3, v1, v2, v3;
	if (ra+rb+2.0*dr<r) return 0;
	if (ra+dr+r<rb) return 0;
	if (rb+dr+r<ra) return 0;
	u1=(ra+rb)/r;
	v1=(ra-rb)/r;
	u2=(ra+rb+2.0*dr)/r;
	v2=(ra-rb)/r;
	u3=(ra+rb+2.0*dr)/r;
	v3=(ra-rb)/r;
	integral1=FreedOverlap(r, u1, u2, v1, v2);
	integral2=FreedOverlap(r, u2, u3, v2, v3);
	//return integral1*integral2*4.0/(r*r*r*pi);
	return integral1;
}

#endif
