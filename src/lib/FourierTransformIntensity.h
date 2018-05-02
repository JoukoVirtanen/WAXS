#ifndef _included_FourierTransformIntensity
#define _included_FourierTransformIntensity

void FourierTransformIntensity(IntensityStruct &i, vector<Real> &PR, Real bin)
{
        int points=i.s.size();
        int GreatestDistance=PR.size();
        Real r=bin, sinc=i.s[1]-i.s[0];
        for (int distance=1;distance<GreatestDistance;distance++)
        {
                for (int n=0;n<points;n++)
                {
                        if (distance==1)
                        {
                                cout <<"PR["<<distance<<"]= "<<PR[distance]
                                <<" i.calc["<<n<<"]= "<<i.calc[n]
                                <<" i.s= "<<i.s[n]<<" r= "<<r<<endl;
                        }
                        PR[distance]+=i.calc[n]*i.s[n]*sin(i.s[n]*r)*sinc;
                }
                PR[distance]=r*PR[distance]*bin/(2.0*pi);
                r+=bin;
        }
}

#endif

