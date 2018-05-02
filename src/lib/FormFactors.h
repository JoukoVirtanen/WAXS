#ifndef _included_FormFactors
#define _included_FormFactors

void FourGaussianParameters(Real a[][NumAtomTypes], Real b[][NumAtomTypes], Real c[])
{

        //Data obtained from:
        //Dirac-Fock calculations of X-ray scattering factors. Acta Cryst. (1994). A50, 481-497	

        for (int m=0;m<NumAtomTypes;m++)
        {
                c[m]=0;
                for (int n=0;n<5;n++)
                {
                        a[n][m]=0;
                        b[n][m]=0;
                }
        }

        a[0][HYDROGEN]=0.493002;                  
        b[0][HYDROGEN]=10.5109;                  
        a[1][HYDROGEN]=0.322912;
        b[1][HYDROGEN]=26.1257;
        a[2][HYDROGEN]=0.140191;
        b[2][HYDROGEN]=3.14236;
        a[3][HYDROGEN]=0.040810;
        b[3][HYDROGEN]=57.7997;

        a[0][CARBON]=2.6158;
        b[0][CARBON]=11.3632;
        a[1][CARBON]=0.2279;
        b[1][CARBON]=3.0830;
        a[2][CARBON]=1.5983;
        b[2][CARBON]=0.3787;
        a[3][CARBON]=1.5602;
        b[3][CARBON]=49.7088;

        a[0][NITROGEN]=0.4743;
        b[0][NITROGEN]=0.1041;
        a[1][NITROGEN]=2.9082;
        b[1][NITROGEN]=9.1890;
        a[2][NITROGEN]=2.2778;
        b[2][NITROGEN]=27.0869;
        a[3][NITROGEN]=1.3332;
        b[3][NITROGEN]=0.4612;

        a[0][OXYGEN]=3.4716;
        b[0][OXYGEN]=11.9964;
        a[1][OXYGEN]=1.8289;
        b[1][OXYGEN]=4.7941;
        a[2][OXYGEN]=1.7198;
        b[2][OXYGEN]=0.2372;
        a[3][OXYGEN]=0.9790;
        b[3][OXYGEN]=31.9917;

        a[0][SULFUR]=7.1301;
        b[0][SULFUR]=1.4247;
        a[1][SULFUR]=5.0712;
        b[1][SULFUR]=21.7545;
        a[2][SULFUR]=2.0611;
        b[2][SULFUR]=0.0994;
        a[3][SULFUR]=1.7369;
        b[3][SULFUR]=54.2128;

        a[0][IRON]=8.7267;
        b[0][IRON]=3.8443;
        a[1][IRON]=8.0070;
        b[1][IRON]=0.2319;
        a[2][IRON]=6.4440;
        b[2][IRON]=9.4874;
        a[3][IRON]=2.8058;
        b[3][IRON]=65.8373;

        a[0][PHOSPHORUS]=7.1583;
        b[0][PHOSPHORUS]=1.7276;
        a[1][PHOSPHORUS]=3.4832;
        b[1][PHOSPHORUS]=23.7455;
        a[2][PHOSPHORUS]=2.0781;
        b[2][PHOSPHORUS]=0.1140;
        a[3][PHOSPHORUS]=2.2770;
        b[3][PHOSPHORUS]=57.1173;
        /****************************************
          a[0][O_Minus]=3.106934;
          b[0][O_Minus]=19.86808;
          a[1][O_Minus]=3.235142;
          b[1][O_Minus]=6.960252;
          a[2][O_Minus]=1.148886;
          b[2][O_Minus]=0.170043;
          a[3][O_Minus]=0.783981;
          b[3][O_Minus]=65.693509;
          a[4][O_Minus]=0.676953;
          b[4][O_Minus]=0.630757;

          c[O_Minus]=0.046136;
         *****************************************/
        a[0][NA_Plus]=4.4278;
        b[0][NA_Plus]=5.2355;
        a[1][NA_Plus]=2.4274;
        b[1][NA_Plus]=2.2658;
        a[2][NA_Plus]=1.7182;
        b[2][NA_Plus]=0.1272;
        a[3][NA_Plus]=1.4258;
        b[3][NA_Plus]=12.6031;

        for (int m=0;m<NumAtomTypes;m++)
        {
                for (int n=0;n<5;n++)
                {
                        b[n][m]/=16.0*pi*pi;
                }
        }
}

void FiveGaussianParameters(Real a[][NumAtomTypes], Real b[][NumAtomTypes], Real c[])
{
        //Calculates scattering factors for atoms.

        for (int m=0;m<NumAtomTypes;m++)
        {
                c[m]=0;
                for (int n=0;n<5;n++)
                {
                        a[n][m]=0;
                        b[n][m]=0;
                }
        }

        a[0][HYDROGEN]=0.493002;                     //Data obtained from:
        b[0][HYDROGEN]=10.5109;                      //Waasmaier, D., and Kirfel, A. 1995. New Analytical Scattering-Factor Functions for Free Atoms and Ions. Acta Cryst.A51: 416-431
        a[1][HYDROGEN]=0.322912;
        b[1][HYDROGEN]=26.1257;
        a[2][HYDROGEN]=0.140191;
        b[2][HYDROGEN]=3.14236;
        a[3][HYDROGEN]=0.040810;
        b[3][HYDROGEN]=57.7997;
        a[4][HYDROGEN]=0;
        b[4][HYDROGEN]=0;

        c[HYDROGEN]=0.003038;

        a[0][CARBON]=2.657506;
        b[0][CARBON]=14.780758;
        a[1][CARBON]=1.078079;
        b[1][CARBON]=0.7767775;
        a[2][CARBON]=1.490909;
        b[2][CARBON]=42.086843;
        a[3][CARBON]=-4.24107;
        b[3][CARBON]=-0.000294;
        a[4][CARBON]=0.713791;
        b[4][CARBON]=0.239535;

        c[CARBON]=4.297983;

        a[0][NITROGEN]=11.893780;
        b[0][NITROGEN]=0.000158;
        a[1][NITROGEN]=3.277479;
        b[1][NITROGEN]=10.232723;
        a[2][NITROGEN]=1.858092;
        b[2][NITROGEN]=30.34469;
        a[3][NITROGEN]=0.858927;
        b[3][NITROGEN]=0.656065;
        a[4][NITROGEN]=0.912985;
        b[4][NITROGEN]=0.217287;

        c[NITROGEN]=-11.804902;

        a[0][OXYGEN]=2.960427;
        b[0][OXYGEN]=14.182259;
        a[1][OXYGEN]=2.508818;
        b[1][OXYGEN]=5.936858;
        a[2][OXYGEN]=0.637853;
        b[2][OXYGEN]=0.112726;
        a[3][OXYGEN]=0.722838;
        b[3][OXYGEN]=34.958481;
        a[4][OXYGEN]=1.142756;
        b[4][OXYGEN]=0.39024;

        c[OXYGEN]=0.027014;

        a[0][SULFUR]=6.372157;
        b[0][SULFUR]=1.514347;
        a[1][SULFUR]=5.154568;
        b[1][SULFUR]=22.092528;
        a[2][SULFUR]=1.473732;
        b[2][SULFUR]=0.061373;
        a[3][SULFUR]=1.635073;
        b[3][SULFUR]=55.445176;
        a[4][SULFUR]=1.209372;
        b[4][SULFUR]=0.646925;

        c[SULFUR]=0.154722;

        a[0][IRON]=12.311098;
        b[0][IRON]=5.009415;
        a[1][IRON]=1.876623;
        b[1][IRON]=0.014461;
        a[2][IRON]=3.066177;
        b[2][IRON]=18.743041;
        a[3][IRON]=2.070451;
        b[3][IRON]=82.767874;
        a[4][IRON]=6.975185;
        b[4][IRON]=0.346506;

        c[IRON]=-0.304931;

        a[0][PHOSPHORUS]=1.950541;
        b[0][PHOSPHORUS]=0.908139;
        a[1][PHOSPHORUS]=4.14930;
        b[1][PHOSPHORUS]=27.044953;
        a[2][PHOSPHORUS]=1.494560;
        b[2][PHOSPHORUS]=0.071280;
        a[3][PHOSPHORUS]=1.522042;
        b[3][PHOSPHORUS]=67.520190;
        a[4][PHOSPHORUS]=5.729711;
        b[4][PHOSPHORUS]=1.981173;

        c[PHOSPHORUS]=0.155233;

        a[0][O_Minus]=3.106934;
        b[0][O_Minus]=19.86808;
        a[1][O_Minus]=3.235142;
        b[1][O_Minus]=6.960252;
        a[2][O_Minus]=1.148886;
        b[2][O_Minus]=0.170043;
        a[3][O_Minus]=0.783981;
        b[3][O_Minus]=65.693509;
        a[4][O_Minus]=0.676953;
        b[4][O_Minus]=0.630757;

        c[O_Minus]=0.046136;

        a[0][NA_Plus]=3.148690;
        b[0][NA_Plus]=2.594987;
        a[1][NA_Plus]=4.073989;
        b[1][NA_Plus]=6.046925;
        a[2][NA_Plus]=0.767888;
        b[2][NA_Plus]=0.070139;
        a[3][NA_Plus]=0.995612;
        b[3][NA_Plus]=14.122657;
        a[4][NA_Plus]=0.968249;
        b[4][NA_Plus]=0.217037;

        c[NA_Plus]=0.0445300;
        for (int m=0;m<NumAtomTypes;m++)
        {
                for (int n=0;n<5;n++)
                {
                        b[n][m]/=16.0*pi*pi;
                }
        }
}

#endif
