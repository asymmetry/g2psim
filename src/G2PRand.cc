// This file defines a namespace G2PRand.
// It is the rand of the whole simulation package.
// It uses TRandom2 as uniform random number generator. TRandom2 is fastest.
//
// History:
//   Jan 2013, C. Gu, First public version.
//   Feb 2013, C. Gu, Change it to a namespace
//

#include "TROOT.h"
#include "TRandom2.h"

#include "G2PRand.hh"

TRandom2 RG(0);

namespace G2PRand
{
    int Integer(int max)
    {
        return RG.Integer(max);
    }

    double Uniform()
    {
        return RG.Uniform();
    }

    double Uniform(double low, double high)
    {
        return low+(high-low)*RG.Uniform();
    }

    double Gaus(double mean = 0.0, double sigma = 1.0)
    {
        return RG.Gaus(mean, sigma);
    }

    double Linear(double a = 1.0, double c = 0.0 , double low = 0.0, double high = 1.0)
    {
        // Return random number in [low,high] following a*x+c prob density

        double x,y;

        double ylow = a*low+c;
        double yhigh = a*high+c;

        if ((ylow<0)||(yhigh<0)) return 0;

        do {
            x = low+(high-low)*Uniform();
            y = ylow+(yhigh-ylow)*Uniform();
        } while ((a*x+c)<y);

        return x;
    }

    double Rand1D(pf_Func1D f, double low, double high, double ylow, double yhigh)
    {
        double x,y;

        if ((ylow<0)||(yhigh<0)) return 0;

        do {
            x = low+(high-low)*Uniform();
            y = ylow+(yhigh-ylow)*Uniform();
        } while (f(x)<y);

        return x;
    }
} // end of namespace G2PRand
