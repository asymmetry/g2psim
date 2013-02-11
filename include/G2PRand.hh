// This file defines a namespace G2PRand.
// It is the rand of the whole simulation package.
// It uses TRandom2 as uniform random number generator. TRandom2 is fastest.
//
// History:
//   Jan 2013, C. Gu, First public version.
//   Feb 2013, C. Gu, Change it to a namespace
//

#ifndef G2P_RAND_H
#define G2P_RAND_H

namespace G2PRand
{    
    typedef double (*pf_Func1D)(double);
    
    int Integer(int max);
    double Uniform();
    double Uniform(double low, double high);
    double Gaus(double mean, double sigma);
    double Linear(double a, double c, double low, double high);
    double Rand1D(pf_Func1D f, double low, double high, double ylow, double yhigh);
}

#endif
