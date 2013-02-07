// This file defines a class G2PRand.
// It is the rand of the whole simulation package.
// It uses TRandom2 as uniform random number generator. TRandom2 is fastest.
// All member functions are virtual so they can be overwrite.
//
// History:
//   Jan 2013, C. Gu, First public version.
//

#include <cstdlib>
#include <cmath>

#include "TROOT.h"
#include "TObject.h"
#include "TRandom2.h"

#include "G2PRand.hh"

ClassImp(G2PRand);

G2PRand::G2PRand()
{
    pRG = new TRandom2(0);
}

G2PRand::~G2PRand()
{
    delete[] pRG;
    pRG = NULL;
}

int G2PRand::Integer(int max)
{
    return pRG->Integer(max);
}

double G2PRand::Uniform()
{
    return pRG->Uniform();
}

double G2PRand::Uniform(double low, double high)
{
    return low+(high-low)*pRG->Uniform();
}

double G2PRand::Gaus(double mean = 0.0, double sigma = 1.0)
{
    return pRG->Gaus(mean, sigma);
}

double G2PRand::Linear(double a = 1.0, double c = 0.0 , double low = 0.0, double high = 1.0)
{
    // Return random number in [low,high] following a*x+c prob density
    double x,y;

    double ylow = a*low + c;
    double yhigh = a*high + c;

    if ((ylow<0)||(yhigh<0)) return 0;

    do {
        x = low + (high-low)*Uniform();
        y = ylow + (yhigh-ylow)*Uniform();
    } while ((a*x+c)<y);

    return x;
}

double G2PRand::Func(pf_Func1D f, double low, double high, double ylow, double yhigh)
{
    double x,y;

    do {
        x = low + (high-low)*Uniform();
        y = ylow + (yhigh-ylow)*Uniform();
    } while (f(x)<y);

    return x;
}


