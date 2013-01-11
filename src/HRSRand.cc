#include <cstdlib>
#include <cmath>

#include "TROOT.h"
#include "TRandom2.h"

#include "HRSRand.hh"

HRSRand::HRSRand()
{
    iRG = new TRandom2(0);
}

HRSRand::~HRSRand()
{
    delete[] iRG;
    iRG = NULL;
}

double HRSRand::Uniform()
{
    return iRG->Uniform();
}

double HRSRand::Uniform(double low, double high)
{
    return low+(high-low)*iRG->Uniform();
}

double HRSRand::Gaus(double mean = 0.0, double sigma = 1.0)
{
    return iRG->Gaus(mean, sigma);
}

double HRSRand::Linear(double a = 1.0, double c = 0.0 , double low = 0.0, double high = 1.0)
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

double HRSRand::Func(func f, double low, double high, double ylow, double yhigh)
{
    double x,y;

    do {
        x = low + (high-low)*Uniform();
        y = ylow + (yhigh-ylow)*Uniform();
    } while (f(x)<y);

    return x;
}


