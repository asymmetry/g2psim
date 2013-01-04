#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "TRandom2.h"

#include "Rand.hh"

TRandom2 iRG(0);

double fRand()
{
    return iRG.Uniform();
}

double fGausRand(double mean = 0.0, double sigma = 1.0)
{
    return iRG.Gaus(mean, sigma);
}

double fLinearRand(double a = 1.0, double c = 0.0 , double low = 0.0, double high = 1.0)
{
    // Return random number in [low,high] following a*x+c prob density
    double x,y;

    double ylow = a*low + c;
    double yhigh = a*high + c;

    if ((ylow<0)||(yhigh<0)) return 0;

    do {
        x = low + (high-low)*fRand();
        y = ylow + (yhigh-ylow)*fRand();
    } while ((a*x+c)<y);

    return x;
}
