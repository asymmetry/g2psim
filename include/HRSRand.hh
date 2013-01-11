#ifndef HRS_RAND_H
#define HRS_RAND_H

#include "TRandom2.h"

typedef double (*func)(double);

class HRSRand
{
public:
    HRSRand();
    ~HRSRand();
    
    double Uniform();
    double Uniform(double low, double high);
    double Gaus(double mean, double sigma);
    double Linear(double a, double c, double low, double high);
    double Func(func f, double low, double high, double ylow, double yhigh);

private:
    TRandom2 *iRG;
};

#endif
