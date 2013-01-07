#ifndef HRS_RAND_H
#define HRS_RAND_H

typedef double (*func)(double);

double fRand();
double fGausRand(double mean, double sigma);
double fLinearRand(double a, double c, double low, double high);
double fFuncRand(func f, double low, double high, double ylow, double yhigh)

#endif
