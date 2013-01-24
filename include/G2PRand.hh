// This file defined a class G2PRand.
// It is the rand of the whole simulation package.
// It use TRandom2 as Uniform random number generator. TRandom2 is fastest.
// All member functions are virtual so they can be overwrite.
//
// History:
// By C. Gu, Jan 12, 2013, First public version.
//

#ifndef G2P_RAND_H
#define G2P_RAND_H

#include "TROOT.h"
#include "TObject.h"
#include "TRandom2.h"

class G2PRand : public TObject
{
public:
    G2PRand();
    ~G2PRand();
    
    typedef double (*pf_Func1D)(double);

    virtual int Integer(int max);
    virtual double Uniform();
    virtual double Uniform(double low, double high);
    virtual double Gaus(double mean, double sigma);
    virtual double Linear(double a, double c, double low, double high);
    virtual double Func(pf_Func1D f, double low, double high, double ylow, double yhigh);

private:
    TRandom2 *pRG;
    ClassDef(G2PRand,1);
};

#endif
