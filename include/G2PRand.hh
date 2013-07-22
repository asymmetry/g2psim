// -*- C++ -*-

/* class G2PRand
 * This file defines a class G2PRand.
 * It is the random number generator of the simulation package.
 * It uses TRandom2 class as uniform random number generator. TRandom2 is fastest in ROOT random number generators.
 */

// History:
//   Jan 2013, C. Gu, First public version.
//   Feb 2013, C. Gu, Change it to a class
//

#ifndef G2P_RAND_H
#define G2P_RAND_H

#include "TObject.h"
#include "TRandom2.h"

class G2PRand : public TObject {
public:
    static G2PRand* GetInstance();

    typedef double (*pfFunc1D_)(double);

    void SetSeed(int n) {
        pRG->SetSeed(n);
    }

    int Integer(int max);
    double Uniform();
    double Uniform(double low, double high);
    double Gaus(double mean, double sigma);
    double Linear(double a, double c, double low, double high);
    double Rand1D(pfFunc1D_ f, double low, double high, double ylow, double yhigh);

private:
    G2PRand(); // Only allow one instance
    ~G2PRand();

    TRandom2* pRG;

    static G2PRand* pG2PRand;

    ClassDef(G2PRand, 1)
};

#endif
