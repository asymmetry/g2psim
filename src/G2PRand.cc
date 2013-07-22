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

#include <cstdlib>
#include <cmath>

#include "TROOT.h"
#include "TObject.h"
#include "TError.h"
#include "TRandom2.h"

#include "G2PRand.hh"

G2PRand* G2PRand::pG2PRand = NULL;

G2PRand::G2PRand() {
    pRG = new TRandom2(0);
}

G2PRand::~G2PRand() {
    delete pRG;
}

G2PRand* G2PRand::GetInstance() {
    if (!pG2PRand) {
        static G2PRand instance;
        pG2PRand = &instance;
    }
    return pG2PRand;
}

int G2PRand::Integer(int max) {
    return pRG->Integer(max);
}

double G2PRand::Uniform() {
    return pRG->Uniform();
}

double G2PRand::Uniform(double low, double high) {
    return low + (high - low) * pRG->Uniform();
}

double G2PRand::Gaus(double mean = 0.0, double sigma = 1.0) {
    return pRG->Gaus(mean, sigma);
}

double G2PRand::Linear(double a = 1.0, double c = 0.0, double low = 0.0, double high = 1.0) {
    // Return random number in [low,high] following a*x+c prob density
    double x, y;

    double ylow = a * low + c;
    double yhigh = a * high + c;

    if ((ylow < 0) || (yhigh < 0)) return 0;

    do {
        x = low + (high - low) * Uniform();
        y = ylow + (yhigh - ylow) * Uniform();
    } while ((a * x + c) < y);

    return x;
}

double G2PRand::Rand1D(pfFunc1D_ f, double low, double high, double ylow, double yhigh) {
    double x, y;

    if ((ylow < 0) || (yhigh < 0)) return 0;

    do {
        x = low + (high - low) * Uniform();
        y = ylow + (yhigh - ylow) * Uniform();
    } while (f(x) < y);

    return x;
}

ClassImp(G2PRand)
