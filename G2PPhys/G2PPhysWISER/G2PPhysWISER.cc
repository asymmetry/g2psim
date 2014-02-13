// -*- C++ -*-

/* class G2PPhysWISER
 * Class for WISER model.
 * Unit is ub/MeV-sr.
 * Photoproduction of pion/nucleons in DIS region.
 */

// History:
//   Apr 2013, C. Gu, First public version.
//

#include <cstdio>
#include <vector>
#include <cmath>

#include "G2PPhysBase.hh"

#include "G2PPhysWISER.hh"

using namespace std;

extern "C" {
void wiser_(int* Z, int* N, int* PART, double* Ei, double* Ep, double* ang, double* radlen, double* xs);
}

static const double kDEG = 3.14159265358979323846 / 180.0;
static const double kMEV = 1.0e-3;

static double WISER(int Z, int A, int PART, double Ei, double Ef, double theta, double radlen)
{
    double XS;
    int N;
    N = A - Z;
    wiser_(&Z, &N, &PART, &Ei, &Ef, &theta, &radlen, &XS);
    return XS;
}

G2PPhysWISER::G2PPhysWISER() :
fRadLen(0.0)
{
    // Nothing to do
}

G2PPhysWISER::~G2PPhysWISER()
{
    // Nothing to do
}

void G2PPhysWISER::SetPars(double* array, int n)
{
    G2PPhysBase::SetPars(array, n);

    switch (n) {
    case 0:
        break;
    case 1:
        fRadLen = fPars[0];
        break;
    default:
        printf("Error: G2PPhysWISER::SetPars(): Invalid number of pars.\n");
        break;
    }
}

double G2PPhysWISER::GetXS(double Ei, double Pf, double theta)
{
    if (fabs(fRadLen) < 1e-8) {
        printf("Error: G2PPhysWISER::GetXS(): Must give radiation lenth.\n");
        return -1;
    }
    switch (fPID) {
    case 211: // pi+
        return WISER(fZ, fA, 1, Ei, Pf, theta, fRadLen);
        break;
    case -211: // pi-
        return WISER(fZ, fA, 2, Ei, Pf, theta, fRadLen);
        break;
    case 321: // k+
        return WISER(fZ, fA, 3, Ei, Pf, theta, fRadLen);
        break;
    case -321: // k-
        return WISER(fZ, fA, 4, Ei, Pf, theta, fRadLen);
        break;
    case 2212: // p
        return WISER(fZ, fA, 5, Ei, Pf, theta, fRadLen);
        break;
    default:
        return -1;
        break;
    }
}
