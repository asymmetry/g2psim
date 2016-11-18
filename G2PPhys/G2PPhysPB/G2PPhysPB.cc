// -*- C++ -*-

/* class G2PPhysPB
 * Class for P. Bosted model.
 * Unit is ub/MeV-sr.
 * Valid for all W<3 GeV and all Q2<10 GeV2.
 *
 * Radiative correction parameters:
 * [1] Tb: total radiative length before scattering in radiation length;
 * [2] Ta: total radiative length after scattering in radiation length.
 */

// History:
//   Mar 2013, C. Gu, First public version.
//

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>

#include "G2PPhysBase.hh"

#include "G2PPhysPB.hh"

using namespace std;

extern "C" {
    void bosted_(double *Z, double *A, double *Ei, double *Ep, double *ang, double *xs, double *Tb, double *Ta);
}

static const double kDEG = 3.14159265358979323846 / 180.0;
static const double kMEV = 1.0e-3;

static double PBosted(int Z, int A, double Ei, double Ef, double theta, double Tb, double Ta)
{
    double NZ, NA;
    NZ = Z;
    NA = A;
    double XS;
    bosted_(&NZ, &NA, &Ei, &Ef, &theta, &XS, &Tb, &Ta);
    return XS;
}

G2PPhysPB::G2PPhysPB() : fTb(0.0), fTa(0.0)
{
    // Nothing to do
}

G2PPhysPB::~G2PPhysPB()
{
    // Nothing to do
}

void G2PPhysPB::SetPar(int id, double value)
{
    switch (id) {
    case 1:
        fTb = value;
        break;

    case 2:
        fTa = value;
        break;

    default:
        printf("Error: G2PPhysPB::SetPars(): Invalid parameter id.\n");
        break;
    }
}

double G2PPhysPB::GetXS(double Ei, double Ef, double theta)
{
    if (fPID == 11) return PBosted(fZ, fA, Ei, Ef, theta, fTb, fTa);
    else return -1;
}
