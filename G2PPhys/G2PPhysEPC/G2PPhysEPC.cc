// -*- C++ -*-

/* class G2PPhysEPC
 * Class for EPC model.
 * Unit is ub/MeV-sr.
 * predict (e,N) cross sections to within a factor of 2 for an incident electron in the energy range 0.5-5 GeV and for nucleon kinetic energies greater than 50 MeV.
 */

// History:
//   Feb 2014, C. Gu, First public version.
//

#include <cstdio>
#include <vector>
#include <cmath>

#include "G2PPhysBase.hh"

#include "G2PPhysEPC.hh"

using namespace std;

extern "C" {
void epc_(int* PART, int* Z, int* N, double* Ei, double* Pf, double* ang, double* xs);
}

static const double kDEG = 3.14159265358979323846 / 180.0;
static const double kMEV = 1.0e-3;

static double EPC(int Z, int A, int PART, double Ei, double Pf, double theta)
{
    double XS;
    int N;
    N = A - Z;
    Ei = Ei / kMEV;
    Pf = Pf / kMEV;
    theta = theta / kDEG;
    epc_(&PART, &Z, &N, &Ei, &Pf, &theta, &XS);
    return XS;
}

G2PPhysEPC::G2PPhysEPC()
{
    // Nothing to do
}

G2PPhysEPC::~G2PPhysEPC()
{
    // Nothing to do
}

void G2PPhysEPC::SetPars(double* array, int n)
{
    G2PPhysBase::SetPars(array, n);

    switch (n) {
    case 0:
        break;
    default:
        printf("Error: G2PPhysEPC::SetPars(): Invalid number of pars.\n");
        break;
    }
}

double G2PPhysEPC::GetXS(double Ei, double Pf, double theta)
{
    switch (fPID) {
    case 211: // pi+
    case -211: // pi-
    case 111: // pi0
    case 2212: // p
    case 2112: // n
        return EPC(fPID, fZ, fA, Ei, Pf, theta);
        break;
    default:
        return -1;
        break;
    }
}
