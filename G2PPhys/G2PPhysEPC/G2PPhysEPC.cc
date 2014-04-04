// -*- C++ -*-

/* class G2PPhysEPC
 * Class for EPC model.
 * Unit is ub/MeV-sr.
 * The cross section is calculated per nuclei. (Notice the difference with WISER model)
 * This model considers 2 situations: single pion production and multiple pion production.
 * 
 * How to set parameters:
 * If set 1 parameters with SetPars(pars,1), then pars[0]=0 means to calculate
 *   single pion production only, default is multiple pion production; 
 * Other uses will be considered as invalid.
 */

// History:
//   Feb 2014, C. Gu, First public version.
//   Apr 2014, C. Gu, Updated with multiple pion production
//

#include <cstdio>
#include <vector>
#include <cmath>

#include "G2PPhysBase.hh"

#include "G2PPhysEPC.hh"

using namespace std;

extern "C" {
void epc_(int* PART, int* Z, int* N, double* Ei, double* Pf, double* ang, double* xs, int* mpi);
}

static const double kDEG = 3.14159265358979323846 / 180.0;
static const double kMEV = 1.0e-3;

static double EPC(int PART, int Z, int A, double Ei, double Pf, double theta, int mpi)
{
    double XS;
    int N;
    N = A - Z;
    Ei = Ei / kMEV;
    Pf = Pf / kMEV;
    theta = theta / kDEG;
    epc_(&PART, &Z, &N, &Ei, &Pf, &theta, &XS, &mpi);
    return XS;
}

G2PPhysEPC::G2PPhysEPC() : fMPI(1)
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
    case 1:
        if (fPars[0] > 0.5) fMPI = 1;
        else fMPI = 0;
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
        return EPC(fPID, fZ, fA, Ei, Pf, theta, fMPI);
        break;
    default:
        return -1;
        break;
    }
}
