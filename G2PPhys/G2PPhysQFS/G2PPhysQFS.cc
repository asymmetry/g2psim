// -*- C++ -*-

/* class G2PPhysQFS
 * Class for QFS model.
 * Unit is ub/MeV-sr.
 * Predict (e,e') cross sections to within 20% for an incident electron in the energy range 0.5-5 GeV and for energy losses greater than 50 MeV.
 * 
 * Radiative correction parameters:
 * Tb: total radiative length before scattering in radiation length;
 * Ta: total radiative length after scattering in radiation length;
 *
 * QFS model parameters:
 * EPS: separation energy in MeV;
 * EPSD: delta separation energy in MeV;
 * FP: Fermi momentum in MeV/c;
 */

// History:
//   Mar 2013, C. Gu, First public version.
//

#include <cstdio>
#include <vector>

#include "G2PPhysBase.hh"

#include "G2PPhysQFS.hh"

using namespace std;

extern "C" {
void qfs_(int* tgt_Z, int* tgt_A, double* Ei, double* Ep, double* ang, double* xs, double* EPS, double* EPSD, double* FP, double* Tb, double* Ta);
}

static const double kDEG = 3.14159265358979323846 / 180.0;
static const double kMEV = 1.0e-3;

static double QFS(int Z, int A, double Ei, double Ef, double theta, double EPS, double EPSD, double FP, double Tb, double Ta)
{
    Ei = Ei / kMEV;
    Ef = Ef / kMEV;
    theta = theta / kDEG;
    double XS;
    qfs_(&Z, &A, &Ei, &Ef, &theta, &XS, &EPS, &EPSD, &FP, &Tb, &Ta);
    return XS;
}

G2PPhysQFS::G2PPhysQFS() :
fEPS(10.0), fEPSD(-10.0), fFP(220.0), fTb(0.0), fTa(0.0)
{
    // Nothing to do
}

G2PPhysQFS::~G2PPhysQFS()
{
    // Nothing to do
}

void G2PPhysQFS::SetPars(double* array, int n)
{
    G2PPhysBase::SetPars(array, n);

    switch (n) {
    case 0:
        break;
    case 2:
        fTb = fPars[0];
        fTa = fPars[1];
        break;
    case 3:
        fEPS = fPars[0];
        fEPSD = fPars[1];
        fFP = fPars[2];
        break;
    case 5:
        fEPS = fPars[0];
        fEPSD = fPars[1];
        fFP = fPars[2];
        fTb = fPars[3];
        fTa = fPars[4];
        break;
    default:
        printf("Error: G2PPhysQFS::SetPars(): Invalid number of pars.\n");
        break;
    }
}

double G2PPhysQFS::GetXS(double Ei, double Ef, double theta)
{
    if (fPID == 11) return QFS(fZ, fA, Ei, Ef, theta, fEPS, fEPSD, fFP, fTb, fTa);
    else return -1;
}


