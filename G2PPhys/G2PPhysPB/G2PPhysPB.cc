#include <cstdio>
#include <vector>

#include "G2PPhysBase.hh"

#include "G2PPhysPB.hh"

using namespace std;

extern "C" 
{
    // void bosted_(int *tgt_Z, int *tgt_A, double *Ei, double *Ep, double *ang, double *EPS, double *EPSD, double *FP, double *Tb, double *Ta, double *xs);
    double bostedxs_(double *Z, double* A, double *E, double *Ep, double* theta);
}

static const double kDEG = 3.14159265358979323846/180.0;
static const double kMEV = 1.0e-3;

static double PBosted(int Z, int A, double Ei, double Ef, double theta, double EPS, double EPSD, double FP, double Tb, double Ta)
{
    double NZ, NA;
    NZ = Z;
    NA = A;
    return bostedxs_(&NZ, &NA, &Ei, &Ef, &theta);
}

G2PPhysPB::G2PPhysPB() :
    fEPS(10.0), fEPSD(-10.0), fFP(220.0), fTb(0.0), fTa(0.0)
{
    // Nothing to do
}

G2PPhysPB::~G2PPhysPB()
{
    // Nothing to do
}

void G2PPhysPB::SetPars(double* array, int n)
{
    G2PPhysBase::SetPars(array, n);

    switch (n) {
    case 0:
        break;
    case 2:
        fTb = fPars[0]; fTa = fPars[1];
        break;
    case 3:
        fEPS = fPars[0]; fEPSD = fPars[1]; fFP = fPars[2];
        break;
    case 5:
        fEPS = fPars[0]; fEPSD = fPars[1]; fFP = fPars[2];
        fTb = fPars[3]; fTa = fPars[4];
        break;
    default:
        printf("Error: G2PPhysPB::SetPars(): Invalid number of pars.\n");
        break;
    }
}

double G2PPhysPB::GetXS(double Ei, double Ef, double theta)
{
    if (iPID==11) return PBosted(iZ, iA, Ei, Ef, theta, fEPS, fEPSD, fFP, fTb, fTa);
    else return -1;
}
