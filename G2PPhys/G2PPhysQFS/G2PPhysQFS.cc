#include <cstdio>
#include <vector>

#include "G2PPhysBase.hh"

#include "G2PPhysQFS.hh"

using namespace std;

extern "C" 
{
    void qfs_(int *tgt_Z, int *tgt_A, double *Ei, double *Ep, double *ang, double *EPS, double *EPSD, double *FP, double *Tb, double *Ta, double *xs);
}

static const double kDEG = 3.14159265358979323846/180.0;
static const double kMEV = 1.0e-3;

static double QFS(int Z, int A, double Ei, double Ef, double theta, double EPS, double EPSD, double FP, double Tb, double Ta)
{
    Ei = Ei/kMEV;
    Ef = Ef/kMEV;
    theta = theta/kDEG;
    double XS;
    qfs_(&Z, &A, &Ei, &Ef, &theta, &EPS, &EPSD, &FP, &Tb, &Ta, &XS);
    return XS;
}

G2PPhysQFS::G2PPhysQFS()
    :fEPS(10.0), fEPSD(-10.0), fFP(220.0), fTb(0.0), fTa(0.0)
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
        printf("Error: G2PPhysQFS::SetPars(): Invalid number of pars.\n");
        break;
    }
}

double G2PPhysQFS::GetXS(double Ei, double Ef, double theta)
{
    if (iPID==11) return QFS(iZ, iA, Ei, Ef, theta, fEPS, fEPSD, fFP, fTb, fTa);
    else return -1;
}


