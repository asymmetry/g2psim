#include "TROOT.h"
#include "TMath.h"

#include "G2PPModsQFS.hh"

extern "C" 
{
    void qfs_(int *tgt_Z, int *tgt_A, double *Ei, double *Ep, double *ang, double *EPS, double *EPSD, double *FP, double *Tb, double *Ta, double *xs);
}

const double kDEG = TMath::Pi()/180.0;
const double kMEV = 1.0e-3;

G2PPModsQFS::G2PPModsQFS()
    :fEPS(10.0), fEPSD(-10.0), fFP(220.0), fTb(0.0), fTa(0.0)
{
    // Nothing to do
}

G2PPModsQFS::~G2PPModsQFS()
{
    // Nothing to do
}

void G2PPModsQFS::operator() (int Z, int A, double Eb, double Ef, double theta, double *XS)
{
    Eb=Eb/kMEV;
    Ef=Ef/kMEV;
    theta=theta/kDEG;
    qfs_(&Z, &A, &Eb, &Ef, &theta, &fEPS, &fEPSD, &fFP, &fTb, &fTa, XS);
}

void G2PPModsQFS::operator() (int Z, int A, double Eb, double Ef, double theta, double EPS, double EPSD, double FP, double *XS)
{
    Eb=Eb/kMEV;
    Ef=Ef/kMEV;
    theta=theta/kDEG;
    qfs_(&Z, &A, &Eb, &Ef, &theta, &EPS, &EPSD, &FP, &fTb, &fTa, XS);
}

void G2PPModsQFS::operator() (int Z, int A, double Eb, double Ef, double theta, double Tb, double Ta, double *XS)
{
    Eb=Eb/kMEV;
    Ef=Ef/kMEV;
    theta=theta/kDEG;
    qfs_(&Z, &A, &Eb, &Ef, &theta, &fEPS, &fEPSD, &fFP, &Tb, &Ta, XS);
}

void G2PPModsQFS::operator() (int Z, int A, double Eb, double Ef, double theta, double EPS, double EPSD, double FP, double Tb, double Ta, double *XS)
{
    Eb=Eb/kMEV;
    Ef=Ef/kMEV;
    theta=theta/kDEG;
    qfs_(&Z, &A, &Eb, &Ef, &theta, &EPS, &EPSD, &FP, &Tb, &Ta, XS);
}
