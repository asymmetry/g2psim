#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "TROOT.h"
#include "TObject.h"
#include "TError.h"

#include "G2PAppsBase.hh"
#include "G2PGlobals.hh"
#include "G2PDrift.hh"
#include "G2PRand.hh"
#include "G2PRunBase.hh"

#include "G2PBPM.hh"

G2PBPM* G2PBPM::pG2PBPM = NULL;

G2PBPM::G2PBPM() :
    fBeamEnergy(2.254), fBPMZ(0.0), fBPMAZ(-956.0e-3),
    fBPMBZ(-690.0e-3), fBPMAX(0.0), fBPMAY(0.0),
    fBPMBX(0.0), fBPMBY(0.0), fBPMARes(0.3e-3), fBPMBRes(0.3e-3),
    pDrift(NULL)
{
    if (pG2PBPM) {
        Error("G2PBPM()", "Only one instance of G2PBPM allowed.");
        MakeZombie();
        return;
    }
    pG2PBPM = this;
}

G2PBPM::~G2PBPM()
{
    if (pG2PBPM==this) pG2PBPM = NULL;
}

G2PAppsBase::EStatus G2PBPM::Init()
{
    static const char* const here = "Init()";

    if (G2PAppsBase::Init()) return fStatus;

    fBeamEnergy = gG2PRun->GetBeamEnergy();

    if (SetBPM()) {
        Error(here, "Cannot initialize.");
        return (fStatus = kINITERROR);
    }

    return (fStatus = kOK);
}

void G2PBPM::GetBPMValue(const double* V5beam_lab, double* V5bpm_lab)
{
    static const char* const here = "GetBPMValue()";

    double x[3] = { V5beam_lab[0], V5beam_lab[2], V5beam_lab[4] };
    double p[3] = { fBeamEnergy*sin(V5beam_lab[1])*cos(V5beam_lab[3]),
                    fBeamEnergy*sin(V5beam_lab[1])*sin(V5beam_lab[3]),
                    fBeamEnergy*cos(V5beam_lab[1]) };

    pDrift->Drift(x, p, fBPMAZ, 10.0, x, p);
    double res = fBPMBRes/(fBPMBZ-fBPMAZ);
    // here theta phi are special coords defined by Pengjia
    // theta is dy/dz in hall coords
    // phi is dx/dz in hall coords
    double theta = pRand->Gaus(atan2(p[1], p[2]), res);
    double phi = pRand->Gaus(atan2(p[0], p[2]), res);
    p[1] = p[2]*tan(theta);
    p[0] = p[2]*tan(phi);
    double normF = fBeamEnergy/sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
    p[0] *= normF; p[1] *= normF; p[2] *= normF;
    x[0] = pRand->Gaus(x[0], fBPMARes);
    x[1] = pRand->Gaus(x[1], fBPMARes);
    pDrift->Drift(x, p, 0.0, 10.0, x, p);
    V5bpm_lab[0] = x[0];
    V5bpm_lab[1] = atan2(p[1], p[2]);
    V5bpm_lab[2] = x[1];
    V5bpm_lab[3] = atan2(p[0], p[2]);
    V5bpm_lab[4] = x[2];

    if (fDebug>1) Info(here, "%10.3e %10.3e %10.3e %10.3e %10.3e", V5bpm_lab[0], V5bpm_lab[1], V5bpm_lab[2], V5bpm_lab[3], V5bpm_lab[4]);
}

void G2PBPM::TransBPM2Lab(const double* V5_bpm, double* V5_lab)
{
    V5_lab[0] = V5_bpm[0];
    V5_lab[2] = V5_bpm[2];
    V5_lab[4] = V5_bpm[4];

    double p[3] = { tan(V5_bpm[3]), tan(V5_bpm[1]), 1.0 };
    double pp = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
    V5_lab[1] = acos(1.0/pp);
    V5_lab[3] = atan2(p[1], p[0]);
}

int G2PBPM::RegisterModel()
{
    pDrift = G2PDrift::GetInstance();
    if (!pDrift) {
        pDrift = new G2PDrift();
        gG2PApps->Add(pDrift);
    }
    fApps->Add(pDrift);

    return 0;
}

int G2PBPM::SetBPM()
{
    static const char* const here = "SetBPM()";

    fBPMAZ = -956e-3;
    fBPMBZ = -690e-3;
    fBPMZ = (fBPMAZ+fBPMBZ)/2.0;

    if (fDebug>0) Info(here, "Done!");

    return 0;
}

ClassImp(G2PBPM)
