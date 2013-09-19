// -*- C++ -*-

/* class G2PFwdProc
 * It simulates the movement of the scatted particles in the spectrometers.
 * G2PDrift, G2PHRS and G2PSieve are used in this class.
 * Input variables: fV5tg_tr, fV5react_lab (register in gG2PVars).
 */

// History:
//   Apr 2013, C. Gu, First public version.
//

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"

#include "G2PAppBase.hh"
#include "G2PAppList.hh"
#include "G2PDrift.hh"
#include "G2PGlobals.hh"
#include "G2PHRS.hh"
#include "G2PProcBase.hh"
#include "G2PRand.hh"
#include "G2PSieve.hh"
#include "G2PVar.hh"
#include "G2PVarDef.hh"
#include "G2PVarList.hh"

#include "G2PFwdProc.hh"

using namespace std;

G2PFwdProc* G2PFwdProc::pG2PFwdProc = NULL;

G2PFwdProc::G2PFwdProc() :
fHRSAngle(0.0), fHRSMomentum(0.0), fSieveOn(false), pDrift(NULL), pHRS(NULL), pSieve(NULL) {
    if (pG2PFwdProc) {
        Error("G2PFwdProc()", "Only one instance of G2PFwdProc allowed.");
        MakeZombie();
        return;
    }
    pG2PFwdProc = this;

    fPriority = 3;

    Clear();
}

G2PFwdProc::~G2PFwdProc() {
    if (pG2PFwdProc == this) pG2PFwdProc = NULL;
}

int G2PFwdProc::Init() {
    //static const char* const here = "Init()";

    if (G2PProcBase::Init() != 0) return fStatus;

    pDrift = static_cast<G2PDrift*> (gG2PApps->Find("G2PDrift"));
    if (!pDrift) {
        pDrift = new G2PDrift();
        gG2PApps->Add(pDrift);
    }

    pSieve = static_cast<G2PSieve*> (gG2PApps->Find("G2PSieve"));
    if (!pSieve) {
        pSieve = new G2PSieve();
        gG2PApps->Add(pSieve);
    }

    pHRS = static_cast<G2PHRS*> (gG2PApps->Find("G2PHRS"));
    if (!pHRS) return (fStatus == kINITERROR);

    return (fStatus = kOK);
}

int G2PFwdProc::Begin() {
    //static const char* const here = "Begin()";

    if (G2PProcBase::Begin() != 0) return fStatus;

    return (fStatus = kOK);
}

int G2PFwdProc::Process() {
    static const char* const here = "Process()";

    double V5react_lab[5], V5tg_tr[5];

    V5react_lab[0] = gG2PVars->FindSuffix("react.l_x")->GetValue();
    V5react_lab[1] = gG2PVars->FindSuffix("react.l_t")->GetValue();
    V5react_lab[2] = gG2PVars->FindSuffix("react.l_y")->GetValue();
    V5react_lab[3] = gG2PVars->FindSuffix("react.l_p")->GetValue();
    V5react_lab[4] = gG2PVars->FindSuffix("react.l_z")->GetValue();

    V5tg_tr[0] = gG2PVars->FindSuffix("tp.x")->GetValue();
    V5tg_tr[1] = gG2PVars->FindSuffix("tp.t")->GetValue();
    V5tg_tr[2] = gG2PVars->FindSuffix("tp.y")->GetValue();
    V5tg_tr[3] = gG2PVars->FindSuffix("tp.p")->GetValue();
    V5tg_tr[4] = gG2PVars->FindSuffix("tp.d")->GetValue();

    double x[3] = {V5react_lab[0], V5react_lab[2], V5react_lab[4]};
    double pp = fHRSMomentum * (1 + V5tg_tr[4]);
    double p[3] = {pp * sin(V5react_lab[1]) * cos(V5react_lab[3]),
                   pp * sin(V5react_lab[1]) * sin(V5react_lab[3]),
                   pp * cos(V5react_lab[1])};
    // Local dump front face
    pDrift->Drift(x, p, 640.0e-3, 10.0, x, p);
    if ((fabs(x[0]) < 46.0e-3) || (fabs(x[0]) > 87.0e-3)) return -1;
    if ((x[1]<-43.0e-3) || (x[1] > 50.0e-3)) return -1;
    // Local dump back face
    pDrift->Drift(x, p, 790.0e-3, 10.0, x, p);
    if ((fabs(x[0]) < 58.0e-3) || (fabs(x[0]) > 106.0e-3)) return -1;
    if ((x[1]<-53.0e-3) || (x[1] > 58.0e-3)) return -1;

    pDrift->Drift(V5tg_tr, fHRSMomentum, 0.0, fHRSAngle, pSieve->GetZ(), 10.0, fV5sieve_tr);

    if (fDebug > 1) {
        Info(here, "sieve_tr  : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5sieve_tr[0], fV5sieve_tr[1], fV5sieve_tr[2], fV5sieve_tr[3], fV5sieve_tr[4]);
    }

    if (fSieveOn) {
        if (!pSieve->CanPass(fV5sieve_tr)) return -1;
    }

    Project(fV5sieve_tr[0], fV5sieve_tr[2], pSieve->GetZ(), 0.0, fV5sieve_tr[1], fV5sieve_tr[3], fV5projtg_tr[0], fV5projtg_tr[2]);
    fV5projtg_tr[1] = fV5sieve_tr[1];
    fV5projtg_tr[3] = fV5sieve_tr[3];
    fV5projtg_tr[4] = fV5sieve_tr[4];

    if (fDebug > 1) {
        Info(here, "projtg_tr : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5projtg_tr[0], fV5projtg_tr[1], fV5projtg_tr[2], fV5projtg_tr[3], fV5projtg_tr[4]);
    }

    if (!pHRS->Forward(fV5projtg_tr, fV5fp_tr)) return -1;
    ApplyVDCRes(fV5fp_tr);
    TRCS2FCS(fV5fp_tr, fHRSAngle, fV5fp_rot);

    if (fDebug > 1) {
        Info(here, "fp_tr     : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5fp_tr[0], fV5fp_tr[1], fV5fp_tr[2], fV5fp_tr[3], fV5fp_tr[4]);
    }

    return 0;
}

void G2PFwdProc::Clear() {
    memset(fV5sieve_tr, 0, sizeof (fV5sieve_tr));
    memset(fV5projtg_tr, 0, sizeof (fV5projtg_tr));
    memset(fV5fp_tr, 0, sizeof (fV5fp_tr));
    memset(fV5fp_rot, 0, sizeof (fV5fp_rot));
}

void G2PFwdProc::ApplyVDCRes(double* V5fp) {
    double WireChamberResX = 0.0013; //m;
    double WireChamberResY = 0.0013; //m;
    double WireChamberResT = 0.0003; //rad;
    double WireChamberResP = 0.0003; //rad;

    V5fp[0] = pRand->Gaus(V5fp[0], WireChamberResX);
    V5fp[1] = pRand->Gaus(V5fp[1], WireChamberResT);
    V5fp[2] = pRand->Gaus(V5fp[2], WireChamberResY);
    V5fp[3] = pRand->Gaus(V5fp[3], WireChamberResP);
}

int G2PFwdProc::Configure(EMode mode) {
    if (mode == kREAD || mode == kTWOWAY) {
        if (fIsInit) return 0;
        else fIsInit = true;
    }

    ConfDef confs[] = {
        {"run.hrs.angle", "HRS Angle", kDOUBLE, &fHRSAngle},
        {"run.hrs.p0", "HRS Momentum", kDOUBLE, &fHRSMomentum},
        {"run.sieve.on", "Sieve On", kBOOL, &fSieveOn},
        {0}
    };

    return ConfigureFromList(confs, mode);
}

int G2PFwdProc::DefineVariables(EMode mode) {
    if (mode == kDEFINE && fIsSetup) return 0;
    fIsSetup = (mode == kDEFINE);

    VarDef vars[] = {
        {"sieve.x", "Sieve X", kDOUBLE, &fV5sieve_tr[0]},
        {"sieve.t", "Sieve T", kDOUBLE, &fV5sieve_tr[1]},
        {"sieve.y", "Sieve Y", kDOUBLE, &fV5sieve_tr[2]},
        {"sieve.p", "Sieve P", kDOUBLE, &fV5sieve_tr[3]},
        {"tp.proj.x", "Project to target plane X", kDOUBLE, &fV5projtg_tr[0]},
        {"tp.proj.t", "Project to target plane T", kDOUBLE, &fV5projtg_tr[1]},
        {"tp.proj.y", "Project to target plane Y", kDOUBLE, &fV5projtg_tr[2]},
        {"tp.proj.p", "Project to target plane P", kDOUBLE, &fV5projtg_tr[3]},
        {"fp.x", "Focus plane X", kDOUBLE, &fV5fp_tr[0]},
        {"fp.t", "Focus plane T", kDOUBLE, &fV5fp_tr[1]},
        {"fp.y", "Focus plane Y", kDOUBLE, &fV5fp_tr[2]},
        {"fp.p", "Focus plane P", kDOUBLE, &fV5fp_tr[3]},
        {"fp.r_x", "Focus plane X (rot)", kDOUBLE, &fV5fp_rot[0]},
        {"fp.r_t", "Focus plane T (rot)", kDOUBLE, &fV5fp_rot[1]},
        {"fp.r_y", "Focus plane Y (rot)", kDOUBLE, &fV5fp_rot[2]},
        {"fp.r_p", "Focus plane P (rot)", kDOUBLE, &fV5fp_rot[3]},
        {0}
    };

    return DefineVarsFromList(vars, mode);
}

void G2PFwdProc::MakePrefix() {
    const char* basename = "fwd";

    G2PAppBase::MakePrefix(basename);
}

ClassImp(G2PFwdProc)
