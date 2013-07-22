// -*- C++ -*-

/* class G2PFwdProc
 * This file defines a class G2PFwdProc.
 * It simulates the transport of the scatted particles in the spectrometers.
 * G2PDBRec, G2PDrift, G2PHRSTrans are used in this class.
 * Input variables: fV5tg_tr, fV5react_lab (register in G2PRun).
 */

// History:
//   Apr 2013, C. Gu, First public version.
//

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <map>

#include "TROOT.h"
#include "TObject.h"
#include "TError.h"

#include "G2PAppBase.hh"
#include "G2PDBRec.hh"
#include "G2PDrift.hh"
#include "G2PGlobals.hh"
#include "G2PHRSTrans.hh"
#include "G2PProcBase.hh"
#include "G2PRunBase.hh"
#include "G2PSieve.hh"
#include "G2PVarDef.hh"
#include "G2PVarList.hh"

#include "G2PFwdProc.hh"

G2PFwdProc::G2PFwdProc() :
fHRSAngle(0.0), fHRSMomentum(0.0), pDrift(NULL), pHRS(NULL), pDBRec(NULL) {
    mName["fV5tg_tr"] = fV5tg_tr;
    mLength["fV5tg_tr"] = 5;
    mName["fV5react_lab"] = fV5react_lab;
    mLength["fV5react_lab"] = 5;
    mName["fV5sieve_tr"] = fV5sieve_tr;
    mLength["fV5sieve_tr"] = 5;
    mName["fV5projtg_tr"] = fV5projtg_tr;
    mLength["fV5projtg_tr"] = 5;
    mName["fV5fp_tr"] = fV5fp_tr;
    mLength["fV5fp_tr"] = 5;
    mName["fV5fp_rot"] = fV5fp_rot;
    mLength["fV5fp_rot"] = 5;

    fAppsList.push_back("G2PDBRec");
    fAppsList.push_back("G2PDrift");
    fAppsList.push_back("G2PHRSTrans");

    Clear();
}

G2PFwdProc::~G2PFwdProc() {
    // Nothing to do
}

int G2PFwdProc::Init() {
    //static const char* const here = "Init()";

    if (G2PProcBase::Init() != 0) return fStatus;

    pDrift = G2PDrift::GetInstance();
    pHRS = G2PHRSTrans::GetInstance();
    pDBRec = G2PDBRec::GetInstance();

    fApps->Add(pDrift);
    fApps->Add(pHRS);
    fApps->Add(pDBRec);

    return (fStatus = kOK);
}

int G2PFwdProc::Begin() {
    //static const char* const here = "Begin()";

    if (G2PProcBase::Begin() != 0) return fStatus;

    fHRSAngle = gG2PRun->GetHRSAngle();
    fHRSMomentum = gG2PRun->GetHRSMomentum();

    SetSieve(fHRSAngle);

    return (fStatus = kOK);
}

int G2PFwdProc::Process() {
    static const char* const here = "Process()";

    double x[3] = {fV5react_lab[0], fV5react_lab[2], fV5react_lab[4]};
    double pp = fHRSMomentum * (1 + fV5tg_tr[4]);
    double p[3] = {pp * sin(fV5react_lab[1]) * cos(fV5react_lab[3]),
                   pp * sin(fV5react_lab[1]) * sin(fV5react_lab[3]),
                   pp * cos(fV5react_lab[1])};
    // Local dump front face
    pDrift->Drift(x, p, 640.0e-3, 10.0, x, p);
    if ((fabs(x[0]) < 46.0e-3) || (fabs(x[0]) > 87.0e-3)) return -1;
    if ((x[1]<-43.0e-3) || (x[1] > 50.0e-3)) return -1;
    // Local dump back face
    pDrift->Drift(x, p, 790.0e-3, 10.0, x, p);
    if ((fabs(x[0]) < 58.0e-3) || (fabs(x[0]) > 106.0e-3)) return -1;
    if ((x[1]<-53.0e-3) || (x[1] > 58.0e-3)) return -1;

    pDrift->Drift(fV5tg_tr, fHRSMomentum, 0.0, fHRSAngle, fSieve.fZ, 10.0, fV5sieve_tr);

    if (fDebug > 1) {
        Info(here, "sieve_tr  : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5sieve_tr[0], fV5sieve_tr[1], fV5sieve_tr[2], fV5sieve_tr[3], fV5sieve_tr[4]);
    }

    Project(fV5sieve_tr[0], fV5sieve_tr[2], fSieve.fZ, 0.0, fV5sieve_tr[1], fV5sieve_tr[3], fV5projtg_tr[0], fV5projtg_tr[2]);
    fV5projtg_tr[1] = fV5sieve_tr[1];
    fV5projtg_tr[3] = fV5sieve_tr[3];
    fV5projtg_tr[4] = fV5sieve_tr[4];

    if (fDebug > 1) {
        Info(here, "projtg_tr : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5projtg_tr[0], fV5projtg_tr[1], fV5projtg_tr[2], fV5projtg_tr[3], fV5projtg_tr[4]);
    }

    if (!pHRS->Forward(fV5projtg_tr, fV5fp_tr)) return -1;
    ApplyVDCRes(fV5fp_tr);
    pDBRec->TransTr2Rot(fV5fp_tr, fV5fp_rot);

    if (fDebug > 1) {
        Info(here, "fp_tr     : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5fp_tr[0], fV5fp_tr[1], fV5fp_tr[2], fV5fp_tr[3], fV5fp_tr[4]);
    }

    return 0;
}

void G2PFwdProc::Clear() {
    memset(fV5tg_tr, 0, sizeof (fV5tg_tr));
    memset(fV5react_lab, 0, sizeof (fV5react_lab));
    memset(fV5sieve_tr, 0, sizeof (fV5sieve_tr));
    memset(fV5projtg_tr, 0, sizeof (fV5projtg_tr));
    memset(fV5fp_tr, 0, sizeof (fV5fp_tr));
    memset(fV5fp_rot, 0, sizeof (fV5fp_rot));
}

int G2PFwdProc::DefineVariables(EMode mode) {
    if (mode == kDefine && bIsSetup) return 0;
    bIsSetup = (mode == kDefine);

    VarDef vars[] = {
        {"sieve.x", "Sieve X", kDouble, &fV5sieve_tr[0]},
        {"sieve.t", "Sieve T", kDouble, &fV5sieve_tr[1]},
        {"sieve.y", "Sieve Y", kDouble, &fV5sieve_tr[2]},
        {"sieve.p", "Sieve P", kDouble, &fV5sieve_tr[3]},
        {"projtg.x", "Project to target plane X", kDouble, &fV5projtg_tr[0]},
        {"projtg.t", "Project to target plane T", kDouble, &fV5projtg_tr[1]},
        {"projtg.y", "Project to target plane Y", kDouble, &fV5projtg_tr[2]},
        {"projtg.p", "Project to target plane P", kDouble, &fV5projtg_tr[3]},
        {"focus.x", "Focus plane X", kDouble, &fV5fp_tr[0]},
        {"focus.t", "Focus plane T", kDouble, &fV5fp_tr[1]},
        {"focus.y", "Focus plane Y", kDouble, &fV5fp_tr[2]},
        {"focus.p", "Focus plane P", kDouble, &fV5fp_tr[3]},
        {"focus.r_x", "Focus plane X (rot)", kDouble, &fV5fp_rot[0]},
        {"focus.r_t", "Focus plane T (rot)", kDouble, &fV5fp_rot[1]},
        {"focus.r_y", "Focus plane Y (rot)", kDouble, &fV5fp_rot[2]},
        {"focus.r_p", "Focus plane P (rot)", kDouble, &fV5fp_rot[3]},
        {0}
    };

    return DefineVarsFromList(vars, mode);
}

void G2PFwdProc::MakePrefix() {
    const char* basename = "fwd";

    G2PAppBase::MakePrefix(basename);
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

ClassImp(G2PFwdProc)
