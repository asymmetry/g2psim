#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <map>

#include "TROOT.h"
#include "TObject.h"
#include "TError.h"

#include "G2PAppBase.hh"
#include "G2PDrift.hh"
#include "G2PGlobals.hh"
#include "G2PHRSTrans.hh"
#include "G2PProcBase.hh"
#include "G2PRecUseDB.hh"
#include "G2PRunBase.hh"
#include "G2PSieve.hh"
#include "G2PVarDef.hh"
#include "G2PVarList.hh"

#include "G2PFwdProc.hh"

G2PFwdProc::G2PFwdProc() :
    fHRSAngle(0.0), fHRSMomentum(0.0),
    pDrift(NULL), pHRS(NULL), pRecDB(NULL)
{
    Clear();
}

G2PFwdProc::~G2PFwdProc()
{
    // Nothing to do
}

int G2PFwdProc::Init()
{
    static const char* const here = "Init()";

    if (G2PProcBase::Init()!=0) return fStatus;

    pDrift = G2PDrift::GetInstance();
    if (!pDrift) {
        Error(here, "Cannot initialize, no G2PDrift found.");
        return (fStatus = kINITERROR);
    }

    pHRS = G2PHRSTrans::GetInstance();
    if (!pHRS) {
        Error(here, "Cannot initialize, no G2PHRSTrans found.");
        return (fStatus = kINITERROR);
    }

    pRecDB = G2PRecUseDB::GetInstance();
    if (!pRecDB) {
        Error(here, "Cannot initialize, no G2PRecUseDB found.");
        return (fStatus = kINITERROR);
    }

    fApps->Add(pDrift);
    fApps->Add(pHRS);
    fApps->Add(pRecDB);

    return (fStatus = kOK);
}

int G2PFwdProc::Begin()
{
    //static const char* const here = "Begin()";

    if (G2PProcBase::Begin()!=0) return fStatus;

    fHRSAngle = gG2PRun->GetHRSAngle();
    fHRSMomentum = gG2PRun->GetHRSMomentum();

    mName["fV5tg_tr"] = fV5tg_tr; mLength["fV5tg_tr"] = 5;
    mName["fV5sieve_tr"] = fV5sieve_tr; mLength["fV5sieve_tr"] = 5;
    mName["fV5projtg_tr"] = fV5projtg_tr; mLength["fV5projtg_tr"] = 5;
    mName["fV5fp_tr"] = fV5fp_tr; mLength["fV5fp_tr"] = 5;
    mName["fV5fp_rot"] = fV5fp_rot; mLength["fV5fp_rot"] = 5;

    SetSieve(fHRSAngle);

    return (fStatus = kOK);
}

int G2PFwdProc::Process()
{
    static const char* const here = "Process()";

    //double V5[5];

    pDrift->Drift(fV5tg_tr, fHRSMomentum, 0.0, fHRSAngle, fSieve.fZ, 10.0, fV5sieve_tr);

    if (fDebug>1) {
        Info(here, "sieve_tr  : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5sieve_tr[0], fV5sieve_tr[1], fV5sieve_tr[2], fV5sieve_tr[3], fV5sieve_tr[4]);
    }

    Project(fV5sieve_tr[0], fV5sieve_tr[2], fSieve.fZ, 0.0, fV5sieve_tr[1], fV5sieve_tr[3], fV5projtg_tr[0], fV5projtg_tr[2]);
    fV5projtg_tr[1] = fV5sieve_tr[1];
    fV5projtg_tr[3] = fV5sieve_tr[3];
    fV5projtg_tr[4] = fV5sieve_tr[4];

    if (fDebug>1) {
        Info(here, "projtg_tr : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5projtg_tr[0], fV5projtg_tr[1], fV5projtg_tr[2], fV5projtg_tr[3], fV5projtg_tr[4]);
    }

    bIsGood = pHRS->Forward(fV5projtg_tr, fV5fp_tr);
    ApplyVDCRes(fV5fp_tr);
    pRecDB->TransTr2Rot(fV5fp_tr, fV5fp_rot);

    if (fDebug>1) {
        Info(here, "fp_tr     : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5fp_tr[0], fV5fp_tr[1], fV5fp_tr[2], fV5fp_tr[3], fV5fp_tr[4]);
    }

    return 0;
}

void G2PFwdProc::Clear()
{
    G2PProcBase::Clear();

    memset(fV5tg_tr, 0, sizeof(fV5tg_tr));
    memset(fV5sieve_tr, 0, sizeof(fV5sieve_tr));
    memset(fV5projtg_tr, 0, sizeof(fV5projtg_tr));
    memset(fV5fp_tr, 0, sizeof(fV5fp_tr));
    memset(fV5fp_rot, 0, sizeof(fV5fp_rot));
}

int G2PFwdProc::DefineVariables(EMode mode)
{
    if (mode==kDefine&&bIsSetup) return 0;
    bIsSetup = (mode==kDefine);

    VarDef vars[] = {
        { "isgood",    "Reach focus plane", kBool, &bIsGood },
        { "sieve.x",   "Sieve X", kDouble, &fV5sieve_tr[0] },
        { "sieve.t",   "Sieve T", kDouble, &fV5sieve_tr[1] },
        { "sieve.y",   "Sieve Y", kDouble, &fV5sieve_tr[2] },
        { "sieve.p",   "Sieve P", kDouble, &fV5sieve_tr[3] },
        { "projtg.x",  "Project to target plane X", kDouble, &fV5projtg_tr[0] },
        { "projtg.t",  "Project to target plane T", kDouble, &fV5projtg_tr[1] },
        { "projtg.y",  "Project to target plane Y", kDouble, &fV5projtg_tr[2] },
        { "projtg.p",  "Project to target plane P", kDouble, &fV5projtg_tr[3] },
        { "focus.x",   "Focus plane X", kDouble, &fV5fp_tr[0] },
        { "focus.t",   "Focus plane T", kDouble, &fV5fp_tr[1] },
        { "focus.y",   "Focus plane Y", kDouble, &fV5fp_tr[2] },
        { "focus.p",   "Focus plane P", kDouble, &fV5fp_tr[3] },
        { "focus.r_x", "Focus plane X (rot)", kDouble, &fV5fp_rot[0] },
        { "focus.r_t", "Focus plane T (rot)", kDouble, &fV5fp_rot[1] },
        { "focus.r_y", "Focus plane Y (rot)", kDouble, &fV5fp_rot[2] },
        { "focus.r_p", "Focus plane P (rot)", kDouble, &fV5fp_rot[3] },
        { 0 }
    };

    return DefineVarsFromList(vars, mode);
}

void G2PFwdProc::MakePrefix()
{
    const char* basename = "fwd";

    G2PAppBase::MakePrefix(basename);
}

void G2PFwdProc::ApplyVDCRes(double* V5fp)
{
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
