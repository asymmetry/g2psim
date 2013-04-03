#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <map>

#include "TROOT.h"
#include "TObject.h"
#include "TError.h"

#include "G2PAppBase.hh"
#include "G2PBPM.hh"
#include "G2PDrift.hh"
#include "G2PGlobals.hh"
#include "G2PGun.hh"
#include "G2PProcBase.hh"
#include "G2PRunBase.hh"
#include "G2PVarDef.hh"
#include "G2PVarList.hh"

#include "G2PGunProc.hh"

G2PGunProc::G2PGunProc() :
    fHRSAngle(0.0), fHRSMomentum(0.0),
    pBPM(NULL), pDrift(NULL), pGun(NULL)
{
    mName["fV5beam_lab"] = fV5beam_lab; mLength["fV5beam_lab"] = 5;
    mName["fV5react_tr"] = fV5react_tr; mLength["fV5react_tr"] = 5;
    mName["fV5react_lab"] = fV5react_lab; mLength["fV5react_lab"] = 5;
    mName["fV5bpm_bpm"] = fV5bpm_bpm; mLength["fV5bpm_bpm"] = 5;
    mName["fV5bpm_lab"] = fV5bpm_lab; mLength["fV5bpm_lab"] = 5;
    mName["fV5tg_tr"] = fV5tg_tr; mLength["fV5tg_tr"] = 5;

    fAppsList.push_back("G2PBPM");
    fAppsList.push_back("G2PDrift");
    fAppsList.push_back("G2PGun");

    Clear();
}

G2PGunProc::~G2PGunProc()
{
    // Nothing to do
}

int G2PGunProc::Init()
{
    //static const char* const here = "Init()";

    if (G2PProcBase::Init()!=0) return fStatus;

    pBPM = G2PBPM::GetInstance();
    pDrift = G2PDrift::GetInstance();
    pGun = G2PGun::GetInstance();

    fApps->Add(pBPM);
    fApps->Add(pDrift);
    fApps->Add(pGun);
    
    return (fStatus = kOK);
}

int G2PGunProc::Begin()
{
    //static const char* const here = "Begin()";

    if (G2PProcBase::Begin()!=0) return fStatus;

    fHRSAngle = gG2PRun->GetHRSAngle();
    fHRSMomentum = gG2PRun->GetHRSMomentum();

    return (fStatus = kOK);
}

int G2PGunProc::Process()
{
    static const char* const here = "Process()";

    double V5[5];

    pGun->Shoot(fV5beam_lab, fV5react_tr);
    ArrayCopy(fV5react_lab, fV5beam_lab, 5);
    TCS2HCS(fV5react_tr[1], fV5react_tr[3], fHRSAngle, fV5react_lab[1], fV5react_lab[3]);

    if (fDebug>1) {
        Info(here, "beam_lab  : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5beam_lab[0], fV5beam_lab[1], fV5beam_lab[2], fV5beam_lab[3], fV5beam_lab[4]);
    }

    pBPM->GetBPMValue(fV5beam_lab, fV5bpm_bpm);
    pBPM->TransBPM2Lab(fV5bpm_bpm, fV5bpm_lab);

    if (fDebug>1) {
        Info(here, "bpm_bpm   : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5bpm_bpm[0], fV5bpm_bpm[1], fV5bpm_bpm[2], fV5bpm_bpm[3], fV5bpm_bpm[4]);
    }

    HCS2TCS(fV5beam_lab[0], fV5beam_lab[2], fV5beam_lab[4], fHRSAngle, V5[0], V5[2], V5[4]);
    pDrift->Drift(fV5react_tr, fHRSMomentum, V5[4], fHRSAngle, 0.0, 10.0, fV5tg_tr);

    if (fDebug>1) {
        Info(here, "tg_tr     : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5tg_tr[0], fV5tg_tr[1], fV5tg_tr[2], fV5tg_tr[3], fV5tg_tr[4]);
    }

    return 0;
}

void G2PGunProc::Clear()
{
    memset(fV5beam_lab, 0, sizeof(fV5beam_lab));
    memset(fV5react_tr, 0, sizeof(fV5react_tr));
    memset(fV5react_lab, 0, sizeof(fV5react_lab));
    memset(fV5bpm_bpm, 0, sizeof(fV5bpm_bpm));
    memset(fV5bpm_lab, 0, sizeof(fV5bpm_lab));
    memset(fV5tg_tr, 0, sizeof(fV5tg_tr));
}

int G2PGunProc::DefineVariables(EMode mode)
{
    if (mode==kDefine&&bIsSetup) return 0;
    bIsSetup = (mode==kDefine);

    VarDef vars[] = {
        { "react.l_x", "React point X (lab)", kDouble, &fV5beam_lab[0] },
        { "beam.l_t",  "Beam T (lab)",        kDouble, &fV5beam_lab[1] },
        { "react.l_y", "React point Y (lab)", kDouble, &fV5beam_lab[2] },
        { "beam.l_p",  "Beam P (lab)",        kDouble, &fV5beam_lab[3] },
        { "react.l_z", "React point Z (lab)", kDouble, &fV5beam_lab[4] },
        { "react.x",   "React point X",       kDouble, &fV5react_tr[0] },
        { "react.t",   "Scatted particle T",  kDouble, &fV5react_tr[1] },
        { "react.y",   "React point Y",       kDouble, &fV5react_tr[2] },
        { "react.p",   "Scatted particle P",  kDouble, &fV5react_tr[3] },
        { "react.d",   "Delta (init)",        kDouble, &fV5react_tr[4] },
        { "react.l_t", "Scatted particle T (lab)", kDouble, &fV5react_lab[1] },
        { "react.l_p", "Scatted particle P (lab)", kDouble, &fV5react_lab[3] },
        { "bpm.t",     "BPM T",           kDouble, &fV5bpm_bpm[1] },
        { "bpm.p",     "BPM P",           kDouble, &fV5bpm_bpm[3] },
        { "bpm.l_x",   "BPM X (lab)",     kDouble, &fV5bpm_lab[0] },
        { "bpm.l_t",   "BPM T (lab)",     kDouble, &fV5bpm_lab[1] },
        { "bpm.l_y",   "BPM Y (lab)",     kDouble, &fV5bpm_lab[2] },
        { "bpm.l_p",   "BPM P (lab)",     kDouble, &fV5bpm_lab[3] },
        { "bpm.l_z",   "BPM Z (lab)",     kDouble, &fV5bpm_lab[4] },
        { "target.x",  "Target plane X",  kDouble, &fV5tg_tr[0] },
        { "target.t",  "Target plane T",  kDouble, &fV5tg_tr[1] },
        { "target.y",  "Target plane Y",  kDouble, &fV5tg_tr[2] },
        { "target.p",  "Target plane P",  kDouble, &fV5tg_tr[3] },
        { 0 }
    };

    return DefineVarsFromList(vars, mode);
}

void G2PGunProc::MakePrefix()
{
    const char* basename = "gun";

    G2PAppBase::MakePrefix(basename);
}

ClassImp(G2PGunProc)
