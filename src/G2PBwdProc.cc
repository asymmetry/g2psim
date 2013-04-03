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
#include "G2PDBRec.hh"
#include "G2PDrift.hh"
#include "G2PField.hh"
#include "G2PGlobals.hh"
#include "G2PHRSTrans.hh"
#include "G2PProcBase.hh"
#include "G2PRunBase.hh"
#include "G2PSieve.hh"
#include "G2PVarDef.hh"
#include "G2PVarList.hh"

#include "G2PBwdProc.hh"

G2PBwdProc::G2PBwdProc() :
    fBeamEnergy(0.0), fHRSAngle(0.0), fHRSMomentum(0.0), fFieldRatio(0.0), 
    pDrift(NULL), pHRS(NULL), pDBRec(NULL)
{
    mName["fV5bpm_bpm"] = fV5bpm_bpm; mLength["fV5bpm_bpm"] = 5;
    mName["fV5projtg_tr"] = fV5projtg_tr; mLength["fV5projtg_tr"] = 5;
    mName["fV5fp_tr"] = fV5fp_tr; mLength["fV5fp_tr"] = 5;
    mName["fV5rectg_tr"] = fV5rectg_tr; mLength["fV5rectg_tr"] = 5;
    mName["fV5recsiv_tr"] = fV5recsiv_tr; mLength["fV5recsiv_tr"] = 5;
    mName["fV5rec_tr"] = fV5rec_tr; mLength["fV5rec_tr"] = 5;
    mName["fV5rec_lab"] = fV5rec_lab; mLength["fV5rec_lab"] = 5;

    fAppsList.push_back("G2PBPM");
    fAppsList.push_back("G2PDBRec");
    fAppsList.push_back("G2PDrift");
    fAppsList.push_back("G2PHRSTrans");

    Clear();
}

G2PBwdProc::~G2PBwdProc()
{
    // Nothing to do
}

int G2PBwdProc::Init()
{
    //static const char* const here = "Init()";

    if (G2PProcBase::Init()!=0) return fStatus;

    pBPM = G2PBPM::GetInstance();
    pDrift = G2PDrift::GetInstance();
    pHRS = G2PHRSTrans::GetInstance();
    pDBRec = G2PDBRec::GetInstance();

    fApps->Add(pBPM);
    fApps->Add(pDrift);
    fApps->Add(pHRS);
    fApps->Add(pDBRec);

    return (fStatus = kOK);
}

int G2PBwdProc::Begin()
{
    //static const char* const here = "Begin()";

    if (G2PProcBase::Begin()!=0) return fStatus;

    fBeamEnergy = gG2PRun->GetBeamEnergy();
    fHRSAngle = gG2PRun->GetHRSAngle();
    fHRSMomentum = gG2PRun->GetHRSMomentum();
    G2PField* field = G2PField::GetInstance();
    if (field) fFieldRatio = field->GetRatio();
    else fFieldRatio = 0.0;

    SetSieve(fHRSAngle);

    return (fStatus = kOK);
}

int G2PBwdProc::Process()
{
    static const char* const here = "Process()";

    double V5[5];

    pBPM->TransBPM2Lab(fV5bpm_bpm, V5);
    double V5bpm_tr[5], z_tr;
    HCS2TCS(V5[0], V5[2], V5[4], fHRSAngle, V5bpm_tr[0], V5bpm_tr[2], z_tr);
    HCS2TCS(V5[1], V5[3], fHRSAngle, V5bpm_tr[1], V5bpm_tr[3]);
    V5bpm_tr[4] = 0.0;
    pDrift->Drift(V5bpm_tr, fBeamEnergy, z_tr, fHRSAngle, 0.0, 10.0, V5bpm_tr);
    double fEffXbeam = GetEffBPM(V5bpm_tr[0], fV5fp_tr);

    if (fDebug>1) {
        Info(here, "xeff_tr   : %10.3e", fEffXbeam);
    }

    fV5fp_tr[4] = fEffXbeam;
    if(!pHRS->Backward(fV5fp_tr, fV5rectg_tr)) return -1;

    if (fDebug>1) {
        Info(here, "rectg_tr  : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5rectg_tr[0], fV5rectg_tr[1], fV5rectg_tr[2], fV5rectg_tr[3], fV5rectg_tr[4]);
    }

    Project(fV5rectg_tr[0], fV5rectg_tr[2], 0.0, fSieve.fZ, fV5rectg_tr[1], fV5rectg_tr[3], fV5recsiv_tr[0], fV5recsiv_tr[2]);
    fV5recsiv_tr[1] = fV5rectg_tr[1];
    fV5recsiv_tr[3] = fV5rectg_tr[3];
    fV5recsiv_tr[4] = fV5rectg_tr[4];

    if (fDebug>1) {
        Info(here, "recsiv_tr : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5recsiv_tr[0], fV5recsiv_tr[1], fV5recsiv_tr[2], fV5recsiv_tr[3], fV5recsiv_tr[4]);
    }

    pDrift->Drift(fV5recsiv_tr, fHRSMomentum, fSieve.fZ, fHRSAngle, 0.0, 10.0, fV5rec_tr);
    TCS2HCS(fV5rec_tr[0], fV5rec_tr[2], 0.0, fHRSAngle, fV5rec_lab[0], fV5rec_lab[2], fV5rec_lab[4]);
    TCS2HCS(fV5rec_tr[1], fV5rec_tr[3], fHRSAngle, fV5rec_lab[1], fV5rec_lab[3]);

    if (fDebug>1) {
        Info(here, "rec_tr    : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5rec_tr[0], fV5rec_tr[1], fV5rec_tr[2], fV5rec_tr[3], fV5rec_tr[4]);
    }

    return 0;
}

void G2PBwdProc::Clear()
{
    memset(fV5bpm_bpm, 0, sizeof(fV5bpm_bpm));
    memset(fV5projtg_tr, 0, sizeof(fV5projtg_tr));
    memset(fV5fp_tr, 0, sizeof(fV5fp_tr));
    memset(fV5rectg_tr, 0, sizeof(fV5rectg_tr));
    memset(fV5recsiv_tr, 0, sizeof(fV5recsiv_tr));
    memset(fV5rec_tr, 0, sizeof(fV5rec_tr));
    memset(fV5rec_lab, 0, sizeof(fV5rec_lab));
}

int G2PBwdProc::DefineVariables(EMode mode)
{
    if (mode==kDefine&&bIsSetup) return 0;
    bIsSetup = (mode==kDefine);

    VarDef vars[] = {
        { "rectg.x", "SNAKE rec to target plane X", kDouble, &fV5rectg_tr[0] },
        { "rectg.t", "SNAKE rec to target plane T", kDouble, &fV5rectg_tr[1] },
        { "rectg.y", "SNAKE rec to target plane Y", kDouble, &fV5rectg_tr[2] },
        { "rectg.p", "SNAKE rec to target plane P", kDouble, &fV5rectg_tr[3] },
        { "recsieve.x", "Project to sieve X", kDouble, &fV5recsiv_tr[0] },
        { "recsieve.t", "Project to sieve T", kDouble, &fV5recsiv_tr[1] },
        { "recsieve.y", "Project to sieve Y", kDouble, &fV5recsiv_tr[2] },
        { "recsieve.p", "Project to sieve P", kDouble, &fV5recsiv_tr[3] },
        { "rec.x",   "Rec target plane X", kDouble, &fV5rec_tr[0] },
        { "rec.t",   "Rec target plane T", kDouble, &fV5rec_tr[1] },
        { "rec.y",   "Rec target plane Y", kDouble, &fV5rec_tr[2] },
        { "rec.p",   "Rec target plane P", kDouble, &fV5rec_tr[3] },
        { "rec.d",   "Delta (rec)",        kDouble, &fV5rec_tr[4] },
        { "rec.l_x", "Rec target plane X (lab)", kDouble, &fV5rec_lab[0] },
        { "rec.l_t", "Rec target plane T (lab)", kDouble, &fV5rec_lab[1] },
        { "rec.l_y", "Rec target plane Y (lab)", kDouble, &fV5rec_lab[2] },
        { "rec.l_p", "Rec target plane P (lab)", kDouble, &fV5rec_lab[3] },
        { "rec.l_z", "Rec target plane Z (lab)", kDouble, &fV5rec_lab[4] },
        { 0 }
    };

    return DefineVarsFromList(vars, mode);
}

void G2PBwdProc::MakePrefix()
{
    const char* basename = "bwd";

    G2PAppBase::MakePrefix(basename);
}

double G2PBwdProc::GetEffBPM(double xbpm_tr, const double* V5fp)
{
    // (Xbpm_tr-Xtg_tr) vs Z
    // ([0]+[1]*x)
    // Fitting result of (Xbpm_tr-Xtg_tr) vs Z @ 2.5T
    // p0                        =    -0.011838   +/-   0.00132798 
    // p1                        =      49.856    +/-   0.163115

    // (Xtg_tr-Xtgproj_tr) vs P
    // ([0]+[1]/x)
    // Fitting result of (Xtg_tr-Xtgproj_tr) vs P @ 2.5T
    // p0                        =    0.0183611   +/-   0.0105237   
    // p1                        =      3.14345   +/-   0.0105453
    // Fitting result of (Xtg_tr-Xtgproj_tr) vs P @ 5.0T
    // p0                        =      0.14139   +/-   0.018683    
    // p1                        =      6.11766   +/-   0.0187211

    double xbpm_tr_eff = xbpm_tr;

    if (fabs(fFieldRatio)<1e-8) return xbpm_tr_eff;

    double V5_fp[5] = { V5fp[0], V5fp[1], V5fp[2], V5fp[3], V5fp[4] };
    double V5_tg[5] = { 0, 0, 0, 0, 0 };

    double p = fHRSMomentum;
    double ratio = fFieldRatio;

    if (ratio<0.75) xbpm_tr_eff -= (3.14345/p + 0.0183611)/1000*ratio/0.5;
    else xbpm_tr_eff -= (6.11766/p + 0.14139)/1000*ratio/1.0;

    V5_fp[4] = xbpm_tr_eff;
    pHRS->Backward(V5_fp, V5_tg);

    p = (1+V5_tg[4])*fHRSMomentum;
    xbpm_tr_eff = xbpm_tr;
    if (ratio<0.75) xbpm_tr_eff -= (3.14345/p + 0.0183611)/1000*ratio/0.5;
    else xbpm_tr_eff -= (6.11766/p + 0.14139)/1000*ratio/1.0;

    V5_fp[4] = xbpm_tr_eff;
    pHRS->Backward(V5_fp, V5_tg);

    p = (1+V5_tg[4])*fHRSMomentum;
    xbpm_tr_eff = xbpm_tr;
    if (ratio<0.75) xbpm_tr_eff -= (3.14345/p + 0.0183611)/1000*ratio/0.5;
    else xbpm_tr_eff -= (6.11766/p + 0.14139)/1000*ratio/1.0;

    return xbpm_tr_eff;
}

ClassImp(G2PBwdProc)
