#include <cstdlib>
#include <cmath>

#include "TROOT.h"
#include "TObject.h"
#include "TError.h"
#include "TTree.h"

#include "G2PBPM.hh"
#include "G2PDrift.hh"
#include "G2PFieldBase.hh"
#include "G2PGlobals.hh"
#include "G2PGunBase.hh"
#include "G2PRand.hh"
#include "G2PRunBase.hh"
#include "G2PHRSTrans.hh"
#include "G2PPhys.hh"
#include "G2PRecUseDB.hh"

#include "G2PRun.hh"

static const double kDEG = 3.14159265358979323846/180.0;

G2PRun::G2PRun() :
    bUseEffBPM(true), pfRun(NULL)
{   
    Clear();
}

G2PRun::~G2PRun()
{
    // Nothing to do
}

G2PAppsBase::EStatus G2PRun::Init()
{
    // static const char* const here = "Init()";

    if (G2PRunBase::Init()) return fStatus;

    if (pGun->UseData()) pfRun = &G2PRun::RunData;
    else pfRun = &G2PRun::RunSim;

    SetSieve(fHRSAngle);

    return (fStatus = kOK);
}

void G2PRun::Clear()
{
    bIsGood = false;

    memset(fV5beam_lab, 0, sizeof(fV5beam_lab));
    memset(fV5react_tr, 0, sizeof(fV5react_tr));
    memset(fV5tg_tr, 0, sizeof(fV5tg_tr));
    memset(fV5tg_lab, 0, sizeof(fV5tg_lab));

    memset(fV5fpdata_tr, 0, sizeof(fV5fpdata_tr));
    memset(fV5fpdata_rot, 0, sizeof(fV5fpdata_rot));

    memset(fV5bpm_lab, 0, sizeof(fV5bpm_lab));
    memset(fV5bpm_tr, 0, sizeof(fV5bpm_tr));

    memset(fV5fp_tr, 0, sizeof(fV5fp_tr));
    memset(fV5fp_rot, 0, sizeof(fV5fp_rot));

    memset(fV5rec_tr, 0, sizeof(fV5rec_tr));
    memset(fV5rec_lab, 0, sizeof(fV5rec_lab));
    memset(fV5recdb_tr, 0, sizeof(fV5recdb_tr));
    memset(fV5recdb_lab, 0, sizeof(fV5recdb_lab));

    memset(fV5sieve_tr, 0, sizeof(fV5sieve_tr));
    memset(fV5tgproj_tr, 0, sizeof(fV5tgproj_tr));
    memset(fV5rectg_tr, 0, sizeof(fV5rectg_tr));
    memset(fV5recsieve_tr, 0, sizeof(fV5recsieve_tr));

    fRealXbeam = 0.0; fEffXbeam = 0.0;

    fThetainit = 0.0; fThetarec = 0.0;
    fXSinit = 1.0; fXSrec = 1.0;
}

int G2PRun::DefineVariables(TTree *t)
{
    t->Branch("IsGood", &bIsGood, "IsGood/O");

    t->Branch("Xbeam_lab", &fV5beam_lab[0], "Xbeam_lab/D");
    t->Branch("Tbeam_lab", &fV5beam_lab[1], "Tbeam_lab/D");
    t->Branch("Ybeam_lab", &fV5beam_lab[2], "Ybeam_lab/D");
    t->Branch("Pbeam_lab", &fV5beam_lab[3], "Pbeam_lab/D");
    t->Branch("Zbeam_lab", &fV5beam_lab[4], "Zbeam_lab/D");

    t->Branch("Xbpm_tr", &fV5bpm_tr[0], "Xbpm_tr/D");
    t->Branch("Tbpm_tr", &fV5bpm_tr[1], "Tbpm_tr/D");
    t->Branch("Ybpm_tr", &fV5bpm_tr[2], "Ybpm_tr/D");
    t->Branch("Pbpm_tr", &fV5bpm_tr[3], "Pbpm_tr/D");

    t->Branch("Xbpm_lab", &fV5bpm_lab[0], "Xbpm_lab/D");
    t->Branch("Tbpm_lab", &fV5bpm_lab[1], "Tbpm_lab/D");
    t->Branch("Ybpm_lab", &fV5bpm_lab[2], "Ybpm_lab/D");
    t->Branch("Pbpm_lab", &fV5bpm_lab[3], "Pbpm_lab/D");
    t->Branch("Zbpm_lab", &fV5bpm_lab[4], "Zbpm_lab/D");

    t->Branch("Xtg_tr", &fV5tg_tr[0], "Xtg_tr/D");
    t->Branch("Ttg_tr", &fV5tg_tr[1], "Ttg_tr/D");
    t->Branch("Ytg_tr", &fV5tg_tr[2], "Ytg_tr/D");
    t->Branch("Ptg_tr", &fV5tg_tr[3], "Ptg_tr/D");
    t->Branch("Delta", &fV5tg_tr[4], "Delta/D");

    t->Branch("Xreact_tr", &fV5react_tr[0], "Xreact_tr/D");
    t->Branch("Treact_tr", &fV5react_tr[1], "Treact_tr/D");
    t->Branch("Yreact_tr", &fV5react_tr[2], "Yreact_tr/D");
    t->Branch("Preact_tr", &fV5react_tr[3], "Preact_tr/D");

    t->Branch("Xtg_lab", &fV5tg_lab[0], "Xtg_lab/D");
    t->Branch("Ttg_lab", &fV5tg_lab[1], "Ttg_lab/D");
    t->Branch("Ytg_lab", &fV5tg_lab[2], "Ytg_lab/D");
    t->Branch("Ptg_lab", &fV5tg_lab[3], "Ptg_lab/D");

    t->Branch("Xfpdata_tr", &fV5fpdata_tr[0], "Xfpdata_tr/D");
    t->Branch("Tfpdata_tr", &fV5fpdata_tr[1], "Tfpdata_tr/D");
    t->Branch("Yfpdata_tr", &fV5fpdata_tr[2], "Yfpdata_tr/D");
    t->Branch("Pfpdata_tr", &fV5fpdata_tr[3], "Pfpdata_tr/D");

    t->Branch("Xfpdata_rot", &fV5fpdata_rot[0], "Xfpdata_rot/D");
    t->Branch("Tfpdata_rot", &fV5fpdata_rot[1], "Tfpdata_rot/D");
    t->Branch("Yfpdata_rot", &fV5fpdata_rot[2], "Yfpdata_rot/D");
    t->Branch("Pfpdata_rot", &fV5fpdata_rot[3], "Pfpdata_rot/D");
    
    t->Branch("Xfp_tr", &fV5fp_tr[0], "Xfp_tr/D");
    t->Branch("Tfp_tr", &fV5fp_tr[1], "Tfp_tr/D");
    t->Branch("Yfp_tr", &fV5fp_tr[2], "Yfp_tr/D");
    t->Branch("Pfp_tr", &fV5fp_tr[3], "Pfp_tr/D");

    t->Branch("Xfp_rot", &fV5fp_rot[0], "Xfp_rot/D");
    t->Branch("Tfp_rot", &fV5fp_rot[1], "Tfp_rot/D");
    t->Branch("Yfp_rot", &fV5fp_rot[2], "Yfp_rot/D");
    t->Branch("Pfp_rot", &fV5fp_rot[3], "Pfp_rot/D");

    t->Branch("Xrec_tr", &fV5rec_tr[0], "Xrec_tr/D");
    t->Branch("Trec_tr", &fV5rec_tr[1], "Trec_tr/D");
    t->Branch("Yrec_tr", &fV5rec_tr[2], "Yrec_tr/D");
    t->Branch("Prec_tr", &fV5rec_tr[3], "Prec_tr/D");
    t->Branch("Deltarec", &fV5rec_tr[4], "Deltarec/D");

    t->Branch("Xrec_lab", &fV5rec_lab[0], "Xrec_lab/D");
    t->Branch("Trec_lab", &fV5rec_lab[1], "Trec_lab/D");
    t->Branch("Yrec_lab", &fV5rec_lab[2], "Yrec_lab/D");
    t->Branch("Prec_lab", &fV5rec_lab[3], "Prec_lab/D");

    t->Branch("Xrecdb_tr", &fV5recdb_tr[0], "Xrecdb_tr/D");
    t->Branch("Trecdb_tr", &fV5recdb_tr[1], "Trecdb_tr/D");
    t->Branch("Yrecdb_tr", &fV5recdb_tr[2], "Yrecdb_tr/D");
    t->Branch("Precdb_tr", &fV5recdb_tr[3], "Precdb_tr/D");
    t->Branch("Deltarecdb", &fV5recdb_tr[4], "Deltarecdb/D");

    t->Branch("Xrecdb_lab", &fV5recdb_lab[0], "Xrecdb_lab/D");
    t->Branch("Trecdb_lab", &fV5recdb_lab[1], "Trecdb_lab/D");
    t->Branch("Yrecdb_lab", &fV5recdb_lab[2], "Yrecdb_lab/D");
    t->Branch("Precdb_lab", &fV5recdb_lab[3], "Precdb_lab/D");

    t->Branch("Xsieve_tr", &fV5sieve_tr[0], "Xsieve_tr/D");
    t->Branch("Tsieve_tr", &fV5sieve_tr[1], "Tsieve_tr/D");
    t->Branch("Ysieve_tr", &fV5sieve_tr[2], "Ysieve_tr/D");
    t->Branch("Psieve_tr", &fV5sieve_tr[3], "Psieve_tr/D");

    t->Branch("Xtgproj_tr", &fV5tgproj_tr[0], "Xtgproj_tr/D");
    t->Branch("Ttgproj_tr", &fV5tgproj_tr[1], "Ttgproj_tr/D");
    t->Branch("Ytgproj_tr", &fV5tgproj_tr[2], "Ytgproj_tr/D");
    t->Branch("Ptgproj_tr", &fV5tgproj_tr[3], "Ptgproj_tr/D");

    t->Branch("Xrectg_tr", &fV5rectg_tr[0], "Xrectg_tr/D");
    t->Branch("Trectg_tr", &fV5rectg_tr[1], "Trectg_tr/D");
    t->Branch("Yrectg_tr", &fV5rectg_tr[2], "Yrectg_tr/D");
    t->Branch("Prectg_tr", &fV5rectg_tr[3], "Prectg_tr/D");

    t->Branch("Xrecsieve_tr", &fV5recsieve_tr[0], "Xrecsieve_tr/D");
    t->Branch("Trecsieve_tr", &fV5recsieve_tr[1], "Trecsieve_tr/D");
    t->Branch("Yrecsieve_tr", &fV5recsieve_tr[2], "Yrecsieve_tr/D");
    t->Branch("Precsieve_tr", &fV5recsieve_tr[3], "Precsieve_tr/D");

    t->Branch("Xeff_tr", &fEffXbeam, "Xeff_tr/D");
    t->Branch("Xreal_tr", &fRealXbeam, "Xreal_tr/D");

    t->Branch("Thetainit", &fThetainit, "Thetainit/D");
    t->Branch("Thetarec", &fThetarec, "Thetarec/D");
    t->Branch("XSinit", &fXSinit, "XSinit/D");
    t->Branch("XSrec", &fXSrec, "XSrec/D");

    return 0;
}

int G2PRun::RunSim()
{
    static const char* const here = "Run()";

    double x_tr, y_tr, z_tr;

    if (!pGun->Shoot(fV5beam_lab, fV5react_tr)) return -1;

    fXSinit = CalXS(fV5beam_lab, fV5react_tr, fThetainit);

    pBPM->GetBPMValue(fV5beam_lab, fV5bpm_lab);
    double V5[5];
    pBPM->TransBPM2Lab(fV5bpm_lab, V5);
    HCS2TCS(V5[0], V5[2], V5[4], fHRSAngle, fV5bpm_tr[0], fV5bpm_tr[2], z_tr);
    HCS2TCS(V5[1], V5[3], fHRSAngle, fV5bpm_tr[1], fV5bpm_tr[3]);
    fV5bpm_tr[4] = 0.0;
    pDrift->Drift(fV5bpm_tr, fBeamEnergy, z_tr, fHRSAngle, 0.0, 10.0, fV5bpm_tr);
    if (fDebug>1) {
        Info(here, "bpm_tr  : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5bpm_tr[0], fV5bpm_tr[1], fV5bpm_tr[2], fV5bpm_tr[3], fV5bpm_tr[4]);
    }

    HCS2TCS(fV5beam_lab[0], fV5beam_lab[2], fV5beam_lab[4], fHRSAngle, x_tr, y_tr, z_tr);
    pDrift->Drift(fV5react_tr, fHRSMomentum, z_tr, fHRSAngle, 0.0, 10.0, fV5tg_tr);
    if (fDebug>1) {
        Info(here, "tg_tr   : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5tg_tr[0], fV5tg_tr[1], fV5tg_tr[2], fV5tg_tr[3], fV5tg_tr[4]);
    }
    
    pDrift->Drift(fV5tg_tr, fHRSMomentum, 0.0, fHRSAngle, fSieve.fZ, 10.0, fV5sieve_tr);
    if (fDebug>1) {
        Info(here, "sieve_tr: %10.3e %10.3e %10.3e %10.3e %10.3e", fV5sieve_tr[0], fV5sieve_tr[1], fV5sieve_tr[2], fV5sieve_tr[3], fV5sieve_tr[4]);
    }

    Project(fV5sieve_tr[0], fV5sieve_tr[2], fSieve.fZ, 0.0, fV5sieve_tr[1], fV5sieve_tr[3], fV5tgproj_tr[0], fV5tgproj_tr[2]);
    fV5tgproj_tr[1] = fV5sieve_tr[1];
    fV5tgproj_tr[3] = fV5sieve_tr[3];
    fV5tgproj_tr[4] = fV5sieve_tr[4];
    fRealXbeam = fV5tgproj_tr[0];
    if (fDebug>1) {
        Info(here, "proj_tr : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5tgproj_tr[0], fV5tgproj_tr[1], fV5tgproj_tr[2], fV5tgproj_tr[3], fV5tgproj_tr[4]);
    }

    bIsGood = pHRS->Forward(fV5tgproj_tr, fV5fp_tr);
    pRecUseDB->TransTr2Rot(fV5fp_tr, fV5fp_rot);

    if (!bIsGood) return 0;

    SmearVDC(fV5fp_tr);

    if (bUseEffBPM) {
        fEffXbeam = GetEffBPM(fV5bpm_tr[0], fV5fp_tr);
        if (fDebug>1) {
            Info(here, "bpm_eff : %10.3e %10.3e", fRealXbeam, fEffXbeam);
        }
        fV5fp_tr[4] = fEffXbeam;
        fV5fp_rot[4] = fEffXbeam;
    }
    else {
        fV5fp_tr[4] = fRealXbeam;
        fV5fp_rot[4] = fRealXbeam;
    }

    pRecUseDB->CalcTargetCoords(fV5fp_rot, fV5rectg_tr);
    Project(fV5rectg_tr[0], fV5rectg_tr[2], 0.0, fSieve.fZ, fV5rectg_tr[1], fV5rectg_tr[3], fV5recsieve_tr[0], fV5recsieve_tr[2]);
    fV5recsieve_tr[1] = fV5rectg_tr[1];
    fV5recsieve_tr[3] = fV5rectg_tr[3];
    fV5recsieve_tr[4] = fV5rectg_tr[4];
    pDrift->Drift(fV5recsieve_tr, fHRSMomentum, fSieve.fZ, fHRSAngle, 0.0, 10.0, fV5recdb_tr);

    bIsGood &= pHRS->Backward(fV5fp_tr, fV5rectg_tr);
    if (fDebug>1) {
        Info(here, "recj_tr : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5rectg_tr[0], fV5rectg_tr[1], fV5rectg_tr[2], fV5rectg_tr[3], fV5rectg_tr[4]);
    }

    Project(fV5rectg_tr[0], fV5rectg_tr[2], 0.0, fSieve.fZ, fV5rectg_tr[1], fV5rectg_tr[3], fV5recsieve_tr[0], fV5recsieve_tr[2]);
    fV5recsieve_tr[1] = fV5rectg_tr[1];
    fV5recsieve_tr[3] = fV5rectg_tr[3];
    fV5recsieve_tr[4] = fV5rectg_tr[4];
    if (fDebug>1) {
        Info(here, "rcsiv_tr: %10.3e %10.3e %10.3e %10.3e %10.3e", fV5recsieve_tr[0], fV5recsieve_tr[1], fV5recsieve_tr[2], fV5recsieve_tr[3], fV5recsieve_tr[4]);
    }

    pDrift->Drift(fV5recsieve_tr, fHRSMomentum, fSieve.fZ, fHRSAngle, 0.0, 10.0, fV5rec_tr);
    if (fDebug>1) {
        Info(here, "rec_tr  : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5rec_tr[0], fV5rec_tr[1], fV5rec_tr[2], fV5rec_tr[3], fV5rec_tr[4]);
    }

    fXSrec = CalXS(V5, fV5rec_tr, fThetarec);

    return 0;
}

int G2PRun::RunData()
{
    // static const char* const here = "Run()";

    if (!pGun->Shoot(fV5bpm_lab, fV5tg_tr, fV5fpdata_tr)) return -1;

    fV5fpdata_tr[4] = -fV5bpm_lab[2];

    pRecUseDB->TransTr2Rot(fV5fpdata_tr, fV5fpdata_rot);
    pRecUseDB->CalcTargetCoords(fV5fpdata_rot, fV5recdb_tr);

    fV5recdb_tr[0] = -fV5bpm_lab[2];
    bIsGood = pHRS->Forward(fV5recdb_tr, fV5fp_tr);
    pRecUseDB->TransTr2Rot(fV5fp_tr, fV5fp_rot);

    bIsGood &= pHRS->Backward(fV5fpdata_tr, fV5rec_tr);

    return 0;
}

void G2PRun::SmearVDC(double* V5_fp)
{
    double WireChamberResX = 0.0013; //m;
    double WireChamberResY = 0.0013; //m;
    double WireChamberResT = 0.0003; //rad;
    double WireChamberResP = 0.0003; //rad;

    V5_fp[0] = pRand->Gaus(V5_fp[0], WireChamberResX);
    V5_fp[1] = pRand->Gaus(V5_fp[1], WireChamberResT);
    V5_fp[2] = pRand->Gaus(V5_fp[2], WireChamberResY);
    V5_fp[3] = pRand->Gaus(V5_fp[3], WireChamberResP);
}

double G2PRun::GetEffBPM(double xbpm_tr, const double* V5fp)
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

    if (!pDrift->GetField()) return xbpm_tr_eff;

    double V5_fp[5] = { V5fp[0], V5fp[1], V5fp[2], V5fp[3], V5fp[4] };
    double V5_tg[5] = { 0, 0, 0, 0, 0 };

    double ratio = pDrift->GetField()->GetRatio();

    double p = fHRSMomentum;
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

double G2PRun::CalXS(const double* V5lab, const double* V5tr, double& scatangle)
{
    double Eb[3] = { sin(V5lab[1])*cos(V5lab[3]),
                     sin(V5lab[1])*sin(V5lab[3]),
                     cos(V5lab[1]) };

    double theta, phi;
    TCS2HCS(V5tr[1], V5tr[3], fHRSAngle, theta, phi);

    double Ef[3] = { sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta) };

    scatangle = acos(Eb[0]*Ef[0]+Eb[1]*Ef[1]+Eb[2]*Ef[2]);

    double Ebval = fBeamEnergy;
    double Efval = (1+V5tr[4])*fHRSMomentum;

    return pPhys->GetXS(Ebval, Efval, scatangle);
}

// double G2PRun::DriftPath()
// {   
//     if (bUseField) {
//         double x[3] = { fV5bpm_lab[0], fV5bpm_lab[2], fV5bpm_lab[4] };
//         double p[3] = { fBeamEnergy*sin(fV5bpm_lab[1])*cos(fV5bpm_lab[3]),
//                         fBeamEnergy*sin(fV5bpm_lab[1])*sin(fV5bpm_lab[3]),
//                         fBeamEnergy*cos(fV5bpm_lab[1]) };

//         double xrec[3];
//         HRSTransTCSNHCS::X_TCS2HCS(fV5rec_tr[0], fV5rec_tr[2], 0.0, fHRSAngle, xrec[0], xrec[1], xrec[2]);
//         double theta,phi;
//         HRSTransTCSNHCS::P_TCS2HCS(fV5rec_tr[1], fV5rec_tr[3], fHRSAngle, theta, phi);
//         double prec[3] = { fHRSMomentum*(1+fV5rec_tr[4])*sin(theta)*cos(phi),
//                            fHRSMomentum*(1+fV5rec_tr[4])*sin(theta)*sin(phi),
//                            fHRSMomentum*(1+fV5rec_tr[4])*cos(theta) };
//         pDrift->Drift(xrec, prec, 0, 10.0, xrec, prec);

//         double z = 0.0;
//         double zmin = 0.0;
//         double distmin = sqrt((x[0]-xrec[0])*(x[0]-xrec[0])+(x[1]-xrec[1])*(x[1]-xrec[1]));
//         for (int i = 1; i<150; i++) {
//             z = 0.1e-3*i;
//             pDrift->Drift(x, p, z, 10.0, x, p);
//             pDrift->Drift(xrec, prec, z, 10.0, xrec, prec);
//             //printf("%e\t%e\t%e\t%e\t%e\n", z, x[0], x[1], xrec[0], xrec[1]);
//             double distance = sqrt((x[0]-xrec[0])*(x[0]-xrec[0])+(x[1]-xrec[1])*(x[1]-xrec[1]));
//             if (distance<distmin) {
//                 zmin = z;
//                 distmin = distance;
//             }
//         }

//         // for (int i = 20; i<200; i++) {
//         //     z = 1e-3*i;
//         //     pDrift->Drift(xrec, prec, z, 10.0, xrec, prec);
//         //     printf("%e\t%e\t%e\t%e\t%e\n", z, 1000.0, 1000.0, xrec[0], xrec[1]);
//         // }

//         for (int i = 1; i<150; i++) {
//             z = -0.1e-3*i;
//             pDrift->Drift(x, p, z, 10.0, x, p);
//             pDrift->Drift(xrec, prec, z, 10.0, xrec, prec);
//             //printf("%e\t%e\t%e\t%e\t%e\n", z, x[0], x[1], xrec[0], xrec[1]);
//             double distance = sqrt((x[0]-xrec[0])*(x[0]-xrec[0])+(x[1]-xrec[1])*(x[1]-xrec[1]));
//             if (distance<distmin) {
//                 zmin = z;
//                 distmin = distance;
//             }
//         }

//         // for (int i = -20; i>-200; i--) {
//         //     z = 1e-3*i;
//         //     pDrift->Drift(x, p, z, 10.0, x, p);
//         //     printf("%e\t%e\t%e\t%e\t%e\n", z, x[0], x[1], 1000.0, 1000.0);
//         // }

//         pDrift->Drift(x, p, zmin, 10.0, x, p);
//         pDrift->Drift(xrec, prec, zmin, 10.0, xrec, prec);

//         double x_tr,y_tr,z_tr;
//         HRSTransTCSNHCS::X_HCS2TCS(xrec[0], xrec[1], xrec[2], fHRSAngle, x_tr, y_tr, z_tr);
//         return (z_tr);
//     }
//     else return 0;
// }

ClassImp(G2PRun)
