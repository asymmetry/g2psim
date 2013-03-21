#include <cstdlib>
#include <cmath>
#include <vector>

#include "TROOT.h"
#include "TObject.h"
#include "TError.h"

#include "G2PBwdProc.hh"
#include "G2PFwdProc.hh"
#include "G2PGlobals.hh"
#include "G2PGunProc.hh"
#include "G2PRunBase.hh"
#include "G2PPhyProc.hh"
#include "G2PProcBase.hh"

#include "G2PRun.hh"

using namespace std;

static const double kDEG = 3.14159265358979323846/180.0;

G2PRun::G2PRun()
{
    // Nothing to do
}

G2PRun::~G2PRun()
{
    // Nothing to do
}

int G2PRun::Init()
{
    static const char* const here = "Init()";

    // Notice: the order is important
    fProcs->Add(new G2PGunProc());
    fProcs->Add(new G2PFwdProc());
    fProcs->Add(new G2PBwdProc());
    fProcs->Add(new G2PPhyProc());

    vector<const char*> gunreqs; gunreqs.clear();
    fProcReqs.push_back(gunreqs);

    vector<const char*> fwdreqs; fwdreqs.clear();
    fwdreqs.push_back("fV5tg_tr");
    fProcReqs.push_back(fwdreqs);

    vector<const char*> bwdreqs; bwdreqs.clear();
    bwdreqs.push_back("fV5bpm_bpm");
    bwdreqs.push_back("fV5projtg_tr");
    bwdreqs.push_back("fV5fp_tr");
    fProcReqs.push_back(bwdreqs);

    vector<const char*> phyreqs; phyreqs.clear();
    phyreqs.push_back("fV5beam_lab");
    phyreqs.push_back("fV5bpm_lab");
    phyreqs.push_back("fV5react_tr");
    phyreqs.push_back("fV5rec_tr");
    fProcReqs.push_back(phyreqs);

    if (G2PRunBase::Init()!=0) return fStatus;

    return (fStatus = kOK);
}

// int G2PRun::RunSim()
// {
//     static const char* const here = "Run()";

//     double x_tr, y_tr, z_tr;

//     if (!pGun->Shoot(fV5beam_lab, fV5react_tr)) return -1;

//     fXSinit = CalXS(fV5beam_lab, fV5react_tr, fThetainit);

//     pBPM->GetBPMValue(fV5beam_lab, fV5bpm_lab);
//     double V5[5];
//     pBPM->TransBPM2Lab(fV5bpm_lab, V5);
//     HCS2TCS(V5[0], V5[2], V5[4], fHRSAngle, fV5bpm_tr[0], fV5bpm_tr[2], z_tr);
//     HCS2TCS(V5[1], V5[3], fHRSAngle, fV5bpm_tr[1], fV5bpm_tr[3]);
//     fV5bpm_tr[4] = 0.0;
//     pDrift->Drift(fV5bpm_tr, fBeamEnergy, z_tr, fHRSAngle, 0.0, 10.0, fV5bpm_tr);
//     if (fDebug>1) {
//         Info(here, "bpm_tr  : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5bpm_tr[0], fV5bpm_tr[1], fV5bpm_tr[2], fV5bpm_tr[3], fV5bpm_tr[4]);
//     }

//     HCS2TCS(fV5beam_lab[0], fV5beam_lab[2], fV5beam_lab[4], fHRSAngle, x_tr, y_tr, z_tr);
//     pDrift->Drift(fV5react_tr, fHRSMomentum, z_tr, fHRSAngle, 0.0, 10.0, fV5tg_tr);
//     if (fDebug>1) {
//         Info(here, "tg_tr   : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5tg_tr[0], fV5tg_tr[1], fV5tg_tr[2], fV5tg_tr[3], fV5tg_tr[4]);
//     }
    
//     pDrift->Drift(fV5tg_tr, fHRSMomentum, 0.0, fHRSAngle, fSieve.fZ, 10.0, fV5sieve_tr);
//     if (fDebug>1) {
//         Info(here, "sieve_tr: %10.3e %10.3e %10.3e %10.3e %10.3e", fV5sieve_tr[0], fV5sieve_tr[1], fV5sieve_tr[2], fV5sieve_tr[3], fV5sieve_tr[4]);
//     }

//     Project(fV5sieve_tr[0], fV5sieve_tr[2], fSieve.fZ, 0.0, fV5sieve_tr[1], fV5sieve_tr[3], fV5tgproj_tr[0], fV5tgproj_tr[2]);
//     fV5tgproj_tr[1] = fV5sieve_tr[1];
//     fV5tgproj_tr[3] = fV5sieve_tr[3];
//     fV5tgproj_tr[4] = fV5sieve_tr[4];
//     fRealXbeam = fV5tgproj_tr[0];
//     if (fDebug>1) {
//         Info(here, "proj_tr : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5tgproj_tr[0], fV5tgproj_tr[1], fV5tgproj_tr[2], fV5tgproj_tr[3], fV5tgproj_tr[4]);
//     }

//     bIsGood = pHRS->Forward(fV5tgproj_tr, fV5fp_tr);
//     pRecUseDB->TransTr2Rot(fV5fp_tr, fV5fp_rot);

//     if (!bIsGood) return 0;

//     SmearVDC(fV5fp_tr);

//     if (bUseEffBPM) {
//         fEffXbeam = GetEffBPM(fV5bpm_tr[0], fV5fp_tr);
//         if (fDebug>1) {
//             Info(here, "bpm_eff : %10.3e %10.3e", fRealXbeam, fEffXbeam);
//         }
//         fV5fp_tr[4] = fEffXbeam;
//         fV5fp_rot[4] = fEffXbeam;
//     }
//     else {
//         fV5fp_tr[4] = fRealXbeam;
//         fV5fp_rot[4] = fRealXbeam;
//     }

//     pRecUseDB->CalcTargetCoords(fV5fp_rot, fV5rectg_tr);
//     Project(fV5rectg_tr[0], fV5rectg_tr[2], 0.0, fSieve.fZ, fV5rectg_tr[1], fV5rectg_tr[3], fV5recsieve_tr[0], fV5recsieve_tr[2]);
//     fV5recsieve_tr[1] = fV5rectg_tr[1];
//     fV5recsieve_tr[3] = fV5rectg_tr[3];
//     fV5recsieve_tr[4] = fV5rectg_tr[4];
//     pDrift->Drift(fV5recsieve_tr, fHRSMomentum, fSieve.fZ, fHRSAngle, 0.0, 10.0, fV5recdb_tr);

//     bIsGood &= pHRS->Backward(fV5fp_tr, fV5rectg_tr);
//     if (fDebug>1) {
//         Info(here, "recj_tr : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5rectg_tr[0], fV5rectg_tr[1], fV5rectg_tr[2], fV5rectg_tr[3], fV5rectg_tr[4]);
//     }

//     Project(fV5rectg_tr[0], fV5rectg_tr[2], 0.0, fSieve.fZ, fV5rectg_tr[1], fV5rectg_tr[3], fV5recsieve_tr[0], fV5recsieve_tr[2]);
//     fV5recsieve_tr[1] = fV5rectg_tr[1];
//     fV5recsieve_tr[3] = fV5rectg_tr[3];
//     fV5recsieve_tr[4] = fV5rectg_tr[4];
//     if (fDebug>1) {
//         Info(here, "rcsiv_tr: %10.3e %10.3e %10.3e %10.3e %10.3e", fV5recsieve_tr[0], fV5recsieve_tr[1], fV5recsieve_tr[2], fV5recsieve_tr[3], fV5recsieve_tr[4]);
//     }

//     pDrift->Drift(fV5recsieve_tr, fHRSMomentum, fSieve.fZ, fHRSAngle, 0.0, 10.0, fV5rec_tr);
//     if (fDebug>1) {
//         Info(here, "rec_tr  : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5rec_tr[0], fV5rec_tr[1], fV5rec_tr[2], fV5rec_tr[3], fV5rec_tr[4]);
//     }

//     fXSrec = CalXS(V5, fV5rec_tr, fThetarec);

//     return 0;
// }

// int G2PRun::RunData()
// {
//     // static const char* const here = "Run()";

//     if (!pGun->Shoot(fV5bpm_lab, fV5tg_tr, fV5fpdata_tr)) return -1;

//     fV5fpdata_tr[4] = -fV5bpm_lab[2];

//     pRecUseDB->TransTr2Rot(fV5fpdata_tr, fV5fpdata_rot);
//     pRecUseDB->CalcTargetCoords(fV5fpdata_rot, fV5recdb_tr);

//     fV5recdb_tr[0] = -fV5bpm_lab[2];
//     bIsGood = pHRS->Forward(fV5recdb_tr, fV5fp_tr);
//     pRecUseDB->TransTr2Rot(fV5fp_tr, fV5fp_rot);

//     bIsGood &= pHRS->Backward(fV5fpdata_tr, fV5rec_tr);

//     return 0;
// }

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
