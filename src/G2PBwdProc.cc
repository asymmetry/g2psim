// -*- C++ -*-

/* class G2PBwdProc
 * It simulates the reconstruction of g2p kinematics.
 * G2PDrift, G2PHRS and G2PSieve are used in this class.
 * Input variables: fV5bpm_lab, fV5fp_tr (register in gG2PVars).
 */

// History:
//   Apr 2013, C. Gu, First public version.
//

#include <cstdlib>
#include <cstdio>
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
#include "G2PSieve.hh"
#include "G2PVar.hh"
#include "G2PVarDef.hh"
#include "G2PVarList.hh"

#include "G2PBwdProc.hh"

#define USE_DEFAULT_ENDZ 1

using namespace std;

G2PBwdProc* G2PBwdProc::pG2PBwdProc = NULL;

G2PBwdProc::G2PBwdProc() :
fBeamEnergy(0.0), fHRSAngle(0.0), fHRSMomentum(0.0), fFieldRatio(0.0), pDrift(NULL), pHRS(NULL), pSieve(NULL) {
    if (pG2PBwdProc) {
        Error("G2PBwdProc()", "Only one instance of G2PBwdProc allowed.");
        MakeZombie();
        return;
    }
    pG2PBwdProc = this;

    fPriority = 5;
    Clear();
}

G2PBwdProc::~G2PBwdProc() {
    if (pG2PBwdProc == this) pG2PBwdProc = NULL;
}

int G2PBwdProc::Init() {
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

int G2PBwdProc::Begin() {
    //static const char* const here = "Begin()";

    if (G2PProcBase::Begin() != 0) return fStatus;

    return (fStatus = kOK);
}

int G2PBwdProc::Process() {
    static const char* const here = "Process()";

    if (fDebug > 2) Info(here, " ");

    double V5bpm_lab[5], V5fp_tr[5];

    V5bpm_lab[0] = gG2PVars->FindSuffix("bpm.l_x")->GetValue();
    V5bpm_lab[1] = gG2PVars->FindSuffix("bpm.l_t")->GetValue();
    V5bpm_lab[2] = gG2PVars->FindSuffix("bpm.l_y")->GetValue();
    V5bpm_lab[3] = gG2PVars->FindSuffix("bpm.l_p")->GetValue();
    V5bpm_lab[4] = gG2PVars->FindSuffix("bpm.l_z")->GetValue();

    V5fp_tr[0] = gG2PVars->FindSuffix("fp.x")->GetValue();
    V5fp_tr[1] = gG2PVars->FindSuffix("fp.t")->GetValue();
    V5fp_tr[2] = gG2PVars->FindSuffix("fp.y")->GetValue();
    V5fp_tr[3] = gG2PVars->FindSuffix("fp.p")->GetValue();

    V5fp_tr[4] = GetEffBPM(V5bpm_lab, V5fp_tr);

    if (!pHRS->Backward(V5fp_tr, fV5tpsnake_tr)) return -1;

    if (fDebug > 1) {
        Info(here, "tpsnake_tr: %10.3e %10.3e %10.3e %10.3e %10.3e", fV5tpsnake_tr[0], fV5tpsnake_tr[1], fV5tpsnake_tr[2], fV5tpsnake_tr[3], fV5tpsnake_tr[4]);
    }

    Project(fV5tpsnake_tr[0], fV5tpsnake_tr[2], 0.0, pSieve->GetZ(), fV5tpsnake_tr[1], fV5tpsnake_tr[3], fV5sieveproj_tr[0], fV5sieveproj_tr[2]);
    fV5sieveproj_tr[1] = fV5tpsnake_tr[1];
    fV5sieveproj_tr[3] = fV5tpsnake_tr[3];
    fV5sieveproj_tr[4] = fV5tpsnake_tr[4];

    if (fDebug > 1) {
        Info(here, "sivproj_tr: %10.3e %10.3e %10.3e %10.3e %10.3e", fV5sieveproj_tr[0], fV5sieveproj_tr[1], fV5sieveproj_tr[2], fV5sieveproj_tr[3], fV5sieveproj_tr[4]);
    }

    pDrift->Drift(fV5sieveproj_tr, fHRSMomentum, pSieve->GetZ(), fHRSAngle, 0.0, 10.0, fV5tprec_tr);
    TCS2HCS(fV5tprec_tr[0], fV5tprec_tr[2], 0.0, fHRSAngle, fV5tprec_lab[0], fV5tprec_lab[2], fV5tprec_lab[4]);
    TCS2HCS(fV5tprec_tr[1], fV5tprec_tr[3], fHRSAngle, fV5tprec_lab[1], fV5tprec_lab[3]);

    // FIXME : should use variables to control this
#ifndef USE_DEFAULT_ENDZ
    double x[3] = {fV5tprec_lab[0], fV5tprec_lab[2], fV5tprec_lab[4]};
    double p[3] = {fHRSMomentum * (1 + fV5tprec_tr[4]) * sin(fV5tprec_lab[1]) * cos(fV5tprec_lab[3]),
                   fHRSMomentum * (1 + fV5tprec_tr[4]) * sin(fV5tprec_lab[1]) * sin(fV5tprec_lab[3]),
                   fHRSMomentum * (1 + fV5tprec_tr[4]) * cos(fV5tprec_lab[1])};
    pDrift->Drift(x, p, -13.6271e-3, 10.0, x, p);
    double tempd;
    fV5tprec_lab[0] = x[0];
    fV5tprec_lab[1] = acos(p[2] / (fHRSMomentum * (1 + fV5tprec_tr[4])));
    fV5tprec_lab[2] = x[1];
    fV5tprec_lab[3] = atan2(p[1], p[0]);
    fV5tprec_lab[4] = x[2];
    HCS2TCS(x[0], x[1], x[2], fHRSAngle, fV5tprec_tr[0], fV5tprec_tr[2], tempd);
    HCS2TCS(fV5tprec_lab[1], fV5tprec_lab[2], fHRSAngle, fV5tprec_tr[1], fV5tprec_tr[3]);
#endif

    if (fDebug > 1) {
        Info(here, "tprec_tr  : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5tprec_tr[0], fV5tprec_tr[1], fV5tprec_tr[2], fV5tprec_tr[3], fV5tprec_tr[4]);
    }

    return 0;
}

void G2PBwdProc::Clear() {
    memset(fV5tpsnake_tr, 0, sizeof (fV5tpsnake_tr));
    memset(fV5sieveproj_tr, 0, sizeof (fV5sieveproj_tr));
    memset(fV5tprec_tr, 0, sizeof (fV5tprec_tr));
    memset(fV5tprec_lab, 0, sizeof (fV5tprec_lab));
}

double G2PBwdProc::GetEffBPM(const double* V5bpm_lab, const double* V5fp) {
    static const char* const here = "GetEffBPM()";

    // Fit result:
    //
    // (Xbpm_tr-Xtg_tr) vs Z
    // ([0]+[1]*x)
    // Fitting result of (Xbpm_tr-Xtg_tr) vs Z @ 2.5T
    // p0                        =    -0.011838   +/-   0.00132798
    // p1                        =      49.856    +/-   0.163115
    //
    // (Xtg_tr-Xtgproj_tr) vs P
    // ([0]+[1]/x)
    // Fitting result of (Xtg_tr-Xtgproj_tr) vs P @ 2.5T
    // p0                        =    0.0183611   +/-   0.0105237
    // p1                        =      3.14345   +/-   0.0105453
    // Fitting result of (Xtg_tr-Xtgproj_tr) vs P @ 5.0T
    // p0                        =      0.14139   +/-   0.018683
    // p1                        =      6.11766   +/-   0.0187211
    //

    int save = pDrift->GetDebugLevel();
    if (fDebug <= 3) pDrift->SetDebugLevel(0);
    double V5bpm_tr[5], z_tr;
    HCS2TCS(V5bpm_lab[0], V5bpm_lab[2], V5bpm_lab[4], fHRSAngle, V5bpm_tr[0], V5bpm_tr[2], z_tr);
    HCS2TCS(V5bpm_lab[1], V5bpm_lab[3], fHRSAngle, V5bpm_tr[1], V5bpm_tr[3]);
    V5bpm_tr[4] = 0.0;
    pDrift->Drift(V5bpm_tr, fBeamEnergy, z_tr, fHRSAngle, 0.0, 10.0, V5bpm_tr);

    double xbpm_tr = V5bpm_tr[0];
    double xbpm_tr_eff = xbpm_tr;

    if (fabs(fFieldRatio) > 1e-8) {

        double V5_fp[5] = {V5fp[0], V5fp[1], V5fp[2], V5fp[3], V5fp[4]};
        double V5_tg[5] = {0, 0, 0, 0, 0};

        double p = fHRSMomentum;
        double ratio = fFieldRatio;

        if (ratio < 0.75) xbpm_tr_eff -= (3.14345 / p + 0.0183611) / 1000 * ratio / 0.5;
        else xbpm_tr_eff -= (6.11766 / p + 0.14139) / 1000 * ratio / 1.0;

        V5_fp[4] = xbpm_tr_eff;
        pHRS->Backward(V5_fp, V5_tg);

        p = (1 + V5_tg[4]) * fHRSMomentum;
        xbpm_tr_eff = xbpm_tr;
        if (ratio < 0.75) xbpm_tr_eff -= (3.14345 / p + 0.0183611) / 1000 * ratio / 0.5;
        else xbpm_tr_eff -= (6.11766 / p + 0.14139) / 1000 * ratio / 1.0;

        V5_fp[4] = xbpm_tr_eff;
        pHRS->Backward(V5_fp, V5_tg);

        p = (1 + V5_tg[4]) * fHRSMomentum;
        xbpm_tr_eff = xbpm_tr;
        if (ratio < 0.75) xbpm_tr_eff -= (3.14345 / p + 0.0183611) / 1000 * ratio / 0.5;
        else xbpm_tr_eff -= (6.11766 / p + 0.14139) / 1000 * ratio / 1.0;
    }

    pDrift->SetDebugLevel(save);

    if (fDebug > 2) {
        Info(here, "%10.3e", xbpm_tr_eff);
    }

    return xbpm_tr_eff;
}

int G2PBwdProc::Configure(EMode mode) {
    if (mode == kREAD || mode == kTWOWAY) {
        if (fIsInit) return 0;
        else fIsInit = true;
    }

    ConfDef confs[] = {
        {"run.e0", "Beam Energy", kDOUBLE, &fBeamEnergy},
        {"field.ratio", "Field Ratio", kDOUBLE, &fFieldRatio},
        {"run.hrs.angle", "HRS Angle", kDOUBLE, &fHRSAngle},
        {"run.hrs.p0", "HRS Momentum", kDOUBLE, &fHRSMomentum},
        {0}
    };

    return ConfigureFromList(confs, mode);
}

int G2PBwdProc::DefineVariables(EMode mode) {
    if (mode == kDEFINE && fIsSetup) return 0;
    fIsSetup = (mode == kDEFINE);

    VarDef vars[] = {
        {"tp.snake.x", "SNAKE rec to Target Plane X", kDOUBLE, &fV5tpsnake_tr[0]},
        {"tp.snake.t", "SNAKE rec to Target Plane T", kDOUBLE, &fV5tpsnake_tr[1]},
        {"tp.snake.y", "SNAKE rec to Target Plane Y", kDOUBLE, &fV5tpsnake_tr[2]},
        {"tp.snake.p", "SNAKE rec to Target Plane P", kDOUBLE, &fV5tpsnake_tr[3]},
        {"tp.snake.d", "SNAKE rec to Target Plane D", kDOUBLE, &fV5tpsnake_tr[4]},
        {"sieve.proj.x", "Project to Sieve X", kDOUBLE, &fV5sieveproj_tr[0]},
        {"sieve.proj.t", "Project to Sieve T", kDOUBLE, &fV5sieveproj_tr[1]},
        {"sieve.proj.y", "Project to Sieve Y", kDOUBLE, &fV5sieveproj_tr[2]},
        {"sieve.proj.p", "Project to Sieve P", kDOUBLE, &fV5sieveproj_tr[3]},
        {"sieve.proj.d", "Project to Sieve D", kDOUBLE, &fV5sieveproj_tr[3]},
        {"tp.rec.x", "Rec X", kDOUBLE, &fV5tprec_tr[0]},
        {"tp.rec.t", "Rec T", kDOUBLE, &fV5tprec_tr[1]},
        {"tp.rec.y", "Rec Y", kDOUBLE, &fV5tprec_tr[2]},
        {"tp.rec.p", "Rec P", kDOUBLE, &fV5tprec_tr[3]},
        {"tp.rec.d", "Rec D", kDOUBLE, &fV5tprec_tr[4]},
        {"tp.rec.l_x", "Rec X (lab)", kDOUBLE, &fV5tprec_lab[0]},
        {"tp.rec.l_t", "Rec T (lab)", kDOUBLE, &fV5tprec_lab[1]},
        {"tp.rec.l_y", "Rec Y (lab)", kDOUBLE, &fV5tprec_lab[2]},
        {"tp.rec.l_p", "Rec P (lab)", kDOUBLE, &fV5tprec_lab[3]},
        {"tp.rec.l_z", "Rec Z (lab)", kDOUBLE, &fV5tprec_lab[4]},
        {0}
    };

    return DefineVarsFromList(vars, mode);
}

void G2PBwdProc::MakePrefix() {
    const char* basename = "bwd";

    G2PAppBase::MakePrefix(basename);
}

ClassImp(G2PBwdProc)

//double G2PRun::DriftPath() {
//    if (bUseField) {
//        double x[3] = {fV5bpm_lab[0], fV5bpm_lab[2], fV5bpm_lab[4]};
//        double p[3] = {fBeamEnergy * sin(fV5bpm_lab[1]) * cos(fV5bpm_lab[3]), fBeamEnergy * sin(fV5bpm_lab[1]) * sin(fV5bpm_lab[3]), fBeamEnergy * cos(fV5bpm_lab[1])};
//
//        double xrec[3];
//        HRSTransTCSNHCS::X_TCS2HCS(fV5tprec_tr[0], fV5tprec_tr[2], 0.0, fHRSAngle, xrec[0], xrec[1], xrec[2]);
//        double theta, phi;
//        HRSTransTCSNHCS::P_TCS2HCS(fV5tprec_tr[1], fV5tprec_tr[3], fHRSAngle, theta, phi);
//        double prec[3] = {fHRSMomentum * (1 + fV5tprec_tr[4]) * sin(theta) * cos(phi), fHRSMomentum * (1 + fV5tprec_tr[4]) * sin(theta) * sin(phi), fHRSMomentum * (1 + fV5tprec_tr[4]) * cos(theta)};
//        pDrift->Drift(xrec, prec, 0, 10.0, xrec, prec);
//
//        double z = 0.0;
//        double zmin = 0.0;
//        double distmin = sqrt((x[0] - xrec[0])*(x[0] - xrec[0])+(x[1] - xrec[1])*(x[1] - xrec[1]));
//        for (int i = 1; i < 150; i++) {
//            z = 0.1e-3 * i;
//            pDrift->Drift(x, p, z, 10.0, x, p);
//            pDrift->Drift(xrec, prec, z, 10.0, xrec, prec);
//            //printf("%e\t%e\t%e\t%e\t%e\n", z, x[0], x[1], xrec[0], xrec[1]);
//            double distance = sqrt((x[0] - xrec[0])*(x[0] - xrec[0])+(x[1] - xrec[1])*(x[1] - xrec[1]));
//            if (distance < distmin) {
//                zmin = z;
//                distmin = distance;
//            }
//        }
//
//        // for (int i = 20; i<200; i++) {
//        //     z = 1e-3*i;
//        //     pDrift->Drift(xrec, prec, z, 10.0, xrec, prec);
//        //     printf("%e\t%e\t%e\t%e\t%e\n", z, 1000.0, 1000.0, xrec[0], xrec[1]);
//        // }
//
//        for (int i = 1; i < 150; i++) {
//            z = -0.1e-3 * i;
//            pDrift->Drift(x, p, z, 10.0, x, p);
//            pDrift->Drift(xrec, prec, z, 10.0, xrec, prec);
//            //printf("%e\t%e\t%e\t%e\t%e\n", z, x[0], x[1], xrec[0], xrec[1]);
//            double distance = sqrt((x[0] - xrec[0])*(x[0] - xrec[0])+(x[1] - xrec[1])*(x[1] - xrec[1]));
//            if (distance < distmin) {
//                zmin = z;
//                distmin = distance;
//            }
//        }
//
//        // for (int i = -20; i>-200; i--) {
//        //     z = 1e-3*i;
//        //     pDrift->Drift(x, p, z, 10.0, x, p);
//        //     printf("%e\t%e\t%e\t%e\t%e\n", z, x[0], x[1], 1000.0, 1000.0);
//        // }
//
//        pDrift->Drift(x, p, zmin, 10.0, x, p);
//        pDrift->Drift(xrec, prec, zmin, 10.0, xrec, prec);
//
//        double x_tr, y_tr, z_tr;
//        HRSTransTCSNHCS::X_HCS2TCS(xrec[0], xrec[1], xrec[2], fHRSAngle, x_tr, y_tr, z_tr);
//        return (z_tr);
//    }
//    else return 0;
//}
