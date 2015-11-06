// -*- C++ -*-

/* class G2PRec
 * The real class to do the reconstruction.
 */

// History:
//   Apr 2013, C. Gu, First public version.
//   Feb 2014, C. Gu, Modified for G2PRec.
//   Jan 2015, C. Gu, Move this file back to the g2psim package.
//

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"

#include "G2PAppBase.hh"
#include "G2PAppList.hh"
#include "G2PGlobals.hh"
#include "G2PProcBase.hh"
#include "G2PSieve.hh"
#include "G2PVar.hh"
#include "G2PVarDef.hh"
#include "G2PVarList.hh"

#include "G2PRec.hh"

G2PRec *G2PRec::pG2PRec = NULL;

G2PRec::G2PRec() : fE0(0.0), fFieldRatio(0.0), pSieve(NULL)
{
    if (pG2PRec) {
        Error("G2PRec()", "Only one instance of G2PRec allowed.");
        MakeZombie();
        return;
    }

    pG2PRec = this;

    memset(fFitPars, 0, sizeof(fFitPars));
    memset(fCorT, 0, sizeof(fCorT));
    memset(fCorP, 0, sizeof(fCorP));
    memset(fCorD, 0, sizeof(fCorD));

    fPriority = 1;

    Clear();
}

G2PRec::~G2PRec()
{
    if (pG2PRec == this)
        pG2PRec = NULL;
}

int G2PRec::Begin()
{
    //static const char *const here = "Begin()";

    if (G2PProcBase::Begin() != 0)
        return (fStatus = kBEGINERROR);

    pSieve = static_cast<G2PSieve *>(gG2PApps->Find("G2PSieve"));

    return (fStatus = kOK);
}

int G2PRec::Process()
{
    //static const char *const here = "Process()";

    if (gG2PVars->FindSuffix("bpm.b_x") && gG2PVars->FindSuffix("tp.mat.x")) {
        fV5bpm_bpm[0] = gG2PVars->FindSuffix("bpm.b_x")->GetValue();
        fV5bpm_bpm[1] = gG2PVars->FindSuffix("bpm.b_t")->GetValue();
        fV5bpm_bpm[2] = gG2PVars->FindSuffix("bpm.b_y")->GetValue();
        fV5bpm_bpm[3] = gG2PVars->FindSuffix("bpm.b_p")->GetValue();
        fV5bpm_bpm[4] = gG2PVars->FindSuffix("bpm.b_z")->GetValue();

        fV5tpmat_tr[0] = gG2PVars->FindSuffix("tp.mat.x")->GetValue();
        fV5tpmat_tr[1] = gG2PVars->FindSuffix("tp.mat.t")->GetValue();
        fV5tpmat_tr[2] = gG2PVars->FindSuffix("tp.mat.y")->GetValue();
        fV5tpmat_tr[3] = gG2PVars->FindSuffix("tp.mat.p")->GetValue();
        fV5tpmat_tr[4] = gG2PVars->FindSuffix("tp.mat.d")->GetValue();
    } else
        return -1;

    return Process(fV5bpm_bpm, fV5tpmat_tr, fV5tprec_tr, fV5tprec_lab);
}

int G2PRec::Process(const double *V5bpm_bpm, const double *V5tpmat_tr, double *V5rec_tr, double *V5rec_lab)
{
    static const char *const here = "Process()";

    if (fDebug > 1)
        Info(here, "bpm_bpm   : %10.3e %10.3e %10.3e %10.3e %10.3e", V5bpm_bpm[0], V5bpm_bpm[1], V5bpm_bpm[2], V5bpm_bpm[3], V5bpm_bpm[4]);

    double recz_lab = V5bpm_bpm[4];
    double V5bpm_tr[5] = {0};
    double bpmz_tr = 0.0;

    double bpm_temp[5];
    BPM2HCS(V5bpm_bpm, bpm_temp);
    HCS2TCS(bpm_temp, V5bpm_tr, bpmz_tr);

    V5bpm_tr[4] = fE0 / fHRSMomentum - 1;

    int save = fDebug;
    fDebug = 0;

    if (bpmz_tr < 0.0)
        Drift("forward", V5bpm_tr, bpmz_tr, 0.0, V5bpm_tr); // Drift to target plane (z_tr = 0)
    else
        Drift("backward", V5bpm_tr, bpmz_tr, 0.0, V5bpm_tr);

    fDebug = save;

    if (fDebug > 1) {
        Info(here, "bpm_tr    : %10.3e %10.3e %10.3e %10.3e %10.3e", V5bpm_tr[0], V5bpm_tr[1], V5bpm_tr[2], V5bpm_tr[3], V5bpm_tr[4]);
        Info(here, "tpmat_tr  : %10.3e %10.3e %10.3e %10.3e %10.3e", V5tpmat_tr[0], V5tpmat_tr[1], V5tpmat_tr[2], V5tpmat_tr[3], V5tpmat_tr[4]);
    }

    // first iteration
    GetEffBPM(V5tpmat_tr, V5bpm_tr, bpm_temp);
    Correct(bpm_temp, V5tpmat_tr, fV5tpcorr_tr);

    if (fDebug > 1)
        Info(here, "tpcorr_tr : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5tpcorr_tr[0], fV5tpcorr_tr[1], fV5tpcorr_tr[2], fV5tpcorr_tr[3], fV5tpcorr_tr[4]);

    // second iteration
    GetEffBPM(fV5tpcorr_tr, V5bpm_tr, bpm_temp);
    Correct(bpm_temp, V5tpmat_tr, fV5tpcorr_tr);

    if (fDebug > 1)
        Info(here, "tpcorr_tr : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5tpcorr_tr[0], fV5tpcorr_tr[1], fV5tpcorr_tr[2], fV5tpcorr_tr[3], fV5tpcorr_tr[4]);

    // third iteration
    GetEffBPM(fV5tpcorr_tr, V5bpm_tr, bpm_temp);
    Correct(bpm_temp, V5tpmat_tr, fV5tpcorr_tr);

    GetEffBPM(fV5tpcorr_tr, V5bpm_tr, bpm_temp);
    fV5tpcorr_tr[0] = bpm_temp[0];
    fV5tpcorr_tr[2] = bpm_temp[2];

    if (fDebug > 1)
        Info(here, "tpcorr_tr : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5tpcorr_tr[0], fV5tpcorr_tr[1], fV5tpcorr_tr[2], fV5tpcorr_tr[3], fV5tpcorr_tr[4]);

    Project(fV5tpcorr_tr, 0.0, pSieve->GetZ(), fV5sieveproj_tr);

    if (fDebug > 1)
        Info(here, "sivproj_tr: %10.3e %10.3e %10.3e %10.3e %10.3e", fV5sieveproj_tr[0], fV5sieveproj_tr[1], fV5sieveproj_tr[2], fV5sieveproj_tr[3], fV5sieveproj_tr[4]);

    double l = Drift("backward", fV5sieveproj_tr, pSieve->GetZ(), 0.0, fV5tprec_tr);

    if (l > 2.0) {
        for (int i = 0; i < 5; i++) {
            fV5tprec_tr[i] = 1.e+38;
            fV5tprec_lab[i] = 1.e+38;
        }

        return -1;
    }

    TCS2HCS(fV5tprec_tr, 0.0, fV5tprec_lab);

    if (fDebug > 1)
        Info(here, "tprec_tr  : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5tprec_tr[0], fV5tprec_tr[1], fV5tprec_tr[2], fV5tprec_tr[3], fV5tprec_tr[4]);

    if (fabs(recz_lab) > 1.0e-5) {
        double z_tr;

        if (fV5tprec_lab[4] < recz_lab)
            Drift("forward", fV5tprec_tr, 0.0, recz_lab, V5rec_tr, z_tr);
        else
            Drift("backward", fV5tprec_tr, 0.0, recz_lab, V5rec_tr, z_tr);

        TCS2HCS(V5rec_tr, z_tr, V5rec_lab);
    }

    if (fDebug > 1) {
        Info(here, "rec_tr    : %10.3e %10.3e %10.3e %10.3e %10.3e", V5rec_tr[0], V5rec_tr[1], V5rec_tr[2], V5rec_tr[3], V5rec_tr[4]);
        Info(here, "rec_lab   : %10.3e %10.3e %10.3e %10.3e %10.3e", V5rec_lab[0], V5rec_lab[1], V5rec_lab[2], V5rec_lab[3], V5rec_lab[4]);
    }

    return 0;
}

void G2PRec::Clear(Option_t * /*option*/)
{
    memset(fV5bpm_bpm, 0, sizeof(fV5bpm_bpm));
    memset(fV5tpmat_tr, 0, sizeof(fV5tpmat_tr));
    memset(fV5tpcorr_tr, 0, sizeof(fV5tpcorr_tr));
    memset(fV5sieveproj_tr, 0, sizeof(fV5sieveproj_tr));
    memset(fV5tprec_tr, 0, sizeof(fV5tprec_tr));
    memset(fV5tprec_lab, 0, sizeof(fV5tprec_lab));
    memset(fV5rec_tr, 0, sizeof(fV5rec_tr));
    memset(fV5rec_lab, 0, sizeof(fV5rec_lab));
}

void G2PRec::GetEffBPM(const double *V5tp_tr, const double *V5bpm_tr, double *V5bpmeff_tr)
{
    static const char *const here = "GetEffBPM()";

    V5bpmeff_tr[0] = V5bpm_tr[0];
    V5bpmeff_tr[1] = V5bpm_tr[1];
    V5bpmeff_tr[2] = V5bpm_tr[2];
    V5bpmeff_tr[3] = V5bpm_tr[3];
    V5bpmeff_tr[4] = V5bpm_tr[4];

    if (fFieldRatio > 1e-5) {
        // Fit:
        // (Xbpm_tr-Xeffbpm_tr) vs P
        // ([0]+[1]/x)
        double p = (1 + V5tp_tr[4]) * fHRSMomentum;
        V5bpmeff_tr[0] = V5bpm_tr[0] - (fFitPars[0][0] + (fFitPars[0][1] + fFitPars[0][2] * V5bpm_tr[2]) / p) / 1000;
        V5bpmeff_tr[2] = V5bpm_tr[2] - (fFitPars[1][0] + (fFitPars[1][1] + fFitPars[1][2] * V5bpm_tr[0]) / p) / 1000;
    }

    if (fDebug > 2)
        Info(here, "effbpm_tr : %10.3e %10.3e", V5bpmeff_tr[0], V5bpmeff_tr[2]);
}

void G2PRec::Correct(const double *V5bpm_tr, const double *V5tp_tr, double *V5corr_tr)
{
    V5corr_tr[0] = V5tp_tr[0];
    V5corr_tr[1] = V5tp_tr[1] + fCorT[0] + fCorT[1] * V5bpm_tr[0] + fCorT[2] * V5bpm_tr[2];
    V5corr_tr[2] = V5tp_tr[2];
    V5corr_tr[3] = V5tp_tr[3] + fCorP[0] + fCorP[1] * V5bpm_tr[0] + fCorP[2] * V5bpm_tr[2];
    V5corr_tr[4] = V5tp_tr[4] + fCorD[0] + fCorD[1] * V5bpm_tr[0] + fCorD[2] * V5bpm_tr[2];
}

int G2PRec::Configure(EMode mode)
{
    if ((mode == kREAD || mode == kTWOWAY) && fConfigured)
        return 0;

    if (G2PProcBase::Configure(mode) != 0)
        return -1;

    ConfDef confs[] = {
        {"run.e0", "Beam Energy", kDOUBLE, &fE0},
        {"field.ratio", "Field Ratio", kDOUBLE, &fFieldRatio},
        {"x.p0", "Effective X p0", kDOUBLE, &fFitPars[0][0]},
        {"x.p1", "Effective X p1", kDOUBLE, &fFitPars[0][1]},
        {"x.p2", "Effective X p2", kDOUBLE, &fFitPars[0][2]},
        {"y.p0", "Effective Y p0", kDOUBLE, &fFitPars[1][0]},
        {"y.p1", "Effective Y p1", kDOUBLE, &fFitPars[1][1]},
        {"y.p2", "Effective Y p2", kDOUBLE, &fFitPars[1][2]},
        {"t.p0", "T Const Correction ", kDOUBLE, &fCorT[0]},
        {"t.px", "T Correction vs Target X", kDOUBLE, &fCorT[1]},
        {"t.py", "T Correction vs Target Y", kDOUBLE, &fCorT[2]},
        {"p.p0", "P Const Correction", kDOUBLE, &fCorP[0]},
        {"p.px", "P Correction vs Target X", kDOUBLE, &fCorP[1]},
        {"p.py", "P Correction vs Target Y", kDOUBLE, &fCorP[2]},
        {"d.p0", "D Const Correction", kDOUBLE, &fCorD[0]},
        {"d.px", "D Correction vs Target X", kDOUBLE, &fCorD[1]},
        {"d.py", "D Correction vs Target Y", kDOUBLE, &fCorD[2]},
        {0}
    };

    return ConfigureFromList(confs, mode);
}

int G2PRec::DefineVariables(EMode mode)
{
    if (mode == kDEFINE && fDefined)
        return 0;

    if (G2PProcBase::DefineVariables(mode) != 0)
        return -1;

    VarDef vars[] = {
        {"tp.corr.x", "After Correction X", kDOUBLE, &fV5tpcorr_tr[0]},
        {"tp.corr.t", "After Correction T", kDOUBLE, &fV5tpcorr_tr[1]},
        {"tp.corr.y", "After Correction Y", kDOUBLE, &fV5tpcorr_tr[2]},
        {"tp.corr.p", "After Correction P", kDOUBLE, &fV5tpcorr_tr[3]},
        {"tp.corr.d", "After Correction D", kDOUBLE, &fV5tpcorr_tr[4]},
        {"sieve.proj.x", "Project to Sieve X", kDOUBLE, &fV5sieveproj_tr[0]},
        {"sieve.proj.t", "Project to Sieve T", kDOUBLE, &fV5sieveproj_tr[1]},
        {"sieve.proj.y", "Project to Sieve Y", kDOUBLE, &fV5sieveproj_tr[2]},
        {"sieve.proj.p", "Project to Sieve P", kDOUBLE, &fV5sieveproj_tr[3]},
        {"sieve.proj.d", "Project to Sieve D", kDOUBLE, &fV5sieveproj_tr[4]},
        {"tp.rec.x", "Rec to Target Plane X", kDOUBLE, &fV5tprec_tr[0]},
        {"tp.rec.t", "Rec to Target Plane T", kDOUBLE, &fV5tprec_tr[1]},
        {"tp.rec.y", "Rec to Target Plane Y", kDOUBLE, &fV5tprec_tr[2]},
        {"tp.rec.p", "Rec to Target Plane P", kDOUBLE, &fV5tprec_tr[3]},
        {"tp.rec.d", "Rec to Target Plane D", kDOUBLE, &fV5tprec_tr[4]},
        {"x", "Rec X", kDOUBLE, &fV5rec_tr[0]},
        {"t", "Rec T", kDOUBLE, &fV5rec_tr[1]},
        {"y", "Rec Y", kDOUBLE, &fV5rec_tr[2]},
        {"p", "Rec P", kDOUBLE, &fV5rec_tr[3]},
        {"d", "Rec D", kDOUBLE, &fV5rec_tr[4]},
        {"l_x", "Rec X (lab)", kDOUBLE, &fV5rec_lab[0]},
        {"l_t", "Rec T (lab)", kDOUBLE, &fV5rec_lab[1]},
        {"l_y", "Rec Y (lab)", kDOUBLE, &fV5rec_lab[2]},
        {"l_p", "Rec P (lab)", kDOUBLE, &fV5rec_lab[3]},
        {"l_z", "Rec Z (lab)", kDOUBLE, &fV5rec_lab[4]},
        {0}
    };

    return DefineVarsFromList(vars, mode);
}

void G2PRec::MakePrefix()
{
    const char *base = "rec";

    G2PAppBase::MakePrefix(base);
}

ClassImp(G2PRec)
