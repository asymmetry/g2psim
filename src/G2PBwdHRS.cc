// -*- C++ -*-

/* class G2PBwdHRS
 * It simulates the reconstruction of g2p kinematics.
 * G2PDrift, G2PHRS and G2PGeoSieve are used in this class.
 * Input variables: fV5bpm_lab, fV5fp_tr (register in gG2PVars).
 */

// History:
//   Apr 2013, C. Gu, First public version.
//   Apr 2014, C. Gu, New effective bpm fitting.
//

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <map>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"

#include "HRSTransBase.hh"
#include "G2PTrans400016/G2PTrans400016.hh"
#include "G2PTrans400016OLD/G2PTrans400016OLD.hh"
#include "G2PTrans484816/G2PTrans484816.hh"
#include "G2PTrans484816OLD/G2PTrans484816OLD.hh"

#include "G2PAppBase.hh"
#include "G2PAppList.hh"
#include "G2PGlobals.hh"
#include "G2PProcBase.hh"
#include "G2PSieve.hh"
#include "G2PVar.hh"
#include "G2PVarDef.hh"
#include "G2PVarList.hh"

#include "G2PBwdHRS.hh"

using namespace std;

G2PBwdHRS *G2PBwdHRS::pG2PBwdHRS = NULL;

G2PBwdHRS::G2PBwdHRS()
{
    //Only for ROOT I/O
}

G2PBwdHRS::G2PBwdHRS(const char *name) : fFieldRatio(0.0), fSetting(1), frecz_lab(0.0), pSieve(NULL), pModel(NULL)
{
    if (pG2PBwdHRS) {
        Error("G2PBwdHRS()", "Only one instance of G2PBwdHRS allowed.");
        MakeZombie();
        return;
    }

    pG2PBwdHRS = this;

    map<string, int> model_map;
    model_map["484816"] = 1;
    model_map["403216"] = 2;
    model_map["400016"] = 3;
    model_map["484816OLD"] = 11;
    model_map["400016OLD"] = 21;

    fSetting = model_map[name];
    fConfigIsSet.insert((unsigned long) &fSetting);

    memset(fFitPars, 0, sizeof(fFitPars));

    fPriority = 5;

    Clear();
}

G2PBwdHRS::~G2PBwdHRS()
{
    delete pModel;

    if (pG2PBwdHRS == this)
        pG2PBwdHRS = NULL;
}

int G2PBwdHRS::Begin()
{
    static const char *const here = "Begin()";

    if (G2PProcBase::Begin() != 0)
        return (fStatus = kBEGINERROR);

    pSieve = static_cast<G2PSieve *>(gG2PApps->Find("G2PSieve"));

    switch (fSetting) {
    case 1:
        pModel = new G2PTrans484816();
        break;

    case 2:
        pModel = new G2PTrans484816(); // FIXME: should be 403216 here
        break;

    case 3:
        pModel = new G2PTrans400016();
        break;

    case 11:
        pModel = new G2PTrans484816OLD();
        break;

    case 21:
        pModel = new G2PTrans400016OLD();
        break;

    default:
        Error(here, "Invalid setting.");
        return (fStatus = kBEGINERROR);
        break;
    }

    return (fStatus = kOK);
}

int G2PBwdHRS::Process()
{
    static const char *const here = "Process()";

    if (fDebug > 2)
        Info(here, " ");

    if (gG2PVars->FindSuffix("bpm.x") && gG2PVars->FindSuffix("fp.x")) {
        fV5bpm_tr[0] = gG2PVars->FindSuffix("bpm.x")->GetValue();
        fV5bpm_tr[1] = gG2PVars->FindSuffix("bpm.t")->GetValue();
        fV5bpm_tr[2] = gG2PVars->FindSuffix("bpm.y")->GetValue();
        fV5bpm_tr[3] = gG2PVars->FindSuffix("bpm.p")->GetValue();
        fV5bpm_tr[4] = gG2PVars->FindSuffix("bpm.z")->GetValue();

        fV5fp_tr[0] = gG2PVars->FindSuffix("fp.x")->GetValue();
        fV5fp_tr[1] = gG2PVars->FindSuffix("fp.t")->GetValue();
        fV5fp_tr[2] = gG2PVars->FindSuffix("fp.y")->GetValue();
        fV5fp_tr[3] = gG2PVars->FindSuffix("fp.p")->GetValue();
    } else
        return -1;

    fV5fp_tr[4] = fV5bpm_tr[0];

    if (!Backward(fV5fp_tr, fV5tpsnake_tr))
        return -1;

    double delta_old = 1e38;

    while (fabs(fV5tpsnake_tr[4] - delta_old) > 1.e-3) {
        delta_old = fV5tpsnake_tr[4];
        fV5fp_tr[4] = GetEffBPM(0); // Get effective bpm x

        if (!Backward(fV5fp_tr, fV5tpsnake_tr))
            return -1;
    }

    fV5tpsnake_tr[2] = GetEffBPM(1); // Get effective bpm y

    if (fDebug > 1)
        Info(here, "tpsnake_tr: %10.3e %10.3e %10.3e %10.3e %10.3e", fV5tpsnake_tr[0], fV5tpsnake_tr[1], fV5tpsnake_tr[2], fV5tpsnake_tr[3], fV5tpsnake_tr[4]);

    Project(fV5tpsnake_tr, 0.0, pSieve->GetZ(), fV5sieveproj_tr);

    if (fDebug > 1)
        Info(here, "sivproj_tr: %10.3e %10.3e %10.3e %10.3e %10.3e", fV5sieveproj_tr[0], fV5sieveproj_tr[1], fV5sieveproj_tr[2], fV5sieveproj_tr[3], fV5sieveproj_tr[4]);

    Drift("backward", fV5sieveproj_tr, pSieve->GetZ(), 0.0, fV5tprec_tr);
    TCS2HCS(fV5tprec_tr, 0.0, fV5tprec_lab);

    if (fabs(frecz_lab) > 1.0e-5) {
        double z_tr;

        if (fV5tprec_lab[4] < frecz_lab)
            Drift("forward", fV5tprec_tr, 0.0, frecz_lab, fV5tprec_tr, z_tr);
        else
            Drift("backward", fV5tprec_tr, 0.0, frecz_lab, fV5tprec_tr, z_tr);

        TCS2HCS(fV5tprec_tr, z_tr, fV5tprec_lab);
    }

    if (fDebug > 1)
        Info(here, "tprec_tr  : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5tprec_tr[0], fV5tprec_tr[1], fV5tprec_tr[2], fV5tprec_tr[3], fV5tprec_tr[4]);

    return 0;
}

void G2PBwdHRS::Clear(Option_t *opt)
{
    memset(fV5bpm_tr, 0, sizeof(fV5bpm_tr));
    memset(fV5fp_tr, 0, sizeof(fV5fp_tr));
    memset(fV5tpsnake_tr, 0, sizeof(fV5tpsnake_tr));
    memset(fV5sieveproj_tr, 0, sizeof(fV5sieveproj_tr));
    memset(fV5tprec_tr, 0, sizeof(fV5tprec_tr));
    memset(fV5tprec_lab, 0, sizeof(fV5tprec_lab));

    G2PProcBase::Clear(opt);
}

void G2PBwdHRS::SetParsX(const double *pars)
{
    fFitPars[0][0] = pars[0];
    fFitPars[0][1] = pars[1];
    fFitPars[0][2] = pars[2];

    fConfigIsSet.insert((unsigned long) &fFitPars[0][0]);
    fConfigIsSet.insert((unsigned long) &fFitPars[0][1]);
    fConfigIsSet.insert((unsigned long) &fFitPars[0][2]);
}

void G2PBwdHRS::SetParsY(const double *pars)
{
    fFitPars[1][0] = pars[0];
    fFitPars[1][1] = pars[1];
    fFitPars[1][2] = pars[2];

    fConfigIsSet.insert((unsigned long) &fFitPars[1][0]);
    fConfigIsSet.insert((unsigned long) &fFitPars[1][1]);
    fConfigIsSet.insert((unsigned long) &fFitPars[1][2]);
}

void G2PBwdHRS::SetRecZ(double z)
{
    frecz_lab = z;

    fConfigIsSet.insert((unsigned long) &frecz_lab);
}

double G2PBwdHRS::GetEffBPM(int axis)
{
    static const char *const here = "GetEffBPM()";

    double xbpm_tr = fV5bpm_tr[0];
    double ybpm_tr = fV5bpm_tr[2];

    double effbpm_tr;

    if (axis == 0)
        effbpm_tr = xbpm_tr;
    else if (axis == 1)
        effbpm_tr = ybpm_tr;
    else
        return 1e38;

    if (fFieldRatio > 1e-5) {
        // Fit:
        // (Xbpm_tr-Xeffbpm_tr) vs P
        // ([0]+[1]/x)
        double p = (1 + fV5tpsnake_tr[4]) * fHRSMomentum;

        if (axis == 0)
            effbpm_tr = xbpm_tr - (fFitPars[0][0] + (fFitPars[0][1] + fFitPars[0][2] * ybpm_tr) / p) / 1000;
        else if (axis == 1)
            effbpm_tr = ybpm_tr - (fFitPars[1][0] + (fFitPars[1][1] + fFitPars[1][2] * xbpm_tr) / p) / 1000;
    }

    if (fDebug > 2)
        Info(here, "effbpm_tr :%10.3e", effbpm_tr);

    return effbpm_tr;
}

bool G2PBwdHRS::Backward(const double *V5fp_tr, double *V5tp_tr)
{
    static const char *const here = "Backward()";

    // Definition of variables
    // V5fp_tr = {x_fp, theta_fp, y_fp, phi_fp, x_tp};
    // V5tp_tr = {x_tp, theta_tp, y_tp, phi_tp, delta@tp};
    // x_tp does not change

    double V5[5];

    V5[0] = V5fp_tr[0];
    V5[1] = tan(V5fp_tr[1]);
    V5[2] = V5fp_tr[2];
    V5[3] = tan(V5fp_tr[3]);
    V5[4] = V5fp_tr[4];

    if (fHRSAngle > 0) {
        pModel->ReconLeftHRS(V5);
        //pModel->CoordsCorrection(fModelAngle-fHRSAngle, V5);
    } else {
        pModel->ReconRightHRS(V5);
        //pModel->CoordsCorrection(-fModelAngle-fHRSAngle, V5);
    }

    V5tp_tr[0] = V5[0];
    V5tp_tr[1] = atan(V5[1]);
    V5tp_tr[2] = V5[2];
    V5tp_tr[3] = atan(V5[3]);
    V5tp_tr[4] = V5[4];

    bool isgood = false;

    if (V5tp_tr[4] < 1.0)
        isgood = true;

    if (fDebug > 2)
        Info(here, "%10.3e %10.3e %10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e %10.3e %10.3e", V5fp_tr[0], V5fp_tr[1], V5fp_tr[2], V5fp_tr[3], V5fp_tr[4], V5tp_tr[0], V5tp_tr[1], V5tp_tr[2], V5tp_tr[3], V5tp_tr[4]);

    return isgood;
}

int G2PBwdHRS::Configure(EMode mode)
{
    if ((mode == kREAD || mode == kTWOWAY) && fConfigured)
        return 0;

    if (G2PProcBase::Configure(mode) != 0)
        return -1;

    ConfDef confs[] = {
        {"field.ratio", "Field Ratio", kDOUBLE, &fFieldRatio},
        {"fit.x.p0", "Effective X p0", kDOUBLE, &fFitPars[0][0]},
        {"fit.x.p1", "Effective X p1", kDOUBLE, &fFitPars[0][1]},
        {"fit.x.p2", "Effective X p2", kDOUBLE, &fFitPars[0][2]},
        {"fit.y.p0", "Effective Y p0", kDOUBLE, &fFitPars[1][0]},
        {"fit.y.p1", "Effective Y p1", kDOUBLE, &fFitPars[1][1]},
        {"fit.y.p2", "Effective Y p2", kDOUBLE, &fFitPars[1][2]},
        {"l_z", "Rec Z (lab)", kDOUBLE, &frecz_lab},
        {0}
    };

    return ConfigureFromList(confs, mode);
}

int G2PBwdHRS::DefineVariables(EMode mode)
{
    if (mode == kDEFINE && fDefined)
        return 0;

    if (G2PProcBase::DefineVariables(mode) != 0)
        return -1;

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

void G2PBwdHRS::MakePrefix()
{
    const char *base = "bwd";

    G2PAppBase::MakePrefix(base);
}

ClassImp(G2PBwdHRS)
