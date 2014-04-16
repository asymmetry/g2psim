// -*- C++ -*-

/* class G2PHRS
 * Interface class of HRSTrans package.
 * It provides HRS transport functions.
 * G2PProcBase classes will call Forward() to get focus plane kinematics to focus plane,
 * will call Backward() to get target plane kinematics.
 *
 * The definition of variables and the list of available models can be found in the comments in the body.
 */

// History:
//   Jan 2013, C. Gu, First public version.
//   Feb 2013, C. Gu, Add correction function.
//   Sep 2013, M. Huang, Add G2PTrans484816R15 module
//   April 2014, M. Huang, complete G2PTrans484816R15 module, and add in G2PTrans400016R15 module

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <map>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"

#include "HRSTransBase.hh"
#include "G2PTrans400016/G2PTrans400016.hh"
#include "G2PTrans484816/G2PTrans484816.hh"
#include "G2PTrans484816R00/G2PTrans484816R00.hh"
#include "G2PTrans484816R15/G2PTrans484816R15.hh" 
#include "G2PTrans400016R15/G2PTrans400016R15.hh"

#include "G2PAppBase.hh"
#include "G2PGlobals.hh"

#include "G2PHRS.hh"

using namespace std;

static const double kDEG = 3.14159265358979323846 / 180.0;

G2PHRS* G2PHRS::pG2PHRS = NULL;

G2PHRS::G2PHRS() : pModel(NULL) {
    // Only for ROOT I/O
}

G2PHRS::G2PHRS(const char* name) :
fSetting(1), fHRSAngle(5.767 * kDEG), fModelAngle(5.767 * kDEG), pModel(NULL) {
    if (pG2PHRS) {
        Error("G2PHRS()", "Only one instance of G2PHRS allowed.");
        MakeZombie();
        return;
    }
    pG2PHRS = this;

    map<string, int> model_map;
    model_map["484816"] = 1;
    model_map["400016"] = 3;
    model_map["484816R00"] = 11;
    model_map["484816R15"] = 12;
    model_map["400016R15"] = 13;

    fSetting = model_map[name];
    fConfigIsSet.insert((unsigned long) &fSetting);
}

G2PHRS::G2PHRS(int setting) :
fSetting(setting), fHRSAngle(5.767 * kDEG), fModelAngle(5.767 * kDEG), pModel(NULL) {
    if (pG2PHRS) {
        Error("G2PHRS()", "Only one instance of G2PHRS allowed.");
        MakeZombie();
        return;
    }
    pG2PHRS = this;

    fConfigIsSet.insert((unsigned long) &fSetting);
}

G2PHRS::~G2PHRS() {
    if (pG2PHRS == this) pG2PHRS = NULL;
}

////////////////////////////////////////////////////////////////////////////////
// Transport particles through HRS using SNAKE model
// Use iModelIndex to identify which SNAKE model to be used
// 1: 484816 with shim, 5.65 deg, 3 cm raster, by JJL
// 2: 403216 with shim, 5.65 deg, SNAKE Model not ready yet
// 3: 400016 with shim, 5.65 deg, 3 cm raster, by Min
// Index > 10 means test
// 11: 484816 with shim, 5.76 deg, no raster, by Min
// 12: 484816 with shim, 5.785 deg, 3cm raster, by Min
// 13: 400016 with shim, 5.77 deg, 3 cm raster, by Min
// May add more HRS packages later
////////////////////////////////////////////////////////////////////////////////

int G2PHRS::Begin() {
    static const char* const here = "Begin()";

    if (G2PAppBase::Begin() != 0) return fStatus;

    switch (fSetting) {
    case 1:
        pModel = new G2PTrans484816();
        break;
    case 3:
        pModel = new G2PTrans400016();
        break;
    case 11:
        pModel = new G2PTrans484816R00();
        break;
    case 12:
        pModel = new G2PTrans484816R15();
        break;
    case 13:
        pModel = new G2PTrans400016R15();
        break;
    default:
        Error(here, "Cannot initialize, invalid setting.");
        return (fStatus = kINITERROR);
        break;
    }

    fModelAngle = pModel->GetAngle();

    if (fDebug > 0) Info(here, "Model angle is %10.3e.", fModelAngle / kDEG);

    return (fStatus = kOK);
}

bool G2PHRS::Forward(const double* V5_tg, double* V5_fp) {
    static const char* const here = "Forward()";

    // Definition of variables
    // V5_tg = {x_tg, theta_tg, y_tg, phi_tg, delta@tg};
    // V5_fp = {x_fg, theta_fp, y_fp, phi_fp, delta@tg};
    // delta does not change
    double V5[5];

    V5[0] = V5_tg[0];
    V5[1] = tan(V5_tg[1]);
    V5[2] = V5_tg[2];
    V5[3] = tan(V5_tg[3]);
    V5[4] = V5_tg[4];

    bool isgood = false;

    if (fHRSAngle > 0) {
        //pModel->CoordsCorrection(fHRSAngle-fModelAngle, V5);
        isgood = pModel->TransLeftHRS(V5);
        //pModel->FPCorrLeft(V5_tg, V5);
    }
    else {
        //pModel->CoordsCorrection(fHRSAngle+fModelAngle, V5);
        isgood = pModel->TransRightHRS(V5);
        //pModel->FPCorrRight(V5_tg, V5);
    }

    V5_fp[0] = V5[0];
    V5_fp[1] = atan(V5[1]);
    V5_fp[2] = V5[2];
    V5_fp[3] = atan(V5[3]);
    V5_fp[4] = V5[4];

    if (fDebug > 2) {
        Info(here, "%10.3e %10.3e %10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e %10.3e %10.3e", V5_tg[0], V5_tg[1], V5_tg[2], V5_tg[3], V5_tg[4], V5_fp[0], V5_fp[1], V5_fp[2], V5_fp[3], V5_fp[4]);
    }

    return isgood;
}

bool G2PHRS::Backward(const double* V5_fp, double* V5_tg) {
    static const char* const here = "Backward()";

    // Definition of variables
    // V5_fp = {x_fp, theta_fp, y_fp, phi_fp, x_tg};
    // V5_tg = {x_tg, theta_tg, y_tg, phi_tg, delta@tg};
    // x_tg does not change

    double V5[5];

    V5[0] = V5_fp[0];
    V5[1] = tan(V5_fp[1]);
    V5[2] = V5_fp[2];
    V5[3] = tan(V5_fp[3]);
    V5[4] = V5_fp[4];

    if (fHRSAngle > 0) {
        pModel->ReconLeftHRS(V5);
        //pModel->CoordsCorrection(fModelAngle-fHRSAngle, V5);
    }
    else {
        pModel->ReconRightHRS(V5);
        //pModel->CoordsCorrection(-fModelAngle-fHRSAngle, V5);
    }

    V5_tg[0] = V5[0];
    V5_tg[1] = atan(V5[1]);
    V5_tg[2] = V5[2];
    V5_tg[3] = atan(V5[3]);
    V5_tg[4] = V5[4];

    bool isgood = false;

    if (V5_tg[4] < 1.0) isgood = true;

    if (fDebug > 2) {
        Info(here, "%10.3e %10.3e %10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e %10.3e %10.3e", V5_fp[0], V5_fp[1], V5_fp[2], V5_fp[3], V5_fp[4], V5_tg[0], V5_tg[1], V5_tg[2], V5_tg[3], V5_tg[4]);
    }

    return isgood;
}

int G2PHRS::Configure(EMode mode) {
    if (mode == kREAD || mode == kTWOWAY) {
        if (fIsInit) return 0;
        else fIsInit = true;
    }

    ConfDef confs[] = {
        {"run.hrs.angle", "Beam Energy", kDOUBLE, &fHRSAngle},
        {"model.id", "Setting ID", kINT, &fSetting},
        {0}
    };

    return ConfigureFromList(confs, mode);
}

void G2PHRS::MakePrefix() {
    const char* base = "hrs";

    G2PAppBase::MakePrefix(base);
}

ClassImp(G2PHRS)
