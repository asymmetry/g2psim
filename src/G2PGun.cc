// -*- C++ -*-

/* class G2PGun
 * This file defines a class G2PGun.
 * It is the base class of g2p event generator classes.
 * G2PProcBase classes will call Shoot() to get reaction point kinematics.
 * Shoot() is a pure virtual method so each derived class has its own implement.
 * G2PDrift is used in this class.
 */

// History:
//   Mar 2013, C. Gu, First public version.
//

#include <cstdlib>
#include <cmath>

#include "TROOT.h"
#include "TObject.h"
#include "TError.h"

#include "G2PAppBase.hh"
#include "G2PDrift.hh"
#include "G2PField.hh"
#include "G2PGlobals.hh"
#include "G2PRunBase.hh"

#include "G2PGun.hh"

static const double kDEG = 3.14159265358979323846 / 180.0;

G2PGun* G2PGun::pG2PGun = NULL;

G2PGun::G2PGun() :
fHRSAngle(5.767 * kDEG), fHRSMomentum(2.251), fBeamEnergy(2.254), fFieldRatio(0.0), fBeamX_lab(0.0), fBeamY_lab(0.0), fBeamTiltAngle(0.0), fBeamR_lab(0.0), fReactZLow_lab(0.0), fReactZHigh_lab(0.0), fTargetThLow_tr(0.0), fTargetThHigh_tr(0.0), fTargetPhLow_tr(0.0), fTargetPhHigh_tr(0.0), fDeltaLow(0.0), fDeltaHigh(0.0), fSigmaPos_lab(0.0), fSigmaAng_lab(0.0), fSigmaAng_tr(0.0), fSigmaDelta(0.0), pDrift(NULL) {
    if (pG2PGun) {
        Error("G2PGun()", "Only one instance of G2PGun allowed.");
        MakeZombie();
        return;
    }
    pG2PGun = this;
}

G2PGun::~G2PGun() {
    if (pG2PGun == this) pG2PGun = NULL;
}

int G2PGun::Init() {
    //static const char* const here = "Init()";

    if (G2PAppBase::Init() != 0) return fStatus;

    pDrift = G2PDrift::GetInstance();
    if (pDrift == NULL) {
        pDrift = new G2PDrift();
        gG2PApps->Add(pDrift);
    }

    return (fStatus = kOK);
}

int G2PGun::Begin() {
    //static const char* const here = "Begin()";

    if (G2PAppBase::Begin() != 0) return fStatus;

    fHRSAngle = gG2PRun->GetHRSAngle();
    fHRSMomentum = gG2PRun->GetHRSMomentum();
    fBeamEnergy = gG2PRun->GetBeamEnergy();

    G2PField* field = G2PField::GetInstance();
    if (field == NULL) fFieldRatio = 0;
    else fFieldRatio = field->GetRatio();

    SetTiltAngle();

    return (fStatus = kOK);
}

void G2PGun::SetTiltAngle() {
    static const char* const here = "SetTiltAngle()";

    if (fabs(fFieldRatio - 0.5) < 1e-8) {
        if (fabs(fBeamEnergy - 2.254) < 0.2) fBeamTiltAngle = 3.31 * kDEG;
        else if (fabs(fBeamEnergy - 1.706) < 0.2) fBeamTiltAngle = 4.03 * kDEG;
        else if (fabs(fBeamEnergy - 1.159) < 0.2) fBeamTiltAngle = 5.97 * kDEG;
        else fBeamTiltAngle = 0.0;
    }
    else if (fabs(fFieldRatio - 1.0) < 1e-8) {
        if (fabs(fBeamEnergy - 2.254) < 0.2) fBeamTiltAngle = 0.0;
        else if (fabs(fBeamEnergy - 3.355) < 0.2) fBeamTiltAngle = 0.0;
        else fBeamTiltAngle = 0.0;
    }
    else fBeamTiltAngle = 0.0;

    if (fDebug > 0) Info(here, "Beam tilt angle is %10.3e deg.", fBeamTiltAngle / kDEG);
}

void G2PGun::GetReactPoint(double x, double y, double z, double* V5) {
    static const char* const here = "GetReactPoint()";

    double xb[3] = {x, y, 0.0};
    double pb[3] = {0.0, fBeamEnergy * sin(fBeamTiltAngle), fBeamEnergy * cos(fBeamTiltAngle)};

    pDrift->Drift(xb, pb, z, 10.0, xb, pb);
    V5[0] = xb[0];
    V5[1] = acos(pb[2] / fBeamEnergy);
    V5[2] = xb[1];
    V5[3] = atan2(pb[1], pb[0]);
    V5[4] = xb[2];

    if (fDebug > 2) Info(here, "%10.3e %10.3e %10.3e %10.3e %10.3e", V5[0], V5[1], V5[2], V5[3], V5[4]);
}

ClassImp(G2PGun)
