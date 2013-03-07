// This file defines a class G2PGunBase.
// This class is the base class of G2PGun.
// G2PRun class will call Shoot() to get kinematic variables. It is a virtual
//+function so each derived class will have its own method.
//
// History:
//   Mar 2013, C. Gu, First public version.
//

#include <cstdlib>
#include <cmath>

#include "TROOT.h"
#include "TObject.h"
#include "TError.h"

#include "G2PAppsBase.hh"
#include "G2PDrift.hh"
#include "G2PFieldBase.hh"
#include "G2PGlobals.hh"
#include "G2PRunBase.hh"

#include "G2PGunBase.hh"

static const double kDEG = 3.14159265358979323846/180.0;

G2PGunBase* G2PGunBase::pG2PGunBase = NULL;

G2PGunBase::G2PGunBase() :
    fHRSAngle(5.767*kDEG), fHRSMomentum(2.251), fBeamEnergy(2.254),
    fBeamX_lab(0.0), fBeamY_lab(0.0), fBeamTiltAngle(0.0),
    fBeamR_lab(0.0), fReactZLow_lab(0.0), fReactZHigh_lab(0.0),
    fTargetThLow_tr(0.0), fTargetThHigh_tr(0.0), fTargetPhLow_tr(0.0),
    fTargetPhHigh_tr(0.0), fDeltaLow(0.0), fDeltaHigh(0.0),
    fSigmaPos_lab(0.0), fSigmaAng_lab(0.0), fSigmaAng_tr(0.0),
    fSigmaDelta(0.0), pDrift(NULL)
{
    if (pG2PGunBase) {
        Error("G2PGunBase()", "Only one instance of G2PGunBase allowed.");
        MakeZombie();
        return;
    }
    pG2PGunBase = this;
}

G2PGunBase::~G2PGunBase()
{
    if (pG2PGunBase==this) pG2PGunBase = NULL;
}

G2PAppsBase::EStatus G2PGunBase::Init()
{
    static const char* const here = "Init()";

    fHRSAngle = gG2PRun->GetHRSAngle();
    fHRSMomentum = gG2PRun->GetHRSMomentum();
    fBeamEnergy = gG2PRun->GetBeamEnergy();

    if (SetTiltAngle()) {
        Error(here, "Cannot initialize.");
        return (fStatus = kINITERROR);
    }

    return (fStatus = kOK);
}

int G2PGunBase::SetTiltAngle()
{
    static const char* const here = "SetTiltAngle()";

    if (pDrift->GetField()) {
        if (fabs(pDrift->GetField()->GetRatio()-0.5)<1e-8) {
            if (fabs(fBeamEnergy-2.254)<0.2) fBeamTiltAngle = 3.31*kDEG;
            else if (fabs(fBeamEnergy-1.706)<0.2) fBeamTiltAngle = 4.03*kDEG;
            else if (fabs(fBeamEnergy-1.159)<0.2) fBeamTiltAngle = 5.97*kDEG;
            else fBeamTiltAngle = 0.0;
        }
        else if (fabs(pDrift->GetField()->GetRatio()-1.0)<1e-8) {
            if (fabs(fBeamEnergy-2.254)<0.2) fBeamTiltAngle = 0.0; 
            else if (fabs(fBeamEnergy-3.355)<0.2) fBeamTiltAngle = 0.0;
            else fBeamTiltAngle = 0.0;
        }
        else {
            fBeamTiltAngle = 0.0;
        }
    }
    else {
        fBeamTiltAngle = 0.0;
    }

    if (fDebug>0) Info(here, "Beam tile angle is %10.3e deg", fBeamTiltAngle/kDEG);

    return 0;
}

int G2PGunBase::RegisterModel()
{
    pDrift = G2PDrift::GetInstance();
    if (!pDrift) {
        pDrift = new G2PDrift();
        gG2PApps->Add(pDrift);
    }
    fApps->Add(pDrift);

    return 0;
}

void G2PGunBase::GetReactPoint(double x, double y, double z, double* V5)
{
    static const char* const here = "GetReactPoint()";

    double xb[3] = { x, y, 0.0 };
    double pb[3] = { 0.0, fBeamEnergy*sin(fBeamTiltAngle), fBeamEnergy*cos(fBeamTiltAngle) };

    pDrift->Drift(xb, pb, z, 10.0, xb, pb);
    V5[0] = xb[0];
    V5[1] = acos(pb[2]/fBeamEnergy);
    V5[2] = xb[1];
    V5[3] = atan2(pb[1], pb[0]);
    V5[4] = xb[2];

    if (fDebug>2) Info(here, "%10.3e %10.3e %10.3e %10.3e %10.3e", V5[0], V5[1], V5[2], V5[3], V5[4]);
}

ClassImp(G2PGunBase)
