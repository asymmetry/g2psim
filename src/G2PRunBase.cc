#include <cstdlib>
#include <cmath>

#include "TROOT.h"
#include "TObject.h"
#include "TError.h"

#include "G2PAppsBase.hh"
#include "G2PBPM.hh"
#include "G2PDrift.hh"
#include "G2PGlobals.hh"
#include "G2PGunBase.hh"
#include "G2PPointGun.hh"
#include "G2PHRSTrans.hh"
#include "G2PPhys.hh"
#include "G2PRand.hh"
#include "G2PRecUseDB.hh"

#include "G2PRunBase.hh"

static const double e = 1.60217656535e-19;
static const double kDEG = 3.14159265358979323846/180.0;

G2PRunBase* G2PRunBase::pG2PRunBase = NULL;

G2PRunBase::G2PRunBase() :
    fHRSAngle(5.767*kDEG), fHRSMomentum(2.251),
    fBeamEnergy(2.254), iTargetZ(1), iTargetA(1), fTargetMass(0.0), 
    fEnergyLoss(0.0), iParticlePID(11), fParticleM0(0.511e-3),
    fParticleQ(-1*e)
{
    if (pG2PRunBase) {
        Error("G2PRunBase()", "Only one instance of G2PRunBase allowed.");
        MakeZombie();
        return;
    }
    pG2PRunBase = this;
}

G2PRunBase::~G2PRunBase()
{
    if (pG2PRunBase==this) pG2PRunBase = NULL;
}

G2PAppsBase::EStatus G2PRunBase::Init()
{
    // static const char* const here = "Init()";

    if (G2PAppsBase::Init()) return fStatus;

    return (fStatus = kOK);
}

int G2PRunBase::RegisterModel()
{
    pBPM = G2PBPM::GetInstance();
    if (!pBPM) {
        pBPM = new G2PBPM();
        gG2PApps->Add(pBPM);
    }

    pDrift = G2PDrift::GetInstance();
    if (!pDrift) {
        pDrift = new G2PDrift();
        gG2PApps->Add(pDrift);
    }

    pGun = G2PGunBase::GetInstance();
    if (!pGun) {
        pGun = new G2PPointGun();
        gG2PApps->Add(pGun);
    }

    pHRS = G2PHRSTrans::GetInstance();
    if (!pHRS) {
        pHRS = new G2PHRSTrans("484816");
        gG2PApps->Add(pHRS);
    }

    pPhys = G2PPhys::GetInstance();
    if (!pPhys) {
        pPhys = new G2PPhys("qfs");
        gG2PApps->Add(pPhys);
    }

    pRecUseDB = G2PRecUseDB::GetInstance();
    if (!pRecUseDB) {
        pRecUseDB = new G2PRecUseDB();
        gG2PApps->Add(pRecUseDB);
    }

    fApps->Add(pBPM);
    fApps->Add(pDrift);
    fApps->Add(pGun);
    fApps->Add(pHRS);
    fApps->Add(pPhys);
    fApps->Add(pRecUseDB);

    return 0;   
}

void G2PRunBase::Project(double x, double y, double z, double z_out, double t, double p, double &xout, double &yout)
{
    xout = x + (z_out-z)*tan(t);
    yout = y + (z_out-z)*tan(p);
}

ClassImp(G2PRunBase)
