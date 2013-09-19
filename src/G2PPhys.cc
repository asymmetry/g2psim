// -*- C++ -*-

/* class G2PPhys
 * It will calculate cross sections at reaction point.
 *
 * Meaning of parameters:
 * fPID: incident particle ID, following the PDG definition:
 *       2212 for p        ;   2112 for n     ;   211 for pi+   ;
 *       -211 for pi-      ;   111  for pi0   ;   11  for e-    ;
 *       22   for photon   ;
 * fZ, fA: proton and mass number of the nucleus.
 *
 * Radiative correction parameters:
 * Tb: total radiative length before scattering in radiation length;
 * Ta: total radiative length after scattering in radiation length;
 * If they are not provided, they will be set to 0, which means not taking target radiative length into account;
 *
 * QFS model parameters:
 * EPS: separation energy in MeV;
 * EPSD: delta separation energy in MeV;
 * FP - Fermi momentum in MeV/c;
 * If they are not provided, they will be set to the recommended value (EPS=10, EPSD=-10, FP=220);
 */

// History:
//   Mar 2013, C. Gu, First public version.
//   Apr 2013, C. Gu, Add P. Bosted's model.
//   May 2013, C. Gu, Add L. Cardman's C12 elastic model.
//

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <map>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"

#include "G2PPhysBase.hh"
#include "G2PPhysEl/G2PPhysEl.hh"
#include "G2PPhysPB/G2PPhysPB.hh"
#include "G2PPhysQFS/G2PPhysQFS.hh"
#include "G2PPhysWISER/G2PPhysWISER.hh"

#include "G2PAppBase.hh"
#include "G2PAppList.hh"
#include "G2PGlobals.hh"
#include "G2PProcBase.hh"
#include "G2PVar.hh"
#include "G2PVarDef.hh"
#include "G2PVarList.hh"

#include "G2PPhys.hh"

using namespace std;

static const double kDEG = 3.14159265358979323846 / 180.0;

G2PPhys* G2PPhys::pG2PPhys = NULL;

G2PPhys::G2PPhys() : pModel(NULL) {
    // Only for ROOT I/O
}

G2PPhys::G2PPhys(const char *model) :
fSetting(1), fPID(11), fZ(1), fA(1), fTargetMass(0.0), fPars(NULL), fNPars(0), fHRSAngle(5.767 * kDEG), fHRSMomentum(2.251), fBeamEnergy(2.254), pModel(NULL) {
    if (pG2PPhys) {
        Error("G2PPhys()", "Only one instance of G2PPhys allowed.");
        MakeZombie();
        return;
    }
    pG2PPhys = this;

    fPriority = 7;
    map<string, int> model_map;
    model_map["elastic"] = 1;
    model_map["pbosted"] = 11;
    model_map["qfs"] = 12;
    model_map["wiser"] = 21;

    fSetting = model_map[model];
}

G2PPhys::~G2PPhys() {
    if (pG2PPhys == this) pG2PPhys = NULL;
}

int G2PPhys::Begin() {
    static const char* const here = "Begin()";

    if (G2PProcBase::Begin() != 0) return fStatus;

    switch (fSetting) {
    case 1:
        pModel = new G2PPhysEl();
        break;
    case 11:
        pModel = new G2PPhysPB();
        break;
    case 12:
        pModel = new G2PPhysQFS();
        break;
    case 21:
        pModel = new G2PPhysWISER();
        break;
    default:
        Error(here, "Cannot initialize, invalid setting.");
        return (fStatus = kINITERROR);
        break;
    }

    pModel->SetTarget(fZ, fA);
    pModel->SetParticle(fPID);
    pModel->SetPars(fPars, fNPars);

    return (fStatus = kOK);
}

int G2PPhys::Process() {
    static const char* const here = "Process()";

    double V51[5], V52[5];

    V51[0] = gG2PVars->FindSuffix("gun.beam.l_x")->GetValue();
    V51[1] = gG2PVars->FindSuffix("gun.beam.l_t")->GetValue();
    V51[2] = gG2PVars->FindSuffix("gun.beam.l_y")->GetValue();
    V51[3] = gG2PVars->FindSuffix("gun.beam.l_p")->GetValue();
    V51[4] = gG2PVars->FindSuffix("gun.beam.l_z")->GetValue();

    V52[0] = gG2PVars->FindSuffix("gun.react.x")->GetValue();
    V52[1] = gG2PVars->FindSuffix("gun.react.t")->GetValue();
    V52[2] = gG2PVars->FindSuffix("gun.react.y")->GetValue();
    V52[3] = gG2PVars->FindSuffix("gun.react.p")->GetValue();
    V52[4] = gG2PVars->FindSuffix("gun.react.d")->GetValue();

    fXSreact = CalXS(V51, V52, fTHreact);

    V51[0] = gG2PVars->FindSuffix("bpm.l_x")->GetValue();
    V51[1] = gG2PVars->FindSuffix("bpm.l_t")->GetValue();
    V51[2] = gG2PVars->FindSuffix("bpm.l_y")->GetValue();
    V51[3] = gG2PVars->FindSuffix("bpm.l_p")->GetValue();
    V51[4] = gG2PVars->FindSuffix("bpm.l_z")->GetValue();

    V52[0] = gG2PVars->FindSuffix("tp.rec.x")->GetValue();
    V52[1] = gG2PVars->FindSuffix("tp.rec.t")->GetValue();
    V52[2] = gG2PVars->FindSuffix("tp.rec.y")->GetValue();
    V52[3] = gG2PVars->FindSuffix("tp.rec.p")->GetValue();
    V52[4] = gG2PVars->FindSuffix("tp.rec.d")->GetValue();

    fXSrec = CalXS(V51, V52, fTHrec);

    if (fDebug > 1) {
        Info(here, "phys_react: %10.3e %10.3e", fTHreact / kDEG, fXSreact);
        Info(here, "phys_rec  : %10.3e %10.3e", fTHrec / kDEG, fXSrec);
    }

    return 0;
}

void G2PPhys::Clear() {
    fTHreact = 0.0;
    fXSreact = 0.0;
    fTHrec = 0.0;
    fXSrec = 0.0;
}

double G2PPhys::CalXS(const double* V5lab, const double* V5tr, double& scatangle) {
    double Eb[3] = {sin(V5lab[1]) * cos(V5lab[3]), sin(V5lab[1]) * sin(V5lab[3]), cos(V5lab[1])};

    double theta, phi;
    TCS2HCS(V5tr[1], V5tr[3], fHRSAngle, theta, phi);

    double Ef[3] = {sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)};

    scatangle = acos(Eb[0] * Ef[0] + Eb[1] * Ef[1] + Eb[2] * Ef[2]);

    double Ebval = fBeamEnergy;
    double Efval = (1 + V5tr[4]) * fHRSMomentum;

    return pModel->GetXS(Ebval, Efval, scatangle);
}

int G2PPhys::Configure(EMode mode) {
    if (mode == kREAD || mode == kTWOWAY) {
        if (fIsInit) return 0;
        else fIsInit = true;
    }

    ConfDef confs[] = {
        {"run.debuglevel", "Global Debug Level", kINT, &fDebug},
        {"run.particle.id", "Particle ID", kINT, &fPID},
        {"run.target.z", "Target Z", kINT, &fZ},
        {"run.target.a", "Target A", kINT, &fA},
        {"run.target.mass", "Target Mass", kDOUBLE, &fTargetMass},
        {"run.hrs.angle", "HRS Angle", kDOUBLE, &fHRSAngle},
        {"run.hrs.p0", "HRS Momentum", kDOUBLE, &fHRSMomentum},
        {"run.e0", "Beam Energy", kDOUBLE, &fBeamEnergy},
        {0}
    };

    return ConfigureFromList(confs, mode);
}

int G2PPhys::DefineVariables(EMode mode) {
    if (mode == kDEFINE && fIsSetup) return 0;
    fIsSetup = (mode == kDEFINE);

    VarDef vars[] = {
        {"react.angle", "Real scattering angle", kDOUBLE, &fTHreact},
        {"react.xs", "Cross section with real kins", kDOUBLE, &fXSreact},
        {"rec.angle", "Rec scattering angle", kDOUBLE, &fTHrec},
        {"rec.xs", "Cross section with rec kins", kDOUBLE, &fXSrec},
        {0}
    };

    return DefineVarsFromList(vars, mode);
}

void G2PPhys::MakePrefix() {
    const char* base = "phys";

    G2PAppBase::MakePrefix(base);
}
