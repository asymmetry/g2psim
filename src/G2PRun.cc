// -*- C++ -*-

/* class G2PRun
 * Run manager for g2p simulation.
 * Parse the configuration file and store all run parameters.
 * It will allocate G2PFwdProc and G2PBwdProc automatically.
 */

// History:
//   Jan 2013, C. Gu, First public version.
//   May 2013, C. Gu, Add G2PProcBase classes, G2PRunBase is more general.
//   Sep 2013, C. Gu, Combine G2PRunBase and G2PRun to be the new run manager.
//

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <map>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TList.h"

#include "G2PAppBase.hh"
#include "G2PAppList.hh"
#include "G2PFwdProc.hh"
#include "G2PBwdProc.hh"
#include "G2PGlobals.hh"
#include "G2PVarDef.hh"

#include "G2PRun.hh"

using namespace std;

static const double e = 1.60217656535e-19;
static const double kDEG = 3.14159265358979323846 / 180.0;

G2PRun* G2PRun::pG2PRun = NULL;

G2PRun::G2PRun() : fConfigFile(NULL) {
    // Constructor

    if (pG2PRun) {
        Error("G2PRun()", "Only one instance of G2PRun allowed.");
        MakeZombie();
        return;
    }
    pG2PRun = this;

    fConfig.clear();
    fConfig["run.debuglevel"] = 0;
    fConfig["run.hrs.angle"] = 5.767 * kDEG;
    fConfig["run.hrs.p0"] = 2.251;
    fConfig["run.particle.id"] = 11;
    fConfig["run.particle.mass"] = 0.51099892811e-3;
    fConfig["run.particle.charge"] = -1 * e;
    fConfig["run.e0"] = 2.254;
    fConfig["run.target.z"] = 1;
    fConfig["run.target.a"] = 1;
    fConfig["run.target.mass"] = 1.008 * 0.931494028;
    fConfig["field.ratio"] = 0.0;

    gG2PApps->Add(new G2PFwdProc());
    gG2PApps->Add(new G2PBwdProc());

    gG2PRun = this;
}

G2PRun::~G2PRun() {
    // Destructor

    Clear();

    if (pG2PRun == this) pG2PRun = NULL;
    gG2PRun = NULL;
}

int G2PRun::Init() {
    static const char* const here = "Init()";

    Info(here, "Initializing ...");

    return 0;
}

int G2PRun::Begin() {
    // Default does nothing

    return 0;
}

int G2PRun::End() {
    // Default does nothing

    return 0;
}

void G2PRun::Clear() {
    fConfig.clear();

    return;
}

int G2PRun::GetConfig(const ConfDef* item, const char* prefix) {
    // Get value of item

    static const char* const here = "GetConfig()";

    if (!item) {
        Error(here, "Bad variable.");
        return 0;
    }

    if (item->var) {
        double dval;
        string keystr(prefix);
        keystr.append(item->name);
        string key(item->name);
        if (fConfig.count(key) > 0) {
            dval = fConfig[key];
        }
        else if (fConfig.count(keystr) > 0) {
            dval = fConfig[keystr];
        }
        else {
            return 0;
        }

        switch (item->type) {
        case kBOOL:
            *((bool*) item->var) = (bool) dval;
            break;
        case kCHAR:
            *((char*) item->var) = (char) dval;
            break;
        case kINT:
            *((int*) item->var) = (int) dval;
            break;
        case kSHORT:
            *((short*) item->var) = (short) dval;
            break;
        case kLONG:
            *((long*) item->var) = (long) dval;
            break;
        case kFLOAT:
            *((float*) item->var) = (float) dval;
            break;
        case kDOUBLE:
            *((double*) item->var) = dval;
            break;
        }
    }

    return 1;
}

int G2PRun::SetConfig(const ConfDef* item, const char* prefix) {
    // Get value of item

    static const char* const here = "SetConfig()";

    if (!item) {
        Error(here, "Bad variable.");
        return 0;
    }

    if (item->var) {
        double dval;
        switch (item->type) {
        case kBOOL:
            dval = (double) (*((bool*) item->var));
            break;
        case kCHAR:
            dval = (double) (*((char*) item->var));
            break;
        case kINT:
            dval = (double) (*((int*) item->var));
            break;
        case kSHORT:
            dval = (double) (*((short*) item->var));
            break;
        case kLONG:
            dval = (double) (*((long*) item->var));
            break;
        case kFLOAT:
            dval = (double) (*((float*) item->var));
            break;
        case kDOUBLE:
            dval = *((double*) item->var);
            break;
        }
        string keystr(prefix);
        keystr.append(item->name);
        string key(item->name);
        if (fConfig.count(keystr) > 0) {
            fConfig[keystr] = dval;
        }
        else if (fConfig.count(key) > 0) {
            fConfig[key] = dval;
        }
        else {
            fConfig[keystr] = dval;
        }
    }

    return 1;
}

int G2PRun::GetDebugLevel() {
    return (int) (fConfig["run.debuglevel"]);
}

double G2PRun::GetHRSAngle() {
    return fConfig["run.hrs.angle"];
}

double G2PRun::GetHRSMomentum() {
    return fConfig["run.hrs.p0"];
}

double G2PRun::GetBeamEnergy() {
    return fConfig["run.e0"];
}

int G2PRun::GetParticleID() {
    return (int) (fConfig["run.particle.id"]);
}

double G2PRun::GetParticleMass() {
    return fConfig["run.particle.mass"];
}

double G2PRun::GetParticleCharge() {
    return fConfig["run.particle.charge"];
}

int G2PRun::GetTargetZ() {
    return (int) (fConfig["run.target.z"]);
}

int G2PRun::GetTargetA() {
    return (int) (fConfig["run.target.a"]);
}

double G2PRun::GetTargetMass() {
    return fConfig["run.target.mass"];
}

double G2PRun::GetFieldRatio() {
    return fConfig["field.ratio"];
}

//double G2PRun::GetEnergyLoss() {
//    return fEnergyLoss;
//}

void G2PRun::SetConfigFile(const char* file) {
    fConfigFile = file;
}

void G2PRun::SetDebugLevel(int n) {
    fConfig["run.debuglevel"] = (double) n;
}

void G2PRun::SetSeed(unsigned n) {
    G2PAppBase::SetSeed(n);
}

void G2PRun::SetHRSAngle(double angle) {
    fConfig["run.hrs.angle"] = angle;
}

void G2PRun::SetHRSMomentum(double P0) {
    fConfig["run.hrs.p0"] = P0;
}

void G2PRun::SetParticleID(int pid) {
    fConfig["run.particle.id"] = (double) pid;
}

void G2PRun::SetParticleMass(double M0) {
    fConfig["run.particle.mass"] = M0;
}

void G2PRun::SetParticleCharge(double Q) {
    fConfig["run.particle.charge"] = Q;
}

void G2PRun::SetBeamEnergy(double E) {
    fConfig["run.e0"] = E;
}

void G2PRun::SetTarget(int Z, int A) {
    fConfig["run.target.z"] = (double) Z;
    fConfig["run.target.a"] = (double) A;
}

void G2PRun::SetTargetMass(double M) {
    fConfig["run.target.mass"] = M;
}

void G2PRun::SetFieldRatio(double ratio) {
    fConfig["field.ratio"] = ratio;
}

void G2PRun::SetSieve() {
    fConfig["run.sieve.on"] = (double) (true);
}

//void G2PRun::SetEnergyLoss(double E) {
//    fEnergyLoss = E;
//}

ClassImp(G2PRun)
