// -*- C++ -*-

/* class G2PRunBase
 * This file defines a class G2PRunBase.
 * It is the base class of g2p run classes.
 * For a particular simulation, it should be derived as a new class.
 * Processes and their input variables should be registered in the new class's Init() method.
 * G2PRun class is an example.
 * G2PSim class will call Process() to process any of registered G2PProcBase classes.
 */

// History:
//   Jan 2013, C. Gu, First public version.
//   May 2013, C. Gu, Add G2PProcBase classes, G2PRunBase is more general.
//

#include <cstdlib>
#include <cmath>
#include <vector>

#include "TROOT.h"
#include "TObject.h"
#include "TError.h"
#include "TList.h"

#include "G2PAppBase.hh"
#include "G2PGlobals.hh"
#include "G2PProcBase.hh"

#include "G2PRunBase.hh"

static const double e = 1.60217656535e-19;
static const double kDEG = 3.14159265358979323846 / 180.0;

G2PRunBase* G2PRunBase::pG2PRunBase = NULL;

G2PRunBase::G2PRunBase() :
fStatus(kNOTINIT), fBeamEnergy(2.254), fHRSAngle(5.767 * kDEG), fHRSMomentum(2.251), iTargetZ(1), iTargetA(1), fTargetMass(0.0), fEnergyLoss(0.0), iParticlePID(11), fParticleM0(0.51099892811e-3), fParticleQ(-1 * e), fProcs(NULL) {
    if (pG2PRunBase) {
        Error("G2PRunBase()", "Only one instance of G2PRunBase allowed.");
        MakeZombie();
        return;
    }
    pG2PRunBase = this;

    fProcs = new TList();
    fProcReqs.clear();
}

G2PRunBase::~G2PRunBase() {
    if (pG2PRunBase == this) pG2PRunBase = NULL;

    TIter next(fProcs);
    while (G2PProcBase * aobj = static_cast<G2PProcBase*> (next())) {
        fProcs->Remove(aobj);
        aobj->Delete();
    }
    delete fProcs;
}

int G2PRunBase::Init() {
    static const char* const here = "Init()";

    Info(here, "Initializing ...");

    EStatus status = kOK;
    TIter next(fProcs);
    while (G2PProcBase * aobj = static_cast<G2PProcBase*> (next())) {
        aobj->SetDebug(fDebug);
        if (aobj->Init() != 0) status = kINITERROR;
    }

    return (fStatus = kOK);
}

int G2PRunBase::Begin() {
    static const char* const here = "Begin()";

    Info(here, "Beginning ...");

    EStatus status = kOK;
    TIter next(fProcs);
    while (G2PProcBase * aobj = static_cast<G2PProcBase*> (next())) {
        if (aobj->Begin() != 0) status = kERROR;
    }

    return (fStatus = status);
}

int G2PRunBase::Process() {
    static const char* const here = "Process()";

    double Vtemp[100];

    TIter prociter(fProcs);
    vector<vector<const char*> >::iterator procreqs = fProcReqs.begin();
    while (G2PProcBase * next = static_cast<G2PProcBase*> (prociter())) {
        if (procreqs == fProcReqs.end()) break;
        for (vector<const char*>::iterator name = procreqs->begin(); name != procreqs->end(); name++) {
            TIter tempiter(fProcs);
            bool found = false;
            while (G2PProcBase * ori = static_cast<G2PProcBase*> (tempiter())) {
                if (ori->GetValue(*name, Vtemp) == 0) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                Error(here, "Critical error, check request vars list.");
                return -1;
            }

            next->SetValue(*name, Vtemp);
        }

        if (next->Process() != 0) return -1;

        procreqs++;
    }

    if ((procreqs != fProcReqs.end()) || (prociter() != NULL)) {
        Error(here, "Critical error, check process list.");
        return -1;
    }

    return 0;
}

void G2PRunBase::Clear() {
    TIter next(fProcs);
    while (G2PProcBase * aobj = static_cast<G2PProcBase*> (next())) {
        aobj->Clear();
    }
}

ClassImp(G2PRunBase)
