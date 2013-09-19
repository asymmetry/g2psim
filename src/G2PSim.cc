// -*- C++ -*-

/* class G2PSim
 * This file defines a class G2PSim.
 * It is the main class of this simulation package.
 */

// History:
//   Jan 2013, C. Gu, First public version.
//   Sep 2013, C. Gu, Rewrite the structure of the simulation.
//

#include <cstdlib>
#include <cstdio>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TClass.h"
#include "TFile.h"

#include "G2PAppBase.hh"
#include "G2PAppList.hh"
#include "G2POutput.hh"
#include "G2PProcBase.hh"
#include "G2PRun.hh"
#include "G2PVarDef.hh"
#include "G2PVarList.hh"

#include "G2PSim.hh"

using namespace std;

G2PAppList* gG2PApps = new G2PAppList();
G2PRun* gG2PRun = NULL;
G2PVarList* gG2PVars = new G2PVarList();

G2PSim* G2PSim::pG2PSim = NULL;

G2PSim::G2PSim() :
fFile(NULL), fOutFile(NULL), fN(50000), fIndex(1), fDebug(1), fIsGood(false), pOutput(NULL), fApps(NULL), fProcs(NULL), pRun(NULL) {
    if (pG2PSim) {
        Error("G2PSim::G2PSim()", "Only one instance of G2PSim allowed.");
        MakeZombie();
        return;
    }
    pG2PSim = this;
}

G2PSim::~G2PSim() {
    if (pG2PSim == this) pG2PSim = NULL;

    TIter next(fApps);
    while (G2PAppBase * aobj = static_cast<G2PAppBase*> (next())) {
        fApps->Remove(aobj);
        aobj->Delete();
    }

    if (pOutput) delete pOutput;
}

int G2PSim::Run() {
    static const char* const here = "Run()";

    if (Init() != 0) {
        Error(here, "Cannot initialize, program will stop.");
        return -1;
    }

    if (Begin() != 0) {
        Error(here, "Critical error, program will stop.");
        return -1;
    }

    while (fIndex <= fN) {
        fIsGood = false;
        if (fDebug > 1) Info(here, "Processing event %d ...", fIndex);
        if (Process() == 0) fIsGood = true;
        pOutput->Process();
        if ((fIndex % 100 == 0)&&(fDebug == 0)) Info(here, "%d events processed ...", fIndex);
        fIndex++;
    }

    End();

    return 0;
}

int G2PSim::Init() {
    static const char* const here = "Init()";

    fApps = gG2PApps;
    pRun = gG2PRun;

    if (pRun == NULL) Error((here), "No run manager found.");

    fDebug = pRun->GetDebugLevel();

    if (fDebug > 0) Info(here, "Initialize tools ......");

    if (pRun->Init() != 0) return -1;

    TIter next(fApps);
    while (G2PAppBase * aobj = static_cast<G2PAppBase*> (next())) {
        if (aobj->IsZombie()) {
            fApps->Remove(aobj);
            continue;
        }
        if (aobj->Init() != 0) return -1;
    }

    fProcs = fApps->FindList("G2PProcBase");

    fFile = new TFile(fOutFile, "RECREATE");
    fFile->cd();

    gG2PVars->DefineByType("event", "Event number", &fIndex, kINT);
    gG2PVars->DefineByType("isgood", "Good event", &fIsGood, kBOOL);

    pOutput = new G2POutput();
    if (pOutput->Init() != 0) return -1;

    if (fDebug > 0) Info(here, "Initialize done!");

    return 0;
}

int G2PSim::Begin() {
    static const char* const here = "Begin()";

    if (fDebug > 0) Info(here, "Start run ......");

    if (pRun->Begin() != 0) return -1;

    TIter next(fApps);
    while (G2PAppBase * aobj = static_cast<G2PAppBase*> (next())) {
        if (aobj->Begin() != 0) return -1;
    }

    if (fDebug > 0) Info(here, "Ready to go!");

    return 0;
}

int G2PSim::End() {
    static const char* const here = "End()";

    if (fDebug > 0) Info(here, "Cleaning ......");

    pOutput->End();
    fFile->Close();

    if (fDebug > 0) Info(here, "Run finished!");

    return 0;
}

int G2PSim::Process() {
    TIter ne(fProcs);
    while (G2PProcBase * pobj = static_cast<G2PProcBase*> (ne())) pobj->SetStage(G2PProcBase::kWAIT);

    int step = 0;
    int status = 0;
    while (!IsAllDone(fProcs)) {
        G2PAppList* list = fProcs->FindList(step);
        TIter next(list);
        while (G2PProcBase * pobj = static_cast<G2PProcBase*> (next())) {
            if (pobj->Process() != 0) status = -1;
            pobj->SetStage(G2PProcBase::kDONE);
        }
        step++;
    }

    return status;
}

bool G2PSim::IsAllDone(G2PAppList* procs) {
    static const char* const g2pproc = "G2PProcBase";

    bool done = true;
    TIter next(procs);
    while (TObject * obj = next()) {
        if (obj->IsA()->InheritsFrom(g2pproc)) {
            G2PProcBase* pobj = static_cast<G2PProcBase*> (obj);
            if (pobj->GetStage() == G2PProcBase::kWAIT) done = false;
        }
        else {
            Error("Process()", "Processing list contains non G2PProcBase classes.");
            continue;
        }
    }

    return done;
}

ClassImp(G2PSim)
