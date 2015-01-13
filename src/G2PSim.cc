// -*- C++ -*-

/* class G2PSim
 * It is the main class to do the simulation.
 */

// History:
//   Jan 2013, C. Gu, First public version.
//   Sep 2013, C. Gu, Rewrite the structure of the simulation.
//   Jan 2015, C. Gu, Remove function Init().
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

G2PAppList *gG2PApps = new G2PAppList();
G2PAppList *gG2PGeos = new G2PAppList();
G2PRun *gG2PRun = NULL;
G2PVarList *gG2PVars = new G2PVarList();

G2PSim *G2PSim::pG2PSim = NULL;

G2PSim::G2PSim() : fFile(NULL), fOutFile(NULL), fN(50000), fIndex(1), fDebug(1), fIsGood(false), pOutput(NULL), fApps(NULL), fProcs(NULL), pRun(NULL)
{
    if (pG2PSim) {
        Error("G2PSim::G2PSim()", "Only one instance of G2PSim allowed.");
        MakeZombie();
        return;
    }

    pG2PSim = this;
}

G2PSim::~G2PSim()
{
    TIter next(fApps);

    while (G2PAppBase *aobj = static_cast<G2PAppBase *>(next())) {
        fApps->Remove(aobj);
        aobj->Delete();
    }

    delete gG2PApps;
    delete gG2PGeos;
    delete gG2PRun;
    delete gG2PVars;

    if (pOutput)
        delete pOutput;

    if (pG2PSim == this)
        pG2PSim = NULL;
}

int G2PSim::Run()
{
    static const char *const here = "Run()";

    if (Begin() != 0) {
        Error(here, "Critical error, program will stop.");
        return -1;
    }

    while (fIndex <= fN) {
        fIsGood = false;

        if (fDebug > 1)
            Info(here, "Processing event %d ...", fIndex);

        if (Process() == 0)
            fIsGood = true;

        pOutput->Process();

        if ((fIndex % 1000 == 0) && (fDebug <= 1))
            Info(here, "%d events processed ...", fIndex);

        fIndex++;
    }

    End();

    return 0;
}

void G2PSim::SetNEvent(int n)
{
    fN = n;
}

void G2PSim::SetOutFile(const char *name)
{
    fOutFile = name;
}

int G2PSim::Begin()
{
    static const char *const here = "Begin()";

    // Set run manager
    pRun = gG2PRun;

    if (pRun == NULL)
        Error((here), "No run manager found.");
    else if (pRun->Begin() != 0)
        return -1;

    ConfDef debug = {"run.debuglevel", "Debug Level", kINT, &fDebug};
    fDebug = pRun->GetConfig(&debug, "");

    if (fDebug > 0)
        Info(here, "Starting run ......");

    // Set tools
    fApps = gG2PApps;
    TIter next(fApps);

    while (G2PAppBase *aobj = static_cast<G2PAppBase *>(next())) {
        if (aobj->IsZombie()) {
            fApps->Remove(aobj);
            continue;
        }

        if (!aobj->IsInit()) {
            aobj->SetDebugLevel(fDebug);

            if (aobj->Begin() != 0)
                return -1;
        }
    }

    fProcs = fApps->FindList("G2PProcBase");

    if (fDebug > 0)
        pRun->Print();

    // Set output manager
    gG2PVars->DefineByType("event", "Event number", &fIndex, kINT);
    gG2PVars->DefineByType("isgood", "Good event", &fIsGood, kBOOL);

    fFile = new TFile(fOutFile, "RECREATE");
    fFile->cd();

    pOutput = new G2POutput();

    if (pOutput->Begin() != 0)
        return -1;

    if (fDebug > 0)
        Info(here, "Ready to go!");

    return 0;
}

int G2PSim::End()
{
    static const char *const here = "End()";

    TIter next(fApps);

    while (G2PAppBase *aobj = static_cast<G2PAppBase *>(next()))
        aobj->End();

    pRun->End();

    if (fDebug > 0)
        Info(here, "Cleaning ......");

    pOutput->End();
    fFile->Close();

    if (fDebug > 0)
        Info(here, "Run finished!");

    return 0;
}

int G2PSim::Process()
{
    TIter next(fProcs);

    while (G2PProcBase *pobj = static_cast<G2PProcBase *>(next())) {
        pobj->Clear();
        pobj->SetStage(G2PProcBase::kREADY);
    }

    int step = 0;
    int status = 0;

    while (!IsAllDone(fProcs)) {
        G2PAppList *list = fProcs->FindList(step);
        TIter it(list);

        while (G2PProcBase *pobj = static_cast<G2PProcBase *>(it())) {
            if (pobj->Process() != 0) {
                status = -1;
                pobj->SetStage(G2PProcBase::kSTOP);
                break;
            } else
                pobj->SetStage(G2PProcBase::kDONE);
        }

        step++;
    }

    return status;
}

bool G2PSim::IsAllDone(G2PAppList *procs)
{
    static const char *const g2pproc = "G2PProcBase";

    bool done = true;
    TIter next(procs);

    while (TObject *obj = next()) {
        if (obj->IsA()->InheritsFrom(g2pproc)) {
            G2PProcBase *pobj = static_cast<G2PProcBase *>(obj);

            if (pobj->GetStage() == G2PProcBase::kREADY)
                done = false;

            if (pobj->GetStage() == G2PProcBase::kSTOP)
                return true;
        } else {
            Error("Process()", "Processing list contains non G2PProcBase classes.");
            continue;
        }
    }

    return done;
}

ClassImp(G2PSim)
