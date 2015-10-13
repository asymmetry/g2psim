// -*- C++ -*-

/* class G2PSim
 * It is the main class to do the simulation.
 */

// History:
//   Jan 2013, C. Gu, First public version.
//   Sep 2013, C. Gu, Rewrite the structure of the simulation.
//   Jan 2015, C. Gu, Remove function Init().
//   Oct 2015, C. Gu, Modify it as a G2PProcBase class.
//

#include <cstdlib>
#include <cstdio>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TClass.h"

#include "G2PAppBase.hh"
#include "G2PAppList.hh"
#include "G2PGlobals.hh"
#include "G2POutput.hh"
#include "G2PProcBase.hh"
#include "G2PRand.hh"
#include "G2PRun.hh"
#include "G2PVarDef.hh"
#include "G2PVarList.hh"

#include "G2PSim.hh"

using namespace std;

G2PSim *G2PSim::pG2PSim = NULL;

G2PSim::G2PSim()
{
    // Only for ROOT I/O
}

G2PSim::G2PSim(const char *filename) : fN(50000), fIndex(1), fDebug(1), fIsGood(false), pOutput(NULL), fProcs(NULL)
{
    if (pG2PSim) {
        Error("G2PSim::G2PSim()", "Only one instance of G2PSim allowed.");
        MakeZombie();
        return;
    }

    pG2PSim = this;

    pOutput = new G2POutput(filename);
}

G2PSim::~G2PSim()
{
    if (pOutput) {
        delete pOutput;
        pOutput = NULL;
    }

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

int G2PSim::Begin()
{
    static const char *const here = "Begin()";

    // Set run manager
    if (gG2PRun == NULL) {
        Error(here, "No run manager found.");
        return (fStatus = kBEGINERROR);
    } else if (gG2PRun->Begin() != 0)
        return (fStatus = kBEGINERROR);

    if (fDebug > 0)
        Info(here, "Starting run ......");

    if (G2PProcBase::Begin() != 0)
        return (fStatus = kBEGINERROR);

    // Set tools
    TIter next(gG2PApps);

    while (G2PAppBase *aobj = static_cast<G2PAppBase *>(next())) {
        if (aobj->IsZombie()) {
            gG2PApps->Remove(aobj);
            continue;
        }

        if (!aobj->IsInit())
            if (aobj->Begin() != 0)
                return (fStatus = kBEGINERROR);
    }

    fProcs = gG2PApps->FindList("G2PProcBase");

    if (fDebug > 0)
        gG2PRun->Print();

    if (pOutput->Begin() != 0)
        return (fStatus = kBEGINERROR);

    if (fDebug > 0)
        Info(here, "Ready to go!");

    return (fStatus = kOK);
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

int G2PSim::End()
{
    static const char *const here = "End()";

    if (fDebug > 0)
        Info(here, "Cleaning ......");

    TIter next(gG2PApps);

    while (G2PAppBase *aobj = static_cast<G2PAppBase *>(next()))
        aobj->End();

    pOutput->End();

    gG2PRun->End();

    if (fDebug > 0)
        Info(here, "Run finished!");

    return 0;
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

int G2PSim::Configure(EMode mode)
{
    if ((mode == kREAD || mode == kTWOWAY) && fConfigured)
        return 0;

    ConfDef confs[] = {
        {"run.n", "Event Amount", kINT, &fN},
        {"run.debuglevel", "Debug Level", kINT, &fDebug},
        {0}
    };

    return ConfigureFromList(confs, mode);
}

int G2PSim::DefineVariables(EMode mode)
{
    if (mode == kDEFINE && fDefined)
        return 0;

    VarDef vars[] = {
        {"event", "Event ID", kINT, &fIndex},
        {"isgood", "Good Event", kBOOL, &fIsGood},
        {0}
    };

    return DefineVarsFromList("", vars, mode);
}

void G2PSim::MakePrefix()
{
    const char *base = "sim";

    G2PAppBase::MakePrefix(base);
}

ClassImp(G2PSim)
