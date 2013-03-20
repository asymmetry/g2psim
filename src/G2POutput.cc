#include "TROOT.h"
#include "TObject.h"
#include "TError.h"
#include "TTree.h"

#include "G2PGlobals.hh"
#include "G2PVar.hh"
#include "G2PVarDef.hh"
#include "G2PVarList.hh"

#include "G2POutput.hh"

G2POutput::G2POutput() :
    fNVar(0), fVar(NULL), fTree(NULL)
{
    fVName.clear();
    fVariables.clear();
}

G2POutput::~G2POutput()
{
    if (fTree) delete fTree;
    if (fVar) delete [] fVar;
}

int G2POutput::Init()
{
    static const char* const here = "Init()";

    if (!gG2PVars) {
        Error(here, "Cannot initialize, no global vars.");
        return -1;
    }

    fTree = new TTree("T", "G2PSim Output DST");

    fNVar = 0;
    TIter next(gG2PVars);
    while (G2PVar* pvar = static_cast<G2PVar*>(next())) {
        fNVar++;
        fVName.push_back(pvar->GetName());
    }

    fVar = new double[fNVar];

    for (int i = 0; i < fNVar; i++)
        fTree->Branch(fVName[i], &fVar[i], Form("%s/D", fVName[i]));

    if (Attach()!=0) return -1;

    return 0;
}

int G2POutput::Process()
{
    G2PVar *pvar;
    for (int i = 0; i<fNVar; i++) {
        pvar = fVariables[i];
        if (pvar) fVar[i] = pvar->GetValue();
    }

    if (fTree!=0) fTree->Fill();

    return 0;
}

int G2POutput::End()
{
    if (fTree!=0) fTree->Write();

    return 0;
}

int G2POutput::Attach()
{
    static const char* const here = "Attach()";

    if (!gG2PVars) return -1;

    G2PVar* pvar;
    fVariables.resize(fNVar);

    for (int i = 0; i<fNVar; i++) {
        pvar = gG2PVars->Find(fVName[i]);
        if (pvar) {
            fVariables[i] = pvar;
        }
        else {
            Error(here, "Global variable %s does not exist.", fVName[i]);
        }
    }

    return 0;
}

ClassImp(G2POutput)
