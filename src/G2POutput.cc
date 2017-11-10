// -*- C++ -*-

/* class G2POutput
 * This file defines a class G2POutput.
 * It generate the output root file of the simulation package.
 * G2PProcBase classes define output variables to gG2PVars in method DefineVariables().
 * This class will read gG2PVars and allocate the tree.
 * G2PSim will call Process() to fill a event into the tree.
 */

// History:
//   Apr 2013, C. Gu, First public version.
//

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <sstream>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TFile.h"
#include "TMacro.h"
#include "TTree.h"

#include "G2PGlobals.hh"
#include "G2PRun.hh"
#include "G2PVar.hh"
#include "G2PVarDef.hh"
#include "G2PVarList.hh"

#include "G2POutput.hh"

using namespace std;

G2POutput::G2POutput()
{
    // Only for ROOT I/O
}

G2POutput::G2POutput(const char *filename) : fNVar(0), fVar(NULL), fConfig(NULL)
{
    fFile = new TFile(filename, "RECREATE");

    fVName.clear();
    fVariables.clear();
}

G2POutput::~G2POutput()
{
    if (fFile) {
        delete fFile;
        fFile = NULL;
    }
}

int G2POutput::Begin()
{
    static const char *const here = "Begin()";

    if (!gG2PVars) {
        Error(here, "Cannot initialize, no global vars.");
        return -1;
    }

    fTree = new TTree("T", "G2PSim Output Tree");

    fNVar = 0;
    TIter next(gG2PVars);

    while (G2PVar *pvar = static_cast<G2PVar *>(next())) {
        fNVar++;
        fVName.push_back(pvar->GetName());
    }

    fVar = new double[fNVar];

    for (int i = 0; i < fNVar; i++)
        fTree->Branch(fVName[i], &fVar[i], Form("%s/D", fVName[i]));

    if (Attach() != 0)
        return -1;

    return 0;
}

int G2POutput::Process()
{
    G2PVar *pvar;

    for (int i = 0; i < fNVar; i++) {
        pvar = fVariables[i];

        if (pvar)
            fVar[i] = pvar->GetValue();
    }

    fTree->Fill();

    return 0;
}

int G2POutput::End()
{
    if (fVar) {
        delete[] fVar;
        fVar = NULL;
    }

    fTree->Write();
    delete fTree;
    fTree = NULL;

    fConfig = new TMacro("config");

    ConfDef *conf = NULL;
    int n = gG2PRun->GetConfigList(conf);

    fConfig->AddLine("{");

    for (int i = 0; i < n; i++) {
        ostringstream ostr; // String stream has a better output format
        ostr << conf[i].name << " = " << *((double *)conf[i].var);
        fConfig->AddLine(Form("cout << \"%s\" << endl;", ostr.str().c_str()));
        delete ((double *) conf[i].var);
    }

    fConfig->AddLine("}");

    delete[] conf;
    conf = NULL;

    fConfig->Write();
    delete fConfig;
    fConfig = NULL;

    fFile->Purge();
    fFile->Close();

    return 0;
}

int G2POutput::Attach()
{
    static const char *const here = "Attach()";

    if (!gG2PVars)
        return -1;

    G2PVar *pvar;
    fVariables.resize(fNVar);

    for (int i = 0; i < fNVar; i++) {
        pvar = gG2PVars->Find(fVName[i]);

        if (pvar)
            fVariables[i] = pvar;
        else
            Error(here, "Global variable %s does not exist.", fVName[i]);
    }

    return 0;
}

ClassImp(G2POutput)
