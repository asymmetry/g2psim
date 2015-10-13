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

#ifndef G2P_OUTPUT_H
#define G2P_OUTPUT_H

#include <vector>

#include "TObject.h"

using namespace std;

class TFile;
class TMacro;
class TTree;
class G2PVar;

class G2POutput : public TObject
{
public:
    G2POutput(const char *filename);
    ~G2POutput();

    int Begin();
    int Process();
    int End();

protected:
    G2POutput(); // Only for ROOT I/O
    int Attach();

    TFile *fFile;

    int fNVar;
    double *fVar;
    vector<const char *> fVName;
    vector<G2PVar *> fVariables;

    TTree *fTree;
    TMacro *fConfig;

private:
    ClassDef(G2POutput, 1)
};

#endif
