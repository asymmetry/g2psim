// -*- C++ -*-

/* class G2PSim
 * It is the main class to do the simulation.
 */

// History:
//   Jan 2013, C. Gu, First public version.
//   Sep 2013, C. Gu, Rewrite the structure of the simulation.
//

#ifndef G2P_SIM_H
#define G2P_SIM_H

#include "TObject.h"

class TFile;

class G2PAppList;
class G2POutput;
class G2PRun;

class G2PSim : public TObject
{
public:
    G2PSim();
    virtual ~G2PSim();

    virtual int Run();

    // Gets

    // Sets
    void SetNEvent(int n);
    void SetOutFile(const char *name);

protected:
    virtual int Begin();
    virtual int End();

    int Process();
    bool IsAllDone(G2PAppList *procs);

    TFile *fFile;
    const char *fOutFile;

    int fN;
    int fIndex;

    int fDebug;

    bool fIsGood;

    G2POutput *pOutput;

    G2PAppList *fApps;
    G2PAppList *fProcs;
    G2PRun *pRun;

private:
    static G2PSim *pG2PSim;

    ClassDef(G2PSim, 0)
};

#endif
