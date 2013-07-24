// -*- C++ -*-

/* class G2PSim
 * This file defines a class G2PSim.
 * It is the main class of this simulation package.
 */

// History:
//   Jan 2013, C. Gu, First public version.
//

#ifndef G2P_SIM_H
#define G2P_SIM_H

#include "TObject.h"

class TFile;
class TList;
class G2PRunBase;
class G2POutput;

class G2PSim : public TObject {
public:
    G2PSim();
    ~G2PSim();

    void SetNEvent(int n) {
        nEvent = n;
    }
    void SetSeed(int n);

    void SetDebug(int n) {
        fDebug = n;
    }

    void SetOutFile(const char* name) {
        pOutFile = name;
    }

    void SetRun(G2PRunBase* run) {
        pRun = run;
    }

    int Init();
    int Begin();
    int End();

    void Run();

    void Run(int n) {
        nEvent = n;
        Run();
    }

    static G2PSim* GetInstance() {
        return pG2PSim;
    }

protected:
    int fDebug;

    TFile* fFile;
    const char* pOutFile;

    int nEvent;
    int nCounter;

    bool bIsGood;

    TList* fApps;
    G2PRunBase* pRun;

    G2POutput* pOutput;

private:
    static G2PSim* pG2PSim;

    ClassDef(G2PSim, 1)
};

#endif
