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

#ifndef G2P_SIM_H
#define G2P_SIM_H

#include "G2PProcBase.hh"

class G2PAppList;
class G2POutput;

class G2PSim : public G2PProcBase
{
public:
    G2PSim(const char *filename);
    virtual ~G2PSim();

    int Run();

    virtual int Begin();
    virtual int Process();
    virtual int End();

protected:
    G2PSim(); // Only for ROOT I/O

    bool IsAllDone(G2PAppList *procs);

    virtual int Configure(EMode mode = kTWOWAY);
    virtual int DefineVariables(EMode mode = kDEFINE);
    virtual void MakePrefix();

    int fN;
    int fIndex;
    int fDebug;
    bool fIsGood;

    G2POutput *pOutput;
    G2PAppList *fProcs;

private:
    static G2PSim *pG2PSim;

    ClassDef(G2PSim, 1)
};

#endif
