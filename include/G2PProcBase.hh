// -*- C++ -*-

/* class G2PProcBase
 * Abstract base class for g2p simulation processes.
 * It provides fundamental functions like variable registration.
 * No instance allowed for this class.
 * Derived class must set its own internal variables and register them to the global variable list.
 */

// History:
//   Apr 2013, C. Gu, First public version.
//   Sep 2013, C. Gu, Move global variable functions from G2PAppBase to here.
//

#ifndef G2P_PROCBASE_H
#define G2P_PROCBASE_H

#include "G2PAppBase.hh"

using namespace std;

class G2PProcBase : public G2PAppBase {
public:
    virtual ~G2PProcBase();

    enum EStage {
        kWAIT = 0, kDONE
    };

    virtual int Init();
    virtual int Begin();
    virtual int Process() = 0;
    virtual int End();
    virtual void Clear();

    // Gets
    EStage GetStage();

    // Sets
    void SetStage(EStage stage);

protected:
    G2PProcBase(); // No instance allowed for this class

    int ArrayCopy(double* out, const double* in, int length);

    virtual int Configure(EMode mode = kTWOWAY) = 0;

    // Global variable functions
    virtual int DefineVariables(EMode mode = kDEFINE) = 0;
    int DefineVarsFromList(const VarDef* list, EMode mode = kDEFINE) const;
    virtual int RemoveVariables();

    virtual void MakePrefix() = 0;

    EStage fStage;

private:
    ClassDef(G2PProcBase, 1)
};

// inline functions

inline G2PProcBase::EStage G2PProcBase::GetStage() {
    return fStage;
}

inline void G2PProcBase::SetStage(EStage stage) {
    fStage = stage;
}

#endif
