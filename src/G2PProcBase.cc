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

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"

#include "G2PAppBase.hh"
#include "G2PAppList.hh"
#include "G2PGlobals.hh"
#include "G2PVarDef.hh"
#include "G2PVarList.hh"

#include "G2PProcBase.hh"

using namespace std;

G2PProcBase::G2PProcBase() : fStage(kREADY)
{
    // Nothing to do
}

G2PProcBase::~G2PProcBase()
{
    // Nothing to do
}

int G2PProcBase::Init()
{
    //static const char* const here = "Init()";

    if (G2PAppBase::Init() != 0)
        return fStatus;

    EStatus status = kOK;

    if (DefineVariables(kDEFINE))
        status = kINITERROR;

    return (fStatus = status);
}

G2PProcBase::EStage G2PProcBase::GetStage()
{
    return fStage;
}

void G2PProcBase::SetStage(EStage stage)
{
    fStage = stage;
}

int G2PProcBase::ArrayCopy(double *out, const double *in, int length)
{
    if (!(out && in))
        return -1;

    for (int i = 0; i < length; i++)
        out[i] = in[i];

    return 0;
}

int G2PProcBase::DefineVarsFromList(const VarDef *list, EMode mode) const
{
    // Add or delete global variables in "list" to the global list

    static const char *const here = "DefineVarsFromList()";

    if (!gG2PVars) {
        Warning(here, "No global variable list found.");
        return (mode == kDEFINE ? kINITERROR : kOK);
    }

    if (mode == kDEFINE)
        gG2PVars->DefineVariables(list, fPrefix);
    else if (mode == kDELETE) {
        const VarDef *item;

        while ((item = list++) && item->name)
            gG2PVars->RemoveName(Form("%s%s", fPrefix, item->name));
    } else
        return kINITERROR;

    return kOK;
}

int G2PProcBase::RemoveVariables()
{
    // Default method for removing global variables

    return DefineVariables(kDELETE);
}

ClassImp(G2PProcBase)
