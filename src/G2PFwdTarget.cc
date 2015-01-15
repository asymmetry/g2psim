// -*- C++ -*-

/* class G2PFwdTarget
 * It simulates the movement of the scatted particles only in the target field without energy loss.
 * Input variables: fV5tp_tr, fV5react_lab (register in gG2PVars).
 */

// History:
//   May 2014, C. Gu, First public version.
//

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <map>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TString.h"

#include "G2PAppBase.hh"
#include "G2PAppList.hh"
#include "G2PGlobals.hh"
#include "G2PProcBase.hh"
#include "G2PSieve.hh"
#include "G2PVarDef.hh"
#include "G2PVarList.hh"

#include "G2PFwdTarget.hh"

using namespace std;

G2PFwdTarget *G2PFwdTarget::pG2PFwdTarget = NULL;

G2PFwdTarget::G2PFwdTarget() : fSieveOn(false), fHoleID(-1), pSieve(NULL)
{
    if (pG2PFwdTarget) {
        Error("G2PFwdTarget()", "Only one instance of G2PFwdTarget allowed.");
        MakeZombie();
        return;
    }

    pG2PFwdTarget = this;

    fPriority = 3;

    Clear();
}

G2PFwdTarget::~G2PFwdTarget()
{
    if (pG2PFwdTarget == this)
        pG2PFwdTarget = NULL;
}

int G2PFwdTarget::Begin()
{
    //static const char* const here = "Begin()";

    if (G2PProcBase::Begin() != 0)
        return (fStatus = kBEGINERROR);

    pSieve = static_cast<G2PSieve *>(gG2PApps->Find("G2PSieve"));

    return (fStatus = kOK);
}

int G2PFwdTarget::Process()
{
    static const char *const here = "Process()";

    if (fDebug > 2)
        Info(here, " ");

    fV5react_lab[0] = gG2PVars->FindSuffix("react.l_x")->GetValue();
    fV5react_lab[1] = gG2PVars->FindSuffix("react.l_t")->GetValue();
    fV5react_lab[2] = gG2PVars->FindSuffix("react.l_y")->GetValue();
    fV5react_lab[3] = gG2PVars->FindSuffix("react.l_p")->GetValue();
    fV5react_lab[4] = gG2PVars->FindSuffix("react.l_z")->GetValue();

    fV5react_tr[0] = gG2PVars->FindSuffix("react.x")->GetValue();
    fV5react_tr[1] = gG2PVars->FindSuffix("react.t")->GetValue();
    fV5react_tr[2] = gG2PVars->FindSuffix("react.y")->GetValue();
    fV5react_tr[3] = gG2PVars->FindSuffix("react.p")->GetValue();
    fV5react_tr[4] = gG2PVars->FindSuffix("react.d")->GetValue();

    freactz_tr = gG2PVars->FindSuffix("react.z")->GetValue();

    double pp = fHRSMomentum * (1 + fV5react_tr[4]);
    double x[3] = {fV5react_lab[0], fV5react_lab[2], fV5react_lab[4]};
    double p[3] = {pp * sin(fV5react_lab[1]) *cos(fV5react_lab[3]), pp * sin(fV5react_lab[1]) *sin(fV5react_lab[3]), pp * cos(fV5react_lab[1])};

    // Local dump front face
    Drift("forward", x, p, 640.0e-3, x, p); // along z direction in lab coordinate

    if ((fabs(x[0]) < 46.0e-3) || (fabs(x[0]) > 87.0e-3))
        return -1;
    else if ((x[1] < -43.0e-3) || (x[1] > 50.0e-3))
        return -1;

    // Local dump back face
    Drift("forward", x, p, 790.0e-3, x, p);

    if ((fabs(x[0]) < 58.0e-3) || (fabs(x[0]) > 106.0e-3))
        return -1;
    else if ((x[1] < -53.0e-3) || (x[1] > 58.0e-3))
        return -1;

    Drift("forward", fV5react_tr, freactz_tr, pSieve->GetZ(), fV5sieve_tr);

    if (fDebug > 1)
        Info(here, "sieve_tr  : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5sieve_tr[0], fV5sieve_tr[1], fV5sieve_tr[2], fV5sieve_tr[3], fV5sieve_tr[4]);

    if (fSieveOn) {
        if (pSieve->CanPass(fV5sieve_tr, fHoleID))
            return -1;
    }

    Project(fV5sieve_tr, pSieve->GetZ(), 0.0, fV5tpproj_tr);

    if (fDebug > 1)
        Info(here, "tpproj_tr : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5tpproj_tr[0], fV5tpproj_tr[1], fV5tpproj_tr[2], fV5tpproj_tr[3], fV5tpproj_tr[4]);

    return 0;
}

void G2PFwdTarget::Clear(Option_t *opt)
{
    freactz_tr = 0;

    memset(fV5react_lab, 0, sizeof(fV5react_lab));
    memset(fV5react_tr, 0, sizeof(fV5react_tr));
    memset(fV5sieve_tr, 0, sizeof(fV5sieve_tr));
    memset(fV5tpproj_tr, 0, sizeof(fV5tpproj_tr));

    G2PProcBase::Clear(opt);
}

void G2PFwdTarget::SetSieve(const char *opt)
{
    TString str(opt);

    if (str == "in") {
        fSieveOn = true;
        fConfigIsSet.insert((unsigned long) &fSieveOn);
    } else {
        fSieveOn = false;
        fConfigIsSet.insert((unsigned long) &fSieveOn);
    }
}

int G2PFwdTarget::DefineVariables(EMode mode)
{
    if (mode == kDEFINE && fDefined)
        return 0;

    if (G2PProcBase::DefineVariables(mode) != 0)
        return -1;

    VarDef vars[] = {
        {"id", "Hole ID", kINT, &fHoleID},
        {"sieve.x", "Sieve X", kDOUBLE, &fV5sieve_tr[0]},
        {"sieve.t", "Sieve T", kDOUBLE, &fV5sieve_tr[1]},
        {"sieve.y", "Sieve Y", kDOUBLE, &fV5sieve_tr[2]},
        {"sieve.p", "Sieve P", kDOUBLE, &fV5sieve_tr[3]},
        {"sieve.d", "Sieve D", kDOUBLE, &fV5sieve_tr[4]},
        {"tp.proj.x", "Project to Target Plane X", kDOUBLE, &fV5tpproj_tr[0]},
        {"tp.proj.t", "Project to Target Plane T", kDOUBLE, &fV5tpproj_tr[1]},
        {"tp.proj.y", "Project to Target Plane Y", kDOUBLE, &fV5tpproj_tr[2]},
        {"tp.proj.p", "Project to Target Plane P", kDOUBLE, &fV5tpproj_tr[3]},
        {"tp.proj.d", "Project to Target Plane D", kDOUBLE, &fV5tpproj_tr[4]},
        {0}
    };

    return DefineVarsFromList(vars, mode);
}

void G2PFwdTarget::MakePrefix()
{
    const char *base = "fwd";

    G2PAppBase::MakePrefix(base);
}

ClassImp(G2PFwdTarget)
