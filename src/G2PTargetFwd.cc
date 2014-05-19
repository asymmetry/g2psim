// -*- C++ -*-

/* class G2PTargetFwd
 * It simulates the movement of the scatted particles only in the target field without any cuts and energy loss.
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
#include "G2PDrift.hh"
#include "G2PGlobals.hh"
#include "G2PProcBase.hh"
#include "G2PSieve.hh"
#include "G2PVar.hh"
#include "G2PVarDef.hh"
#include "G2PVarList.hh"

#include "G2PTargetFwd.hh"

using namespace std;

G2PTargetFwd* G2PTargetFwd::pG2PTargetFwd = NULL;

G2PTargetFwd::G2PTargetFwd() :
fHRSAngle(0.0), fHRSMomentum(0.0), pDrift(NULL), pSieve(NULL)
{
    if (pG2PTargetFwd) {
        Error("G2PTargetFwd()", "Only one instance of G2PTargetFwd allowed.");
        MakeZombie();
        return;
    }
    pG2PTargetFwd = this;

    fPriority = 3;

    Clear();
}

G2PTargetFwd::~G2PTargetFwd()
{
    if (pG2PTargetFwd == this) pG2PTargetFwd = NULL;
}

int G2PTargetFwd::Init()
{
    //static const char* const here = "Init()";

    if (G2PProcBase::Init() != 0) return fStatus;

    pDrift = static_cast<G2PDrift*> (gG2PApps->Find("G2PDrift"));
    if (!pDrift) {
        pDrift = new G2PDrift();
        gG2PApps->Add(pDrift);
    }

    pSieve = static_cast<G2PSieve*> (gG2PApps->Find("G2PSieve"));
    if (!pSieve) {
        pSieve = new G2PSieve();
        gG2PApps->Add(pSieve);
    }

    return (fStatus = kOK);
}

int G2PTargetFwd::Begin()
{
    //static const char* const here = "Begin()";

    if (G2PProcBase::Begin() != 0) return fStatus;

    return (fStatus = kOK);
}

int G2PTargetFwd::Process()
{
    static const char* const here = "Process()";

    if (fDebug > 2) Info(here, " ");

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

    double x_tr, y_tr, z_tr;

    HCS2TCS(fV5react_lab[0], fV5react_lab[2], fV5react_lab[4], fHRSAngle, x_tr, y_tr, z_tr);

    double t_lab, p_lab;
    TCS2HCS(fV5react_tr[1], fV5react_tr[3], fHRSAngle, t_lab, p_lab);
    double pp = fHRSMomentum * (1 + fV5react_tr[4]);
    double x[3] = {fV5react_lab[0], fV5react_lab[2], fV5react_lab[4]};
    double p[3] = {pp * sin(t_lab) * cos(p_lab), pp * sin(t_lab) * sin(p_lab), pp * cos(t_lab)};

    // Local dump front face
    pDrift->Drift(x, p, 640.0e-3, x, p); // along z direction in lab coordinate
    if ((fabs(x[0]) < 46.0e-3) || (fabs(x[0]) > 87.0e-3)) return -1;
    if ((x[1] < -43.0e-3) || (x[1] > 50.0e-3)) return -1;

    // Local dump back face
    pDrift->Drift(x, p, 790.0e-3, x, p);
    if ((fabs(x[0]) < 58.0e-3) || (fabs(x[0]) > 106.0e-3)) return -1;
    if ((x[1] < -53.0e-3) || (x[1] > 58.0e-3)) return -1;

    pDrift->Drift(fV5react_tr, z_tr, fHRSMomentum, fHRSAngle, pSieve->GetZ(), fV5sieve_tr);

    if (fDebug > 1) {
        Info(here, "sieve_tr  : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5sieve_tr[0], fV5sieve_tr[1], fV5sieve_tr[2], fV5sieve_tr[3], fV5sieve_tr[4]);
    }

    Project(fV5sieve_tr[0], fV5sieve_tr[2], pSieve->GetZ(), 0.0, fV5sieve_tr[1], fV5sieve_tr[3], fV5tpproj_tr[0], fV5tpproj_tr[2]);
    fV5tpproj_tr[1] = fV5sieve_tr[1];
    fV5tpproj_tr[3] = fV5sieve_tr[3];
    fV5tpproj_tr[4] = fV5sieve_tr[4];

    if (fDebug > 1) {
        Info(here, "tpproj_tr : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5tpproj_tr[0], fV5tpproj_tr[1], fV5tpproj_tr[2], fV5tpproj_tr[3], fV5tpproj_tr[4]);
    }

    return 0;
}

void G2PTargetFwd::Clear(Option_t * option)
{
    memset(fV5react_lab, 0, sizeof (fV5react_lab));
    memset(fV5react_tr, 0, sizeof (fV5react_tr));
    memset(fV5sieve_tr, 0, sizeof (fV5sieve_tr));
    memset(fV5tpproj_tr, 0, sizeof (fV5tpproj_tr));

    G2PProcBase::Clear(option);
}

int G2PTargetFwd::Configure(EMode mode)
{
    if (mode == kREAD || mode == kTWOWAY) {
        if (fIsInit) return 0;
        else fIsInit = true;
    }

    ConfDef confs[] = {
        {"run.hrs.angle", "HRS Angle", kDOUBLE, &fHRSAngle},
        {"run.hrs.p0", "HRS Momentum", kDOUBLE, &fHRSMomentum},
        {0}
    };

    return ConfigureFromList(confs, mode);
}

int G2PTargetFwd::DefineVariables(EMode mode)
{
    if (mode == kDEFINE && fIsSetup) return 0;
    fIsSetup = (mode == kDEFINE);

    VarDef vars[] = {
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

void G2PTargetFwd::MakePrefix()
{
    const char* base = "fwd";

    G2PAppBase::MakePrefix(base);
}

ClassImp(G2PTargetFwd)
