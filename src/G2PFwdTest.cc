// -*- C++ -*-

/* class G2PFwdTest
 * It simulates the movement of the scatted particles only in the target field without energy loss.
 * Input variables: fV5tp_tr, fV5react_lab (register in gG2PVars).
 */

// History:
//   May 2014, C. Gu, First public version.
//

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <map>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"

#include "G2PAppBase.hh"
#include "G2PAppList.hh"
#include "G2PGlobals.hh"
#include "G2PProcBase.hh"
#include "G2PSieve.hh"
#include "G2PVar.hh"
#include "G2PVarDef.hh"
#include "G2PVarList.hh"

#include "G2PFwdTest.hh"

using namespace std;

G2PFwdTest *G2PFwdTest::pG2PFwdTest = NULL;

G2PFwdTest::G2PFwdTest() : fSieveOn(false), fHoleID(-1), pSieve(NULL)
{
    if (pG2PFwdTest) {
        Error("G2PFwdTest()", "Only one instance of G2PFwdTest allowed.");
        MakeZombie();
        return;
    }

    pG2PFwdTest = this;

    fPriority = 3;

    Clear();
}

G2PFwdTest::~G2PFwdTest()
{
    if (pG2PFwdTest == this)
        pG2PFwdTest = NULL;
}

int G2PFwdTest::Begin()
{
    //static const char* const here = "Begin()";

    if (G2PProcBase::Begin() != 0)
        return (fStatus = kBEGINERROR);

    pSieve = static_cast<G2PSieve *>(gG2PApps->Find("G2PSieve"));

    return (fStatus = kOK);
}

int G2PFwdTest::Process()
{
    static const char *const here = "Process()";

    if (fDebug > 2)
        Info(here, " ");

    if (gG2PVars->FindSuffix("react.x") && gG2PVars->FindSuffix("react.z")) {
        fV5react_tr[0] = gG2PVars->FindSuffix("react.x")->GetValue();
        fV5react_tr[1] = gG2PVars->FindSuffix("react.t")->GetValue();
        fV5react_tr[2] = gG2PVars->FindSuffix("react.y")->GetValue();
        fV5react_tr[3] = gG2PVars->FindSuffix("react.p")->GetValue();
        fV5react_tr[4] = gG2PVars->FindSuffix("react.d")->GetValue();

        freactz_tr = gG2PVars->FindSuffix("react.z")->GetValue();
    }

    double V5troj[5] = {fV5react_tr[0], fV5react_tr[1], fV5react_tr[2], fV5react_tr[3], fV5react_tr[4]};
    double ztroj = freactz_tr;

    // Local dump front face
    double x[3];

    Drift("forward", V5troj, ztroj, 640.0e-3, V5troj, ztroj);
    TCS2HCS(V5troj[0], V5troj[2], ztroj, x[0], x[1], x[2]);

    if ((fabs(x[0]) < 46.0e-3) || (fabs(x[0]) > 87.0e-3) || (x[1] < -43.0e-3) || (x[1] > 50.0e-3))
        return -1;

    // Local dump back face
    Drift("forward", V5troj, ztroj, 790.0e-3, V5troj, ztroj);
    TCS2HCS(V5troj[0], V5troj[2], ztroj, x[0], x[1], x[2]);

    if ((fabs(x[0]) < 58.0e-3) || (fabs(x[0]) > 106.0e-3) || (x[1] < -53.0e-3) || (x[1] > 58.0e-3))
        return -1;

    Drift("forward", V5troj, ztroj, pSieve->GetZ(), fV5sieve_tr);

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

void G2PFwdTest::Clear(Option_t *opt)
{
    freactz_tr = 0;

    memset(fV5react_tr, 0, sizeof(fV5react_tr));
    memset(fV5sieve_tr, 0, sizeof(fV5sieve_tr));
    memset(fV5tpproj_tr, 0, sizeof(fV5tpproj_tr));

    G2PProcBase::Clear(opt);
}

void G2PFwdTest::SetSieve(const char *opt)
{
    string str = opt;

    if (str == "in") {
        fSieveOn = true;
        fConfigIsSet.insert((unsigned long) &fSieveOn);
    } else {
        fSieveOn = false;
        fConfigIsSet.insert((unsigned long) &fSieveOn);
    }
}

int G2PFwdTest::Configure(EMode mode)
{
    if ((mode == kREAD || mode == kTWOWAY) && fConfigured)
        return 0;

    if (G2PProcBase::Configure(mode) != 0)
        return -1;

    ConfDef confs[] = {
        {"sieve", "Sieve Slit Switch", kBOOL, &fSieveOn},
        {0}
    };

    return ConfigureFromList(confs, mode);
}

int G2PFwdTest::DefineVariables(EMode mode)
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

void G2PFwdTest::MakePrefix()
{
    const char *base = "fwd";

    G2PAppBase::MakePrefix(base);
}

ClassImp(G2PFwdTest)
