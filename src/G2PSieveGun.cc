// -*- C++ -*-

/* class G2PSieveGun
 * It generates special events which can pass through the holes on sieve slit.
 * Only work for no target field situation.
 */

// History:
//   Mar 2013, C. Gu, First public version.
//

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"

#include "G2PAppList.hh"
#include "G2PGlobals.hh"
#include "G2PGun.hh"
#include "G2PRand.hh"
#include "G2PSieve.hh"
#include "G2PVarDef.hh"

#include "G2PSieveGun.hh"

using namespace std;

G2PSieveGun::G2PSieveGun() : pSieve(NULL)
{
    // Nothing to do
}

G2PSieveGun::~G2PSieveGun()
{
    // Nothing to do
}

int G2PSieveGun::Init()
{
    //static const char* const here = "Init()";

    if (G2PGun::Init() != 0) return fStatus;

    pSieve = static_cast<G2PSieve*> (gG2PApps->Find("G2PSieve"));
    if (!pSieve) {
        pSieve = new G2PSieve();
        gG2PApps->Add(pSieve);
    }

    return (fStatus = kOK);
}

int G2PSieveGun::Shoot(double* V5beam_lab, double* V5react_tr)
{
    static const char* const here = "Shoot()";

    double X_lab, Y_lab;
    if (fBeamR_lab > 1e-5) {
        do {
            X_lab = pRand->Uniform(-fBeamR_lab, fBeamR_lab);
            Y_lab = pRand->Uniform(-fBeamR_lab, fBeamR_lab);
        } while (X_lab * X_lab + Y_lab * Y_lab > fBeamR_lab * fBeamR_lab);
    } else {
        X_lab = 0.0;
        Y_lab = 0.0;
    }

    X_lab += fBeamX_lab;
    Y_lab += fBeamY_lab;
    double Z_lab = pRand->Uniform(fReactZLow_lab, fReactZHigh_lab);
    GetReactPoint(X_lab, Y_lab, Z_lab, V5beam_lab);

    double Xreact_tr, Yreact_tr, Zreact_tr;
    HCS2TCS(V5beam_lab[0], V5beam_lab[2], V5beam_lab[4], fHRSAngle, Xreact_tr, Yreact_tr, Zreact_tr);

    double V3sieve_tr[3];
    double V3pd_tr[3];

    int index = pSieve->GetPos(V3sieve_tr);

    if (fDebug > 2) {
        Info(here, "(%d,%d)", index % pSieve->GetNRow(), index % pSieve->GetNRow());
    }

    V3pd_tr[0] = V3sieve_tr[0] - Xreact_tr;
    V3pd_tr[1] = V3sieve_tr[1] - Yreact_tr;
    V3pd_tr[2] = V3sieve_tr[2] - Zreact_tr;

    double Thetareact_tr = atan(V3pd_tr[0] / V3pd_tr[2]);
    double Phireact_tr = atan(V3pd_tr[1] / V3pd_tr[2]);

    // Calculate delta based on angle
    double V3pd_lab[3];
    TCS2HCS(V3pd_tr[0], V3pd_tr[1], V3pd_tr[2], fHRSAngle, V3pd_lab[0], V3pd_lab[1], V3pd_lab[2]);

    double cosscatangle = V3pd_lab[2] / (sqrt(V3pd_lab[0] * V3pd_lab[0] + V3pd_lab[1] * V3pd_lab[1] + V3pd_lab[2] * V3pd_lab[2]));

    double scatmom = (fTargetMass * fBeamEnergy) / (fTargetMass + fBeamEnergy - fBeamEnergy * cosscatangle);

    double Delta = scatmom / fHRSMomentum - 1;

    V5react_tr[0] = Xreact_tr;
    V5react_tr[1] = Thetareact_tr;
    V5react_tr[2] = Yreact_tr;
    V5react_tr[3] = Phireact_tr;
    V5react_tr[4] = Delta;
    freactZ = Zreact_tr;

    if (fDebug > 2) {
        Info(here, "%10.3e %10.3e %10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e %10.3e %10.3e", V5beam_lab[0], V5beam_lab[1], V5beam_lab[2], V5beam_lab[3], V5beam_lab[4], V5react_tr[0], V5react_tr[1], V5react_tr[2], V5react_tr[3], V5react_tr[4]);
    }

    return 0;
}

ClassImp(G2PSieveGun)
