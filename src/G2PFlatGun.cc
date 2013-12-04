// -*- C++ -*-

/* class G2PFlatGun
 * It generates events in flat distribution.
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

#include "G2PGlobals.hh"
#include "G2PGun.hh"
#include "G2PRand.hh"

#include "G2PFlatGun.hh"

using namespace std;

G2PFlatGun::G2PFlatGun()
{
    // Nothing to do
}

G2PFlatGun::~G2PFlatGun()
{
    // Nothing to do
}

int G2PFlatGun::Shoot(double* V5beam_lab, double* V5react_tr)
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

    V5react_tr[0] = Xreact_tr;
    V5react_tr[1] = pRand->Uniform(fTargetThLow_tr, fTargetThHigh_tr);
    V5react_tr[2] = Yreact_tr;
    V5react_tr[3] = pRand->Uniform(fTargetPhLow_tr, fTargetPhHigh_tr);

    if (fForceElastic) {
        double Pi[3] = {sin(V5beam_lab[1]) * cos(V5beam_lab[3]), sin(V5beam_lab[1]) * sin(V5beam_lab[3]), cos(V5beam_lab[1])};
        double theta, phi;
        TCS2HCS(V5react_tr[1], V5react_tr[3], fHRSAngle, theta, phi);
        double Pf[3] = {sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)};
        double cosscatangle = Pi[0] * Pf[0] + Pi[1] * Pf[1] + Pi[2] * Pf[2];
        double scatmom = (fTargetMass * fBeamEnergy) / (fTargetMass + fBeamEnergy - fBeamEnergy * cosscatangle);
        V5react_tr[4] = scatmom / fHRSMomentum - 1;
    } else {
        V5react_tr[4] = pRand->Uniform(fDeltaLow, fDeltaHigh);
    }

    if (fDebug > 2) {
        Info(here, "%10.3e %10.3e %10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e %10.3e %10.3e", V5beam_lab[0], V5beam_lab[1], V5beam_lab[2], V5beam_lab[3], V5beam_lab[4], V5react_tr[0], V5react_tr[1], V5react_tr[2], V5react_tr[3], V5react_tr[4]);
    }

    return 0;
}

ClassImp(G2PFlatGun)

