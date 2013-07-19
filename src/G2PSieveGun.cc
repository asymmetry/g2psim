// -*- C++ -*-

/* class G2PSieveGun
 * This file defines a class G2PSieveGun.
 * It generates special events which can pass through the holes on sieve slit.
 * G2PProcBase classes will call Shoot() to get reaction point kinematics.
 */

// History:
//   Mar 2013, C. Gu, First public version.
//

#include <cstdlib>
#include <cmath>
#include <vector>

#include "TROOT.h"
#include "TObject.h"
#include "TError.h"

#include "G2PDrift.hh"
#include "G2PGlobals.hh"
#include "G2PGun.hh"
#include "G2PRand.hh"
#include "G2PRunBase.hh"
#include "G2PSieve.hh"

#include "G2PSieveGun.hh"

using namespace std;

static const double kU = 0.93149406121;

G2PSieveGun::G2PSieveGun() :
bUseFast(true), fTargetMass(0.0), fEnergyLoss(0.0), fThreshold(0.0), pfGun(NULL) {
    // Nothing to do
}

G2PSieveGun::~G2PSieveGun() {
    // Nothing to do
}

int G2PSieveGun::Begin() {
    //static const char* const here = "Begin()";

    if (G2PGun::Begin() != 0) return fStatus;

    if (fFieldRatio > 0) bUseFast = false;

    if (bUseFast) pfGun = &G2PSieveGun::ShootFast;
    else pfGun = &G2PSieveGun::ShootNormal;

    fTargetMass = gG2PRun->GetTargetMass() * kU;
    fEnergyLoss = gG2PRun->GetEnergyLoss();

    SetSieve(fHRSAngle);
    double ratio = fSieve.fDLargerHole * fSieve.fDLargerHole / (fSieve.fDHole * fSieve.fDHole);
    fThreshold = (ratio - 1) / ((ratio - 1) * fSieve.nLargerHole + fSieve.nRow * fSieve.nCol);

    return (fStatus = kOK);
}

int G2PSieveGun::ShootNormal(double* V5beam_lab, double* V5react_tr, double* reserved) {
    static const char* const here = "Shoot()";

    bool found = false;

    while (!found) {
        // generate flat distribution
        double X_lab, Y_lab;
        if (fBeamR_lab > 1e-5) {
            do {
                X_lab = pRand->Uniform(-fBeamR_lab, fBeamR_lab);
                Y_lab = pRand->Uniform(-fBeamR_lab, fBeamR_lab);
            } while (X_lab * X_lab + Y_lab * Y_lab > fBeamR_lab * fBeamR_lab);
        }
        else {
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
        V5react_tr[4] = pRand->Uniform(fDeltaLow, fDeltaHigh);

        // drift to sieve slit
        double V5siv_tr[5];
        pDrift->Drift(V5react_tr, fHRSMomentum, Zreact_tr, fHRSAngle, fSieve.fZ, 0, V5siv_tr);

        // check if it passed the sieve slit
        for (int i = 0; i < fSieve.nRow * fSieve.nCol; i++) {
            double dhole = fSieve.fDHole;
            for (int j = 0; j < fSieve.nLargerHole; j++)
                if (i == fSieve.iLargerHole[j]) dhole = fSieve.fDLargerHole;
            int col = i / (fSieve.nRow);
            int row = i % (fSieve.nRow);
            if (fSieve.bOpen[i]) {
                double distance = V5siv_tr[0] - fSieve.fX[row] - fSieve.fXOffset * V5siv_tr[0] - fSieve.fX[row] - fSieve.fXOffset;
                distance += V5siv_tr[2] - fSieve.fY[col] - fSieve.fYOffset * V5siv_tr[2] - fSieve.fY[col] - fSieve.fYOffset;
                distance = sqrt(distance);
                if (distance < dhole / 2.0) {
                    found = true;
                    break;
                }
            }
        }
    }

    if (fDebug > 2) {
        Info(here, "%10.3e %10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e %10.3e %10.3e", V5beam_lab[0], V5beam_lab[1], V5beam_lab[2], V5beam_lab[3], V5react_tr[0], V5react_tr[1], V5react_tr[2], V5react_tr[3], V5react_tr[4]);
    }

    return 0;
}

int G2PSieveGun::ShootFast(double* V5beam_lab, double* V5react_tr, double* reserved) {
    static const char* const here = "Shoot()";

    double dhole;

    int selector;
    do {
        dhole = fSieve.fDHole;
        selector = pRand->Integer(fSieve.nRow * fSieve.nCol);
        double temp = pRand->Uniform();
        for (int i = 0; i < fSieve.nLargerHole; i++)
            if (temp < i * fThreshold) {
                selector = fSieve.iLargerHole[i];
                dhole = fSieve.fDLargerHole;
            }
    } while (fSieve.bOpen[selector] == 0);

    int col = selector / (fSieve.nRow);
    int row = selector % (fSieve.nRow);

    double X_lab, Y_lab;
    if (fBeamR_lab > 1e-5) {
        do {
            X_lab = pRand->Uniform(-fBeamR_lab, fBeamR_lab);
            Y_lab = pRand->Uniform(-fBeamR_lab, fBeamR_lab);
        } while (X_lab * X_lab + Y_lab * Y_lab > fBeamR_lab * fBeamR_lab);
    }
    else {
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

    double Xsiv_tr, Ysiv_tr;
    do {
        Xsiv_tr = pRand->Uniform(-dhole / 2, dhole / 2);
        Ysiv_tr = pRand->Uniform(-dhole / 2, dhole / 2);
    } while (Xsiv_tr * Xsiv_tr + Ysiv_tr * Ysiv_tr > dhole * dhole / 4.0);
    V3sieve_tr[0] = fSieve.fXOffset + fSieve.fX[row] + Xsiv_tr;
    V3sieve_tr[1] = fSieve.fYOffset + fSieve.fY[col] + Ysiv_tr;
    V3sieve_tr[2] = fSieve.fZ;

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

    double Delta = scatmom / fHRSMomentum - 1 - fEnergyLoss / fHRSMomentum;

    V5react_tr[0] = Xreact_tr;
    V5react_tr[1] = Thetareact_tr;
    V5react_tr[2] = Yreact_tr;
    V5react_tr[3] = Phireact_tr;
    V5react_tr[4] = Delta;

    if (fDebug > 2) {
        Info(here, "%10.3e %10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e %10.3e %10.3e", V5beam_lab[0], V5beam_lab[1], V5beam_lab[2], V5beam_lab[3], V5react_tr[0], V5react_tr[1], V5react_tr[2], V5react_tr[3], V5react_tr[4]);
    }

    return 0;
}

ClassImp(G2PSieveGun)
