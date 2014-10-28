// -*- C++ -*-

/* class G2PSieve
 * This file defines a class G2PSieve.
 * It defines geometry of the sieve slit.
 * G2PGun classes will inherit this class to get sieve slit geometry.
 */

// History:
//   Mar 2013, C. Gu, First public version.
//

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"

#include "G2PAppBase.hh"
#include "G2PRand.hh"
#include "G2PVarDef.hh"

#include "G2PSieve.hh"

//#define USE_CENTRAL_HOLE 1

using namespace std;

G2PSieve* G2PSieve::pG2PSieve = NULL;

G2PSieve::G2PSieve() :
fHRSAngle(0.0), fNRow(7), fNCol(7), fZ(799.6e-3), fXOffset(0), fYOffset(0), fNLargerHole(0), fDHole(0), fDLargerHole(0)
{
    // Constructor
    if (pG2PSieve) {
        Error("G2PSieve()", "Only one instance of G2PSieve allowed.");
        MakeZombie();

        return;
    }
    pG2PSieve = this;

    fX.clear();
    fY.clear();
    fLargerHole.clear();
    fIsOpen.clear();
}

G2PSieve::~G2PSieve()
{

    if (pG2PSieve == this) pG2PSieve = NULL;
}

int G2PSieve::Begin()
{
    //static const char* const here = "Begin()";

    if (G2PAppBase::Begin() != 0) return fStatus;

    if (fHRSAngle > 0) { // left arm
        const double kSIEVEX[7] = {-3 * 13.3096e-3, -2 * 13.3096e-3, -1 * 13.3096e-3, 0.0, 1 * 13.3096e-3, 2 * 13.3096e-3, 3 * 13.3096e-3};
        const double kSIEVEY[7] = {3 * 6.1214e-3, 2 * 6.1214e-3, 1 * 6.1214e-3, 0.0, -1 * 4.7752e-3, -2 * 4.7752e-3, -3 * 4.7752e-3};
        const int kLARGERHOLE[2] = {15, 24};

#ifdef USE_CENTRAL_HOLE
        const int kSIEVEOPEN[49] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
#else
        const int kSIEVEOPEN[49] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
#endif

        fNRow = 7;
        fNCol = 7;
        fNLargerHole = 2;

        fX.clear();
        for (int i = 0; i < fNRow; i++) fX.push_back(kSIEVEX[i]);
        fY.clear();
        for (int i = 0; i < fNCol; i++) fY.push_back(kSIEVEY[i]);
        fZ = 799.60e-3;
        fXOffset = 0.0;
        fYOffset = 0.0;

        fLargerHole.clear();
        for (int i = 0; i < fNLargerHole; i++) fLargerHole.push_back(kLARGERHOLE[i]);
        fIsOpen.clear();
        for (int i = 0; i < fNRow * fNCol; i++) fIsOpen.push_back((kSIEVEOPEN[i] == 1) ? true : false);

        fDHole = 1.3970e-3;
        fDLargerHole = 2.6924e-3;
    } else { // right arm
        const double kSIEVEX[7] = {-3 * 13.3096e-3, -2 * 13.3096e-3, -1 * 13.3096e-3, 0.0, 1 * 13.3096e-3, 2 * 13.3096e-3, 3 * 13.3096e-3};
        const double kSIEVEY[7] = {-3 * 6.1214e-3, -2 * 6.1214e-3, -1 * 6.1214e-3, 0.0, 1 * 4.7752e-3, 2 * 4.7752e-3, 3 * 4.7752e-3};
        const int kLARGERHOLE[2] = {15, 24};
        const int kSIEVEOPEN[49] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

        fNRow = 7;
        fNCol = 7;
        fNLargerHole = 2;

        fX.clear();
        for (int i = 0; i < fNRow; i++) fX.push_back(kSIEVEX[i]);
        fY.clear();
        for (int i = 0; i < fNCol; i++) fY.push_back(kSIEVEY[i]);
        fZ = 799.46e-3;
        fXOffset = 0.0;
        fYOffset = 0.0;

        fLargerHole.clear();
        for (int i = 0; i < fNLargerHole; i++) fLargerHole.push_back(kLARGERHOLE[i]);
        fIsOpen.clear();

        for (int i = 0; i < fNRow * fNCol; i++) fIsOpen.push_back((kSIEVEOPEN[i] == 1) ? true : false);

        fDHole = 1.3970e-3;
        fDLargerHole = 2.6924e-3;
    }

    double ratio = fDLargerHole * fDLargerHole / (fDHole * fDHole);
    fThreshold = (ratio - 1) / ((ratio - 1) * fNLargerHole + fNRow * fNCol);

    return kOK;
}

int G2PSieve::GetPos(double* V3)
{
    int selector;
    do {
        selector = pRand->Integer(fNRow * fNCol);
        double temp = pRand->Uniform();
        for (int i = 0; i < fNLargerHole; i++) {
            if (temp < (i + 1) * fThreshold) {
                selector = fLargerHole[i];
                break;
            }
        }
    } while (!fIsOpen[selector]);

    double dhole = fDHole;
    for (int i = 0; i < fNLargerHole; i++) {
        if (selector == fLargerHole[i])
            dhole = fDLargerHole;
    }

    int col = selector / fNRow;
    int row = selector % fNRow;

    double dx, dy;
    do {
        dx = pRand->Uniform(-dhole / 2, dhole / 2);
        dy = pRand->Uniform(-dhole / 2, dhole / 2);
    } while (dx * dx + dy * dy > dhole * dhole / 4.0);
    V3[0] = fXOffset + fX[row] + dx;
    V3[1] = fYOffset + fY[col] + dy;
    V3[2] = fZ;

    return selector;
}

void G2PSieve::GetPos(int index, double* V3)
{
    int col = index / fNRow;
    int row = index % fNRow;

    V3[0] = fXOffset + fX[row];
    V3[1] = fYOffset + fY[col];
    V3[2] = fZ;
}

int G2PSieve::CanPass(double* V5)
{
    int id = -1;
    for (int i = 0; i < fNRow * fNCol; i++) {
        double dhole = fDHole;
        for (int j = 0; j < fNLargerHole; j++)
            if (i == fLargerHole[j]) dhole = fDLargerHole;
        int col = i / fNRow;
        int row = i % fNRow;
        if (fIsOpen[i]) {
            double dx = V5[0] - fX[row] - fXOffset;
            double dy = V5[2] - fY[col] - fYOffset;
            double distance = sqrt(dx * dx + dy * dy);
            if (distance < dhole / 2.0) {
                id = i;
                break;
            }
        }
    }

    return id;
}

double G2PSieve::GetZ()
{
    return fZ;
}

int G2PSieve::GetNRow()
{
    return fNRow;
}

int G2PSieve::GetNCol()
{
    return fNCol;
}

int G2PSieve::Configure(EMode mode)
{
    if (mode == kREAD || mode == kTWOWAY) {
        if (fIsInit) return 0;

        else fIsInit = true;
    }

    ConfDef confs[] = {
        {"run.hrs.angle", "HRS Angle", kDOUBLE, &fHRSAngle},
        {0}
    };

    return ConfigureFromList(confs, mode);
}

void G2PSieve::MakePrefix()
{
    const char* base = "sieve";

    G2PAppBase::MakePrefix(base);
}

ClassImp(G2PSieve)
