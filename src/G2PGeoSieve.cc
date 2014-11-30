// -*- C++ -*-

/* class G2PGeoSieve
 * This file defines a class G2PGeoSieve.
 * It defines a sieve slit.
 * Derived from G2PGeoBase so G2PDrift can use it as boundary.
 * Use transport coordinate system in this class.
 */

// History:
//   Mar 2013, C. Gu, First public version.
//   Nov 2014, C. Gu, Rewrite it with G2PGeoPlane class.
//

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"

#include "G2PGeoPlane.hh"
#include "G2PVarDef.hh"

#include "G2PGeoSieve.hh"

using namespace std;

G2PGeoSieve *G2PGeoSieve::pG2PGeoSieve = NULL;

G2PGeoSieve::G2PGeoSieve() : fHRSAngle(0), fNRow(7), fNCol(7), fXOffset(0), fYOffset(0), fZ(799.6e-3), fNLargerHole(0), fDHole(0), fDLargerHole(0), fThreshold(0)
{
    if (pG2PGeoSieve) {
        Error("G2PGeoSieve()", "Only one instance of G2PGeoSieve allowed.");
        MakeZombie();

        return;
    }

    pG2PGeoSieve = this;

    fX.clear();
    fY.clear();
    fLargerHole.clear();
}

G2PGeoSieve::~G2PGeoSieve()
{
    if (pG2PGeoSieve == this)
        pG2PGeoSieve = NULL;
}

int G2PGeoSieve::Begin()
{
    //static const char* const here = "Begin()";

    if (G2PGeoPlane::Begin() != 0)
        return (fStatus = kINITERROR);

    fNRow = 7;
    fNCol = 7;
    fNLargerHole = 2;

    if (fHRSAngle > 0) { // left arm
        const double kSIEVEX[7] = { -3 * 13.3096e-3, -2 * 13.3096e-3, -1 * 13.3096e-3, 0.0, 1 * 13.3096e-3, 2 * 13.3096e-3, 3 * 13.3096e-3};
        const double kSIEVEY[7] = {3 * 6.1214e-3, 2 * 6.1214e-3, 1 * 6.1214e-3, 0.0, -1 * 4.7752e-3, -2 * 4.7752e-3, -3 * 4.7752e-3};
        const int kLARGERHOLE[2] = {15, 24};

        fX.clear();

        for (int i = 0; i < fNRow; i++)
            fX.push_back(kSIEVEX[i]);

        fY.clear();

        for (int i = 0; i < fNCol; i++)
            fY.push_back(kSIEVEY[i]);

        fZ = 799.60e-3;
        fXOffset = 0.0;
        fYOffset = 0.0;

        fLargerHole.clear();

        for (int i = 0; i < fNLargerHole; i++)
            fLargerHole.push_back(kLARGERHOLE[i]);
    } else { // right arm
        const double kSIEVEX[7] = { -3 * 13.3096e-3, -2 * 13.3096e-3, -1 * 13.3096e-3, 0.0, 1 * 13.3096e-3, 2 * 13.3096e-3, 3 * 13.3096e-3};
        const double kSIEVEY[7] = { -3 * 6.1214e-3, -2 * 6.1214e-3, -1 * 6.1214e-3, 0.0, 1 * 4.7752e-3, 2 * 4.7752e-3, 3 * 4.7752e-3};
        const int kLARGERHOLE[2] = {15, 24};

        fX.clear();

        for (int i = 0; i < fNRow; i++)
            fX.push_back(kSIEVEX[i]);

        fY.clear();

        for (int i = 0; i < fNCol; i++)
            fY.push_back(kSIEVEY[i]);

        fZ = 799.46e-3;
        fXOffset = 0.0;
        fYOffset = 0.0;

        fLargerHole.clear();

        for (int i = 0; i < fNLargerHole; i++)
            fLargerHole.push_back(kLARGERHOLE[i]);
    }

    fDHole = 1.3970e-3;
    fDLargerHole = 2.6924e-3;

    double ratio = fDLargerHole * fDLargerHole / (fDHole * fDHole);
    fThreshold = (ratio - 1) / ((ratio - 1) * fNLargerHole + fNRow * fNCol);

    return kOK;
}

bool G2PGeoSieve::CanPass(double *V3)
{
    int temp;
    double V5[5] = {V3[0], 0, V3[1], 0, V3[2]};

    return CanPass(V5, temp);
}

bool G2PGeoSieve::CanPass(double *V5, int &id)
{
    id = -1;

    for (int i = 0; i < fNRow * fNCol; i++) {
        double dhole = fDHole;

        for (int j = 0; j < fNLargerHole; j++)
            if (i == fLargerHole[j])
                dhole = fDLargerHole;

        int col = i / fNRow;
        int row = i % fNRow;

        double dx = V5[0] - fX[row] - fXOffset;
        double dy = V5[2] - fY[col] - fYOffset;
        double distance = sqrt(dx * dx + dy * dy);

        if (distance < dhole / 2.0) {
            id = i;
            break;
        }
    }

    if (id < 0)
        return false;

    return true;
}

bool G2PGeoSieve::AtBoundary(double *V3)
{
    if (V3[2] > fZ) return true;

    return false;
}

double G2PGeoSieve::GetZ()
{
    return fZ;
}

int G2PGeoSieve::GetNRow()
{
    return fNRow;
}

int G2PGeoSieve::GetNCol()
{
    return fNCol;
}

void G2PGeoSieve::GetPos(int index, double *V3)
{
    int col = index / fNRow;
    int row = index % fNRow;

    V3[0] = fXOffset + fX[row];
    V3[1] = fYOffset + fY[col];
    V3[2] = fZ;
}

int G2PGeoSieve::Configure(EMode mode)
{
    if (mode == kREAD || mode == kTWOWAY) {
        if (fIsInit)
            return 0;
        else
            fIsInit = true;
    }

    ConfDef confs[] = {
        {"run.hrs.angle", "HRS Angle", kDOUBLE, &fHRSAngle},
        {0}
    };

    return ConfigureFromList(confs, mode);
}

void G2PGeoSieve::MakePrefix()
{
    const char *base = "sieve";

    G2PAppBase::MakePrefix(base);
}

ClassImp(G2PGeoSieve)
