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

#include "G2PAppBase.hh"
#include "G2PGeoPlane.hh"
#include "G2PVarDef.hh"

#include "G2PGeoSieve.hh"

using namespace std;

G2PGeoSieve *G2PGeoSieve::pG2PGeoSieve = NULL;

G2PGeoSieve::G2PGeoSieve() : fHRSAngle(0), fNRow(7), fNCol(7), fNLargerHole(0), fDHole(0), fDLargerHole(0), fThreshold(0)
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

    fUseTrans = true;
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

    fX.clear();
    fY.clear();
    fLargerHole.clear();

    if (fHRSAngle > 0) { // left arm
        SetOrigin(0.0, 0.0, 799.60e-3);

        const double kSIEVEX[7] = { -3 * 13.3096e-3, -2 * 13.3096e-3, -1 * 13.3096e-3, 0.0, 1 * 13.3096e-3, 2 * 13.3096e-3, 3 * 13.3096e-3};
        const double kSIEVEY[7] = {3 * 6.1214e-3, 2 * 6.1214e-3, 1 * 6.1214e-3, 0.0, -1 * 4.7752e-3, -2 * 4.7752e-3, -3 * 4.7752e-3};
        const int kLARGERHOLE[2] = {15, 24};

        for (int i = 0; i < fNRow; i++)
            fX.push_back(kSIEVEX[i]);

        for (int i = 0; i < fNCol; i++)
            fY.push_back(kSIEVEY[i]);

        for (int i = 0; i < fNLargerHole; i++)
            fLargerHole.push_back(kLARGERHOLE[i]);
    } else { // right arm
        SetOrigin(0.0, 0.0, 799.46e-3);

        const double kSIEVEX[7] = { -3 * 13.3096e-3, -2 * 13.3096e-3, -1 * 13.3096e-3, 0.0, 1 * 13.3096e-3, 2 * 13.3096e-3, 3 * 13.3096e-3};
        const double kSIEVEY[7] = { -3 * 6.1214e-3, -2 * 6.1214e-3, -1 * 6.1214e-3, 0.0, 1 * 4.7752e-3, 2 * 4.7752e-3, 3 * 4.7752e-3};
        const int kLARGERHOLE[2] = {15, 24};

        for (int i = 0; i < fNRow; i++)
            fX.push_back(kSIEVEX[i]);

        for (int i = 0; i < fNCol; i++)
            fY.push_back(kSIEVEY[i]);

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

        double dx = V5[0] - fX[row] - fOrigin[0];
        double dy = V5[2] - fY[col] - fOrigin[1];
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

int G2PGeoSieve::GetNRow()
{
    return fNRow;
}

int G2PGeoSieve::GetNCol()
{
    return fNCol;
}

double G2PGeoSieve::GetZ()
{
    return fOrigin[2];
}

void G2PGeoSieve::GetPos(int index, double *V3)
{
    int col = index / fNRow;
    int row = index % fNRow;

    V3[0] = fOrigin[0] + fX[row];
    V3[1] = fOrigin[1] + fY[col];
    V3[2] = fOrigin[2];
}

void G2PGeoSieve::MakePrefix()
{
    const char *base = "sieve";

    G2PAppBase::MakePrefix(base);
}

ClassImp(G2PGeoSieve)
