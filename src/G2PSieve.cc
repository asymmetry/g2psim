// -*- C++ -*-

/* class G2PSieve
 * This file defines a class G2PSieve.
 * It defines geometry of the sieve slit.
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

#include "G2PSieve.hh"

using namespace std;

G2PSieve *G2PSieve::pG2PSieve = NULL;

G2PSieve::G2PSieve() : fNRow(7), fNCol(7), fNLargerHole(0), fDHole(0), fDLargerHole(0), fThreshold(0)
{
    if (pG2PSieve) {
        Error("G2PSieve()", "Only one instance of G2PSieve allowed.");
        MakeZombie();

        return;
    }

    pG2PSieve = this;

    fX.clear();
    fY.clear();
    fLargerHole.clear();
}

G2PSieve::~G2PSieve()
{
    if (pG2PSieve == this)
        pG2PSieve = NULL;
}

int G2PSieve::Begin()
{
    //static const char* const here = "Begin()";

    if (G2PAppBase::Begin() != 0)
        return (fStatus = kBEGINERROR);

    fNRow = 7;
    fNCol = 7;
    fNLargerHole = 2;

    fX.clear();
    fY.clear();
    fLargerHole.clear();

    if (fHRSAngle > 0) { // left arm
        fZ = 799.60e-3;

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
        fZ = 799.46e-3;

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

    return (fStatus = kOK);
}

bool G2PSieve::CanPass(double *V5, int &id)
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

void G2PSieve::GetPos(int index, double *V3)
{
    int col = index / fNRow;
    int row = index % fNRow;

    V3[0] = fXOffset + fX[row];
    V3[1] = fYOffset + fY[col];
    V3[2] = fZ;
}

int G2PSieve::GetNRow()
{
    return fNRow;
}

int G2PSieve::GetNCol()
{
    return fNCol;
}

double G2PSieve::GetZ()
{
    return fZ;
}

void G2PSieve::MakePrefix()
{
    const char *base = "sieve";

    G2PAppBase::MakePrefix(base);
}

ClassImp(G2PSieve)
