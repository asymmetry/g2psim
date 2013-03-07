#include <cstdlib>
#include <cmath>
#include <vector>

#include "TROOT.h"
#include "TObject.h"
#include "TError.h"

#include "G2PSieve.hh"

using namespace std;

G2PSieve::G2PSieve()
{
    // Nothing to do
}

G2PSieve::~G2PSieve()
{
    // Nothing to do
}

void G2PSieve::SetSieve(double angle)
{
    // static const char* const here = "SetSieve()";

    if (angle>0) { // left arm
        const double kSIEVEX[7] = { -3*13.3096e-3, -2*13.3096e-3, -1*13.3096e-3, 0.0, 1*13.3096e-3, 2*13.3096e-3, 3*13.3096e-3 };
        const double kSIEVEY[7] = { 3*6.1214e-3, 2*6.1214e-3, 1*6.1214e-3, 0.0, -1*4.7752e-3, -2*4.7752e-3, -3*4.7752e-3 };
        const int kLARGERHOLE[2] = { 15, 24 };
        const int kSIEVEOPEN[49] = { 0, 0, 0, 0, 1, 1, 1,
                                     0, 0, 1, 1, 1, 1, 1,
                                     0, 0, 1, 1, 1, 1, 1,
                                     0, 0, 1, 1, 1, 1, 1,
                                     0, 0, 1, 1, 1, 1, 1,
                                     0, 0, 0, 1, 1, 1, 0,
                                     0, 0, 0, 0, 0, 0, 0 };

        fSieve.nRow = 7;
        fSieve.nCol = 7;
        fSieve.nLargerHole = 2;

        fSieve.fX.clear();
        for (int i = 0; i<fSieve.nRow; i++) fSieve.fX.push_back(kSIEVEX[i]);
        fSieve.fY.clear();
        for (int i = 0; i<fSieve.nCol; i++) fSieve.fX.push_back(kSIEVEY[i]);
        fSieve.fZ = 799.60e-3;
        fSieve.fXOffset = 0.0;
        fSieve.fYOffset = 0.0;

        fSieve.iLargerHole.clear();
        for (int i = 0; i<fSieve.nLargerHole; i++) fSieve.iLargerHole.push_back(kLARGERHOLE[i]);
        fSieve.bOpen.clear();
        for (int i = 0; i<fSieve.nRow*fSieve.nCol; i++) fSieve.bOpen.push_back((kSIEVEOPEN[i]==1)?true:false);

        fSieve.fDHole = 1.3970e-3;
        fSieve.fDLargerHole = 2.6924e-3;
    }
    else { // right arm
        const double kSIEVEX[7] = { -3*13.3096e-3, -2*13.3096e-3, -1*13.3096e-3, 0.0, 1*13.3096e-3, 2*13.3096e-3, 3*13.3096e-3 };
        const double kSIEVEY[7] = { -3*6.1214e-3, -2*6.1214e-3, -1*6.1214e-3, 0.0, 1*4.7752e-3, 2*4.7752e-3, 3*4.7752e-3 };
        const int kLARGERHOLE[2] = { 15, 24 };
        const int kSIEVEOPEN[49] = { 0, 0, 0, 0, 1, 1, 1,
                                     0, 0, 1, 1, 1, 1, 1,
                                     0, 0, 1, 1, 1, 1, 1,
                                     0, 0, 1, 1, 1, 1, 1,
                                     0, 0, 1, 1, 1, 1, 1,
                                     0, 0, 0, 1, 1, 1, 0,
                                     0, 0, 0, 0, 0, 0, 0 };

        fSieve.nRow = 7;
        fSieve.nCol = 7;
        fSieve.nLargerHole = 2;

        fSieve.fX.clear();
        for (int i = 0; i<fSieve.nRow; i++) fSieve.fX.push_back(kSIEVEX[i]);
        fSieve.fY.clear();
        for (int i = 0; i<fSieve.nCol; i++) fSieve.fX.push_back(kSIEVEY[i]);
        fSieve.fZ = 799.46e-3;
        fSieve.fXOffset = 0.0;
        fSieve.fYOffset = 0.0;

        fSieve.iLargerHole.clear();
        for (int i = 0; i<fSieve.nLargerHole; i++) fSieve.iLargerHole.push_back(kLARGERHOLE[i]);
        fSieve.bOpen.clear();
        for (int i = 0; i<fSieve.nRow*fSieve.nCol; i++) fSieve.bOpen.push_back((kSIEVEOPEN[i]==1)?true:false);

        fSieve.fDHole = 1.3970e-3;
        fSieve.fDLargerHole = 2.6924e-3;
    }
}

ClassImp(G2PSieve)
