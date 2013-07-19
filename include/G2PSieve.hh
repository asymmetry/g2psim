// -*- C++ -*-

/* class G2PSieve
 * This file defines a class G2PSieve.
 * It defines geometry of the sieve slit.
 * G2PGun classes will inherit this class to get sieve slit geometry.
 */

// History:
//   Mar 2013, C. Gu, First public version.
//

#ifndef G2P_SIEVE_H
#define G2P_SIEVE_H

#include <vector>

#include "TObject.h"

using namespace std;

class G2PSieve {
public:
    virtual ~G2PSieve();

protected:
    G2PSieve(); // No instance allowed for this class

    virtual void SetSieve(double angle);

    typedef struct {
        int nRow;
        int nCol;
        vector<double> fX;
        vector<double> fY;
        double fZ;
        double fXOffset;
        double fYOffset;
        int nLargerHole;
        vector<int> iLargerHole;
        vector<bool> bOpen;
        double fDHole;
        double fDLargerHole;
    } sSieve;

    sSieve fSieve;

    ClassDef(G2PSieve, 1)
};

#endif
