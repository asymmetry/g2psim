#ifndef G2P_SIEVE_H
#define G2P_SIEVE_H

#include <vector>

#include "TObject.h"

using namespace std;

class G2PSieve
{
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
