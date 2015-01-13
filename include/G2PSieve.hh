// -*- C++ -*-

/* class G2PSieve
 * This file defines a class G2PSieve.
 * It defines geometry of the sieve slit.
 */

// History:
//   Mar 2013, C. Gu, First public version.
//

#ifndef G2P_SIEVE_H
#define G2P_SIEVE_H

#include <vector>

#include "G2PAppBase.hh"

using namespace std;

class G2PSieve : public G2PAppBase
{
public:
    G2PSieve();
    virtual ~G2PSieve();

    virtual int Begin();

    virtual bool CanPass(double *V5, int &id);
    void GetPos(int index, double *V3);

    // Gets
    int GetNRow();
    int GetNCol();
    double GetZ();

    // Sets

protected:
    virtual void MakePrefix();

    int fNRow;
    int fNCol;

    double fXOffset;
    vector<double> fX;
    double fYOffset;
    vector<double> fY;
    double fZ;

    int fNLargerHole;
    vector<int> fLargerHole; // index of larger sieve holes

    double fDHole;
    double fDLargerHole;
    double fThreshold;

private:
    static G2PSieve *pG2PSieve;

    ClassDef(G2PSieve, 1)
};

#endif
