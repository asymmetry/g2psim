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

#include "G2PAppBase.hh"

using namespace std;

class G2PSieve : public G2PAppBase {
public:
    G2PSieve();
    virtual ~G2PSieve();

    virtual int Begin();

    virtual int GetPos(double* V3);
    virtual void GetPos(int index, double* V3);
    virtual int CanPass(double* V5);

    // Gets
    double GetZ();
    int GetNRow();
    int GetNCol();

    // Sets

protected:
    virtual int Configure(EMode mode = kTWOWAY);
    virtual void MakePrefix();

    double fHRSAngle;

    int fNRow;
    int fNCol;
    vector<double> fX;
    vector<double> fY;
    double fZ;
    double fXOffset;
    double fYOffset;
    int fNLargerHole;
    vector<int> fLargerHole;
    vector<bool> fIsOpen;
    double fDHole;
    double fDLargerHole;
    double fThreshold;

private:
    static G2PSieve* pG2PSieve;

    ClassDef(G2PSieve, 1)
};

#endif
