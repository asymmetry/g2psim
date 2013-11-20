// -*- C++ -*-

/* class G2PPhysBase
 * Abstract base class of G2PPhys classes.
 * It provides interface functions.
 */

// History:
//   Mar 2013, C. Gu, First public version.
//

#ifndef G2P_PHYSBASE_H
#define G2P_PHYSBASE_H

#include <vector>

using namespace std;

class G2PPhysBase {
public:
    G2PPhysBase();
    virtual ~G2PPhysBase();

    void SetTarget(int Z, int A);

    void SetTargetMass(double value)
    {
        fTargetMass = value;
    }

    void SetParticle(int pid)
    {
        fPID = pid;
    }

    virtual void SetPars(double* array, int n);

    virtual double GetXS(double Ei, double Ef, double theta) = 0;

protected:
    void SetTargetMass();

    int fZ, fA; // Define Target
    double fTargetMass;

    int fPID; // Define particle

    vector<double> fPars;
};

#endif
