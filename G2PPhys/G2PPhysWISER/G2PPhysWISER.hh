#ifndef G2P_PHYSWISER_H
#define G2P_PHYSWISER_H

#include "G2PPhysBase.hh"

class G2PPhysWISER : public G2PPhysBase {
public:
    G2PPhysWISER();
    ~G2PPhysWISER();

    void SetPars(double* array, int n);

    double GetXS(double Ei, double Pf, double theta);

private:
    double fRadLen;
};

#endif
