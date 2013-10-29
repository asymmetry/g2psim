#ifndef G2P_PHYSPB_H
#define G2P_PHYSPB_H

#include "G2PPhysBase.hh"

class G2PPhysPB : public G2PPhysBase {
public:
    G2PPhysPB();
    ~G2PPhysPB();

    void SetPars(double* array, int n);

    double GetXS(double Ei, double Ef, double theta);

private:
    double fTb, fTa;
};

#endif
