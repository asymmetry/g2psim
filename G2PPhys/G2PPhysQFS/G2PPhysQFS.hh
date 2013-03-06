#ifndef G2P_PHYSQFS_H
#define G2P_PHYSQFS_H

#include "G2PPhysBase.hh"

class G2PPhysQFS : public G2PPhysBase
{
public:
    G2PPhysQFS();
    ~G2PPhysQFS();

    void SetPars(double* array, int n);

    double GetXS(double Ei, double Ef, double theta);

private:
    double fEPS, fEPSD, fFP, fTb, fTa;
}; 

#endif
