#ifndef G2P_PHYSEL_H
#define G2P_PHYSEL_H

#include "G2PPhysBase.hh"

class G2PPhysEl : public G2PPhysBase {
public:
    G2PPhysEl();
    ~G2PPhysEl();

    void SetPars(double* array, int n);

    double GetXS(double Ei, double Ef, double theta);

private:
    int iSetting;

    double GetXS_He4(double Ei, double theta);
    double GetXS_C12(double Ei, double theta);
    double GetXS_All(double Ei, double theta);
};

#endif
