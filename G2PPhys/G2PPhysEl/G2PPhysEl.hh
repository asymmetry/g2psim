// -*- C++ -*-

/* class G2PPhysEl
 * Elastic cross section models.
 * He4: calculated from charge and magnetization densities, original coded by M. Friedman.
 * C12: 1) calculated from fitted charge densities.
 *      2) calculated from form factor.
 * N14: same as He4.
 */

// History:
//   Mar 2013, C. Gu, First public version, with C12 models.
//   Apr 2013, C. Gu, Add L. Cardman's C12 charge densities.
//   Nov 2013, C. Gu, Add He charge and magnetization densities from D. Jager, original coded by M. Friedman.
//

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
    double GetXS_N14(double Ei, double theta);
    double GetXS_All(double Ei, double theta);
};

#endif
