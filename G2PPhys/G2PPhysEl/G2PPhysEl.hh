// -*- C++ -*-

/* class G2PPhysEl
 * Class to calculate elastic cross section.
 * Unit is ub/sr.
 * 
 * Elastic cross section models.
 * * All: Form factors from K. C. Stansfield et al., Phys. Rev. C, 3(1971)1448
 * * 1H : Form factors from J. Arrington, Phys. Rev. C, 69(2004)022201
 * * 4He: Charge and magnetization densities from De Jager, At. Data Nucl. Data Tables, 14(1974)
 * * 12C: Charge distribution from L. S. Cardman et al., Phys. Lett. B, 91(1970)203 
 * * 14N: Charge and magnetization densities from De Jager, At. Data Nucl. Data Tables, 14(1974)
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

    double GetXS_H1(double Ei, double theta);
    double GetXS_He4(double Ei, double theta);
    double GetXS_N14(double Ei, double theta);
    double GetXS_All(double Ei, double theta);
};

#endif
