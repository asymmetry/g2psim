// -*- C++ -*-

/* class G2PPhysEl
 * Class to calculate elastic cross section.
 * Unit is ub/sr.
 *
 * Elastic cross section models.
 * * All: Form factors from K. C. Stansfield et al., Phys. Rev. C, 3(1971)1448
 * * H1 : Form factors from S. Venkat et al., Phys. Rev. C, 83(2011)015203 (global fit, with TPE correction)
 *                          J. Arrington et al., Phys. Rev. C 76(2007)035201 (low Q2, with/without TPE correction)
 * * He4: Charge and magnetization densities from De Jager, At. Data Nucl. Data Tables, 14(1974)
 * * C12: Charge and magnetization densities from De Jager, At. Data Nucl. Data Tables, 14(1974)
 *        Charge distribution from L. S. Cardman et al., Phys. Lett. B, 91(1970)203
 * * N14: Charge and magnetization densities from De Jager, At. Data Nucl. Data Tables, 14(1974)
 *
 * Parameters:
 * [1] fRes: resolution of the elastic peak;
 * [2] fSettings:
 *   H1:
 *     2 means to use 2007 low Q2 fit without TPE correction,
 *     3 means to use 2007 low Q2 fit with TPE correction,
 *     default is to use 2011 global fit with TPE correction;
 *   C12:
 *     2 means to use Larry's fit,
 *     default is to use charge densities from De Jager.
 */

// History:
//   Mar 2013, C. Gu, First public version, with C12 models.
//   Apr 2013, C. Gu, Add L. Cardman's C12 charge densities.
//   Nov 2013, C. Gu, Add He4 and N14 charge and magnetization densities from D. Jager, original coded by M. Friedman.
//   Apr 2014, C. Gu, Add H1 form factors from J. Arrington. (Thanks to M. Cummings)
//   Nov 2016, C. Gu, Rewrite the elastic model with a Gaussian type resolution.
//   Apr 2017, C. Gu, Add C12 charge and magnetization densities from D. Jager.
//

#ifndef G2P_PHYSEL_H
#define G2P_PHYSEL_H

#include "G2PPhysBase.hh"

class G2PPhysEl : public G2PPhysBase
{
public:
    G2PPhysEl();
    ~G2PPhysEl();

    void SetPar(int id, double value);
    double GetXS(double Ei, double Ef, double theta);

private:
    double fRes;
    int fSetting;

    double GetXS_H1(double Ei, double theta);
    double GetXS_He4(double Ei, double theta);
    double GetXS_C12(double Ei, double theta);
    double GetXS_N14(double Ei, double theta);
    double GetXS_All(double Ei, double theta);
};

#endif
