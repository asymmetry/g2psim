// -*- C++ -*-

/* class G2PPhysQFS
 * Class for QFS model.
 * Unit is ub/MeV-sr.
 * Predict (e,e') cross sections to within 20% for an incident electron in the energy range 0.5-5 GeV and for energy losses greater than 50 MeV.
 * 
 * Radiative correction parameters:
 * Tb: total radiative length before scattering in radiation length;
 * Ta: total radiative length after scattering in radiation length;
 *
 * QFS model parameters:
 * EPS: separation energy in MeV;
 * EPSD: delta separation energy in MeV;
 * FP: Fermi momentum in MeV/c;
 */

// History:
//   Mar 2013, C. Gu, First public version.
//

#ifndef G2P_PHYSQFS_H
#define G2P_PHYSQFS_H

#include "G2PPhysBase.hh"

class G2PPhysQFS : public G2PPhysBase {
public:
    G2PPhysQFS();
    ~G2PPhysQFS();

    void SetPars(double* array, int n);

    double GetXS(double Ei, double Ef, double theta);

private:
    double fEPS, fEPSD, fFP;
    double fTb, fTa;
};

#endif
