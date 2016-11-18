// -*- C++ -*-

/* class G2PPhysQFS
 * Class for QFS model.
 * Unit is ub/MeV-sr.
 * Predict (e,e') cross sections to within 20% for an incident electron in the energy range 0.5-5 GeV and for energy losses greater than 50 MeV.
 *
 * Radiative correction parameters:
 * [1] Tb: total radiative length before scattering in radiation length;
 * [2] Ta: total radiative length after scattering in radiation length.
 *
 * QFS model parameters:
 * [3] EPS: separation energy in MeV;
 * [4] EPSD: delta separation energy in MeV;
 * [5] FP: Fermi momentum in MeV/c.
 */

// History:
//   Mar 2013, C. Gu, First public version.
//

#ifndef G2P_PHYSQFS_H
#define G2P_PHYSQFS_H

#include "G2PPhysBase.hh"

class G2PPhysQFS : public G2PPhysBase
{
public:
    G2PPhysQFS();
    ~G2PPhysQFS();

    void SetPar(int id, double value);
    double GetXS(double Ei, double Ef, double theta);

private:
    double fEPS, fEPSD, fFP;
    double fTb, fTa;
};

#endif
