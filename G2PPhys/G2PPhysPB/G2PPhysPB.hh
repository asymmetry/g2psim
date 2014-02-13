// -*- C++ -*-

/* class G2PPhysPB
 * Class for P. Bosted model.
 * Unit is ub/MeV-sr.
 * Valid for all W<3 GeV and all Q2<10 GeV2.
 *  
 * Radiative correction parameters:
 * Tb: total radiative length before scattering in radiation length;
 * Ta: total radiative length after scattering in radiation length;
 */

// History:
//   Mar 2013, C. Gu, First public version.
//

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
