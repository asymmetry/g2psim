// -*- C++ -*-

/* class G2PPhysEPC
 * Class for EPC model.
 * Unit is ub/MeV-sr.
 * predict (e,N) cross sections to within a factor of 2 for an incident electron in the energy range 0.5-5 GeV and for nucleon kinetic energies greater than 50 MeV.
 */

// History:
//   Feb 2014, C. Gu, First public version.
//

#ifndef G2P_PHYSEPC_H
#define G2P_PHYSEPC_H

#include "G2PPhysBase.hh"

class G2PPhysEPC : public G2PPhysBase {
public:
    G2PPhysEPC();
    ~G2PPhysEPC();

    void SetPars(double* array, int n);

    double GetXS(double Ei, double Pf, double theta);
};

#endif
