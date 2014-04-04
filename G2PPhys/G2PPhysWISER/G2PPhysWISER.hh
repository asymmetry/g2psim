// -*- C++ -*-

/* class G2PPhysWISER
 * Class for WISER model.
 * Unit is ub/MeV-sr.
 * The cross-section is calculated per nucleon. (Notice the difference with EPC model)
 * Photoproduction of pion/nucleons in DIS region.
 * 
 * Parameters:
 * fRadLen: radiation length of the target, include both external and internal contribution.
 * 
 * How to set parameters:
 * If set 1 parameters with SetPars(pars,1), then pars[0]->fRadLen;
 * Other uses will be considered as invalid.
 */

// History:
//   Apr 2013, C. Gu, First public version.
//

#ifndef G2P_PHYSWISER_H
#define G2P_PHYSWISER_H

#include "G2PPhysBase.hh"

class G2PPhysWISER : public G2PPhysBase {
public:
    G2PPhysWISER();
    ~G2PPhysWISER();

    void SetPars(double* array, int n);

    double GetXS(double Ei, double Pf, double theta);

private:
    double fRadLen;
};

#endif
