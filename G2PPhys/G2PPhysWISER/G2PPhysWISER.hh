// -*- C++ -*-

/* class G2PPhysWISER
 * Class for WISER model.
 * Unit is ub/MeV-sr.
 * The cross-section is calculated per nucleon. (Notice the difference with EPC model)
 * Photoproduction of pion/nucleons in DIS region.
 *
 * Parameters:
 * [1] fRadLen: radiation length of the target, include both external and internal contribution.
 */

// History:
//   Apr 2013, C. Gu, First public version.
//

#ifndef G2P_PHYSWISER_H
#define G2P_PHYSWISER_H

#include "G2PPhysBase.hh"

class G2PPhysWISER : public G2PPhysBase
{
public:
    G2PPhysWISER();
    ~G2PPhysWISER();

    void SetPar(int id, double value);
    double GetXS(double Ei, double Pf, double theta);

private:
    double fRadLen;
};

#endif
