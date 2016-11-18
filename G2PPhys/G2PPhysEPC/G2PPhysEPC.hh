// -*- C++ -*-

/* class G2PPhysEPC
 * Class for EPC model.
 * Unit is ub/MeV-sr.
 * The cross section is calculated per nuclei. (Notice the difference with WISER model)
 * This model considers 2 situations: single pion production and multiple pion production.
 *
 * Parameters:
 * [1] fMPI: 0 means to calculate single pion production only, default is multiple pion production.
 */

// History:
//   Feb 2014, C. Gu, First public version.
//   Apr 2014, C. Gu, Updated with multiple pion production
//

#ifndef G2P_PHYSEPC_H
#define G2P_PHYSEPC_H

#include "G2PPhysBase.hh"

class G2PPhysEPC : public G2PPhysBase
{
public:
    G2PPhysEPC();
    ~G2PPhysEPC();

    void SetPar(int id, double value);
    double GetXS(double Ei, double Pf, double theta);

private:
    int fMPI;
};

#endif
