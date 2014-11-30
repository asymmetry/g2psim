// -*- C++ -*-

/* class G2PGeoPlane
 * Abstract base class for collimators like sieve or local dump.
 * It defines interface functions like CanPass().
 */

// History:
//   Nov 2014, C. Gu, Add this class for g2p geometries.
//

#ifndef G2P_GEOPLANE_H
#define G2P_GEOPLANE_H

#include "G2PGeoBase.hh"

class G2PGeoPlane : public G2PGeoBase
{
public:
    virtual ~G2PGeoPlane();

    virtual bool CanPass(double *V3) = 0;

    // Gets

    // Sets

protected:
    G2PGeoPlane(); // No instance allowed for this class

private:
    ClassDef(G2PGeoPlane, 1)
};

#endif
