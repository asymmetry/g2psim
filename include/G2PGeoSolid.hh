// -*- C++ -*-

/* class G2PGeoSolid
 * Abstract base class for g2p geometry solids.
 * It defines interface functions like IsInside().
 * AtBoundary() is redirected to IsInside() for solids.
 */

// History:
//   Nov 2014, C. Gu, Add this class for g2p geometry solids.
//

#ifndef G2P_GEOSOLID_H
#define G2P_GEOSOLID_H

#include "G2PGeoBase.hh"

class G2PGeoSolid : public G2PGeoBase
{
public:
    virtual ~G2PGeoSolid();

    bool AtBoundary(double *V3);

    // Gets

    // Sets

protected:
    G2PGeoSolid(); // No instance allowed for this class

    virtual bool IsInside(double *V3) = 0;

private:

    ClassDef(G2PGeoSolid, 1)
};

#endif
