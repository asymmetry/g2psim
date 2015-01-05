// -*- C++ -*-

/* class G2PGeoSolid
 * Abstract base class for g2p geometry solids.
 * It defines interface functions like IsInside().
 * TouchBoundaryGeo() is redirected to IsInside() for solids.
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

    virtual bool IsInside(double *V3) = 0;

    // Gets

    // Sets

protected:
    G2PGeoSolid(); // No instance allowed for this class

    bool TouchBoundaryGeo(double x, double y, double z);

private:
    ClassDef(G2PGeoSolid, 1)
};

#endif
