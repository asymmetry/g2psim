// -*- C++ -*-

/* class G2PGeoTube
 * Defines the geometry of tube.
 * It derives from G2PGeoBase.
 */

// History:
//   Dec 2014, C. Gu, Add this class to g2p geometries.
//

#ifndef G2P_GEOTUBE_H
#define G2P_GEOTUBE_H

#include "G2PGeoSolid.hh"

class G2PGeoTube : public G2PGeoSolid
{
public:
    G2PGeoTube(double rin, double rout, double length);
    virtual ~G2PGeoTube();

    bool IsInside(double *V3);

    // Gets

    // Sets

protected:
    G2PGeoTube(); // Only for ROOT I/O

    double fRin, fRout;
    double fLength;

private:
    ClassDef(G2PGeoTube, 1)
};

#endif
