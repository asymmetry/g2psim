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

#include "G2PGeoBase.hh"

class G2PGeoTube : public G2PGeoBase
{
public:
    G2PGeoTube(double rin, double rout, double length);
    G2PGeoTube(double rin, double rout, double length, double phi0, double dphi);
    virtual ~G2PGeoTube();

    // Gets

    // Sets

protected:
    G2PGeoTube(); // Only for ROOT I/O

    bool IsInside(double x, double y, double z);

    double fRin, fRout;
    double fLength;

    bool fUsePhi;
    double fPhi0, fDPhi;

private:
    ClassDef(G2PGeoTube, 1)
};

#endif
