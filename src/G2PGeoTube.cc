// -*- C++ -*-

/* class G2PGeoTube
 * Defines the geometry of tube.
 * It derives from G2PGeoBase.
 */

// History:
//   Dec 2014, C. Gu, Add this class to g2p geometries.
//

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"

#include "G2PGeoBase.hh"

#include "G2PGeoTube.hh"

static const double kPI = 3.14159265358979323846;

G2PGeoTube::G2PGeoTube()
{
    // Nothing to do
}

G2PGeoTube::G2PGeoTube(double rin, double rout, double length) : fRin(rin), fRout(rout), fLength(length), fUsePhi(false), fPhi0(0.0), fDPhi(0.0)
{
    // Nothing to do
}

G2PGeoTube::G2PGeoTube(double rin, double rout, double length, double phi0, double dphi) : fRin(rin), fRout(rout), fLength(length), fUsePhi(true), fPhi0(phi0), fDPhi(dphi)
{
    // Nothing to do
}

G2PGeoTube::~G2PGeoTube()
{
    // Nothing to do
}

bool G2PGeoTube::IsInside(double x, double y, double z)
{
    double r = sqrt(x * x + y * y);

    if (r < fRin || r > fRout || fabs(z) > (fLength / 2.0))
        return false;

    if (fUsePhi) {
        double phi = atan2(y, x);
        phi = (phi < 0) ? (phi + 2 * kPI) : phi;

        if (fPhi0 >= 0) {
            if ((phi < fPhi0) || (phi > fPhi0 + fDPhi))
                return false;
        } else {
            if ((phi < fPhi0 + 2 * kPI) && (phi > fPhi0 + fDPhi))
                return false;
        }
    }

    return true;
}

ClassImp(G2PGeoTube)
