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

G2PGeoTube::G2PGeoTube()
{
    // Nothing to do
}

G2PGeoTube::G2PGeoTube(double rin, double rout, double length) : fRin(rin), fRout(rout), fLength(length)
{
    // Nothing to do
}

G2PGeoTube::~G2PGeoTube()
{
    // Nothing to do
}

bool G2PGeoTube::IsInside(double *V3)
{
    double r = sqrt(V3[0] * V3[0] + V3[1] * V3[1]);

    if (r < fRin || r > fRout || fabs(V3[2]) > (fLength / 2.0))
        return false;

    return true;
}

ClassImp(G2PGeoTube)
