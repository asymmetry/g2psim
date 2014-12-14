// -*- C++ -*-

/* class G2PGeoSolid
 * Abstract base class for g2p geometry solids.
 * It defines interface functions like IsInside().
 * TouchBoundaryGeo() is redirected to IsInside() for solids.
 */

// History:
//   Nov 2014, C. Gu, Add this class for g2p geometry solids.
//

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"

#include "G2PGeoBase.hh"

#include "G2PGeoSolid.hh"

using namespace std;

G2PGeoSolid::G2PGeoSolid()
{
    // Nothing to do
}

G2PGeoSolid::~G2PGeoSolid()
{
    // Nothing to do
}

bool G2PGeoSolid::TouchBoundaryGeo(double x, double y, double z)
{
    return (!IsInside(x, y, z));
}

ClassImp(G2PGeoSolid)
