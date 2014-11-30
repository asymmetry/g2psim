// -*- C++ -*-

/* class G2PGeoPlane
 * Abstract base class for collimators like sieve or local dump.
 * It defines interface functions like CanPass().
 */

// History:
//   Nov 2014, C. Gu, Add this class for g2p geometries.
//

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"

#include "G2PGeoBase.hh"

#include "G2PGeoPlane.hh"

using namespace std;

G2PGeoPlane::G2PGeoPlane()
{
    // Nothing to do
}

G2PGeoPlane::~G2PGeoPlane()
{
    // Nothing to do
}

ClassImp(G2PGeoPlane)
