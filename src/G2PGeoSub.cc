// -*- C++ -*-

/* class G2PGeoSub
 * Defines the subtraction of geometries.
 * It derives from G2PGeoBase.
 */

// History:
//   Jan 2015, C. Gu, Add this class to g2p geometries.
//

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"

#include "G2PAppList.hh"
#include "G2PGeoBase.hh"

#include "G2PGeoSub.hh"

G2PGeoSub::G2PGeoSub()
{
    // Nothing to do
}

G2PGeoSub::G2PGeoSub(G2PGeoBase *geo) : fMinuend(geo)
{
    fSubGeos = new G2PAppList();
}

G2PGeoSub::~G2PGeoSub()
{
    TIter next(fSubGeos);

    while (G2PAppBase *aobj = static_cast<G2PAppBase *>(next()))
        fSubGeos->Remove(aobj);

    delete fSubGeos;
}

bool G2PGeoSub::IsInside(const double *V3)
{
    bool result = false;

    if (fMinuend->IsInside(V3)) result = true;

    TIter next(fSubGeos);

    while (G2PGeoBase *geo = static_cast<G2PGeoBase *>(next())) {
        if (!geo->IsInside(V3)) result = true;
    }

    return result;
}

void G2PGeoSub::Substract(G2PGeoBase *geo)
{
    fSubGeos->Add(geo);
}

bool G2PGeoSub::IsInside(double x, double y, double z)
{
    // Nothing to do

    return true;
}

ClassImp(G2PGeoSub)
