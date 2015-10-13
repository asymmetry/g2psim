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
    // Only for ROOT I/O
}

G2PGeoSub::G2PGeoSub(G2PGeoBase *geo) : fMinuend(geo)
{
    fSubGeos = new G2PAppList();
}

G2PGeoSub::~G2PGeoSub()
{
    TIter next(fSubGeos);

    while (G2PGeoBase *aobj = static_cast<G2PGeoBase *>(next())) {
        if (aobj != NULL) {
            fSubGeos->Remove(aobj);

            if (aobj->GetNUsed() == 0)
                delete aobj;
            else
                aobj->SetNUsed(aobj->GetNUsed() - 1);
        }
    }

    if (fMinuend) {
        delete fMinuend;
        fMinuend = NULL;
    }

    if (fSubGeos) {
        delete fSubGeos;
        fSubGeos = NULL;
    }
}

int G2PGeoSub::Begin()
{
    if (G2PAppBase::Begin() != 0)
        return (fStatus = kBEGINERROR);

    fMinuend->Begin();

    TIter geo_iter(fSubGeos);

    while (G2PGeoBase *geo = static_cast<G2PGeoBase *>(geo_iter())) {
        if (!geo->IsInit())
            if (geo->Begin() != 0)
                return (fStatus = kBEGINERROR);
    }

    return (fStatus = kOK);
}

int G2PGeoSub::End()
{
    fMinuend->End();

    TIter geo_iter(fSubGeos);

    while (G2PGeoBase *geo = static_cast<G2PGeoBase *>(geo_iter()))
        if (geo->End() != 0)
            return -1;

    return (G2PGeoBase::End());
}

bool G2PGeoSub::IsInside(const double *V3)
{
    bool result = false;

    if (fMinuend->IsInside(V3))
        result = true;

    TIter next(fSubGeos);

    while (G2PGeoBase *geo = static_cast<G2PGeoBase *>(next())) {
        if (geo->IsInside(V3))
            result = false;
    }

    return result;
}

void G2PGeoSub::Subtract(G2PGeoBase *geo)
{
    geo->SetNUsed(geo->GetNUsed() + 1);
    fSubGeos->Add(geo);
}

void G2PGeoSub::SetOrigin(double x, double y, double z)
{
    fMinuend->SetOrigin(x, y, z);
}

void G2PGeoSub::SetEulerAngle(double alpha, double beta, double gamma)
{
    fMinuend->SetEulerAngle(alpha, beta, gamma);
}

bool G2PGeoSub::IsInside(double x, double y, double z)
{
    // Nothing to do

    return true;
}

ClassImp(G2PGeoSub)
