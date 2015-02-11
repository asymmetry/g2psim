// -*- C++ -*-

/* class G2PGeoBase
 * Abstract base class for g2p geometries.
 * It provides fundamental functions like rotation.
 * No instance allowed for this class.
 * The rotation matrix is defined with Euler angle in Z-X'-Z" convention.
 * G2PDrift class uses IsInside() to determine whether stop drifting or not.
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

#include "G2PAppBase.hh"

#include "G2PGeoBase.hh"

G2PGeoBase::G2PGeoBase() : fTranslation(false), fRotation(false), pMaterial(NULL), pfLab2Geo(NULL), pfGeo2Lab(NULL)
{
    memset(fOrigin, 0, sizeof(fOrigin));
    memset(fEulerAngle, 0, sizeof(fEulerAngle));
    memset(fRotationMatrix, 0, sizeof(fRotationMatrix));
}

G2PGeoBase::~G2PGeoBase()
{
    // Nothing to do
}

int G2PGeoBase::Begin()
{
    if (G2PAppBase::Begin() != 0)
        return (fStatus = kBEGINERROR);

    SetGeoPosition();

    if (fTranslation) {
        if (fRotation) {
            pfLab2Geo = &G2PGeoBase::Lab2Geo;
            pfGeo2Lab = &G2PGeoBase::Geo2Lab;
        } else {
            pfLab2Geo = &G2PGeoBase::Lab2GeoNoR;
            pfGeo2Lab = &G2PGeoBase::Geo2LabNoR;
        }
    } else {
        if (fRotation) {
            pfLab2Geo = &G2PGeoBase::Lab2GeoNoT;
            pfGeo2Lab = &G2PGeoBase::Geo2LabNoT;
        } else {
            pfLab2Geo = &G2PGeoBase::Lab2GeoNoTNoR;
            pfGeo2Lab = &G2PGeoBase::Geo2LabNoTNoR;
        }
    }

    return (fStatus = kOK);
}

bool G2PGeoBase::IsInside(const double *V5_tr, double z_tr)
{
    double V3_lab[3];

    TCS2HCS(V5_tr[0], V5_tr[2], z_tr, V3_lab[0], V3_lab[1], V3_lab[2]);

    return IsInside(V3_lab);
}

bool G2PGeoBase::IsInside(const double *V3)
{
    double V3_geo[3];

    (this->*pfLab2Geo)(V3, V3_geo);

    return IsInside(V3_geo[0], V3_geo[1], V3_geo[2]);
}

G2PMaterial *G2PGeoBase::GetMaterial() const
{
    return pMaterial;
}

void G2PGeoBase::SetOrigin(double x, double y, double z)
{
    fOrigin[0] = x;
    fOrigin[1] = y;
    fOrigin[2] = z;

    fConfigIsSet.insert((unsigned long) &fOrigin[0]);
    fConfigIsSet.insert((unsigned long) &fOrigin[1]);
    fConfigIsSet.insert((unsigned long) &fOrigin[2]);
}

void G2PGeoBase::SetEulerAngle(double alpha, double beta, double gamma)
{
    // The Euler angle is defined using Z-X'-Z" convention

    fEulerAngle[0] = alpha;
    fEulerAngle[1] = beta;
    fEulerAngle[2] = gamma;

    fConfigIsSet.insert((unsigned long) &fEulerAngle[0]);
    fConfigIsSet.insert((unsigned long) &fEulerAngle[1]);
    fConfigIsSet.insert((unsigned long) &fEulerAngle[2]);
}

void G2PGeoBase::SetMaterial(G2PMaterial *mat)
{
    pMaterial = mat;
}

void G2PGeoBase::SetGeoPosition()
{
    // Set translation and rotation
    if ((fabs(fOrigin[0]) < 1e-5) && (fabs(fOrigin[1]) < 1e-5) && (fabs(fOrigin[2]) < 1e-5))
        fTranslation = false;
    else
        fTranslation = true;

    // The Euler angle is defined using Z-X'-Z" convention
    if ((fabs(fEulerAngle[0]) < 1e-5) && (fabs(fEulerAngle[1]) < 1e-5) && (fabs(fEulerAngle[2]) < 1e-5))
        fRotation = false;
    else {
        double s1 = sin(fEulerAngle[0]);
        double c1 = cos(fEulerAngle[0]);
        double s2 = sin(fEulerAngle[1]);
        double c2 = cos(fEulerAngle[1]);
        double s3 = sin(fEulerAngle[2]);
        double c3 = cos(fEulerAngle[2]);

        fRotationMatrix[0][0][0] = c1 * c3 - c2 * s1 * s3;
        fRotationMatrix[0][1][0] = -c1 * s3 - c2 * s1 * c3;
        fRotationMatrix[0][2][0] = s2 * s1;
        fRotationMatrix[0][0][1] = s1 * c3 + c2 * c1 * s3;
        fRotationMatrix[0][1][1] = -s1 * s3 + c2 * c1 * c3;
        fRotationMatrix[0][2][1] = -s2 * c1;
        fRotationMatrix[0][0][2] = s2 * s3;
        fRotationMatrix[0][1][2] = s2 * c3;
        fRotationMatrix[0][2][2] = c2;

        // inverse matrix
        // this is also the matrix to rotate the geometry to its direction
        fRotationMatrix[1][0][0] = c1 * c3 - c2 * s1 * s3;
        fRotationMatrix[1][0][1] = -c1 * s3 - c2 * c3 * s1;
        fRotationMatrix[1][0][2] = s1 * s2;
        fRotationMatrix[1][1][0] = c3 * s1 + c1 * c2 * s3;
        fRotationMatrix[1][1][1] = c1 * c2 * c3 - s1 * s3;
        fRotationMatrix[1][1][2] = -c1 * s2;
        fRotationMatrix[1][2][0] = s2 * s3;
        fRotationMatrix[1][2][1] = c3 * s2;
        fRotationMatrix[1][2][2] = c2;

        fRotation = true;
    }
}

void G2PGeoBase::Lab2Geo(const double *V3_lab, double *V3_geo)
{
    double temp[3] = {0};

    Lab2GeoNoR(V3_lab, temp);
    Lab2GeoNoT(temp, V3_geo);
}

void G2PGeoBase::Lab2GeoNoT(const double *V3_lab, double *V3_geo)
{
    V3_geo[0] = fRotationMatrix[0][0][0] * V3_lab[0] + fRotationMatrix[0][0][1] * V3_lab[1] + fRotationMatrix[0][0][2] * V3_lab[2];
    V3_geo[1] = fRotationMatrix[0][1][0] * V3_lab[0] + fRotationMatrix[0][1][1] * V3_lab[1] + fRotationMatrix[0][1][2] * V3_lab[2];
    V3_geo[2] = fRotationMatrix[0][2][0] * V3_lab[0] + fRotationMatrix[0][2][1] * V3_lab[1] + fRotationMatrix[0][2][2] * V3_lab[2];
}

void G2PGeoBase::Lab2GeoNoR(const double *V3_lab, double *V3_geo)
{
    V3_geo[0] = V3_lab[0] - fOrigin[0];
    V3_geo[1] = V3_lab[1] - fOrigin[1];
    V3_geo[2] = V3_lab[2] - fOrigin[2];
}

void G2PGeoBase::Lab2GeoNoTNoR(const double *V3_lab, double *V3_geo)
{
    V3_geo[0] = V3_lab[0];
    V3_geo[1] = V3_lab[1];
    V3_geo[2] = V3_lab[2];
}

void G2PGeoBase::Geo2Lab(const double *V3_geo, double *V3_lab)
{
    double temp[3] = {0};

    Geo2LabNoT(V3_geo, temp);
    Geo2LabNoR(temp, V3_lab);
}

void G2PGeoBase::Geo2LabNoT(const double *V3_geo, double *V3_lab)
{
    V3_lab[0] = fRotationMatrix[1][0][0] * V3_geo[0] + fRotationMatrix[1][0][1] * V3_geo[1] + fRotationMatrix[1][0][2] * V3_geo[2];
    V3_lab[1] = fRotationMatrix[1][1][0] * V3_geo[0] + fRotationMatrix[1][1][1] * V3_geo[1] + fRotationMatrix[1][1][2] * V3_geo[2];
    V3_lab[2] = fRotationMatrix[1][2][0] * V3_geo[0] + fRotationMatrix[1][2][1] * V3_geo[1] + fRotationMatrix[1][2][2] * V3_geo[2];
}

void G2PGeoBase::Geo2LabNoR(const double *V3_geo, double *V3_lab)
{
    V3_lab[0] = V3_geo[0] + fOrigin[0];
    V3_lab[1] = V3_geo[1] + fOrigin[1];
    V3_lab[2] = V3_geo[2] + fOrigin[2];
}

void G2PGeoBase::Geo2LabNoTNoR(const double *V3_geo, double *V3_lab)
{
    V3_lab[0] = V3_geo[0];
    V3_lab[1] = V3_geo[1];
    V3_lab[2] = V3_geo[2];
}

void G2PGeoBase::MakePrefix()
{
    // G2PGeoBase classes do not need prefix

    const char *base = "";

    G2PAppBase::MakePrefix(base);
}

ClassImp(G2PGeoBase)
