// -*- C++ -*-

/* class G2PGeoBase
 * Abstract base class for g2p geometries.
 * It provides fundamental functions like rotation.
 *
 * Use SetEulerAngle() to set the pointing direction of the geometry.
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

using namespace std;

G2PGeoBase::G2PGeoBase() : fRotation(false)
{
    memset(fOrigin, 0, sizeof(fOrigin));
    memset(fEulerAngle, 0, sizeof(fEulerAngle));
    memset(fRotationMatrix, 0, sizeof(fRotationMatrix));
}

G2PGeoBase::~G2PGeoBase()
{
    // Nothing to do
}

int G2PGeoBase::Init()
{
    // Default does nothing

    return (G2PAppBase::Init());
}

int G2PGeoBase::Begin()
{
    if (G2PAppBase::Begin() != 0)
        return fStatus;

    SetRotationMatrix();
    return (G2PAppBase::Begin());
}

int G2PGeoBase::End()
{
    // Default does nothing

    return (G2PAppBase::End());
}

void G2PGeoBase::Clear(Option_t *option)
{
    // Default does nothing

    G2PAppBase::Clear(option);
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

void G2PGeoBase::SetRotationMatrix()
{
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
    double temp[3] = {V3_lab[0], V3_lab[1], V3_lab[2]};

    for (int i = 0; i < 3; i++)
        temp[i] -= fOrigin[i];

    if (fRotation) {
        for (int i = 0; i < 3; i++)
            V3_geo[i] = fRotationMatrix[0][i][0] * temp[0] + fRotationMatrix[0][i][1] * temp[1] + fRotationMatrix[0][i][2] * temp[2];
    } else {
        for (int i = 0; i < 3; i++)
            V3_geo[i] = temp[i];
    }
}

void G2PGeoBase::Geo2Lab(const double *V3_geo, double *V3_lab)
{
    double temp[3];

    if (fRotation) {
        for (int i = 0; i < 3; i++)
            temp[i] = fRotationMatrix[1][i][0] * V3_geo[0] + fRotationMatrix[1][i][1] * V3_geo[1] + fRotationMatrix[1][i][2] * V3_geo[2];
    } else {
        for (int i = 0; i < 3; i++)
            temp[i] = V3_geo[i];
    }

    for (int i = 0; i < 3; i++)
        V3_lab[i] = temp[i] + fOrigin[i];
}

void G2PGeoBase::FieldLab2Geo(const double *V3_lab, double *V3_geo)
{
    if (fRotation) {
        for (int i = 0; i < 3; i++)
            V3_geo[i] = fRotationMatrix[0][i][0] * V3_lab[0] + fRotationMatrix[0][i][1] * V3_lab[1] + fRotationMatrix[0][i][2] * V3_lab[2];
    } else {
        for (int i = 0; i < 3; i++)
            V3_geo[i] = V3_lab[i];
    }
}

void G2PGeoBase::FieldGeo2Lab(const double *V3_geo, double *V3_lab)
{
    if (fRotation) {
        for (int i = 0; i < 3; i++)
            V3_lab[i] = fRotationMatrix[1][i][0] * V3_geo[0] + fRotationMatrix[1][i][1] * V3_geo[1] + fRotationMatrix[1][i][2] * V3_geo[2];
    } else {
        for (int i = 0; i < 3; i++)
            V3_lab[i] = V3_geo[i];
    }
}

ClassImp(G2PGeoBase)
