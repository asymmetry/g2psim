// -*- C++ -*-

/* class G2PGeoBase
 * Abstract base class for g2p geometries.
 * It provides fundamental functions like rotation.
 * No instance allowed for this class.
 * The rotation matrix is defined with Euler angle in Z-X'-Z" convention.
 * G2PDrift class uses AtBoundary() to test the boundary of a geometry and stop drifting.
 */

// History:
//   Nov 2014, C. Gu, Add this class for g2p geometries.
//

#ifndef G2P_GEOBASE_H
#define G2P_GEOBASE_H

#include "G2PAppBase.hh"

using namespace std;

class G2PGeoBase : public G2PAppBase
{
public:
    virtual ~G2PGeoBase();

    virtual int Begin();

    virtual bool AtBoundary(double *V3) = 0;

    // Gets

    // Sets
    void SetOrigin(double x, double y, double z);
    void SetEulerAngle(double alpha, double beta, double gamma);

protected:
    G2PGeoBase(); // No instance allowed for this class

    void SetRotationMatrix();

    void Lab2Geo(const double *V3_lab, double *V3_geo);
    void Geo2Lab(const double *V3_geo, double *V3_lab);

    void FieldLab2Geo(const double *V3_lab, double *V3_geo);
    void FieldGeo2Lab(const double *V3_geo, double *V3_lab);

    double fOrigin[3];

    bool fRotation;
    double fEulerAngle[3];
    double fRotationMatrix[2][3][3];

private:
    ClassDef(G2PGeoBase, 1)
};

#endif
