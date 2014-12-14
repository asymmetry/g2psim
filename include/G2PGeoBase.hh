// -*- C++ -*-

/* class G2PGeoBase
 * Abstract base class for g2p geometries.
 * It provides fundamental functions like rotation.
 * No instance allowed for this class.
 * The rotation matrix is defined with Euler angle in Z-X'-Z" convention.
 * G2PDrift class uses TouchBoundary() to determine whether stop drifting or not.
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

    typedef void (G2PGeoBase::*pfX2X_)(const double *, double *);

    virtual int Begin();

    virtual bool TouchBoundary(const double *V5_tr, double z_tr);
    virtual bool TouchBoundary(const double *V3_lab);
    virtual bool TouchBoundary(double x, double y, double z);

    // Gets
    bool UseTrans();

    // Sets
    void SetOrigin(double x, double y, double z);
    void SetEulerAngle(double alpha, double beta, double gamma);

protected:
    G2PGeoBase(); // No instance allowed for this class

    void SetGeoPosition();

    void Lab2Geo(const double *V3_lab, double *V3_geo);
    void Lab2GeoNoT(const double *V3_lab, double *V3_geo);
    void Lab2GeoNoR(const double *V3_lab, double *V3_geo);
    void Lab2GeoNoTNoR(const double *V3_lab, double *V3_geo);

    void Geo2Lab(const double *V3_geo, double *V3_lab);
    void Geo2LabNoT(const double *V3_geo, double *V3_lab);
    void Geo2LabNoR(const double *V3_geo, double *V3_lab);
    void Geo2LabNoTNoR(const double *V3_geo, double *V3_lab);

    virtual bool TouchBoundaryGeo(double x, double y, double z) = 0; // true boundary is defined here

    bool fUseTrans; // flag indicates whether defined in transport coords

    bool fTranslation;
    double fOrigin[3];

    bool fRotation;
    double fEulerAngle[3];
    double fRotationMatrix[2][3][3];

    pfX2X_ pfLab2Geo;
    pfX2X_ pfGeo2Lab;

private:
    ClassDef(G2PGeoBase, 1)
};

#endif
