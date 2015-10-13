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

#ifndef G2P_GEOBASE_H
#define G2P_GEOBASE_H

#include "G2PAppBase.hh"

class G2PMaterial;

class G2PGeoBase : public G2PAppBase
{
public:
    virtual ~G2PGeoBase();

    typedef void (G2PGeoBase::*pfX2X_)(const double *, double *);

    virtual int Begin();

    virtual bool IsInside(const double *V5_tr, double z_tr);
    virtual bool IsInside(const double *V3);

    // Gets
    virtual G2PMaterial *GetMaterial() const;
    int GetNUsed() const;

    // Sets
    virtual void SetOrigin(double x, double y, double z);
    virtual void SetEulerAngle(double alpha, double beta, double gamma);
    virtual void SetMaterial(G2PMaterial *mat);
    void SetNUsed(int n);

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

    virtual bool IsInside(double x, double y, double z) = 0;

    virtual void MakePrefix();

    int fNUsed;

    bool fTranslation;
    double fOrigin[3];

    bool fRotation;
    double fEulerAngle[3];
    double fRotationMatrix[2][3][3];

    G2PMaterial *pMaterial;

    pfX2X_ pfLab2Geo;
    pfX2X_ pfGeo2Lab;

private:
    ClassDef(G2PGeoBase, 1)
};

#endif
