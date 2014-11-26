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

#ifndef G2P_GEOBASE_H
#define G2P_GEOBASE_H

#include "G2PAppBase.hh"

using namespace std;

class G2PGeoBase : public G2PAppBase
{
public:
    G2PGeoBase();
    virtual ~G2PGeoBase();

    virtual int Init();
    virtual int Begin();
    virtual int End();
    virtual void Clear(Option_t * /*option*/ = "");

    // Gets

    // Sets
    void SetOrigin(double x, double y, double z);
    void SetEulerAngle(double alpha, double beta, double gamma);

protected:
    void SetRotationMatrix();

    void Lab2Geo(const double *V3_lab, double *V3_geo);
    void Geo2Lab(const double *V3_geo, double *V3_lab);

    void FieldLab2Geo(const double *V3_lab, double *V3_geo);
    void FieldGeo2Lab(const double *V3_geo, double *V3_lab);

    virtual int Configure(EMode mode = kTWOWAY) = 0;
    virtual void MakePrefix() = 0;

    double fOrigin[3];

    bool fRotation;
    double fEulerAngle[3];
    double fRotationMatrix[2][3][3];

private:
    ClassDef(G2PGeoBase, 1)
};

#endif
