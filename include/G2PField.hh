// -*- C++ -*-

/* class G2PField
 * This file defines a class G2PField.
 * It is the base class of g2p target field classes.
 * It calculates field strength of a particular point from the field map using Lagrange polynomial interpolation, default is 2nd order.
 * G2PProcBase classes will call GetField() to get field values.
 *
 * Field map may have an angle respect to the lab coordinates.
 * Use SetEularAngle() to set this angle in Eular angles and the program will rotate the field map to correct direction.
 */

// History:
//   Mar 2013, C. Gu, First public version.
//

#ifndef G2P_FIELDBASE_H
#define G2P_FIELDBASE_H

#include <vector>

#include "G2PAppBase.hh"

using namespace std;

class G2PField : public G2PAppBase {
public:
    virtual ~G2PField();

    void SetZRange(double zmin, double zmax) {
        fZMin = zmin;
        fZMax = zmax;
    }

    void SetRRange(double rmin, double rmax) {
        fRMin = rmin;
        fRMax = rmax;
    }

    void SetZStep(double stepz) {
        fZStep = stepz;
    }

    void SetRStep(double stepr) {
        fRStep = stepr;
    }

    void SetOrigin(double x, double y, double z) {
        fOrigin[0] = x;
        fOrigin[1] = y;
        fOrigin[2] = z;
    }
    void SetEulerAngle(double alpha, double beta, double gamma);

    void SetRatio(double ratio) {
        fRatio = ratio;
    }

    virtual int Init();
    virtual int Begin();

    virtual void GetField(const double* x, double* b);

    double GetRatio() {
        return fRatio;
    }

    static G2PField* GetInstance() {
        return pG2PField;
    }

protected:
    G2PField(); // No instance allowed for this class

    virtual int ReadMap();
    virtual int CreateMap();

    int Interpolate(const double* x, double* b, int order);

    void TransLab2Field(const double* x, double* xout);
    void TransField2Lab(const double* b, double* bout);

    void SaveRootFile();

    const char* pMapFileName;

    vector<vector<vector<double> > > fBField;

    double fOrigin[3];

    double fZMin, fZMax;
    double fRMin, fRMax;
    double fZStep, fRStep;
    int nZ, nR;

    bool bRotation;
    double fEulerAngle[3];
    double fRotationMatrix[2][3][3];

    double fRatio;

private:
    static G2PField* pG2PField;

    ClassDef(G2PField, 1)
};

#endif
