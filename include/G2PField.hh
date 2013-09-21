// -*- C++ -*-

/* class G2PField
 * Generate field map for HallB magnets.
 * Calculate field strength of a particular point from the field map using Lagrange polynomial interpolation, default is 2nd order.
 * G2PProcBase classes will call GetField() to get field values.
 *
 * Field map may have an angle respect to the lab coordinates.
 * Use SetEulerAngle() to set this angle and the program will rotate the field map to correct direction.
 */

// History:
//   Mar 2013, C. Gu, First public version.
//   Sep 2013, C. Gu, Put HallB field map into G2PField.
//

#ifndef G2P_FIELD_H
#define G2P_FIELD_H

#include <vector>

#include "G2PAppBase.hh"

using namespace std;

class G2PField : public G2PAppBase {
public:
    G2PField();
    virtual ~G2PField();

    virtual int Init();
    virtual int Begin();

    virtual void GetField(const double* x, double* b);

    // Gets

    // Sets
    void SetOrigin(double x, double y, double z);
    void SetEulerAngle(double alpha, double beta, double gamma);
    void SetZRange(double zmin, double zmax);
    void SetRRange(double rmin, double rmax);
    void SetZStep(double stepz);
    void SetRStep(double stepr);

protected:
    void SetRotationMatrix();

    virtual int ReadMap();
    virtual int CreateMap();

    virtual int Interpolate(const double* x, double* b, int order);

    void TransLab2Field(const double* x, double* xout);
    void TransField2Lab(const double* b, double* bout);

    void SaveRootFile();

    virtual int Configure(EMode mode = kTWOWAY);
    virtual void MakePrefix();

    const char* fMapFile;

    vector<vector<vector<double> > > fBField;

    double fOrigin[3];

    double fZMin, fZMax;
    double fRMin, fRMax;
    double fZStep, fRStep;
    int nZ, nR;

    bool fRotation;
    double fEulerAngle[3];
    double fRotationMatrix[2][3][3];

    double fRatio;

private:
    static G2PField* pG2PField;

    ClassDef(G2PField, 1)
};

#endif
