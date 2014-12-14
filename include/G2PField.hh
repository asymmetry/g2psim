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
//   Nov 2014, C. Gu, Rewrite it with G2PGeoBase.
//

#ifndef G2P_FIELD_H
#define G2P_FIELD_H

//#define DEBUGWITHROOT

#include <vector>

#include "G2PGeoBase.hh"

using namespace std;

class G2PField : public G2PGeoBase
{
public:
    G2PField();
    virtual ~G2PField();

    virtual int Init();
    virtual int Begin();

    virtual void GetField(const double *x, double *b);

    // Gets

    // Sets
    void SetZRange(double zmin, double zmax);
    void SetRRange(double rmin, double rmax);
    void SetZStep(double stepz);
    void SetRStep(double stepr);

protected:
    virtual int ReadMap();
    virtual int CreateMap();

    virtual int Interpolate(const double *x, double *b, int order);

    bool TouchBoundaryGeo(double x, double y, double z);

#ifdef DEBUGWITHROOT
    void SaveRootFile();
#endif

    virtual int Configure(EMode mode = kTWOWAY);
    virtual void MakePrefix();

    const char *fMapFile;

    vector<vector<vector<double> > > fBField;

    double fZMin, fZMax;
    double fRMin, fRMax;
    double fZStep, fRStep;
    int nZ, nR;

    double fRatio;

    pfX2X_ pfGeo2LabField;

private:
    static G2PField *pG2PField;

    ClassDef(G2PField, 1)
};

#endif
