#ifndef G2P_TARGETFIELD_H
#define G2P_TARGETFIELD_H

#include "TROOT.h"
#include "TObject.h"

#include "G2PHallBField.hh"

class G2PTargetField : public TObject
{
public:
    G2PTargetField();
    G2PTargetField(const char* name);
    ~G2PTargetField();

    typedef bool (G2PTargetField::*pf_Field)();

    void SetMapFile(const char* mapfile) { pMapFile = mapfile; }

    void SetZRange(double zmin, double zmax) { fZMin = zmin; fZMax = zmax; }
    void SetRRange(double rmin, double rmax) { fRMin = rmin; fRMax = rmax; }
    void SetStepZ(double stepz) { fStepZ = stepz; }
    void SetStepR(double stepr) { fStepR = stepr; }
    void SetOrigin(double x, double y, double z) { fOrigin[0] = x; fOrigin[1] = y; fOrigin[2] = z; }
    void SetEulerAngle(double alpha, double beta, double gamma);
    void SetRatio(double ratio) { fRatio = ratio; }

    bool IsInit() { return bIsInit; }

    void Init();
    void End();

    void GetField(const double* x, double* b);

private:
    bool ReadMap();
    virtual bool CreateMap() { return (this->*pfFieldSelector)(); }
    bool CreateHallBMap();
    bool CreateUniformMap();

    bool Interpolate(const double* x, double* b, const int order);

    void TransLab2Field(double* x);
    void TransField2Lab(double* b);

    void SaveRootFile();

    bool bIsInit;

    const char* pMapFile;

    double*** fBField;

    double fOrigin[3];

    double fZMin, fZMax;
    double fRMin, fRMax;
    double fStepZ, fStepR;
    int nNumZ, nNumR;

    bool bRotation;
    double fEulerAngle[3];
    double fRotationMatrix[2][3][3];

    double fRatio;

    G2PHallBField GetHallBField; // Use this as a function

    pf_Field pfFieldSelector;
    
    ClassDef(G2PTargetField, 1);
};

#endif
