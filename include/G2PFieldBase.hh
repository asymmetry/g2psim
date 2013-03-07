#ifndef G2P_FIELDBASE_H
#define G2P_FIELDBASE_H

#include <vector>

#include "G2PAppsBase.hh"

using namespace std;

class G2PFieldBase : public G2PAppsBase
{
public:
    virtual ~G2PFieldBase();

    void SetRangeZ(double zmin, double zmax) { fZMin = zmin; fZMax = zmax; }
    void SetRangeR(double rmin, double rmax) { fRMin = rmin; fRMax = rmax; }
    void SetStepZ(double stepz) { fStepZ = stepz; }
    void SetStepR(double stepr) { fStepR = stepr; }
    void SetOrigin(double x, double y, double z) { fOrigin[0] = x; fOrigin[1] = y; fOrigin[2] = z; }
    void SetEulerAngle(double alpha, double beta, double gamma);
    void SetRatio(double ratio) { fRatio = ratio; }

    virtual EStatus Init();
    virtual void Clear() { }

    virtual void GetField(const double* x, double* b);

    double GetRatio() { return fRatio; }

    static G2PFieldBase* GetInstance() { return pG2PFieldBase; }

protected:
    G2PFieldBase();

    virtual bool ReadMap();
    virtual bool CreateMap();

    bool Interpolate(const double* x, double* b, int order);

    void TransLab2Field(const double* x, double* xout);
    void TransField2Lab(const double* b, double* bout);

    void SaveRootFile();

    const char* pMapFileName;

    vector<vector<vector<double> > > fBField;

    double fOrigin[3];

    double fZMin, fZMax;
    double fRMin, fRMax;
    double fStepZ, fStepR;
    int nNumZ, nNumR;

    bool bRotation;
    double fEulerAngle[3];
    double fRotationMatrix[2][3][3];

    double fRatio;

private:
    static G2PFieldBase* pG2PFieldBase;

    ClassDef(G2PFieldBase, 1)
};

#endif
