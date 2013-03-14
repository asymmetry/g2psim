#ifndef G2P_FIELDBASE_H
#define G2P_FIELDBASE_H

#include <vector>

#include "G2PAppBase.hh"

using namespace std;

class G2PFieldBase : public G2PAppBase
{
public:
    virtual ~G2PFieldBase();

    void SetZRange(double zmin, double zmax) { fZMin = zmin; fZMax = zmax; }
    void SetRRange(double rmin, double rmax) { fRMin = rmin; fRMax = rmax; }
    void SetZStep(double stepz) { fZStep = stepz; }
    void SetRStep(double stepr) { fRStep = stepr; }
    void SetOrigin(double x, double y, double z) { fOrigin[0] = x; fOrigin[1] = y; fOrigin[2] = z; }
    void SetEulerAngle(double alpha, double beta, double gamma);
    void SetRatio(double ratio) { fRatio = ratio; }

    virtual int Init();
    virtual int Begin();

    virtual void GetField(const double* x, double* b);

    double GetRatio() { return fRatio; }

    static G2PFieldBase* GetInstance() { return pG2PFieldBase; }

protected:
    G2PFieldBase(); // No instance allowed for this class

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
    static G2PFieldBase* pG2PFieldBase;

    ClassDef(G2PFieldBase, 1)
};

#endif
