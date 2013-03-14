#ifndef G2P_DRIFT_H
#define G2P_DRIFT_H

#include "G2PAppBase.hh"

class G2PFieldBase;

class G2PDrift : public G2PAppBase
{
public:
    G2PDrift();
    ~G2PDrift();

    typedef void (G2PDrift::*pfDriftHCS_)(const double*, const double*, double, double, double*, double*);
    typedef void (G2PDrift::*pfDriftTCS_)(const double*, double, double, double, double, double, double*);

    void SetLimit(double lo, double hi) { fErrLoLimit = lo; fErrHiLimit = hi; }

    int Init();
    int Begin();
    void Clear();

    void Drift(const double* x, const double* p, double zlimit, double llimit, double *xout, double *pout);
    void Drift(const double* x, double p, double z_tr, double angle, double zlimit, double llimit, double* xout);

    static G2PDrift* GetInstance() { return pG2PDrift; }

protected:
    void DriftHCS(const double* x, const double* p, double zlimit, double llimit, double *xout, double *pout);
    void DriftHCSNF(const double* x, const double* p, double zlimit, double llimit, double *xout, double *pout);
    void DriftTCS(const double* x, double p, double z_tr, double angle, double zlimit, double llimit, double* xout);
    void DriftTCSNF(const double* x, double p, double z_tr, double angle, double zlimit, double llimit, double* xout);

    void NystromRK4(const double* x, const double* dxdt, double step, double* xo, double* err);
    double DistChord();
    void ComputeRHS(const double* x, double* dxdt);

    double fM0;
    double fQ, fQsave;
    double fStep, fStepLimit, fErrLoLimit, fErrHiLimit;
    double fVelocity, fVelocity2, fGamma;
    double fCof;
    double fField[3];
    double fIPoint[3];
    double fMPoint[3];
    double fEPoint[3];

    G2PFieldBase* pField;

    pfDriftHCS_ pfDriftHCS;
    pfDriftTCS_ pfDriftTCS;

private:
    static G2PDrift* pG2PDrift;

    ClassDef(G2PDrift, 1)
};

#endif
