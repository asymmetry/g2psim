// This file defines a class G2PDrift.
// This class is a tool class.
// It use Nystrom-Runge-Kutta method to derive the trajectory of a charged
//+particle in static magnetic field.
// G2PProcBase classes will call Drift() to get the end point.
// The subroutine Drift() has 2 prototypes, one for lab corrds and one for
//+transportation coords.
//
// History:
//   Feb 2013, C. Gu, First public version.
//   Feb 2013, C. Gu, Change algorithm to Nystrom-Runge-Kutta method.
//   Mar 2013, C. Gu, Add flexible step length and boundary check.
//

#ifndef G2P_DRIFT_H
#define G2P_DRIFT_H

#include "G2PAppBase.hh"

class G2PField;

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

    G2PField* pField;

    pfDriftHCS_ pfDriftHCS;
    pfDriftTCS_ pfDriftTCS;

private:
    static G2PDrift* pG2PDrift;

    ClassDef(G2PDrift, 1)
};

#endif
