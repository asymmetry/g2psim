// -*- C++ -*-

/* class G2PDrift
 * Use Nystrom-Runge-Kutta method to derive the trajectory of a charged particle in static magnetic field.
 * G2PProcBase classes will call Drift() to get the end point position and momentum of the trajectory.
 * Drift() has 2 prototypes, one for lab coordinates and one for HRS transportation coordinates.
 */

// History:
//   Feb 2013, C. Gu, First public version.
//   Feb 2013, C. Gu, Change algorithm to Nystrom-Runge-Kutta method.
//   Mar 2013, C. Gu, Add flexible step length and boundary check.
//   Oct 2013, J. Liu, Add drift function to stop at a cylinder boundary.
//

#ifndef G2P_DRIFT_H
#define G2P_DRIFT_H

#include "G2PAppBase.hh"

class G2PField;

class G2PDrift : public G2PAppBase {
public:
    G2PDrift();
    virtual ~G2PDrift();

    typedef double (G2PDrift::*pfDriftHCS_)(const double*, const double*, double, double*, double*);
    typedef double (G2PDrift::*pfDriftTCS_)(const double*, double, double, double, double, double*);
    typedef double (G2PDrift::*pfDriftCV_)(const double*, double, double, double, double, double*, double&);
    typedef double (G2PDrift::*pfDriftCL_)(const double*, double, double, double, double, double, double*, double&, int&);

    virtual int Init();
    virtual int Begin();
    virtual void Clear(Option_t* /*option*/ = "");

    virtual double Drift(const double* x, const double* p, double zf, double *xout, double *pout); // HCS
    virtual double Drift(const double* x, double z_tr, double p, double angle, double zf_tr, double* xout); // TCS
    virtual double Drift(const double* x, double z_tr, double p, double angle, double rf_lab, double* xout, double& zout); // CV // J. Liu
    virtual double Drift(const double* x, double z_tr, double p, double angle, double rf_lab, double zf_lab, double* xout, double &zout, int& surf); // CL // J. Liu

    // Gets

    // Sets
    void SetStep(double init, double limit);
    void SetErrLimit(double lo, double hi);

protected:
    double DriftHCS(const double* x, const double* p, double zf, double *xout, double *pout);
    double DriftHCSNF(const double* x, const double* p, double zf, double *xout, double *pout);

    double DriftTCS(const double* x, double z_tr, double p, double angle, double zf_tr, double* xout);
    double DriftTCSNF(const double* x, double z_tr, double p, double angle, double zf_tr, double* xout);

    double DriftCV(const double* x, double z_tr, double p, double angle, double rf_lab, double* xout, double& zout); // J. Liu
    double DriftCVNF(const double* x, double z_tr, double p, double angle, double rf_lab, double* xout, double& zout);

    double DriftCL(const double* x, double z_tr, double p, double angle, double rf_lab, double zf_lab, double* xout, double &zout, int& surf); // J. Liu
    double DriftCLNF(const double* x, double z_tr, double p, double angle, double rf_lab, double zf_lab, double* xout, double &zout, int& surf); // J. Liu

    void NystromRK4(const double* x, const double* dxdt, double step, double* xo, double* err);
    double DistChord();
    void ComputeRHS(const double* x, double* dxdt);

    virtual int Configure(EMode mode = kTWOWAY);
    virtual void MakePrefix();

    double fM0;
    double fQ, fQSave;
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
    pfDriftCV_ pfDriftCV;
    pfDriftCL_ pfDriftCL;

private:
    static G2PDrift* pG2PDrift;

    ClassDef(G2PDrift, 1)
};

#endif
