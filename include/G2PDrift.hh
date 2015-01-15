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
//   Dec 2014, C. Gu, Merge Drift() functions into one. The stop condition is set by Condition class.
//

#ifndef G2P_DRIFT_H
#define G2P_DRIFT_H

#include "G2PAppBase.hh"

class G2PField;
class G2PGeoBase;

class Condition
{
public:
    Condition(double zi, double zf);
    Condition(double zi_tr, double zf_tr, double angle);
    Condition(G2PGeoBase *geo);
    ~Condition();

    bool operator()(const double *x);

private:
    int fType;

    double fZi, fZf;
    bool fSign;
    double fSinAng, fCosAng;

    G2PGeoBase *pGeo;
};

class G2PDrift : public G2PAppBase
{
public:
    G2PDrift();
    virtual ~G2PDrift();

    typedef double (G2PDrift::*pfDriftHCS_)(const double *, const double *, Condition &, double *, double *);

    virtual int Begin();

    virtual double Drift(const char *dir, const double *x, const double *p, Condition &stop, double *xout, double *pout);

    // Gets

    // Sets
    void SetStep(double init, double limit);
    void SetErrLimit(double lo, double hi);

protected:
    double DriftHCS(const double *x, const double *p, Condition &stop, double *xout, double *pout);
    double DriftHCSNF(const double *x, const double *p, Condition &stop, double *xout, double *pout);

    void NystromRK4(const double *x, const double *dxdt, double step, double *xo, double *err);
    double DistChord();
    void ComputeRHS(const double *x, double *dxdt);

    virtual int Configure(EMode mode = kTWOWAY);
    virtual void MakePrefix();

    double fM0;
    double fQ;
    double fFieldRatio;

    double fStep, fStepLimit, fErrLoLimit, fErrHiLimit;
    double fVelocity, fVelocity2, fGamma;
    double fCof;
    double fField[3];
    double fIPoint[3];
    double fMPoint[3];
    double fEPoint[3];

    G2PField *pField;

    pfDriftHCS_ pfDriftHCS;

private:
    static G2PDrift *pG2PDrift;

    ClassDef(G2PDrift, 1)
};

#endif
