#ifndef G2P_DRIFT_H
#define G2P_DRIFT_H

#include "TROOT.h"
#include "TObject.h"

#include "G2PTargetField.hh"

class G2PDrift : public TObject
{
public:
    G2PDrift();
    G2PDrift(G2PTargetField *field);
    ~G2PDrift();

    void SetMass(double value) { fM0 = value; }
    void SetCharge(double value) { fQ = value; }
    void SetStepLength(double value) { fStep = value; }

    void SetHRSAngle(double value) { fHRSAngle = value; }
    void SetHRSMomentum(double value) { fHRSMomentum = value; }

    void SetField(G2PTargetField* field) { pField = field; }

    void Drift(const double* x, const double* p, double zlimit, double llimit, double *xout, double *pout);
    void Drift(const double* x, double z_tr, double zlimit, double llimit, double* xout);
    
private:
    void NystromRK4(const double* x, const double* dxdt, double step, double* xo, double* err);
    double DistChord() const;
    void ComputeRHS(const double* x, double* dxdt);

    double fM0;
    double fQ;
    double fStep;

    double fHRSAngle;
    double fHRSMomentum;
        
    double fField[3];
    double fVelocity, fVelocity2;
    double fGamma;
    double fCof;
    double fInPoint[3];
    double fMidPoint[3];
    double fFinPoint[3];

    G2PTargetField* pField;
    
    ClassDef(G2PDrift, 1);
};

#endif
