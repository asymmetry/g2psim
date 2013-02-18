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

    void SetField(G2PTargetField* field) { pField = field; }

    void Drift(double* x, double* p, double zlimit, double llimit);
    
private:
    void NystromRK4(const double* x, const double* dxdt, double step, double* xo, double* err);
    double DistChord() const;
    void ComputeRHS(const double* x, double* dxdt);

    double fM0;
    double fQ;
    double fStep;
        
    double fField[3];
    double fCof;
    double fVelocity2;
    double fInPoint[3];
    double fMidPoint[3];
    double fFinPoint[3];

    G2PTargetField* pField;
    
    ClassDef(G2PDrift, 1);
};

#endif
