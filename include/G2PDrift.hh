#ifndef G2P_DRIFT_H
#define G2P_DRIFT_H

#include "G2PTargetField.hh"

namespace G2PDrift
{   
    void SetField(G2PTargetField* field);
    void SetMass(double value);
    void SetCharge(double value);
    void SetStep(double value);

    bool HasField();

    void Drift(const double* x, const double* p, double zlimit, double llimit, double *xout, double *pout);
    void Drift(const double* x, double p, double angle, double z_tr, double zlimit, double llimit, double* xout);
};

#endif
