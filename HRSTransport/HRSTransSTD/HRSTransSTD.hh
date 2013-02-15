#ifndef G2P_TRANSNOSEPTA_H
#define G2P_TRANSNOSEPTA_H

#include "HRSTransBase.hh"

class HRSTransSTD : public HRSTransBase
{
public:
    HRSTransSTD();
    ~HRSTransSTD();

    bool TransLeftHRS(double* vector_jjl);
    bool TransRightHRS(double* vector_jjl);
    void ReconLeftHRS(double* vector_jjl);
    void ReconRightHRS(double* vector_jjl);

    double GetAngle() { return cModelAngle; }

private:
    const double cModelAngle;
};

#endif
