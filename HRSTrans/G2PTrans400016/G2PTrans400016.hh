#ifndef G2P_TRANS400016_H
#define G2P_TRANS400016_H

#include "HRSTransBase.hh"

class G2PTrans400016 : public HRSTransBase
{
public:
    G2PTrans400016();
    ~G2PTrans400016();

    bool TransLeftHRS(double* vector_jjl);
    bool TransRightHRS(double* vector_jjl);
    void ReconLeftHRS(double* vector_jjl);
    void ReconRightHRS(double* vector_jjl);

    double GetAngle() { return cModelAngle; }

private:
    const double cModelAngle;
};

#endif
