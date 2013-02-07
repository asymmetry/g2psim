#ifndef G2P_TRANS484816R20_H
#define G2P_TRANS484816R20_H

#include "../G2PTrans.hh"

class G2PTrans484816R20 : public G2PTrans
{
public:
    G2PTrans484816R20();
    ~G2PTrans484816R20();
    
    bool TransLeftHRS(double* vector_jjl);
    bool TransRightHRS(double* vector_jjl);
    void ReconLeftHRS(double* vector_jjl);
    void ReconRightHRS(double* vector_jjl);

    double GetAngle() { return cModelAngle; }

private:
    const double cModelAngle;
};

#endif
