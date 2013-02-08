#ifndef G2P_TRANS484816R00_H
#define G2P_TRANS484816R00_H

#include "../G2PTrans.hh"

class G2PTrans484816R00 : public G2PTrans
{
public:
    G2PTrans484816R00();
    ~G2PTrans484816R00();

    bool TransLeftHRS(double* vector_jjl);
    bool TransRightHRS(double* vector_jjl);
    void ReconLeftHRS(double* vector_jjl);
    void ReconRightHRS(double* vector_jjl);

    double GetAngle() { return cModelAngle; }

    void FPCorrection(const double* V5tg, double* V5fp);

private:
    const double cModelAngle;
};

#endif
