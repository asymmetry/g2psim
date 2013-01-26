#ifndef G2P_TRANS484816_H
#define G2P_TRANS484816_H

#include "../G2PTrans.hh"

class G2PTrans484816 : public G2PTrans
{
public:
    G2PTrans484816();
    ~G2PTrans484816();
    
    bool TransLeftHRS(double* vector_jjl);
    bool TransRightHRS(double* vector_jjl);
    void ReconLeftHRS(double* vector_jjl);
    void ReconRightHRS(double* vector_jjl);
};

#endif
