#ifndef G2P_TRANS_H
#define G2P_TRANS_H

class G2PTrans
{
public:
    G2PTrans();
    ~G2PTrans();

    virtual bool TransLeftHRS(double* v) = 0;
    virtual bool TransRightHRS(double* v) = 0;
    virtual void ReconLeftHRS(double* v) = 0;
    virtual void ReconRightHRS(double* v) = 0;
};

#endif
