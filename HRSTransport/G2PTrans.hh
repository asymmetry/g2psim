#ifndef G2P_TRANS_H
#define G2P_TRANS_H

class G2PTrans
{
public:
    G2PTrans();
    virtual ~G2PTrans();

    void SetHRSAngle(double value) { fHRSAngle = value; } 

    virtual bool TransLeftHRS(double* v) = 0;
    virtual bool TransRightHRS(double* v) = 0;
    virtual void ReconLeftHRS(double* v) = 0;
    virtual void ReconRightHRS(double* v) = 0;

protected:
    // Rotate on X axis in transport coordinate
    // Notice the positive direction is anti-clockwise
    // Z is assumed to be 0
    void RotateX(double angle, double* v);

    double fHRSAngle;
};

#endif
