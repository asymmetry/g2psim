// -*- C++ -*-

/* class G2PTransBase
 * Abstract base class of HRSTrans classes.
 * It provides interface functions.
 */

// History:
//   Sep 2013, C. Gu, First public version.
//

#ifndef HRSTRANS_BASE_H
#define HRSTRANS_BASE_H

class HRSTransBase {
public:
    HRSTransBase();
    virtual ~HRSTransBase();

    virtual int TransLeftHRS(double* v, double* pposx, double* pposy) = 0;
    virtual int TransRightHRS(double* v, double* pposx, double* pposy) = 0;
    virtual void ReconLeftHRS(double* v) = 0;
    virtual void ReconRightHRS(double* v) = 0;

    virtual void FPCorrLeft(const double* V5tg, double* V5fp);
    virtual void FPCorrRight(const double* V5tg, double* V5fp);

    // Rotate on X axis in transport coordinate
    // Notice the positive direction is anti-clockwise
    // Z is assumed to be 0
    void CoordsCorrection(double angle, double* v);

    // Gets
    virtual double GetAngle();

    // Sets

protected:
    double fModelAngle;
};

#endif
