// -*- C++ -*-

/* class G2PGun
 * This file defines a class G2PGun.
 * It is the base class of g2p event generator classes.
 * G2PProcBase classes will call Shoot() to get reaction point kinematics.
 * Shoot() is a pure virtual method so each derived class has its own implement.
 * G2PDrift is used in this class.
 */

// History:
//   Mar 2013, C. Gu, First public version.
//

#ifndef G2P_GUNBASE_H
#define G2P_GUNBASE_H

#include "G2PAppBase.hh"

class G2PDrift;

class G2PGun : public G2PAppBase {
public:
    G2PGun();
    ~G2PGun();

    void SetBeamX(double value) {
        fBeamX_lab = value;
    }

    void SetBeamY(double value) {
        fBeamY_lab = value;
    }

    void SetBeamR(double value) {
        fBeamR_lab = value;
    }

    void SetReactZ(double value) {
        fReactZLow_lab = value;
        fReactZHigh_lab = value;
    }

    void SetReactZRange(double low, double high) {
        fReactZLow_lab = low;
        fReactZHigh_lab = high;
    }

    void SetTargetTh(double value) {
        fTargetThLow_tr = value;
        fTargetThHigh_tr = value;
    }

    void SetTargetThRange(double low, double high) {
        fTargetThLow_tr = low;
        fTargetThHigh_tr = high;
    }

    void SetTargetPh(double value) {
        fTargetPhLow_tr = value;
        fTargetPhHigh_tr = value;
    }

    void SetTargetPhRange(double low, double high) {
        fTargetPhLow_tr = low;
        fTargetPhHigh_tr = high;
    }

    void SetDelta(double value) {
        fDeltaLow = value;
        fDeltaHigh = value;
    }

    void SetDeltaRange(double low, double high) {
        fDeltaLow = low;
        fDeltaHigh = high;
    }

    void SetSigmaPosLab(double value) {
        fSigmaPos_lab = value;
    }

    void SetSigmaAngLab(double value) {
        fSigmaAng_lab = value;
    }

    void SetSigmaAngTr(double value) {
        fSigmaAng_tr = value;
    }

    void SetSigmaDelta(double value) {
        fSigmaDelta = value;
    }

    virtual int Init();
    virtual int Begin();

    virtual int Shoot(double* V51, double* V52, double* V53 = NULL) = 0;

    virtual bool UseData() = 0;

    static G2PGun* GetInstance() {
        return pG2PGun;
    }

protected:
    virtual void SetTiltAngle();
    virtual void GetReactPoint(double x, double y, double z, double* V5);

    double fHRSAngle;
    double fHRSMomentum;
    double fBeamEnergy;
    double fFieldRatio;

    double fBeamX_lab, fBeamY_lab;
    double fBeamTiltAngle;
    double fBeamR_lab;

    double fReactZLow_lab;
    double fReactZHigh_lab;

    double fTargetThLow_tr;
    double fTargetThHigh_tr;
    double fTargetPhLow_tr;
    double fTargetPhHigh_tr;

    double fDeltaLow; // in the unit of delta
    double fDeltaHigh;

    double fSigmaPos_lab;
    double fSigmaAng_lab;
    double fSigmaAng_tr;
    double fSigmaDelta;

    G2PDrift* pDrift;

private:
    static G2PGun* pG2PGun;

    ClassDef(G2PGun, 1)
};

#endif
