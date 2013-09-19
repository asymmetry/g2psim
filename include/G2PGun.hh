// -*- C++ -*-

/* class G2PGun
 * Abstract base class of g2p event generator classes.
 * Use function Shoot() to get reaction point kinematics.
 * Shoot() is a pure virtual method so each derived class should define its own implement.
 */

// History:
//   Mar 2013, C. Gu, First public version.
//   Sep 2013, C. Gu, Rewrite it as a G2PProcBase class.
//

#ifndef G2P_GUNBASE_H
#define G2P_GUNBASE_H

#include "G2PProcBase.hh"

class G2PDrift;

class G2PGun : public G2PProcBase {
public:
    G2PGun();
    virtual ~G2PGun();

    virtual int Init();
    virtual int Begin();
    virtual int Process();
    virtual void Clear();

    // Gets

    // Sets
    void SetBeamPos(double x, double y);
    void SetReactZ(double low, double high);
    void SetRasterSize(double val);
    void SetTargetTh(double low, double high);
    void SetTargetPh(double low, double high);
    void SetDelta(double low, double high);

protected:
    virtual int Shoot(double* V51, double* V52) = 0;

    void SetTiltAngle();
    void GetReactPoint(double x, double y, double z, double* V5);

    virtual int Configure(EMode mode = kTWOWAY);
    virtual int DefineVariables(EMode mode = kDEFINE);
    virtual void MakePrefix();

    double fHRSAngle;
    double fHRSMomentum;
    double fBeamEnergy;
    double fFieldRatio;

    double fBeamX_lab, fBeamY_lab;
    double fBeamR_lab;

    double fReactZLow_lab;
    double fReactZHigh_lab;

    double fTargetThLow_tr;
    double fTargetThHigh_tr;
    double fTargetPhLow_tr;
    double fTargetPhHigh_tr;

    double fDeltaLow; // in the unit of delta
    double fDeltaHigh;

    double fBeamTiltAngle;

    double fV5beam_lab[5];
    double fV5react_tr[5];
    double fV5react_lab[5];

    double fV5tg_tr[5];

    G2PDrift* pDrift;

private:
    static G2PGun* pG2PGun;

    ClassDef(G2PGun, 1)
};

inline void G2PGun::SetBeamPos(double x, double y) {
    fBeamX_lab = x;
    fBeamY_lab = y;

    fConfigIsSet[&fBeamX_lab] = true;
    fConfigIsSet[&fBeamY_lab] = true;
}

inline void G2PGun::SetReactZ(double low, double high) {
    fReactZLow_lab = low;
    fReactZHigh_lab = high;

    fConfigIsSet[&fReactZLow_lab] = true;
    fConfigIsSet[&fReactZHigh_lab] = true;
}

inline void G2PGun::SetRasterSize(double val) {
    fBeamR_lab = val;

    fConfigIsSet[&fBeamR_lab] = true;
}

inline void G2PGun::SetTargetTh(double low, double high) {
    fTargetThLow_tr = low;
    fTargetThHigh_tr = high;

    fConfigIsSet[&fTargetThLow_tr] = true;
    fConfigIsSet[&fTargetThHigh_tr] = true;
}

inline void G2PGun::SetTargetPh(double low, double high) {
    fTargetPhLow_tr = low;
    fTargetPhHigh_tr = high;

    fConfigIsSet[&fTargetPhLow_tr] = true;
    fConfigIsSet[&fTargetPhHigh_tr] = true;
}

inline void G2PGun::SetDelta(double low, double high) {
    fDeltaLow = low;
    fDeltaHigh = high;

    fConfigIsSet[&fDeltaLow] = true;
    fConfigIsSet[&fDeltaHigh] = true;
}

#endif
