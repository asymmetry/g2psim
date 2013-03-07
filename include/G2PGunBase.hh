// This file defines a class G2PGunBase.
// This class is the base class of G2PGun.
// G2PRun class will call Shoot() to get kinematic variables. It is a virtual
//+function so each derived class will have its own method.
//
// History:
//   Mar 2013, C. Gu, First public version.
//

#ifndef G2P_GUNBASE_H
#define G2P_GUNBASE_H

#include "G2PAppsBase.hh"

class G2PDrift;

class G2PGunBase : public G2PAppsBase
{
public:
    G2PGunBase();
    ~G2PGunBase();
    
    void SetBeamX(double value) { fBeamX_lab = value; }
    void SetBeamY(double value) { fBeamY_lab = value; }
    void SetBeamR(double value) { fBeamR_lab = value; }
    
    void SetReactZ(double value) { fReactZLow_lab = value; fReactZHigh_lab = value;}
    void SetReactZRange(double low, double high) { fReactZLow_lab = low; fReactZHigh_lab = high; }

    void SetTargetTh(double value) { fTargetThLow_tr = value; fTargetThHigh_tr = value; }
    void SetTargetThRange(double low, double high) { fTargetThLow_tr = low; fTargetThHigh_tr = high; }
    void SetTargetPh(double value) { fTargetPhLow_tr = value; fTargetPhHigh_tr = value; }
    void SetTargetPhRange(double low, double high) { fTargetPhLow_tr = low; fTargetPhHigh_tr = high; }

    void SetDelta(double value) { fDeltaLow = value; fDeltaHigh = value; }
    void SetDeltaRange(double low, double high) { fDeltaLow = low; fDeltaHigh = high; }
    
    void SetSigmaPosLab(double value) { fSigmaPos_lab = value; }
    void SetSigmaAngLab(double value) { fSigmaAng_lab = value; }
    void SetSigmaAngTr(double value) { fSigmaAng_tr = value; }
    void SetSigmaDelta(double value) { fSigmaDelta = value; }

    virtual EStatus Init();
    virtual void Clear() { }

    virtual bool Shoot(double* V51, double* V52, double* V53 = NULL) = 0;

    virtual bool UseData() = 0;

    virtual int RegisterModel();

    static G2PGunBase* GetInstance() { return pG2PGunBase; }

protected:
    virtual int SetTiltAngle();
    virtual void GetReactPoint(double x, double y, double z, double* V5);

    double fHRSAngle;
    double fHRSMomentum;
    double fBeamEnergy;

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
    static G2PGunBase* pG2PGunBase;

    ClassDef(G2PGunBase, 1)
};

#endif
