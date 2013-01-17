// This file defined a class HRSGun.
// This class will be used in g2pSim class as particle gun.
// It has 5 standard gun: ShootDelta(), ShootGaus(), ShootFlat(), ShootSieve()
//+and ShootData(). They can be chosen when initializing this class.
// g2pSim class will call Shoot() to get kinematic variables. It is a virtual
//+function so you can rewrite it by inheriting this class.
//
// History:
// By C. Gu, Jan 12, 2013, First public version.
//

#ifndef HRS_GUN_H
#define HRS_GUN_H

#include <cstdio>

#include "HRSRand.hh"

class HRSGun
{
public:
    HRSGun();
    HRSGun(const char * dist);
    ~HRSGun();

    typedef bool (HRSGun::*pf_Gun)(double *, double *);
    
    void SetHRSAngle(double value) { fHRSAngle = value; }
    void SetTargetX(double value) { fTargetX_lab = value; }
    void SetTargetY(double value) { fTargetY_lab = value; }
    void SetTargetZ(double value) { fTargetZLow_lab = value; }
    void SetTargetZRange(double low, double high) { fTargetZLow_lab = low; fTargetZHigh_lab = high; }
    void SetTargetR(double value) { fTargetR_lab = value; }

    void SetTheta(double value) { fTargetThLow_tr = value; }
    void SetThetaRange(double low, double high) { fTargetThLow_tr = low; fTargetThHigh_tr = high; }
    void SetPhi(double value) { fTargetPhLow_tr = value; }
    void SetPhiRange(double low, double high) { fTargetPhLow_tr = low; fTargetPhHigh_tr = high; }

    void SetDelta(double value) { fDeltaLow = value; }
    void SetDeltaRange(double low, double high) { fDeltaLow = low; fDeltaHigh = high; }
    
    void SetPositionRes(double value) { fPosRes = value; }
    void SetAngleRes(double value) { fAngleRes = value; }
    void SetDeltaRes(double value) { fDeltaRes = value; }
    
    void SetDataFile(const char *name) { pFileName = name; }

    void SetRand(HRSRand * rand) { pRand = rand; }

    bool IsInit() { return bIsInit; }
    bool IsUsingData() { return bUseData; }

    int GetSetting() { return iSetting; }
    double GetPosResolution() { return fPosRes; }
    double GetAngleResolution() { return fAngleRes; }
    double GetDeltaResolution() { return fDeltaRes; }

    virtual void Init();
    virtual bool Shoot(double *V3bpm, double *V5tg) { return (this->*pfGunSelector)(V3bpm, V5tg); }
    virtual void End();

private:
    void SetGun(int setting);

    bool ShootDelta(double *V3bpm, double *V5tg);
    bool ShootGaus(double *V3bpm, double *V5tg);
    bool ShootFlat(double *V3bpm, double *V5tg);
    bool ShootSieve(double *V3bpm, double *V5tg);
    bool ShootData(double *V3bpm, double *V5tg);

    bool bIsInit;
    
    int iSetting;
    bool bUseData;

    double fHRSAngle;

    double fTargetX_lab;
    double fTargetY_lab;
    double fTargetZLow_lab;
    double fTargetZHigh_lab;
    double fTargetR_lab;

    double fTargetThLow_tr;
    double fTargetThHigh_tr;
    double fTargetPhLow_tr;
    double fTargetPhHigh_tr;

    double fDeltaLow; // in the unit of delta
    double fDeltaHigh;

    double fPosRes;
    double fAngleRes;
    double fDeltaRes;

    FILE *pFilePtr;
    const char *pFileName;

    HRSRand *pRand;
    
    pf_Gun pfGunSelector;
};

#endif
