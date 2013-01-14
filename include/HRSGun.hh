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

    typedef void (HRSGun::*gun_ptr)(double *, double *);

    void SetTargetX(double value) { pTargetX_lab = value; }
    void SetTargetY(double value) { pTargetY_lab = value; }
    void SetTargetZ(double value) { pTargetZLow_lab = value; }
    void SetTargetZRange(double low, double high) { pTargetZLow_lab = low; pTargetZHigh_lab = high; }
    void SetTargetR(double value) { pTargetR_lab = value; }

    void SetTheta(double value) { pTargetThetaLow_tr = value; }
    void SetThetaRange(double low, double high) { pTargetThetaLow_tr = low; pTargetThetaHigh_tr = high; }
    void SetPhi(double value) { pTargetPhiLow_tr = value; }
    void SetPhiRange(double low, double high) { pTargetPhiLow_tr = low; pTargetPhiHigh_tr = high; }

    void SetDelta(double value) { pDeltaLow = value; }
    void SetDeltaRange(double low, double high) { pDeltaLow = low; pDeltaHigh = high; }

    void SetHRSAngle(double value) { pHRSAngle = value; }

    void SetPositionRes(double value) { pPosRes = value; }
    void SetAngleRes(double value) { pAngleRes = value; }
    void SetDeltaRes(double value) { pDeltaRes = value; }
    
    void SetDataFile(const char *name) { pFileName = name; pUseData = true; }

    void SetRand(HRSRand * rand) { pRand = rand; }

    bool IsInit() { return pIsInit(); }
    bool IsUsingData() { return UseData(); }

    int GetSetting() { return pSetting; }
    double GetPosResolution() { return pPosRes; }
    double GetAngleResolution() { return pAngleRes; }
    double GetDeltaResolution() { return pDeltaRes; }

    virtual void Init();
    virtual void Shoot(double *pV3, double *pV5) { (this->*pGunSelector)(pV3, pV5); }
    virtual void End();

private:
    void SetGun(int dist);
    void ShootDelta(double *pV3, double *pV5);
    void ShootGaus(double *pV3, double *pV5);
    void ShootFlat(double *pV3, double *pV5);
    void ShootSieve(double *pV3, double *pV5);
    void ShootData(double *pV3, double *pV5);

    bool pIsInit;
    bool pUseData;

    int pSetting;

    double pV3bpm_lab[3];
    double pV5tg_tr[5];
    double pV5fp_tr[5];

    double pTargetX_lab;
    double pTargetY_lab;
    double pTargetZLow_lab;
    double pTargetZHigh_lab;
    double pTargetR_lab;

    double pTargetThetaLow_tr;
    double pTargetThetaHigh_tr;
    double pTargetPhiLow_tr;
    double pTargetPhiHigh_tr;

    double pDeltaLow; // in the unit of delta
    double pDeltaHigh;

    double pPosRes;
    double pAngleRes;
    double pDeltaRes;

    double pHRSAngle;

    gun_ptr pGunSelector;

    FILE *pFilePtr;
    const char *pFileName;

    HRSRand *pRand;
};

#endif
