#ifndef HRS_GUN_H
#define HRS_GUN_H

#include <cstdio>
#include <cstring>

#include "HRSRand.hh"

class HRSGun
{
public:
    HRSGun();
    HRSGun(const char * dist);
    ~HRSGun();

    typedef void (HRSGun::*fun_ptr)(double *, double *);
    
    void SetTargetZRange(double low, double high) { pTargetZLow_lab = low; pTargetZHigh_lab = high; }
    void SetTargetX(double value) { pTargetX_lab = value; }
    void SetTargetY(double value) { pTargetY_lab = value; }
    void SetTargetZ(double value) { pTargetZLow_lab = value; }
    void SetTargetR(double value) { pTargetR_lab = value; }
    void SetThetaRange(double low, double high) { pMomentThetaLow_tr = low; pMomentThetaHigh_tr = high; }
    void SetPhiRange(double low, double high) { pMomentPhiLow_tr = low; pMomentPhiHigh_tr = high; }
    void SetTheta(double value) { pMomentThetaLow_tr = value; }
    void SetPhi(double value) { pMomentPhiLow_tr = value; }
    void SetDeltaRange(double low, double high) { pDeltaLow = low; pDeltaHigh = high; }
    void SetDelta(double value) { pDeltaLow = value; }
    void SetHRSAngle(double value) { pHRSAngle = value; }

    void SetPositionRes(double value) { pPosRes = value; }
    void SetAngleRes(double value) { pAngleRes = value; }
    void SetDeltaRes(double value) { pDeltaRes = value; }
    
    void SetDataFile(const char *name) { strcpy(pDataFile, name); }

    void GetPV5(double *pV3, double *pV5) { (this->*pGunSelector)(pV3, pV5); }

private:
    void SetGun(int dist);
    void GetPV5Delta(double *pV3, double *pV5);
    void GetPV5Gaus(double *pV3, double *pV5);
    void GetPV5Flat(double *pV3, double *pV5);
    void GetPV5Sieve(double *pV3, double *pV5);
    void GetPV5Data(double *pV3, double *pV5);
    void GetPV5Focus(double *pV3, double *pV5);

    double pV3bpm_lab[3];
    double pV5tg_tr[5];
    double pV5fp_tr[5];

    double pTargetX_lab;
    double pTargetY_lab;
    double pTargetZLow_lab;
    double pTargetZHigh_lab;
    double pTargetR_lab;
    
    double pMomentThetaLow_tr;
    double pMomentThetaHigh_tr;
    double pMomentPhiLow_tr;
    double pMomentPhiHigh_tr;

    double pDeltaLow; // in the unit of delta
    double pDeltaHigh;

    double pPosRes;
    double pAngleRes;
    double pDeltaRes;

    double pHRSAngle;

    fun_ptr pGunSelector;

    FILE *pFilePtr;
    char pDataFile[300];

    HRSRand *pRand;
};

#endif
