#include <cstdio>
#include <cstdlib>
#include <string>
#include <map>

#include "HRSRand.hh"
#include "HRSTransTCSNHCS.hh"

#include "HRSGun.hh"

using namespace std;
using namespace Transform;

HRSGun::HRSGun()
  : pTargetX_lab(0), pTargetY_lab(0), pTargetZLow_lab(0),
    pTargetZHigh_lab(0), pTargetR_lab(0.015), pMomentThetaLow_tr(0),
    pMomentThetaHigh_tr(0), pMomentPhiLow_tr(0), pMomentPhiHigh_tr(0),
    pDeltaLow(0), pDeltaHigh(0), pPosRes(0.0001), pAngleRes(0.01),
    pDeltaRes(0.03), pHRSAngle(0), pGunSelector(NULL), pFilePtr(NULL)
{
    memset(pV3bpm_lab, 0, sizeof(pV3bpm_lab));
    memset(pV5tg_tr, 0, sizeof(pV5tg_tr));
    memset(pV5fp_tr, 0, sizeof(pV5fp_tr));
    memset(pDataFile, 0, sizeof(pDataFile));
    pRand = new HRSRand();
}

HRSGun::HRSGun(const char* dist)
  : pTargetX_lab(0), pTargetY_lab(0), pTargetZLow_lab(0),
    pTargetZHigh_lab(0), pTargetR_lab(0.015), pMomentThetaLow_tr(0),
    pMomentThetaHigh_tr(0), pMomentPhiLow_tr(0), pMomentPhiHigh_tr(0),
    pDeltaLow(0), pDeltaHigh(0), pPosRes(0.0001), pAngleRes(0.01),
    pDeltaRes(0.03), pHRSAngle(0), pGunSelector(NULL), pFilePtr(NULL)
{
    memset(pV3bpm_lab, 0, sizeof(pV3bpm_lab));
    memset(pV5tg_tr, 0, sizeof(pV5tg_tr));
    memset(pV5fp_tr, 0, sizeof(pV5fp_tr));
    memset(pDataFile, 0, sizeof(pDataFile));
    pRand = new HRSRand();

    map<string, int> dist_map;
    dist_map["delta"] = 1;
    dist_map["gaus"] = 2;
    dist_map["flat"] = 3;
    dist_map["sieve"] = 4;
    dist_map["data"] = 5;
    dist_map["focus"] = 6;

    if(dist_map[dist] == 0){
        printf("Unknown gun setting!\n");
        SetGun(1);
    }
    else{
        SetGun(dist_map[dist]);
    }
}

HRSGun::~HRSGun()
{
    delete[] pRand;
}

void HRSGun::SetGun(int dist)
{
    switch (dist) {
    case 1:
        pGunSelector = &HRSGun::GetPV5Delta;
        break;
    case 2:
        pGunSelector = &HRSGun::GetPV5Gaus;
        break;
    case 3:
        pGunSelector = &HRSGun::GetPV5Flat;
        break;
    case 4:
        pGunSelector = &HRSGun::GetPV5Sieve;
        break;
    case 5:
        pGunSelector = &HRSGun::GetPV5Data;
        break;
    case 6:
        pGunSelector = &HRSGun::GetPV5Focus;
        break;
    }
}

void HRSGun::GetPV5Delta(double *pV3, double *pV5)
{
    double Xtg_lab = pTargetX_lab;
    double Ytg_lab = pTargetY_lab;
    double Ztg_lab = pTargetZLow_lab;

    pV3[0] = Xtg_lab;
    pV3[1] = Ytg_lab;
    pV3[2] = Ztg_lab;

    double Xtg_tr, Ytg_tr, Ztg_tr;
    
    X_HCS2TCS(Xtg_lab, Ytg_lab, Ztg_lab, pHRSAngle, Xtg_tr, Ytg_tr, Ztg_tr);

    double Thetatg_tr = pMomentThetaLow_tr;
    double Phitg_tr = pMomentPhiLow_tr;

    Project(Xtg_tr, Ytg_tr, Ztg_tr, -Ztg_tr, Thetatg_tr, Phitg_tr);

    double Deltatg = pDeltaLow;
    
    pV5[0] = Xtg_tr;
    pV5[1] = Thetatg_tr;
    pV5[2] = Ytg_tr;
    pV5[3] = Phitg_tr;
    pV5[4] = Deltatg;
}

void HRSGun::GetPV5Gaus(double *pV3, double *pV5)
{
    double Xtg_lab = pRand->Gaus(pTargetX_lab, pPosRes);
    double Ytg_lab = pRand->Gaus(pTargetY_lab, pPosRes);
    double Ztg_lab = pRand->Gaus(pTargetZLow_lab, pPosRes);

    pV3[0] = Xtg_lab;
    pV3[1] = Ytg_lab;
    pV3[2] = Ztg_lab;

    double Xtg_tr, Ytg_tr, Ztg_tr;
    
    X_HCS2TCS(Xtg_lab, Ytg_lab, Ztg_lab, pHRSAngle, Xtg_tr, Ytg_tr, Ztg_tr);

    double Thetatg_tr = pRand->Gaus(pMomentThetaLow_tr, pAngleRes);
    double Phitg_tr = pRand->Gaus(pMomentPhiLow_tr, pAngleRes);

    Project(Xtg_tr, Ytg_tr, Ztg_tr, -Ztg_tr, Thetatg_tr, Phitg_tr);

    double Deltatg = pRand->Gaus(pDeltaLow, pDeltaRes);
    
    pV5[0] = Xtg_tr;
    pV5[1] = Thetatg_tr;
    pV5[2] = Ytg_tr;
    pV5[3] = Phitg_tr;
    pV5[4] = Deltatg;
}

void HRSGun::GetPV5Flat(double *pV3, double *pV5)
{
    double Xtg_lab, Ytg_lab;
    do {
        Xtg_lab = pRand->Uniform(-pTargetR_lab, pTargetR_lab);
        Ytg_lab = pRand->Uniform(-pTargetR_lab, pTargetR_lab);
    } while (Xtg_lab*Xtg_lab+Ytg_lab*Ytg_lab>pTargetR_lab*pTargetR_lab);
    
    double Ztg_lab = pRand->Uniform(pTargetZLow_lab, pTargetZHigh_lab);

    pV3[0] = Xtg_lab + pTargetX_lab;
    pV3[1] = Ytg_lab + pTargetY_lab;
    pV3[2] = Ztg_lab;

    double Xtg_tr, Ytg_tr, Ztg_tr;
    
    X_HCS2TCS(Xtg_lab, Ytg_lab, Ztg_lab, pHRSAngle, Xtg_tr, Ytg_tr, Ztg_tr);

    double Thetatg_tr = pRand->Uniform(pMomentThetaLow_tr, pMomentThetaHigh_tr);
    double Phitg_tr = pRand->Uniform(pMomentPhiLow_tr, pMomentPhiHigh_tr);

    Project(Xtg_tr, Ytg_tr, Ztg_tr, -Ztg_tr, Thetatg_tr, Phitg_tr);

    double Deltatg = pRand->Uniform(pDeltaLow, pDeltaHigh);

    pV5[0] = Xtg_tr;
    pV5[1] = Thetatg_tr;
    pV5[2] = Ytg_tr;
    pV5[3] = Phitg_tr;
    pV5[4] = Deltatg;
}

void HRSGun::GetPV5Sieve(double *pV3, double *pV5)
{
    double Xtg_lab, Ytg_lab;
    do {
        Xtg_lab = pRand->Uniform(-pTargetR_lab, pTargetR_lab);
        Ytg_lab = pRand->Uniform(-pTargetR_lab, pTargetR_lab);
    } while (Xtg_lab*Xtg_lab+Ytg_lab*Ytg_lab>pTargetR_lab*pTargetR_lab);
    
    double Ztg_lab = pRand->Uniform(pTargetZLow_lab, pTargetZHigh_lab);

    pV3[0] = Xtg_lab + pTargetX_lab;
    pV3[1] = Ytg_lab + pTargetY_lab;
    pV3[2] = Ztg_lab;

    double Xtg_tr, Ytg_tr, Ztg_tr;
    
    X_HCS2TCS(Xtg_lab, Ytg_lab, Ztg_lab, pHRSAngle, Xtg_tr, Ytg_tr, Ztg_tr);

    double Thetatg_tr = pRand->Uniform(pMomentThetaLow_tr, pMomentThetaHigh_tr);
    double Phitg_tr = pRand->Uniform(pMomentPhiLow_tr, pMomentPhiHigh_tr);

    Project(Xtg_tr, Ytg_tr, Ztg_tr, -Ztg_tr, Thetatg_tr, Phitg_tr);

    double Deltatg = pRand->Uniform(pDeltaLow, pDeltaHigh);

    pV5[0] = Xtg_tr;
    pV5[1] = Thetatg_tr;
    pV5[2] = Ytg_tr;
    pV5[3] = Phitg_tr;
    pV5[4] = Deltatg;
}

void HRSGun::GetPV5Data(double *pV3, double *pV5)
{
    int temp;

    if (pFilePtr==NULL) pFilePtr=fopen(pDataFile, "r");
    // index, bpm_x, bpm_y, thetatg_tr, ytg_tr, phitg_tr, deltatg 
    fscanf(pFilePtr, "%d%lf%lf%lf%lf%lf%lf", &temp, &pV3[0], &pV3[1], &pV5[1], &pV5[2], &pV5[3], &pV5[4]);

    pV3[2] = 0;
    pV5[0] = 0;
}

void HRSGun::GetPV5Focus(double *pV3, double *pV5)
{
    int temp;

    if (pFilePtr==NULL) pFilePtr=fopen(pDataFile, "r");
    // index, bpm_x, bpm_y, xfp_tr, thetafp_tr, yfp_tr, phitg_tr 
    fscanf(pFilePtr, "%d%lf%lf%lf%lf%lf%lf", &temp, &pV3[0], &pV3[1], &pV5[0], &pV5[1], &pV5[2], &pV5[3]);

    pV3[2] = 0;
    pV5[4] = 0;
}
