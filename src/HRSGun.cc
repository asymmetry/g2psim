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

#include <cstdio>
#include <cstdlib>
#include <string>
#include <map>

#include "TROOT.h"
#include "TMath.h"

#include "HRSRand.hh"
#include "HRSTransTCSNHCS.hh"

#include "HRSGun.hh"

//#define GUN_DEBUG 1

using namespace std;
using namespace Transform;

const double deg = TMath::Pi()/180.0;

HRSGun::HRSGun()
    :pIsInit(false), pUseData(false), pSetting(1), pTargetX_lab(0),
     pTargetY_lab(0), pTargetZLow_lab(0), pTargetZHigh_lab(0),
     pTargetR_lab(0.015), pTargetThetaLow_tr(0), pTargetThetaHigh_tr(0),
     pTargetPhiLow_tr(0), pTargetPhiHigh_tr(0), pDeltaLow(0), pDeltaHigh(0),
     pPosRes(0.0001), pAngleRes(0.01), pDeltaRes(0.03), pHRSAngle(5.767*deg),
     pGunSelector(NULL), pFilePtr(NULL), pFileName(NULL), pRand(NULL)
{
    memset(pV3bpm_lab, 0, sizeof(pV3bpm_lab));
    memset(pV5tg_tr, 0, sizeof(pV5tg_tr));
    memset(pV5fp_tr, 0, sizeof(pV5fp_tr));

    SetGun(1);
    pSetting=1;
}

HRSGun::HRSGun(const char* dist)
    :pIsInit(false), pUseData(false), pSetting(1), pTargetX_lab(0),
     pTargetY_lab(0), pTargetZLow_lab(0), pTargetZHigh_lab(0),
     pTargetR_lab(0.015), pTargetThetaLow_tr(0), pTargetThetaHigh_tr(0),
     pTargetPhiLow_tr(0), pTargetPhiHigh_tr(0), pDeltaLow(0), pDeltaHigh(0),
     pPosRes(0.0001), pAngleRes(0.01), pDeltaRes(0.03), pHRSAngle(5.767*deg),
     pGunSelector(NULL), pFilePtr(NULL), pFileName(NULL), pRand(NULL)
{
    memset(pV3bpm_lab, 0, sizeof(pV3bpm_lab));
    memset(pV5tg_tr, 0, sizeof(pV5tg_tr));
    memset(pV5fp_tr, 0, sizeof(pV5fp_tr));

    map<string, int> dist_map;
    dist_map["delta"] = 1;
    dist_map["gaus"] = 2;
    dist_map["flat"] = 3;
    dist_map["sieve"] = 4;
    dist_map["data"] = 5;

    if(dist_map[dist] == 0){
        printf("You are using an unknown gun setting!\n");
        SetGun(1);
        pSetting=1;
    }
    else{
        SetGun(dist_map[dist]);
        pSetting=dist_map[dist];
    }

    if (pSetting==5) pUseData = true;
}

HRSGun::~HRSGun()
{
}

void HRSGun::Init()
{
    pIsInit = true;
    if (pFileName!=NULL) {
        if ((pFilePtr=fopen(pFileName, "r"))==NULL) {
            pIsInit = false;
        }
    }
}

void HRSGun::End()
{
    if (pFilePtr!=NULL) {
        fclose(pFilePtr);
    }
    pIsInit = false;
}

void HRSGun::SetGun(int dist)
{
    switch (dist) {
    case 1:
        pGunSelector = &HRSGun::ShootDelta;
        break;
    case 2:
        pGunSelector = &HRSGun::ShootGaus;
        break;
    case 3:
        pGunSelector = &HRSGun::ShootFlat;
        break;
    case 4:
        pGunSelector = &HRSGun::ShootSieve;
        break;
    case 5:
        pGunSelector = &HRSGun::ShootData;
        break;
    }
}

void HRSGun::ShootDelta(double *pV3, double *pV5)
{
    double Xtg_lab = pTargetX_lab;
    double Ytg_lab = pTargetY_lab;
    double Ztg_lab = pTargetZLow_lab;

    pV3[0] = Xtg_lab;
    pV3[1] = Ytg_lab;
    pV3[2] = Ztg_lab;

    double Xtg_tr, Ytg_tr, Ztg_tr;
    
    X_HCS2TCS(Xtg_lab, Ytg_lab, Ztg_lab, pHRSAngle, Xtg_tr, Ytg_tr, Ztg_tr);

    double Thetatg_tr = pTargetThetaLow_tr;
    double Phitg_tr = pTargetPhiLow_tr;

    Project(Xtg_tr, Ytg_tr, Ztg_tr, -Ztg_tr, Thetatg_tr, Phitg_tr);

    double Deltatg = pDeltaLow;
    
    pV5[0] = Xtg_tr;
    pV5[1] = Thetatg_tr;
    pV5[2] = Ytg_tr;
    pV5[3] = Phitg_tr;
    pV5[4] = Deltatg;

#ifdef GUN_DEBUG
    printf("%e\t%e\t%e\t%e\t%e\n", pV5[0], pV5[1], pV5[2], pV5[3], pV5[4]);
#endif
}

void HRSGun::ShootGaus(double *pV3, double *pV5)
{
    double Xtg_lab = pRand->Gaus(pTargetX_lab, pPosRes);
    double Ytg_lab = pRand->Gaus(pTargetY_lab, pPosRes);
    double Ztg_lab = pRand->Gaus(pTargetZLow_lab, pPosRes);

    pV3[0] = Xtg_lab;
    pV3[1] = Ytg_lab;
    pV3[2] = Ztg_lab;

    double Xtg_tr, Ytg_tr, Ztg_tr;
    
    X_HCS2TCS(Xtg_lab, Ytg_lab, Ztg_lab, pHRSAngle, Xtg_tr, Ytg_tr, Ztg_tr);

    double Thetatg_tr = pRand->Gaus(pTargetThetaLow_tr, pAngleRes);
    double Phitg_tr = pRand->Gaus(pTargetPhiLow_tr, pAngleRes);

    Project(Xtg_tr, Ytg_tr, Ztg_tr, -Ztg_tr, Thetatg_tr, Phitg_tr);

    double Deltatg = pRand->Gaus(pDeltaLow, pDeltaRes);
    
    pV5[0] = Xtg_tr;
    pV5[1] = Thetatg_tr;
    pV5[2] = Ytg_tr;
    pV5[3] = Phitg_tr;
    pV5[4] = Deltatg;

#ifdef GUN_DEBUG
    printf("%e\t%e\t%e\t%e\t%e\n", pV5[0], pV5[1], pV5[2], pV5[3], pV5[4]);
#endif
}

void HRSGun::ShootFlat(double *pV3, double *pV5)
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

    double Thetatg_tr = pRand->Uniform(pTargetThetaLow_tr, pTargetThetaHigh_tr);
    double Phitg_tr = pRand->Uniform(pTargetPhiLow_tr, pTargetPhiHigh_tr);

    Project(Xtg_tr, Ytg_tr, Ztg_tr, -Ztg_tr, Thetatg_tr, Phitg_tr);

    double Deltatg = pRand->Uniform(pDeltaLow, pDeltaHigh);

    pV5[0] = Xtg_tr;
    pV5[1] = Thetatg_tr;
    pV5[2] = Ytg_tr;
    pV5[3] = Phitg_tr;
    pV5[4] = Deltatg;

#ifdef GUN_DEBUG
    printf("%e\t%e\t%e\t%e\t%e\n", pV5[0], pV5[1], pV5[2], pV5[3], pV5[4]);
#endif
}

void HRSGun::ShootSieve(double *pV3, double *pV5)
{
    int selector = pRand->Integer(17);
    double Xtg_lab = 0, Ytg_lab = 0;
    double Ztg_lab = 0;

    pV3[0] = Xtg_lab;
    pV3[1] = Ytg_lab;
    pV3[2] = Ztg_lab;

    double Xtg_tr, Ytg_tr, Ztg_tr;
    
    X_HCS2TCS(Xtg_lab, Ytg_lab, Ztg_lab, pHRSAngle, Xtg_tr, Ytg_tr, Ztg_tr);
    
    double Thetatg_tr = 0, Phitg_tr = 0;
    switch (selector) {
    case 0:
    case 1:
    case 2:
    case 3:
    case 4:
        Thetatg_tr = pRand->Gaus(0, 0.0005);
        Phitg_tr = pRand->Gaus(0, 0.0005);
        break;
    case 5:
        Thetatg_tr = pRand->Gaus(0.02, 0.0005);
        Phitg_tr = pRand->Gaus(0, 0.0005);
        break;
    case 6:
        Thetatg_tr = pRand->Gaus(-0.02, 0.0005);
        Phitg_tr = pRand->Gaus(0, 0.0005);
        break;
    case 7:
        Thetatg_tr = pRand->Gaus(0, 0.0005);
        Phitg_tr = pRand->Gaus(0.01, 0.0005);
        break;
    case 8:
    case 9:
    case 10:
    case 11:
    case 12:
        Thetatg_tr = pRand->Gaus(0.02, 0.0005);
        Phitg_tr = pRand->Gaus(0.01, 0.0005);
        break;
    case 13:
        Thetatg_tr = pRand->Gaus(-0.02, 0.0005);
        Phitg_tr = pRand->Gaus(0.01, 0.0005);
        break;
    case 14:
        Thetatg_tr = pRand->Gaus(0, 0.0005);
        Phitg_tr = pRand->Gaus(-0.01, 0.0005);
        break;
    case 15:
        Thetatg_tr = pRand->Gaus(0.02, 0.0005);
        Phitg_tr = pRand->Gaus(-0.01, 0.0005);
        break;
    case 16:
        Thetatg_tr = pRand->Gaus(-0.02, 0.0005);
        Phitg_tr = pRand->Gaus(-0.01, 0.0005);
        break;
    }
    
    Project(Xtg_tr, Ytg_tr, Ztg_tr, -Ztg_tr, Thetatg_tr, Phitg_tr);

    //double Deltatg = pRand->Gaus(0, pDeltaRes);
    double Deltatg = 0.0;

    pV5[0] = Xtg_tr;
    pV5[1] = Thetatg_tr;
    pV5[2] = Ytg_tr;
    pV5[3] = Phitg_tr;
    pV5[4] = Deltatg;

#ifdef GUN_DEBUG
    printf("%e\t%e\t%e\t%e\t%e\n", pV5[0], pV5[1], pV5[2], pV5[3], pV5[4]);
#endif
}

void HRSGun::ShootData(double *pV3, double *pV5)
{
    int temp;

    if (!feof(pFilePtr)) {
        // index, bpm_x, bpm_y, thetatg_tr, ytg_tr, phitg_tr, deltatg
        fscanf(pFilePtr, "%d%lf%lf%lf%lf%lf%lf", &temp, &pV3[0], &pV3[1], &pV5[0], &pV5[1], &pV5[2], &pV5[3]);
    }
    else{
        pIsInit = false;
    }

    pV3[2] = 0;
    pV5[4] = 0;

#ifdef GUN_DEBUG
    printf("%e\t%e\t%e\t%e\t%e\n", pV5[0], pV5[1], pV5[2], pV5[3], pV5[4]);
#endif
}
