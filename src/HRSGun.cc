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

const double cDeg = TMath::Pi()/180.0;

HRSGun::HRSGun()
    :bIsInit(false), iSetting(1), bUseData(false), fHRSAngle(5.767*cDeg),
     fTargetX_lab(0), fTargetY_lab(0), fTargetZLow_lab(0),
     fTargetZHigh_lab(0), fTargetR_lab(0.015), fTargetThLow_tr(0),
     fTargetThHigh_tr(0), fTargetPhLow_tr(0), fTargetPhHigh_tr(0),
     fDeltaLow(0), fDeltaHigh(0), fPosRes(0.0001), fAngleRes(0.001),
     fDeltaRes(0.01), pFilePtr(NULL), pFileName(NULL), pRand(NULL),
     pfGunSelector(NULL)
{
    // Nothing to do
}

HRSGun::HRSGun(const char* dist)
    :bIsInit(false), iSetting(1), bUseData(false), fHRSAngle(5.767*cDeg),
     fTargetX_lab(0), fTargetY_lab(0), fTargetZLow_lab(0),
     fTargetZHigh_lab(0), fTargetR_lab(0.015), fTargetThLow_tr(0),
     fTargetThHigh_tr(0), fTargetPhLow_tr(0), fTargetPhHigh_tr(0),
     fDeltaLow(0), fDeltaHigh(0), fPosRes(0.0001), fAngleRes(0.001),
     fDeltaRes(0.01), pFilePtr(NULL), pFileName(NULL), pRand(NULL),
     pfGunSelector(NULL)
{
    map<string, int> dist_map;
    dist_map["delta"] = 1;
    dist_map["gaus"] = 2;
    dist_map["flat"] = 3;
    dist_map["test"] = 4;
    dist_map["sieve"] = 5;
    dist_map["data"] = 6;

    if(dist_map[dist] == 0){
        printf("Unknown gun setting, set to delta distribution ...\n");
        iSetting=1;
    }
    else{
        iSetting=dist_map[dist];
    }

    if (iSetting==6) bUseData = true;
}

HRSGun::~HRSGun()
{
    // Nothing to do
}

void HRSGun::Init()
{
    bool noerror = true;
    SetGun(iSetting);
    if (bUseData) {
        if ((pFilePtr=fopen(pFileName, "r"))==NULL) noerror = false;
    }
    bIsInit = noerror;
}

void HRSGun::End()
{
    if (bUseData&&(pFilePtr!=NULL)) {
        fclose(pFilePtr);
    }
    bIsInit = false;
}

void HRSGun::SetGun(int setting)
{
    switch (setting) {
    case 1:
        pfGunSelector = &HRSGun::ShootDelta;
        break;
    case 2:
        pfGunSelector = &HRSGun::ShootGaus;
        break;
    case 3:
        pfGunSelector = &HRSGun::ShootFlat;
        break;
    case 4:
    case 5:
        pfGunSelector = &HRSGun::ShootSieve;
        break;
    case 6:
        pfGunSelector = &HRSGun::ShootData;
        break;
    }
}

bool HRSGun::ShootDelta(double *V3bpm, double *V5tg)
{
    V3bpm[0] = fTargetX_lab;
    V3bpm[1] = fTargetY_lab;
    V3bpm[2] = fTargetZLow_lab;

    double Xtg_tr, Ytg_tr, Ztg_tr;
    
    X_HCS2TCS(V3bpm[0], V3bpm[1], V3bpm[2], fHRSAngle, Xtg_tr, Ytg_tr, Ztg_tr);

    V5tg[1] = fTargetThLow_tr;
    V5tg[3] = fTargetPhLow_tr;

    Project(Xtg_tr, Ytg_tr, Ztg_tr, -Ztg_tr, V5tg[1], V5tg[3]);
    
    V5tg[0] = Xtg_tr;
    V5tg[2] = Ytg_tr;
    V5tg[4] = fDeltaLow;

#ifdef GUN_DEBUG
    printf("%e\t%e\t%e\t%e\t%e\n", V5tg[0], V5tg[1], V5tg[2], V5tg[3], V5tg[4]);
#endif

    return true;
}

bool HRSGun::ShootGaus(double *V3bpm, double *V5tg)
{
    V3bpm[0] = pRand->Gaus(fTargetX_lab, fPosRes);
    V3bpm[1] = pRand->Gaus(fTargetY_lab, fPosRes);
    V3bpm[2] = pRand->Gaus(fTargetZLow_lab, fPosRes);

    double Xtg_tr, Ytg_tr, Ztg_tr;
    
    X_HCS2TCS(V3bpm[0], V3bpm[1], V3bpm[2], fHRSAngle, Xtg_tr, Ytg_tr, Ztg_tr);

    V5tg[1] = pRand->Gaus(fTargetThLow_tr, fAngleRes);
    V5tg[3] = pRand->Gaus(fTargetPhLow_tr, fAngleRes);

    Project(Xtg_tr, Ytg_tr, Ztg_tr, -Ztg_tr, V5tg[1], V5tg[3]);
    
    V5tg[0] = Xtg_tr;
    V5tg[2] = Ytg_tr;
    V5tg[4] = pRand->Gaus(fDeltaLow, fDeltaRes);

#ifdef GUN_DEBUG
    printf("%e\t%e\t%e\t%e\t%e\n", V5tg[0], V5tg[1], V5tg[2], V5tg[3], V5tg[4]);
#endif

    return true;
}

bool HRSGun::ShootFlat(double *V3bpm, double *V5tg)
{
    double Xtg_lab, Ytg_lab;
    do {
        Xtg_lab = pRand->Uniform(-fTargetR_lab, fTargetR_lab);
        Ytg_lab = pRand->Uniform(-fTargetR_lab, fTargetR_lab);
    } while (Xtg_lab*Xtg_lab+Ytg_lab*Ytg_lab>fTargetR_lab*fTargetR_lab);
    
    double Ztg_lab = pRand->Uniform(fTargetZLow_lab, fTargetZHigh_lab);

    V3bpm[0] = Xtg_lab + fTargetX_lab;
    V3bpm[1] = Ytg_lab + fTargetY_lab;
    V3bpm[2] = Ztg_lab;

    double Xtg_tr, Ytg_tr, Ztg_tr;
    
    X_HCS2TCS(V3bpm[0], V3bpm[1], V3bpm[2], fHRSAngle, Xtg_tr, Ytg_tr, Ztg_tr);

    V5tg[1] = pRand->Uniform(fTargetThLow_tr, fTargetThHigh_tr);
    V5tg[3] = pRand->Uniform(fTargetPhLow_tr, fTargetPhHigh_tr);

    Project(Xtg_tr, Ytg_tr, Ztg_tr, -Ztg_tr, V5tg[1], V5tg[3]);

    V5tg[0] = Xtg_tr;
    V5tg[2] = Ytg_tr;
    V5tg[4] = pRand->Uniform(fDeltaLow, fDeltaHigh);

#ifdef GUN_DEBUG
    printf("%e\t%e\t%e\t%e\t%e\n", V5tg[0], V5tg[1], V5tg[2], V5tg[3], V5tg[4]);
#endif

    return true;
}

bool HRSGun::ShootSieve(double *V3bpm, double *V5tg)
{
    int selector = pRand->Integer(17);
    
    V3bpm[0] = pRand->Gaus(fTargetX_lab, fPosRes);
    V3bpm[1] = pRand->Gaus(fTargetY_lab, fPosRes);
    V3bpm[2] = pRand->Gaus(fTargetZLow_lab, fPosRes);

    double Xtg_tr, Ytg_tr, Ztg_tr;
    
    X_HCS2TCS(V3bpm[0], V3bpm[1], V3bpm[2], fHRSAngle, Xtg_tr, Ytg_tr, Ztg_tr);
    
    double Thetatg_tr = 0.0;
    double Phitg_tr = 0.0;
    switch (selector) {
    case 0:
    case 1:
    case 2:
    case 3:
    case 4:
        Thetatg_tr = pRand->Gaus(0, fAngleRes);
        Phitg_tr = pRand->Gaus(0, fAngleRes);
        break;
    case 5:
        Thetatg_tr = pRand->Gaus(0.02, fAngleRes);
        Phitg_tr = pRand->Gaus(0, fAngleRes);
        break;
    case 6:
        Thetatg_tr = pRand->Gaus(-0.02, fAngleRes);
        Phitg_tr = pRand->Gaus(0, fAngleRes);
        break;
    case 7:
        Thetatg_tr = pRand->Gaus(0, fAngleRes);
        Phitg_tr = pRand->Gaus(0.01, fAngleRes);
        break;
    case 8:
    case 9:
    case 10:
    case 11:
    case 12:
        Thetatg_tr = pRand->Gaus(0.02, fAngleRes);
        Phitg_tr = pRand->Gaus(0.01, fAngleRes);
        break;
    case 13:
        Thetatg_tr = pRand->Gaus(-0.02, fAngleRes);
        Phitg_tr = pRand->Gaus(0.01, fAngleRes);
        break;
    case 14:
        Thetatg_tr = pRand->Gaus(0, fAngleRes);
        Phitg_tr = pRand->Gaus(-0.01, fAngleRes);
        break;
    case 15:
        Thetatg_tr = pRand->Gaus(0.02, fAngleRes);
        Phitg_tr = pRand->Gaus(-0.01, fAngleRes);
        break;
    case 16:
        Thetatg_tr = pRand->Gaus(-0.02, fAngleRes);
        Phitg_tr = pRand->Gaus(-0.01, fAngleRes);
        break;
    }
    
    Project(Xtg_tr, Ytg_tr, Ztg_tr, -Ztg_tr, Thetatg_tr, Phitg_tr);

    V5tg[0] = Xtg_tr;
    V5tg[1] = Thetatg_tr;
    V5tg[2] = Ytg_tr;
    V5tg[3] = Phitg_tr;
    V5tg[4] = pRand->Gaus(0, fDeltaRes);

#ifdef GUN_DEBUG
    printf("%e\t%e\t%e\t%e\t%e\n", V5tg[0], V5tg[1], V5tg[2], V5tg[3], V5tg[4]);
#endif

    return true;
}

bool HRSGun::ShootData(double *V3bpm, double *V5tg)
{
    int temp;
    bool noerror = true;

    if (!feof(pFilePtr)) {
        // index, bpm_x, bpm_y, thetatg_tr, ytg_tr, phitg_tr, deltatg
        fscanf(pFilePtr, "%d%lf%lf%lf%lf%lf%lf", &temp, &V3bpm[0], &V3bpm[1], &V5tg[0], &V5tg[1], &V5tg[2], &V5tg[3]);
    }
    else{
        noerror = false;
    }

    V3bpm[2] = 0;
    V5tg[4] = 0;

#ifdef GUN_DEBUG
    printf("%e\t%e\t%e\t%e\t%e\n", V5tg[0], V5tg[1], V5tg[2], V5tg[3], V5tg[4]);
#endif

    return noerror;
}
