// This file defined a class G2PGun.
// This class will be used in G2PSim class as particle gun.
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
#include "TObject.h"
#include "TMath.h"

#include "G2PRand.hh"
#include "HRSTransTCSNHCS.hh"

#include "G2PGun.hh"

//#define GUN_DEBUG 1

using namespace std;
using namespace Transform;

const double cDeg = TMath::Pi()/180.0;

ClassImp(G2PGun);

G2PGun::G2PGun()
    :bIsInit(false), iSetting(1), bUseData(false), fHRSAngle(5.767*cDeg),
     fHRSMomentum(2.251), fBeamEnergy(2.254), fTargetX_lab(0), fTargetY_lab(0),
     fTargetZLow_lab(0), fTargetZHigh_lab(0), fTargetR_lab(0.015),
     fTargetThLow_tr(0), fTargetThHigh_tr(0), fTargetPhLow_tr(0),
     fTargetPhHigh_tr(0), fDeltaLow(0), fDeltaHigh(0), fPosRes(0.001),
     fAngleRes(0.001), fDeltaRes(0.0002), pFilePtr(NULL), pFileName(NULL),
     pRand(NULL), pfGunSelector(NULL)
{
    // Nothing to do
}

G2PGun::G2PGun(const char *dist)
    :bIsInit(false), iSetting(1), bUseData(false), fHRSAngle(5.767*cDeg),
     fHRSMomentum(2.251), fBeamEnergy(2.254), fTargetX_lab(0), fTargetY_lab(0),
     fTargetZLow_lab(0), fTargetZHigh_lab(0), fTargetR_lab(0.015),
     fTargetThLow_tr(0), fTargetThHigh_tr(0), fTargetPhLow_tr(0),
     fTargetPhHigh_tr(0), fDeltaLow(0), fDeltaHigh(0), fPosRes(0.001),
     fAngleRes(0.001), fDeltaRes(0.0002), pFilePtr(NULL), pFileName(NULL),
     pRand(NULL), pfGunSelector(NULL)
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

G2PGun::~G2PGun()
{
    // Nothing to do
}

void G2PGun::Init()
{
    bool noerror = true;
    SetGun();
    if (bUseData) {
        if ((pFilePtr=fopen(pFileName, "r"))==NULL) noerror = false;
    }
    bIsInit = noerror;
}

void G2PGun::End()
{
    if (bUseData&&(pFilePtr!=NULL)) {
        fclose(pFilePtr);
    }
    bIsInit = false;
}

void G2PGun::SetGun()
{
    switch (iSetting) {
    case 1:
        pfGunSelector = &G2PGun::ShootDelta;
        break;
    case 2:
        pfGunSelector = &G2PGun::ShootGaus;
        break;
    case 3:
        pfGunSelector = &G2PGun::ShootFlat;
        break;
    case 4:
        pfGunSelector = &G2PGun::ShootTest;
        break;
    case 5:
        pfGunSelector = &G2PGun::ShootSieve;
        break;
    case 6:
        pfGunSelector = &G2PGun::ShootData;
        break;
    }
}

bool G2PGun::ShootDelta(double *V3bpm, double *V5tg)
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
    printf("G2PGun: %e\t%e\t%e\t%e\t%e\n", V5tg[0], V5tg[1], V5tg[2], V5tg[3], V5tg[4]);
#endif

    return true;
}

bool G2PGun::ShootGaus(double *V3bpm, double *V5tg)
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
    printf("G2PGun: %e\t%e\t%e\t%e\t%e\n", V5tg[0], V5tg[1], V5tg[2], V5tg[3], V5tg[4]);
#endif

    return true;
}

bool G2PGun::ShootFlat(double *V3bpm, double *V5tg)
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
    printf("G2PGun: %e\t%e\t%e\t%e\t%e\n", V5tg[0], V5tg[1], V5tg[2], V5tg[3], V5tg[4]);
#endif

    return true;
}

bool G2PGun::ShootTest(double *V3bpm, double *V5tg)
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
    printf("G2PGun: %e\t%e\t%e\t%e\t%e\n", V5tg[0], V5tg[1], V5tg[2], V5tg[3], V5tg[4]);
#endif

    return true;
}

bool G2PGun::ShootSieve(double *V3bpm, double *V5tg)
{
    const int nsieverow = 7;
    const int arearatio = 4.0;
    const double thr = (arearatio-1.0)/(arearatio*2+47.0);

    const double sievex[] = { -3*13.3096e-3, -2*13.3096e-3, -1*13.3096e-3, 0.0, 1*13.3096e-3, 2*13.3096e-3, 3*13.3096e-3 };
    const double sievey[] = { 3*6.1214e-3, 2*6.1214e-3, 1*6.1214e-3, 0.0, -1*4.7752e-3, -2*4.7752e-3, -3*4.7752e-3 };
    const double sievez = 799.60e-3;

    const double sieveoffx = 0.0;
    const double sieveoffy = 0.0;

    const double targetmass = 12.0107*0.931494028;
    const double energyloss = 1.009711e-3+0.501422e-3;

    const int sieveon[7][7] ={
        { 0, 0, 0, 0, 1, 1, 1 } ,
        { 0, 0, 1, 1, 1, 1, 1 } ,
        { 0, 0, 1, 1, 1, 1, 1 } ,
        { 0, 0, 1, 1, 1, 1, 1 } ,
        { 0, 0, 1, 1, 1, 1, 1 } ,
        { 0, 0, 0, 1, 1, 1, 0 } ,
        { 0, 0, 0, 0, 0, 0, 0 }
    };

    // const int sieveon[7][7] ={
    //     { 1, 1, 1, 1, 1, 1, 1 } ,
    //     { 1, 1, 1, 1, 1, 1, 1 } ,
    //     { 1, 1, 1, 1, 1, 1, 1 } ,
    //     { 1, 1, 1, 1, 1, 1, 1 } ,
    //     { 1, 1, 1, 1, 1, 1, 1 } ,
    //     { 1, 1, 1, 1, 1, 1, 1 } ,
    //     { 1, 1, 1, 1, 1, 1, 1 }
    // };

    int selector;
    int col, row;
    do {
        selector = pRand->Integer(49);
        double temp = pRand->Uniform();
        if (temp<thr) selector = 24;
        else if (temp<2*thr) selector = 15;
        col = selector/(nsieverow);
        row = selector%(nsieverow);
    } while (sieveon[row][col]==0);
    
    V3bpm[0] = pRand->Gaus(fTargetX_lab, fPosRes);
    V3bpm[1] = pRand->Gaus(fTargetY_lab, fPosRes);
    V3bpm[2] = pRand->Uniform(fTargetZLow_lab, fTargetZHigh_lab);

    double Xtg_tr, Ytg_tr, Ztg_tr;
    
    X_HCS2TCS(V3bpm[0], V3bpm[1], V3bpm[2], fHRSAngle, Xtg_tr, Ytg_tr, Ztg_tr);
    
    double Thetatg_tr = 0.0;
    double Phitg_tr = 0.0;

    double V3sieve_tr[3];
    double V3pd_tr[3];

    V3sieve_tr[0] = sieveoffx + sievex[row];
    V3sieve_tr[1] = sieveoffy + sievey[col];
    V3sieve_tr[2] = sievez;

    V3pd_tr[0] = V3sieve_tr[0]-Xtg_tr;
    V3pd_tr[1] = V3sieve_tr[1]-Ytg_tr;
    V3pd_tr[2] = V3sieve_tr[2]-Ztg_tr;

    Thetatg_tr = pRand->Gaus(V3pd_tr[0]/V3pd_tr[2],fAngleRes);
    Phitg_tr = pRand->Gaus(V3pd_tr[1]/V3pd_tr[2],fAngleRes);
    
    Project(Xtg_tr, Ytg_tr, Ztg_tr, -Ztg_tr, Thetatg_tr, Phitg_tr);

    // Calculate delta based on angle
    double V3pd_lab[3];
    X_TCS2HCS(V3pd_tr[0], V3pd_tr[1], V3pd_tr[2], fHRSAngle, V3pd_lab[0], V3pd_lab[1], V3pd_lab[2]);

    double cosscatangle = V3pd_lab[2]/(sqrt(V3pd_lab[0]*V3pd_lab[0]+V3pd_lab[1]*V3pd_lab[1]+V3pd_lab[2]*V3pd_lab[2]));

    double scatmom = (targetmass*fBeamEnergy)/(targetmass+fBeamEnergy-fBeamEnergy*cosscatangle);
    double censcatmom = (targetmass*fBeamEnergy)/(targetmass+fBeamEnergy-fBeamEnergy*cos(fHRSAngle));

    double dpkinoffset = (scatmom-censcatmom)/fHRSMomentum;

    double dpkin = censcatmom/fHRSMomentum-1-energyloss/fHRSMomentum;

    V5tg[0] = Xtg_tr;
    V5tg[1] = Thetatg_tr;
    V5tg[2] = Ytg_tr;
    V5tg[3] = Phitg_tr;
    V5tg[4] = dpkin + dpkinoffset;

#ifdef GUN_DEBUG
    printf("G2PGun: %e\t%e\t%e\t%e\t%e\n", V5tg[0], V5tg[1], V5tg[2], V5tg[3], V5tg[4]);
#endif

    return true;
}

bool G2PGun::ShootData(double *V3bpm, double *V5fp)
{
    int temp;
    bool noerror = true;

    if (!feof(pFilePtr)) {
        // index, bpm_x, bpm_y, thetatg_tr, ytg_tr, phitg_tr, deltatg
        fscanf(pFilePtr, "%d%lf%lf%lf%lf%lf%lf", &temp, &V3bpm[0], &V3bpm[1], &V5fp[0], &V5fp[1], &V5fp[2], &V5fp[3]);
    }
    else{
        noerror = false;
    }

    V3bpm[2] = pRand->Uniform(fTargetZLow_lab, fTargetZHigh_lab);
    V5fp[4] = 0;

#ifdef GUN_DEBUG
    printf("G2PGun: %e\t%e\t%e\t%e\t%e\n", V5fp[0], V5fp[1], V5fp[2], V5fp[3], V5fp[4]);
#endif

    return noerror;
}
