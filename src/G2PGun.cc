// This file defines a class G2PGun.
// This class is used in G2PSim class as particle gun.
// It has 6 standard gun: ShootDelta(), ShootGaus(), ShootFlat(), ShootTest(),
//+ShootSieve() and ShootData().
// The active gun is chosen during initializing.
// G2PSim class will call Shoot() to get kinematic variables. It is a virtual
//+function so you can rewrite it by inheriting this class.
//
// History:
//   Jan 2013, C. Gu, First public version.
//   Jan 2013, C. Gu, Add ShootSieve() method.
//

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <vector>
#include <cmath>

#include "TROOT.h"
#include "TObject.h"

#include "G2PRand.hh"
#include "HRSTransTCSNHCS.hh"

#include "G2PGun.hh"

//#define GUN_DEBUG 1

using namespace std;
using namespace HRSTransTCSNHCS;

const double kDEG = 3.14159265358979323846/180.0;

ClassImp(G2PGun);

G2PGun::G2PGun()
    :bIsInit(false), iSetting(1), bUseData(false), fHRSAngle(5.767*kDEG),
     fHRSMomentum(2.251), fBeamEnergy(2.254), fTargetX_lab(0),
     fTargetY_lab(0), fTargetZLow_lab(0), fTargetZHigh_lab(0),
     fTargetR_lab(0.015), fTargetThLow_tr(0), fTargetThHigh_tr(0),
     fTargetPhLow_tr(0), fTargetPhHigh_tr(0), fDeltaLow(0), fDeltaHigh(0),
     fPosRes(0.001), fAngleRes(0.001), fDeltaRes(0.0002), iIndex(0),
     pDrift(NULL), bFieldOn(false), pFileName(NULL), pfGunSelector(NULL)
{
    fData.clear();
}

G2PGun::G2PGun(const char* dist)
    :bIsInit(false), iSetting(1), bUseData(false), fHRSAngle(5.767*kDEG),
     fHRSMomentum(2.251), fBeamEnergy(2.254), fTargetX_lab(0),
     fTargetY_lab(0), fTargetZLow_lab(0), fTargetZHigh_lab(0),
     fTargetR_lab(0.015), fTargetThLow_tr(0), fTargetThHigh_tr(0),
     fTargetPhLow_tr(0), fTargetPhHigh_tr(0), fDeltaLow(0), fDeltaHigh(0),
     fPosRes(0.001), fAngleRes(0.001), fDeltaRes(0.0002), iIndex(0),
     pDrift(NULL), bFieldOn(false), pFileName(NULL), pfGunSelector(NULL)
{
    fData.clear();

    map<string, int> dist_map;
    dist_map["delta"] = 1;
    dist_map["gaus"] = 2;
    dist_map["flat"] = 3;
    dist_map["test"] = 4;
    dist_map["sieve"] = 5;
    dist_map["data"] = 10;

    if(dist_map[dist] == 0){
        printf("Unknown gun setting, set to delta distribution ...\n");
        iSetting=1;
    }
    else{
        iSetting=dist_map[dist];
    }

    if (iSetting==10) bUseData = true;
}

G2PGun::~G2PGun()
{
    // Nothing to do
}

void G2PGun::Init()
{
    bool noerror = true;
    SetGun();
    if (bUseData) noerror = LoadData();
    bIsInit = noerror;
}

void G2PGun::GetFP(double* V5fp)
{
    V5fp[0] = fDataAtIndex.xf;
    V5fp[1] = atan(fDataAtIndex.tf);
    V5fp[2] = fDataAtIndex.yf;
    V5fp[3] = atan(fDataAtIndex.pf);
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
    case 10:
        pfGunSelector = &G2PGun::ShootData;
    }
}

bool G2PGun::ShootDelta(double* V5bpm, double* V5tg)
{
    V5bpm[0] = fTargetX_lab;
    V5bpm[2] = fTargetY_lab;
    V5bpm[4] = fTargetZLow_lab;

    double Xtg_tr, Ytg_tr, Ztg_tr;
    
    X_HCS2TCS(V5bpm[0], V5bpm[2], V5bpm[4], fHRSAngle, Xtg_tr, Ytg_tr, Ztg_tr);

    V5tg[0] = Xtg_tr;
    V5tg[1] = fTargetThLow_tr;
    V5tg[2] = Ytg_tr;
    V5tg[3] = fTargetPhLow_tr;
    V5tg[4] = fDeltaLow;

    if (bFieldOn)
        pDrift->Drift(V5tg, Ztg_tr, 0.0, 10.0, V5tg);
    else
        Project(V5tg[0], V5tg[2], Ztg_tr, -Ztg_tr, V5tg[1], V5tg[3]);

#ifdef GUN_DEBUG
    printf("G2PGun: %e\t%e\t%e\t%e\t%e\n", V5tg[0], V5tg[1], V5tg[2], V5tg[3], V5tg[4]);
#endif

    return true;
}

bool G2PGun::ShootGaus(double* V5bpm, double* V5tg)
{
    V5bpm[0] = G2PRand::Gaus(fTargetX_lab, fPosRes);
    V5bpm[2] = G2PRand::Gaus(fTargetY_lab, fPosRes);
    V5bpm[4] = G2PRand::Gaus(fTargetZLow_lab, fPosRes);

    double Xtg_tr, Ytg_tr, Ztg_tr;
    
    X_HCS2TCS(V5bpm[0], V5bpm[2], V5bpm[4], fHRSAngle, Xtg_tr, Ytg_tr, Ztg_tr);
   
    V5tg[0] = Xtg_tr;
    V5tg[1] = G2PRand::Gaus(fTargetThLow_tr, fAngleRes);
    V5tg[2] = Ytg_tr;
    V5tg[3] = G2PRand::Gaus(fTargetPhLow_tr, fAngleRes);
    V5tg[4] = G2PRand::Gaus(fDeltaLow, fDeltaRes);

    if (bFieldOn)
        pDrift->Drift(V5tg, Ztg_tr, 0.0, 10.0, V5tg);
    else
        Project(V5tg[0], V5tg[2], Ztg_tr, -Ztg_tr, V5tg[1], V5tg[3]);

#ifdef GUN_DEBUG
    printf("G2PGun: %e\t%e\t%e\t%e\t%e\n", V5tg[0], V5tg[1], V5tg[2], V5tg[3], V5tg[4]);
#endif

    return true;
}

bool G2PGun::ShootFlat(double* V5bpm, double* V5tg)
{
    double Xtg_lab, Ytg_lab;
    do {
        Xtg_lab = G2PRand::Uniform(-fTargetR_lab, fTargetR_lab);
        Ytg_lab = G2PRand::Uniform(-fTargetR_lab, fTargetR_lab);
    } while (Xtg_lab*Xtg_lab+Ytg_lab*Ytg_lab>fTargetR_lab*fTargetR_lab);
    
    double Ztg_lab = G2PRand::Uniform(fTargetZLow_lab, fTargetZHigh_lab);

    V5bpm[0] = Xtg_lab + fTargetX_lab;
    V5bpm[2] = Ytg_lab + fTargetY_lab;
    V5bpm[4] = Ztg_lab;

    double Xtg_tr, Ytg_tr, Ztg_tr;
    
    X_HCS2TCS(V5bpm[0], V5bpm[2], V5bpm[4], fHRSAngle, Xtg_tr, Ytg_tr, Ztg_tr);

    V5tg[0] = Xtg_tr;
    V5tg[1] = G2PRand::Uniform(fTargetThLow_tr, fTargetThHigh_tr);
    V5tg[2] = Ytg_tr;
    V5tg[3] = G2PRand::Uniform(fTargetPhLow_tr, fTargetPhHigh_tr);
    V5tg[4] = G2PRand::Uniform(fDeltaLow, fDeltaHigh);

    if (bFieldOn)
        pDrift->Drift(V5tg, Ztg_tr, 0.0, 10.0, V5tg);
    else
        Project(V5tg[0], V5tg[2], Ztg_tr, -Ztg_tr, V5tg[1], V5tg[3]);   

#ifdef GUN_DEBUG
    printf("G2PGun: %e\t%e\t%e\t%e\t%e\n", V5tg[0], V5tg[1], V5tg[2], V5tg[3], V5tg[4]);
#endif

    return true;
}

bool G2PGun::ShootTest(double* V5bpm, double* V5tg)
{
    int selector = G2PRand::Integer(17);
    
    V5bpm[0] = G2PRand::Gaus(fTargetX_lab, fPosRes);
    V5bpm[2] = G2PRand::Gaus(fTargetY_lab, fPosRes);
    V5bpm[4] = G2PRand::Gaus(fTargetZLow_lab, fPosRes);

    double Xtg_tr, Ytg_tr, Ztg_tr;
    
    X_HCS2TCS(V5bpm[0], V5bpm[2], V5bpm[4], fHRSAngle, Xtg_tr, Ytg_tr, Ztg_tr);
    
    double Thetatg_tr = 0.0;
    double Phitg_tr = 0.0;
    switch (selector) {
    case 0:
    case 1:
    case 2:
    case 3:
    case 4:
        Thetatg_tr = G2PRand::Gaus(0, fAngleRes);
        Phitg_tr = G2PRand::Gaus(0, fAngleRes);
        break;
    case 5:
        Thetatg_tr = G2PRand::Gaus(0.02, fAngleRes);
        Phitg_tr = G2PRand::Gaus(0, fAngleRes);
        break;
    case 6:
        Thetatg_tr = G2PRand::Gaus(-0.02, fAngleRes);
        Phitg_tr = G2PRand::Gaus(0, fAngleRes);
        break;
    case 7:
        Thetatg_tr = G2PRand::Gaus(0, fAngleRes);
        Phitg_tr = G2PRand::Gaus(0.01, fAngleRes);
        break;
    case 8:
    case 9:
    case 10:
    case 11:
    case 12:
        Thetatg_tr = G2PRand::Gaus(0.02, fAngleRes);
        Phitg_tr = G2PRand::Gaus(0.01, fAngleRes);
        break;
    case 13:
        Thetatg_tr = G2PRand::Gaus(-0.02, fAngleRes);
        Phitg_tr = G2PRand::Gaus(0.01, fAngleRes);
        break;
    case 14:
        Thetatg_tr = G2PRand::Gaus(0, fAngleRes);
        Phitg_tr = G2PRand::Gaus(-0.01, fAngleRes);
        break;
    case 15:
        Thetatg_tr = G2PRand::Gaus(0.02, fAngleRes);
        Phitg_tr = G2PRand::Gaus(-0.01, fAngleRes);
        break;
    case 16:
        Thetatg_tr = G2PRand::Gaus(-0.02, fAngleRes);
        Phitg_tr = G2PRand::Gaus(-0.01, fAngleRes);
        break;
    }

    V5tg[0] = Xtg_tr;
    V5tg[1] = Thetatg_tr;
    V5tg[2] = Ytg_tr;
    V5tg[3] = Phitg_tr;
    V5tg[4] = G2PRand::Gaus(0, fDeltaRes);

    if (bFieldOn)
        pDrift->Drift(V5tg, Ztg_tr, 0.0, 10.0, V5tg);
    else
        Project(V5tg[0], V5tg[2], Ztg_tr, -Ztg_tr, V5tg[1], V5tg[3]);  

#ifdef GUN_DEBUG
    printf("G2PGun: %e\t%e\t%e\t%e\t%e\n", V5tg[0], V5tg[1], V5tg[2], V5tg[3], V5tg[4]);
#endif

    return true;
}

bool G2PGun::ShootSieve(double* V5bpm, double* V5tg)
{
    const int nsieverow = 7;

    const double sievex[] = { -3*13.3096e-3, -2*13.3096e-3, -1*13.3096e-3, 0.0, 1*13.3096e-3, 2*13.3096e-3, 3*13.3096e-3 };
    const double sievey[] = { 3*6.1214e-3, 2*6.1214e-3, 1*6.1214e-3, 0.0, -1*4.7752e-3, -2*4.7752e-3, -3*4.7752e-3 };
    const double sievez = 799.60e-3;

    const double sieveoffx = 0.0;
    const double sieveoffy = 0.0;

    const double dsmallhole = 1.3970e-3;
    const double dlargehole = 2.6924e-3;
    
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

    const double arearatio = 3.714380165;
    const double thr = (arearatio-1.0)/(arearatio*2+47.0);

    int selector;
    int col, row;
    do {
        selector = G2PRand::Integer(49);
        double temp = G2PRand::Uniform();
        if (temp<thr) selector = 24;
        else if (temp<2*thr) selector = 15;
        col = selector/(nsieverow);
        row = selector%(nsieverow);
    } while (sieveon[row][col]==0);
    
    V5bpm[0] = G2PRand::Gaus(fTargetX_lab, fPosRes);
    V5bpm[2] = G2PRand::Gaus(fTargetY_lab, fPosRes);
    V5bpm[4] = G2PRand::Uniform(fTargetZLow_lab, fTargetZHigh_lab);

    double Xtg_tr, Ytg_tr, Ztg_tr;
    
    X_HCS2TCS(V5bpm[0], V5bpm[2], V5bpm[4], fHRSAngle, Xtg_tr, Ytg_tr, Ztg_tr);

    double V3sieve_tr[3];
    double V3pd_tr[3];

    if ((selector==15)||(selector==24)) {
        V3sieve_tr[0] = sieveoffx+G2PRand::Uniform(sievex[row]-dlargehole/2,sievex[row]+dlargehole/2);
        V3sieve_tr[1] = sieveoffy+G2PRand::Uniform(sievey[col]-dlargehole/2,sievey[col]+dlargehole/2);
    }
    else{
        V3sieve_tr[0] = sieveoffx+G2PRand::Uniform(sievex[row]-dsmallhole/2,sievex[row]+dsmallhole/2);
        V3sieve_tr[1] = sieveoffy+G2PRand::Uniform(sievey[col]-dsmallhole/2,sievey[col]+dsmallhole/2);
    }
    V3sieve_tr[2] = sievez;

    V3pd_tr[0] = V3sieve_tr[0]-Xtg_tr;
    V3pd_tr[1] = V3sieve_tr[1]-Ytg_tr;
    V3pd_tr[2] = V3sieve_tr[2]-Ztg_tr;

    double centheta = atan(V3pd_tr[0]/V3pd_tr[2]);
    double cenphi = atan(V3pd_tr[1]/V3pd_tr[2]);

    double Thetatg_tr = G2PRand::Gaus(centheta, fAngleRes);
    double Phitg_tr = G2PRand::Gaus(cenphi, fAngleRes);

    // Calculate delta based on angle
    double V3pd_lab[3];
    X_TCS2HCS(V3pd_tr[0], V3pd_tr[1], V3pd_tr[2], fHRSAngle, V3pd_lab[0], V3pd_lab[1], V3pd_lab[2]);

    double cosscatangle = V3pd_lab[2]/(sqrt(V3pd_lab[0]*V3pd_lab[0]+V3pd_lab[1]*V3pd_lab[1]+V3pd_lab[2]*V3pd_lab[2]));

    double scatmom = (targetmass*fBeamEnergy)/(targetmass+fBeamEnergy-fBeamEnergy*cosscatangle);

    double Delta = scatmom/fHRSMomentum-1-energyloss/fHRSMomentum;

    V5tg[0] = Xtg_tr;
    V5tg[1] = Thetatg_tr;
    V5tg[2] = Ytg_tr;
    V5tg[3] = Phitg_tr;
    V5tg[4] = Delta;

    if (bFieldOn)
        pDrift->Drift(V5tg, Ztg_tr, 0.0, 10.0, V5tg);
    else
        Project(V5tg[0], V5tg[2], Ztg_tr, -Ztg_tr, V5tg[1], V5tg[3]);  

#ifdef GUN_DEBUG
    printf("G2PGun: %e\t%e\t%e\t%e\t%e\n", V5tg[0], V5tg[1], V5tg[2], V5tg[3], V5tg[4]);
#endif

    return true;
}

bool G2PGun::ShootData(double* V5bpm, double* V5tg)
{
    bool noerror = true;

    if (iIndex>=(int)fData.size()) return false;
    int index = fData[iIndex].ind;
    V5bpm[0] = fData[iIndex].xb;
    V5bpm[1] = fData[iIndex].tb;
    V5bpm[2] = fData[iIndex].yb;
    V5bpm[3] = fData[iIndex].pb;
    V5bpm[4] = fData[iIndex].zb;
    fDataAtIndex = fData[iIndex];

    const int nsieverow = 7;

    const double sievex[] = { -3*13.3096e-3, -2*13.3096e-3, -1*13.3096e-3, 0.0, 1*13.3096e-3, 2*13.3096e-3, 3*13.3096e-3 };
    const double sievey[] = { 3*6.1214e-3, 2*6.1214e-3, 1*6.1214e-3, 0.0, -1*4.7752e-3, -2*4.7752e-3, -3*4.7752e-3 };
    const double sievez = 799.60e-3;

    const double sieveoffx = 0.0;
    const double sieveoffy = 0.0;
    
    const double targetmass = 12.0107*0.931494028;
    const double energyloss = 1.009711e-3+0.501422e-3;

    int col = index/(nsieverow);
    int row = index%(nsieverow);

    double Xtg_tr, Ytg_tr, Ztg_tr;
    X_HCS2TCS(V5bpm[0], V5bpm[2], V5bpm[4], fHRSAngle, Xtg_tr, Ytg_tr, Ztg_tr);

    double V3sieve_tr[3];
    double V3pd_tr[3];

    V3sieve_tr[0] = sieveoffx+sievex[row];
    V3sieve_tr[1] = sieveoffy+sievey[col];
    V3sieve_tr[2] = sievez;

    V3pd_tr[0] = V3sieve_tr[0]-Xtg_tr;
    V3pd_tr[1] = V3sieve_tr[1]-Ytg_tr;
    V3pd_tr[2] = V3sieve_tr[2]-Ztg_tr;

    double Thetatg_tr = atan(V3pd_tr[0]/V3pd_tr[2]);
    double Phitg_tr = atan(V3pd_tr[1]/V3pd_tr[2]);

    // Calculate delta based on angle
    double V3pd_lab[3];
    X_TCS2HCS(V3pd_tr[0], V3pd_tr[1], V3pd_tr[2], fHRSAngle, V3pd_lab[0], V3pd_lab[1], V3pd_lab[2]);

    double cosscatangle = V3pd_lab[2]/(sqrt(V3pd_lab[0]*V3pd_lab[0]+V3pd_lab[1]*V3pd_lab[1]+V3pd_lab[2]*V3pd_lab[2]));

    double scatmom = (targetmass*fBeamEnergy)/(targetmass+fBeamEnergy-fBeamEnergy*cosscatangle);

    double Delta = scatmom/fHRSMomentum-1-energyloss/fHRSMomentum;

    V5tg[0] = Xtg_tr;
    V5tg[1] = Thetatg_tr;
    V5tg[2] = Ytg_tr;
    V5tg[3] = Phitg_tr;
    V5tg[4] = Delta;

    if (bFieldOn)
        pDrift->Drift(V5tg, Ztg_tr, 0.0, 10.0, V5tg);
    else
        Project(V5tg[0], V5tg[2], Ztg_tr, -Ztg_tr, V5tg[1], V5tg[3]);

#ifdef GUN_DEBUG
    printf("G2PGun: %e\t%e\t%e\t%e\t%e\n", V5tg[0], V5tg[1], V5tg[2], V5tg[3], V5tg[4]);
#endif

    iIndex++;
    return noerror;
}

bool G2PGun::LoadData()
{
    FILE *fp;

    if ((fp=fopen(pFileName, "r"))==NULL) return false;

    sData temp;
    fscanf(fp, "%d%lf%lf%lf%lf%lf%lf%lf%lf%lf", &temp.ind, &temp.xb, &temp.tb, &temp.yb, &temp.pb, &temp.zb, &temp.xf, &temp.tf, &temp.yf, &temp.pf);
    while (!feof(fp)) {
        fData.push_back(temp);
        fscanf(fp, "%d%lf%lf%lf%lf%lf%lf%lf%lf%lf", &temp.ind, &temp.xb, &temp.tb, &temp.yb, &temp.pb, &temp.zb, &temp.xf, &temp.tf, &temp.yf, &temp.pf);
    }

    fclose(fp);

    if (!fData.empty()) return true;
    else return false;
}
