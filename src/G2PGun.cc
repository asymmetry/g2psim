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
#include "G2PDrift.hh"

#include "G2PGun.hh"

//#define GUN_DEBUG 1

using namespace std;

const double kDEG = 3.14159265358979323846/180.0;

ClassImp(G2PGun);

G2PGun::G2PGun()
    :bIsInit(false), iSetting(1), bUseData(false), bUseField(false),
     fHRSAngle(5.767*kDEG), fHRSMomentum(2.251), fBeamEnergy(2.254),
     fBeamX_lab(0), fBeamY_lab(0), fBeamTh_lab(0), fBeamPh_lab(0),
     fBeamR(0.015), fReactZLow_lab(0), fReactZHigh_lab(0),
     fTargetThLow_tr(0), fTargetThHigh_tr(0), fTargetPhLow_tr(0),
     fTargetPhHigh_tr(0), fDeltaLow(0), fDeltaHigh(0),
     fSigmaPos_lab(0.001), fSigmaAng_lab(0.001), fSigmaAng_tr(0.001),
     fSigmaDelta(0.0002), fBPMZ_lab(-0.8235), fBPMPosRes(0.2e-3),
     fBPMAngRes(0.4e-3), iIndex(0), pFileName(NULL), pfGunSelector(NULL)
{
    fData.clear();
}

G2PGun::G2PGun(const char* dist)
    :bIsInit(false), iSetting(1), bUseData(false), bUseField(false),
     fHRSAngle(5.767*kDEG), fHRSMomentum(2.251), fBeamEnergy(2.254),
     fBeamX_lab(0), fBeamY_lab(0), fBeamTh_lab(0), fBeamPh_lab(0),
     fBeamR(0.015), fReactZLow_lab(0), fReactZHigh_lab(0),
     fTargetThLow_tr(0), fTargetThHigh_tr(0), fTargetPhLow_tr(0),
     fTargetPhHigh_tr(0), fDeltaLow(0), fDeltaHigh(0),
     fSigmaPos_lab(0.001), fSigmaAng_lab(0.001), fSigmaAng_tr(0.001),
     fSigmaDelta(0.0002),  fBPMZ_lab(-0.8235), fBPMPosRes(0.2e-3),
     fBPMAngRes(0.2e-3), iIndex(0), pFileName(NULL), pfGunSelector(NULL)
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
}

G2PGun::~G2PGun()
{
    // Nothing to do
}

void G2PGun::Init()
{
    bool noerror = true;
    SetGun();
    if (iSetting==10) bUseData = true;
    if (G2PDrift::HasField()) bUseField = true;
    if (bUseData) noerror = LoadData();
    bIsInit = noerror;
}

void G2PGun::GetFP(double* V5fp_tr)
{
    V5fp_tr[0] = fDataAtIndex.xf;
    V5fp_tr[1] = atan(fDataAtIndex.tf);
    V5fp_tr[2] = fDataAtIndex.yf;
    V5fp_tr[3] = atan(fDataAtIndex.pf);
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

bool G2PGun::ShootDelta(double* V5beam_lab, double* V5bpm_lab, double* V5tg_tr)
{
    V5beam_lab[0] = fBeamX_lab;
    V5beam_lab[1] = fBeamTh_lab;
    V5beam_lab[2] = fBeamY_lab;
    V5beam_lab[3] = fBeamPh_lab;
    V5beam_lab[4] = fReactZLow_lab;

    double Xtg_tr, Ytg_tr, Ztg_tr;
    
    HRSTransTCSNHCS::X_HCS2TCS(V5beam_lab[0], V5beam_lab[2], V5beam_lab[4], fHRSAngle, Xtg_tr, Ytg_tr, Ztg_tr);

    V5tg_tr[0] = Xtg_tr;
    V5tg_tr[1] = fTargetThLow_tr;
    V5tg_tr[2] = Ytg_tr;
    V5tg_tr[3] = fTargetPhLow_tr;
    V5tg_tr[4] = fDeltaLow;

    GetBPMValue(V5beam_lab, V5bpm_lab);

    if (bUseField)
        G2PDrift::Drift(V5tg_tr, fHRSMomentum, fHRSAngle, Ztg_tr, 0.0, 10.0, V5tg_tr);
    else
        HRSTransTCSNHCS::Project(V5tg_tr[0], V5tg_tr[2], Ztg_tr, -Ztg_tr, V5tg_tr[1], V5tg_tr[3]);

#ifdef GUN_DEBUG
    printf("G2PGun: %e\t%e\t%e\t%e\t%e\n", V5tg_tr[0], V5tg_tr[1], V5tg_tr[2], V5tg_tr[3], V5tg_tr[4]);
#endif

    return true;
}

bool G2PGun::ShootGaus(double* V5beam_lab, double* V5bpm_lab, double* V5tg_tr)
{
    V5beam_lab[0] = G2PRand::Gaus(fBeamX_lab, fSigmaPos_lab);
    V5beam_lab[1] = G2PRand::Gaus(fBeamTh_lab, fSigmaAng_lab);
    V5beam_lab[2] = G2PRand::Gaus(fBeamY_lab, fSigmaPos_lab);
    V5beam_lab[3] = G2PRand::Gaus(fBeamPh_lab, fSigmaAng_lab);
    V5beam_lab[4] = G2PRand::Gaus(fReactZLow_lab, fSigmaPos_lab);

    double Xtg_tr, Ytg_tr, Ztg_tr;
    
    HRSTransTCSNHCS::X_HCS2TCS(V5beam_lab[0], V5beam_lab[2], V5beam_lab[4], fHRSAngle, Xtg_tr, Ytg_tr, Ztg_tr);
   
    V5tg_tr[0] = Xtg_tr;
    V5tg_tr[1] = G2PRand::Gaus(fTargetThLow_tr, fSigmaAng_tr);
    V5tg_tr[2] = Ytg_tr;
    V5tg_tr[3] = G2PRand::Gaus(fTargetPhLow_tr, fSigmaAng_tr);
    V5tg_tr[4] = G2PRand::Gaus(fDeltaLow, fSigmaDelta);

    GetBPMValue(V5beam_lab, V5bpm_lab);

    if (bUseField)
        G2PDrift::Drift(V5tg_tr, fHRSMomentum, fHRSAngle, Ztg_tr, 0.0, 10.0, V5tg_tr);
    else
        HRSTransTCSNHCS::Project(V5tg_tr[0], V5tg_tr[2], Ztg_tr, -Ztg_tr, V5tg_tr[1], V5tg_tr[3]);

#ifdef GUN_DEBUG
    printf("G2PGun: %e\t%e\t%e\t%e\t%e\n", V5tg_tr[0], V5tg_tr[1], V5tg_tr[2], V5tg_tr[3], V5tg_tr[4]);
#endif

    return true;
}

bool G2PGun::ShootFlat(double* V5beam_lab, double* V5bpm_lab, double* V5tg_tr)
{
    double Xbeam_lab, Ybeam_lab;
    do {
        Xbeam_lab = G2PRand::Uniform(-fBeamR, fBeamR);
        Ybeam_lab = G2PRand::Uniform(-fBeamR, fBeamR);
    } while (Xbeam_lab*Xbeam_lab+Ybeam_lab*Ybeam_lab>fBeamR*fBeamR);
    
    double Zreact_lab = G2PRand::Uniform(fReactZLow_lab, fReactZHigh_lab);

    V5beam_lab[0] = Xbeam_lab + fBeamX_lab;
    V5beam_lab[1] = fBeamTh_lab;
    V5beam_lab[2] = Ybeam_lab + fBeamY_lab;
    V5beam_lab[3] = fBeamPh_lab;
    V5beam_lab[4] = Zreact_lab;

    double Xtg_tr, Ytg_tr, Ztg_tr;
    
    HRSTransTCSNHCS::X_HCS2TCS(V5beam_lab[0], V5beam_lab[2], V5beam_lab[4], fHRSAngle, Xtg_tr, Ytg_tr, Ztg_tr);

    V5tg_tr[0] = Xtg_tr;
    V5tg_tr[1] = G2PRand::Uniform(fTargetThLow_tr, fTargetThHigh_tr);
    V5tg_tr[2] = Ytg_tr;
    V5tg_tr[3] = G2PRand::Uniform(fTargetPhLow_tr, fTargetPhHigh_tr);
    V5tg_tr[4] = G2PRand::Uniform(fDeltaLow, fDeltaHigh);

    GetBPMValue(V5beam_lab, V5bpm_lab);

    if (bUseField)
        G2PDrift::Drift(V5tg_tr, fHRSMomentum, fHRSAngle, Ztg_tr, 0.0, 10.0, V5tg_tr);
    else
        HRSTransTCSNHCS::Project(V5tg_tr[0], V5tg_tr[2], Ztg_tr, -Ztg_tr, V5tg_tr[1], V5tg_tr[3]);

#ifdef GUN_DEBUG
    printf("G2PGun: %e\t%e\t%e\t%e\t%e\n", V5tg_tr[0], V5tg_tr[1], V5tg_tr[2], V5tg_tr[3], V5tg_tr[4]);
#endif

    return true;
}

bool G2PGun::ShootTest(double* V5beam_lab, double* V5bpm_lab, double* V5tg_tr)
{
    int selector = G2PRand::Integer(17);
    
    V5beam_lab[0] = G2PRand::Gaus(fBeamX_lab, fSigmaPos_lab);
    V5beam_lab[1] = fBeamTh_lab;
    V5beam_lab[2] = G2PRand::Gaus(fBeamY_lab, fSigmaPos_lab);
    V5beam_lab[3] = fBeamPh_lab;
    V5beam_lab[4] = G2PRand::Gaus(fReactZLow_lab, fSigmaPos_lab);

    double Xtg_tr, Ytg_tr, Ztg_tr;
    
    HRSTransTCSNHCS::X_HCS2TCS(V5beam_lab[0], V5beam_lab[2], V5beam_lab[4], fHRSAngle, Xtg_tr, Ytg_tr, Ztg_tr);
    
    double Thetatg_tr = 0.0;
    double Phitg_tr = 0.0;
    switch (selector) {
    case 0:
    case 1:
    case 2:
    case 3:
    case 4:
        Thetatg_tr = G2PRand::Gaus(0, fSigmaAng_tr);
        Phitg_tr = G2PRand::Gaus(0, fSigmaAng_tr);
        break;
    case 5:
        Thetatg_tr = G2PRand::Gaus(0.02, fSigmaAng_tr);
        Phitg_tr = G2PRand::Gaus(0, fSigmaAng_tr);
        break;
    case 6:
        Thetatg_tr = G2PRand::Gaus(-0.02, fSigmaAng_tr);
        Phitg_tr = G2PRand::Gaus(0, fSigmaAng_tr);
        break;
    case 7:
        Thetatg_tr = G2PRand::Gaus(0, fSigmaAng_tr);
        Phitg_tr = G2PRand::Gaus(0.01, fSigmaAng_tr);
        break;
    case 8:
    case 9:
    case 10:
    case 11:
    case 12:
        Thetatg_tr = G2PRand::Gaus(0.02, fSigmaAng_tr);
        Phitg_tr = G2PRand::Gaus(0.01, fSigmaAng_tr);
        break;
    case 13:
        Thetatg_tr = G2PRand::Gaus(-0.02, fSigmaAng_tr);
        Phitg_tr = G2PRand::Gaus(0.01, fSigmaAng_tr);
        break;
    case 14:
        Thetatg_tr = G2PRand::Gaus(0, fSigmaAng_tr);
        Phitg_tr = G2PRand::Gaus(-0.01, fSigmaAng_tr);
        break;
    case 15:
        Thetatg_tr = G2PRand::Gaus(0.02, fSigmaAng_tr);
        Phitg_tr = G2PRand::Gaus(-0.01, fSigmaAng_tr);
        break;
    case 16:
        Thetatg_tr = G2PRand::Gaus(-0.02, fSigmaAng_tr);
        Phitg_tr = G2PRand::Gaus(-0.01, fSigmaAng_tr);
        break;
    }

    V5tg_tr[0] = Xtg_tr;
    V5tg_tr[1] = Thetatg_tr;
    V5tg_tr[2] = Ytg_tr;
    V5tg_tr[3] = Phitg_tr;
    V5tg_tr[4] = G2PRand::Gaus(0, fSigmaDelta);

    GetBPMValue(V5beam_lab, V5bpm_lab);

    if (bUseField)
        G2PDrift::Drift(V5tg_tr, fHRSMomentum, fHRSAngle, Ztg_tr, 0.0, 10.0, V5tg_tr);
    else
        HRSTransTCSNHCS::Project(V5tg_tr[0], V5tg_tr[2], Ztg_tr, -Ztg_tr, V5tg_tr[1], V5tg_tr[3]);  

#ifdef GUN_DEBUG
    printf("G2PGun: %e\t%e\t%e\t%e\t%e\n", V5tg_tr[0], V5tg_tr[1], V5tg_tr[2], V5tg_tr[3], V5tg_tr[4]);
#endif

    return true;
}

bool G2PGun::ShootSieve(double* V5beam_lab, double* V5bpm_lab, double* V5tg_tr)
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
    
    V5beam_lab[0] = G2PRand::Gaus(fBeamX_lab, fSigmaPos_lab);
    V5beam_lab[1] = fBeamTh_lab;
    V5beam_lab[2] = G2PRand::Gaus(fBeamY_lab, fSigmaPos_lab);
    V5beam_lab[3] = fBeamPh_lab;
    V5beam_lab[4] = G2PRand::Uniform(fReactZLow_lab, fReactZHigh_lab);

    double Xtg_tr, Ytg_tr, Ztg_tr;
    
    HRSTransTCSNHCS::X_HCS2TCS(V5beam_lab[0], V5beam_lab[2], V5beam_lab[4], fHRSAngle, Xtg_tr, Ytg_tr, Ztg_tr);

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

    double Thetatg_tr = G2PRand::Gaus(centheta, fSigmaAng_lab);
    double Phitg_tr = G2PRand::Gaus(cenphi, fSigmaAng_lab);

    // Calculate delta based on angle
    double V3pd_lab[3];
    HRSTransTCSNHCS::X_TCS2HCS(V3pd_tr[0], V3pd_tr[1], V3pd_tr[2], fHRSAngle, V3pd_lab[0], V3pd_lab[1], V3pd_lab[2]);

    double cosscatangle = V3pd_lab[2]/(sqrt(V3pd_lab[0]*V3pd_lab[0]+V3pd_lab[1]*V3pd_lab[1]+V3pd_lab[2]*V3pd_lab[2]));

    double scatmom = (targetmass*fBeamEnergy)/(targetmass+fBeamEnergy-fBeamEnergy*cosscatangle);

    double Delta = scatmom/fHRSMomentum-1-energyloss/fHRSMomentum;

    V5tg_tr[0] = Xtg_tr;
    V5tg_tr[1] = Thetatg_tr;
    V5tg_tr[2] = Ytg_tr;
    V5tg_tr[3] = Phitg_tr;
    V5tg_tr[4] = Delta;

    GetBPMValue(V5beam_lab, V5bpm_lab);

    if (bUseField)
        G2PDrift::Drift(V5tg_tr, fHRSMomentum, fHRSAngle, Ztg_tr, 0.0, 10.0, V5tg_tr);
    else
        HRSTransTCSNHCS::Project(V5tg_tr[0], V5tg_tr[2], Ztg_tr, -Ztg_tr, V5tg_tr[1], V5tg_tr[3]);

#ifdef GUN_DEBUG
    printf("G2PGun: %e\t%e\t%e\t%e\t%e\n", V5tg_tr[0], V5tg_tr[1], V5tg_tr[2], V5tg_tr[3], V5tg_tr[4]);
#endif

    return true;
}

bool G2PGun::ShootData(double* V5beam_lab, double* V5bpm_lab, double* V5tg_tr)
{
    bool noerror = true;

    if (iIndex>=(int)fData.size()) return false;
    int index = fData[iIndex].ind;
    V5bpm_lab[0] = fData[iIndex].xb;
    V5bpm_lab[1] = fData[iIndex].tb;
    V5bpm_lab[2] = fData[iIndex].yb;
    V5bpm_lab[3] = fData[iIndex].pb;
    V5bpm_lab[4] = fData[iIndex].zb;
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
    HRSTransTCSNHCS::X_HCS2TCS(V5bpm_lab[0], V5bpm_lab[2], V5bpm_lab[4], fHRSAngle, Xtg_tr, Ytg_tr, Ztg_tr);

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
    HRSTransTCSNHCS::X_TCS2HCS(V3pd_tr[0], V3pd_tr[1], V3pd_tr[2], fHRSAngle, V3pd_lab[0], V3pd_lab[1], V3pd_lab[2]);

    double cosscatangle = V3pd_lab[2]/(sqrt(V3pd_lab[0]*V3pd_lab[0]+V3pd_lab[1]*V3pd_lab[1]+V3pd_lab[2]*V3pd_lab[2]));

    double scatmom = (targetmass*fBeamEnergy)/(targetmass+fBeamEnergy-fBeamEnergy*cosscatangle);

    double Delta = scatmom/fHRSMomentum-1-energyloss/fHRSMomentum;

    V5tg_tr[0] = Xtg_tr;
    V5tg_tr[1] = Thetatg_tr;
    V5tg_tr[2] = Ytg_tr;
    V5tg_tr[3] = Phitg_tr;
    V5tg_tr[4] = Delta;

    if (bUseField)
        G2PDrift::Drift(V5tg_tr, fHRSMomentum, fHRSAngle, Ztg_tr, 0.0, 10.0, V5tg_tr);
    else
        HRSTransTCSNHCS::Project(V5tg_tr[0], V5tg_tr[2], Ztg_tr, -Ztg_tr, V5tg_tr[1], V5tg_tr[3]);

#ifdef GUN_DEBUG
    printf("G2PGun: %e\t%e\t%e\t%e\t%e\n", V5tg_tr[0], V5tg_tr[1], V5tg_tr[2], V5tg_tr[3], V5tg_tr[4]);
#endif

    iIndex++;
    return noerror;
}

void G2PGun::GetBPMValue(const double* V5beam_lab, double* V5bpm_lab)
{
    double x[3] = { V5beam_lab[0], V5beam_lab[2], V5beam_lab[4] };
    double p[3] = { fBeamEnergy*sin(V5beam_lab[1])*cos(V5beam_lab[3]),
                    fBeamEnergy*sin(V5beam_lab[1])*sin(V5beam_lab[3]),
                    fBeamEnergy*cos(V5beam_lab[1]) };

    if (bUseField) {
        double z_lab = x[2];
        G2PDrift::Drift(x, p, fBPMZ_lab, 10.0, x, p);
        double theta_tr = G2PRand::Gaus(atan(p[0]/p[2]), fBPMAngRes);
        double phi_tr = G2PRand::Gaus(atan(p[1]/p[2]), fBPMAngRes);
        p[0] = p[2]*tan(theta_tr);
        p[1] = p[2]*tan(phi_tr);
        double normF = fBeamEnergy/sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
        p[0] *= normF;
        p[1] *= normF;
        p[2] *= normF;
        x[0] = G2PRand::Gaus(x[0], fBPMPosRes);
        x[1] = G2PRand::Gaus(x[1], fBPMPosRes);
        G2PDrift::Drift(x, p, z_lab, 10.0, x, p);
        V5bpm_lab[0] = x[0];
        V5bpm_lab[1] = atan(p[1]/p[0]);
        V5bpm_lab[2] = x[1];
        V5bpm_lab[3] = acos(p[2]/fBeamEnergy);
        V5bpm_lab[4] = x[2];
    }
    else {
        double theta_tr = atan(p[0]/p[2]);
        double phi_tr = atan(p[1]/p[2]);
        double z_lab = x[2];
        HRSTransTCSNHCS::Project(x[0], x[1], x[2], fBPMZ_lab-x[2], theta_tr, phi_tr);
        theta_tr = G2PRand::Gaus(theta_tr, fBPMAngRes);
        phi_tr = G2PRand::Gaus(phi_tr, fBPMAngRes);
        x[0] = G2PRand::Gaus(x[0], fBPMPosRes);
        x[1] = G2PRand::Gaus(x[1], fBPMPosRes);
        HRSTransTCSNHCS::Project(x[0], x[1], x[2], z_lab-x[2], theta_tr, phi_tr);
        p[0] = p[2]*tan(theta_tr);
        p[1] = p[2]*tan(phi_tr);
        double normF = fBeamEnergy/sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
        p[0] *= normF;
        p[1] *= normF;
        p[2] *= normF;
        printf("G2PGun: %e\t%e\n", theta_tr, phi_tr);
        V5bpm_lab[0] = x[0];
        V5bpm_lab[1] = acos(p[2]/fBeamEnergy);
        V5bpm_lab[2] = x[1];
        V5bpm_lab[3] = atan(p[1]/p[0]);
        V5bpm_lab[4] = x[2];
    }

#ifdef GUN_DEBUG
    printf("G2PGun: %e\t%e\t%e\t%e\t%e\n", V5bpm_lab[0], V5bpm_lab[1], V5bpm_lab[2], V5bpm_lab[3], V5bpm_lab[4]);
#endif
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
