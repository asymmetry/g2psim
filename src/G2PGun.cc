// This file defines a class G2PGun.
// This class is used in G2PSim class as particle gun.
// The active gun is chosen during initializing.
// G2PSim class will call Shoot() to get kinematic variables. It is a virtual
//+function so you can rewrite it by inheriting this class.
//
// History:
//   Jan 2013, C. Gu, First public version.
//   Jan 2013, C. Gu, Add ShootSieve() method.
//   Feb 2013, C. Gu, Rewrite the sieve and data parts.
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

static const double kDEG = 3.14159265358979323846/180.0;

ClassImp(G2PGun);

G2PGun::G2PGun()
    :bIsInit(false), iSetting(1), bUseData(false), bUseField(false),
     fHRSAngle(5.767*kDEG), fHRSMomentum(2.251), fBeamEnergy(2.254),
     fBeamX_lab(0), fBeamY_lab(0), fBeamTiltAngle(0),
     fBeamR(0.015), fReactZLow_lab(0), fReactZHigh_lab(0),
     fTargetThLow_tr(0), fTargetThHigh_tr(0), fTargetPhLow_tr(0),
     fTargetPhHigh_tr(0), fDeltaLow(0), fDeltaHigh(0),
     fSigmaPos_lab(0.001), fSigmaAng_lab(0.001), fSigmaAng_tr(0.001),
     fSigmaDelta(0.0002), fTargetMass(12.0107*0.931494028),
     fEnergyLoss(1.009711e-3+0.501422e-3), pFileName(NULL),
     pfGunSelector(NULL)
{
    fData.clear();
}

G2PGun::G2PGun(const char* dist)
    :bIsInit(false), iSetting(1), bUseData(false), bUseField(false),
     fHRSAngle(5.767*kDEG), fHRSMomentum(2.251), fBeamEnergy(2.254),
     fBeamX_lab(0), fBeamY_lab(0), fBeamTiltAngle(0),
     fBeamR(0.015), fReactZLow_lab(0), fReactZHigh_lab(0),
     fTargetThLow_tr(0), fTargetThHigh_tr(0), fTargetPhLow_tr(0),
     fTargetPhHigh_tr(0), fDeltaLow(0), fDeltaHigh(0),
     fSigmaPos_lab(0.001), fSigmaAng_lab(0.001), fSigmaAng_tr(0.001),
     fSigmaDelta(0.0002), fTargetMass(12.0107*0.931494028),
     fEnergyLoss(1.009711e-3+0.501422e-3), pFileName(NULL),
     pfGunSelector(NULL)
{
    fData.clear();

    map<string, int> dist_map;
    dist_map["gaus"] = 1;
    dist_map["flat"] = 2;
    dist_map["sieve"] = 3;
    dist_map["fastsieve"] = 4;
    dist_map["data"] = 10;
    dist_map["opticsdata"] = 11;

    if (dist_map[dist] == 0) {
        printf("Unknown gun setting, set to delta distribution ...\n");
        iSetting=1;
    }
    else {
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
    if (iSetting==10||iSetting==11) bUseData = true;
    if (iSetting==3||iSetting==4) SetSieve();
    if (G2PDrift::HasField()) {
        bUseField = true;
        if (iSetting==4) iSetting = 3;
    }
    if (bUseData) noerror = LoadData();
    SetGun();
    SetBeamTiltAngle();
    bIsInit = noerror;
}

void G2PGun::SetGun()
{
    switch (iSetting) {
    case 1:
        pfGunSelector = &G2PGun::ShootGaus;
        break;
    case 2:
        pfGunSelector = &G2PGun::ShootFlat;
        break;
    case 3:
        pfGunSelector = &G2PGun::ShootSieve;
        break;
    case 4:
        pfGunSelector = &G2PGun::ShootSieveFast;
        break;
    case 10:
        pfGunSelector = &G2PGun::ShootData;
        break;
    case 11:
        pfGunSelector = &G2PGun::ShootOpticsData;
    }
}

void G2PGun::SetBeamTiltAngle()
{
    if (bUseField) {
        if (fabs(G2PDrift::GetField()->GetRatio()-0.5)<1e-8) {
            if (fabs(fBeamEnergy-2.254)<0.2) fBeamTiltAngle = 3.31*kDEG;
            else if (fabs(fBeamEnergy-1.706)<0.2) fBeamTiltAngle = 4.03*kDEG;
            else if (fabs(fBeamEnergy-1.159)<0.2) fBeamTiltAngle = 5.97*kDEG;
            else fBeamTiltAngle = 0.0;
        }
        else if (fabs(G2PDrift::GetField()->GetRatio()-1.0)<1e-8) {
            if (fabs(fBeamEnergy-2.254)<0.2) fBeamTiltAngle = 0.0; 
            else if (fabs(fBeamEnergy-3.355)<0.2) fBeamTiltAngle = 0.0;
            else fBeamTiltAngle = 0.0;
        }
        else {
            fBeamTiltAngle = 0.0;
        }
    }
    else {
        fBeamTiltAngle = 0.0;
    }
}

void G2PGun::SetSieve()
{
    if (fHRSAngle>0) { // left arm
        const double kSIEVEX[7] = { -3*13.3096e-3, -2*13.3096e-3, -1*13.3096e-3, 0.0, 1*13.3096e-3, 2*13.3096e-3, 3*13.3096e-3 };
        const double kSIEVEY[7] = { 3*6.1214e-3, 2*6.1214e-3, 1*6.1214e-3, 0.0, -1*4.7752e-3, -2*4.7752e-3, -3*4.7752e-3 };
        const int kLARGERHOLE[2] = { 15, 24 };
        const int kSIEVEOPEN[49] = { 0, 0, 0, 0, 1, 1, 1,
                                     0, 0, 1, 1, 1, 1, 1,
                                     0, 0, 1, 1, 1, 1, 1,
                                     0, 0, 1, 1, 1, 1, 1,
                                     0, 0, 1, 1, 1, 1, 1,
                                     0, 0, 0, 1, 1, 1, 0,
                                     0, 0, 0, 0, 0, 0, 0 };

        fSieve.nRow = 7;
        fSieve.nCol = 7;
        fSieve.nLargerHole = 2;

        fSieve.fX.clear();
        for (int i = 0; i<fSieve.nRow; i++) fSieve.fX.push_back(kSIEVEX[i]);
        fSieve.fY.clear();
        for (int i = 0; i<fSieve.nCol; i++) fSieve.fX.push_back(kSIEVEY[i]);
        fSieve.fZ = 799.60e-3;
        fSieve.fXOffset = 0.0;
        fSieve.fYOffset = 0.0;

        fSieve.iLargerHole.clear();
        for (int i = 0; i<fSieve.nLargerHole; i++) fSieve.iLargerHole.push_back(kLARGERHOLE[i]);
        fSieve.bOpen.clear();
        for (int i = 0; i<fSieve.nRow*fSieve.nCol; i++) fSieve.bOpen.push_back((kSIEVEOPEN[i]==1)?true:false);

        fSieve.fDHole = 1.3970e-3;
        fSieve.fDLargerHole = 2.6924e-3;

        double ratio = fSieve.fDLargerHole*fSieve.fDLargerHole/(fSieve.fDHole*fSieve.fDHole);
        fSieve.fThreshold = (ratio-1)/((ratio-1)*fSieve.nLargerHole+fSieve.nRow*fSieve.nCol);
    }
    else { // right arm
        const double kSIEVEX[7] = { -3*13.3096e-3, -2*13.3096e-3, -1*13.3096e-3, 0.0, 1*13.3096e-3, 2*13.3096e-3, 3*13.3096e-3 };
        const double kSIEVEY[7] = { -3*6.1214e-3, -2*6.1214e-3, -1*6.1214e-3, 0.0, 1*4.7752e-3, 2*4.7752e-3, 3*4.7752e-3 };
        const int kLARGERHOLE[2] = { 15, 24 };
        const int kSIEVEOPEN[49] = { 0, 0, 0, 0, 1, 1, 1,
                                     0, 0, 1, 1, 1, 1, 1,
                                     0, 0, 1, 1, 1, 1, 1,
                                     0, 0, 1, 1, 1, 1, 1,
                                     0, 0, 1, 1, 1, 1, 1,
                                     0, 0, 0, 1, 1, 1, 0,
                                     0, 0, 0, 0, 0, 0, 0 };

        fSieve.nRow = 7;
        fSieve.nCol = 7;
        fSieve.nLargerHole = 2;

        fSieve.fX.clear();
        for (int i = 0; i<fSieve.nRow; i++) fSieve.fX.push_back(kSIEVEX[i]);
        fSieve.fY.clear();
        for (int i = 0; i<fSieve.nCol; i++) fSieve.fX.push_back(kSIEVEY[i]);
        fSieve.fZ = 799.46e-3;
        fSieve.fXOffset = 0.0;
        fSieve.fYOffset = 0.0;

        fSieve.iLargerHole.clear();
        for (int i = 0; i<fSieve.nLargerHole; i++) fSieve.iLargerHole.push_back(kLARGERHOLE[i]);
        fSieve.bOpen.clear();
        for (int i = 0; i<fSieve.nRow*fSieve.nCol; i++) fSieve.bOpen.push_back((kSIEVEOPEN[i]==1)?true:false);

        fSieve.fDHole = 1.3970e-3;
        fSieve.fDLargerHole = 2.6924e-3;

        double ratio = fSieve.fDLargerHole*fSieve.fDLargerHole/(fSieve.fDHole*fSieve.fDHole);
        fSieve.fThreshold = (ratio-1)/((ratio-1)*fSieve.nLargerHole+fSieve.nRow*fSieve.nCol);
    }
}

bool G2PGun::ShootGaus(double* V5beam_lab, double* V5react_tr, double* reserved)
{
    double X_lab = G2PRand::Gaus(fBeamX_lab, fSigmaPos_lab);
    double Y_lab = G2PRand::Gaus(fBeamY_lab, fSigmaPos_lab);
    double Z_lab = G2PRand::Gaus(fReactZLow_lab, fSigmaPos_lab);
    GetReactPoint(X_lab, Y_lab, Z_lab, V5beam_lab);

    double Xreact_tr, Yreact_tr, Zreact_tr;
    HRSTransTCSNHCS::X_HCS2TCS(V5beam_lab[0], V5beam_lab[2], V5beam_lab[4], fHRSAngle, Xreact_tr, Yreact_tr, Zreact_tr);
   
    V5react_tr[0] = Xreact_tr;
    V5react_tr[1] = G2PRand::Gaus(fTargetThLow_tr, fSigmaAng_tr);
    V5react_tr[2] = Yreact_tr;
    V5react_tr[3] = G2PRand::Gaus(fTargetPhLow_tr, fSigmaAng_tr);
    V5react_tr[4] = G2PRand::Gaus(fDeltaLow, fSigmaDelta);

#ifdef GUN_DEBUG
    printf("G2PGun: %e\t%e\t%e\t%e\t%e\n", V5react_tr[0], V5react_tr[1], V5react_tr[2], V5react_tr[3], V5react_tr[4]);
#endif

    return true;
}

bool G2PGun::ShootFlat(double* V5beam_lab, double* V5react_tr, double* reserved)
{
    double X_lab, Y_lab;
    if (fBeamR>1e-5) {
        do {
            X_lab = G2PRand::Uniform(-fBeamR, fBeamR);
            Y_lab = G2PRand::Uniform(-fBeamR, fBeamR);
        } while (X_lab*X_lab+Y_lab*Y_lab>fBeamR*fBeamR);
    }
    else {
        X_lab = 0.0;
        Y_lab = 0.0;
    }

    X_lab+=fBeamX_lab;
    Y_lab+=fBeamY_lab;
    double Z_lab = G2PRand::Uniform(fReactZLow_lab, fReactZHigh_lab);
    GetReactPoint(X_lab, Y_lab, Z_lab, V5beam_lab);

    double Xreact_tr, Yreact_tr, Zreact_tr;
    HRSTransTCSNHCS::X_HCS2TCS(V5beam_lab[0], V5beam_lab[2], V5beam_lab[4], fHRSAngle, Xreact_tr, Yreact_tr, Zreact_tr);

    V5react_tr[0] = Xreact_tr;
    V5react_tr[1] = G2PRand::Uniform(fTargetThLow_tr, fTargetThHigh_tr);
    V5react_tr[2] = Yreact_tr;
    V5react_tr[3] = G2PRand::Uniform(fTargetPhLow_tr, fTargetPhHigh_tr);
    V5react_tr[4] = G2PRand::Uniform(fDeltaLow, fDeltaHigh);

#ifdef GUN_DEBUG
    printf("G2PGun: %e\t%e\t%e\t%e\t%e\n", V5react_tr[0], V5react_tr[1], V5react_tr[2], V5react_tr[3], V5react_tr[4]);
#endif

    return true;
}

bool G2PGun::ShootSieve(double* V5beam_lab, double* V5react_tr, double* reserved)
{
    bool found = false;

    while (!found) {
        // generate flat distribution
        double X_lab, Y_lab;
        if (fBeamR>1e-5) {
            do {
                X_lab = G2PRand::Uniform(-fBeamR, fBeamR);
                Y_lab = G2PRand::Uniform(-fBeamR, fBeamR);
            } while (X_lab*X_lab+Y_lab*Y_lab>fBeamR*fBeamR);
        }
        else {
            X_lab = 0.0;
            Y_lab = 0.0;
        }

        X_lab+=fBeamX_lab;
        Y_lab+=fBeamY_lab;
        double Z_lab = G2PRand::Uniform(fReactZLow_lab, fReactZHigh_lab);
        GetReactPoint(X_lab, Y_lab, Z_lab, V5beam_lab);

        double Xreact_tr, Yreact_tr, Zreact_tr;
        HRSTransTCSNHCS::X_HCS2TCS(V5beam_lab[0], V5beam_lab[2], V5beam_lab[4], fHRSAngle, Xreact_tr, Yreact_tr, Zreact_tr);

        V5react_tr[0] = Xreact_tr;
        V5react_tr[1] = G2PRand::Uniform(fTargetThLow_tr, fTargetThHigh_tr);
        V5react_tr[2] = Yreact_tr;
        V5react_tr[3] = G2PRand::Uniform(fTargetPhLow_tr, fTargetPhHigh_tr);
        V5react_tr[4] = G2PRand::Uniform(fDeltaLow, fDeltaHigh);

        // drift to sieve slit
        double V5siv_tr[5];
        if (bUseField)
            G2PDrift::Drift(V5react_tr, fHRSMomentum, fHRSAngle, Zreact_tr, fSieve.fZ, 10.0, V5siv_tr);
        else
            HRSTransTCSNHCS::Project(V5react_tr[0], V5react_tr[2], Zreact_tr, fSieve.fZ-Zreact_tr, V5react_tr[1], V5react_tr[3], V5siv_tr[0], V5siv_tr[2], V5siv_tr[4]);

        // check if it passed the sieve slit
        for (int i = 0; i<fSieve.nRow*fSieve.nCol; i++) {
            double dhole = fSieve.fDHole;
            for (int j = 0; j<fSieve.nLargerHole; j++)
                if (i==fSieve.iLargerHole[j]) dhole = fSieve.fDLargerHole;
            int col = i/(fSieve.nRow);
            int row = i%(fSieve.nRow);
            if (fSieve.bOpen[i]) {
                double distance = V5siv_tr[0]-fSieve.fX[row]-fSieve.fXOffset*V5siv_tr[0]-fSieve.fX[row]-fSieve.fXOffset;
                distance += V5siv_tr[2]-fSieve.fY[col]-fSieve.fYOffset*V5siv_tr[2]-fSieve.fY[col]-fSieve.fYOffset;
                distance = sqrt(distance);
                if (distance<dhole/2.0) {
                    found = true;
                    break;
                }
            }
        }
    }

#ifdef GUN_DEBUG
    printf("G2PGun: %e\t%e\t%e\t%e\t%e\n", V5react_tr[0], V5react_tr[1], V5react_tr[2], V5react_tr[3], V5react_tr[4]);
#endif

    return true;
}

bool G2PGun::ShootSieveFast(double* V5beam_lab, double* V5react_tr, double* reserved)
{
    double dhole;

    int selector;
    do {
        dhole = fSieve.fDHole;
        selector = G2PRand::Integer(fSieve.nRow*fSieve.nCol);
        double temp = G2PRand::Uniform();
        for (int i = 0; i<fSieve.nLargerHole; i++)
            if (temp<i*fSieve.fThreshold) {
                selector = fSieve.iLargerHole[i];
                dhole = fSieve.fDLargerHole;
            }
    } while (fSieve.bOpen[selector]==0);

    int col = selector/(fSieve.nRow);
    int row = selector%(fSieve.nRow);

    double X_lab, Y_lab;
    if (fBeamR>1e-5) {
        do {
            X_lab = G2PRand::Uniform(-fBeamR, fBeamR);
            Y_lab = G2PRand::Uniform(-fBeamR, fBeamR);
        } while (X_lab*X_lab+Y_lab*Y_lab>fBeamR*fBeamR);
    }
    else {
        X_lab = 0.0;
        Y_lab = 0.0;
    }

    X_lab+=fBeamX_lab;
    Y_lab+=fBeamY_lab;
    double Z_lab = G2PRand::Uniform(fReactZLow_lab, fReactZHigh_lab);
    GetReactPoint(X_lab, Y_lab, Z_lab, V5beam_lab);
        
    double Xreact_tr, Yreact_tr, Zreact_tr;
    HRSTransTCSNHCS::X_HCS2TCS(V5beam_lab[0], V5beam_lab[2], V5beam_lab[4], fHRSAngle, Xreact_tr, Yreact_tr, Zreact_tr);

    double V3sieve_tr[3];
    double V3pd_tr[3];

    double Xsiv_tr, Ysiv_tr;
    do {
        Xsiv_tr = G2PRand::Uniform(-dhole/2, dhole/2);
        Ysiv_tr = G2PRand::Uniform(-dhole/2, dhole/2);
    } while (Xsiv_tr*Xsiv_tr+Ysiv_tr*Ysiv_tr>dhole*dhole/4.0);
    V3sieve_tr[0] = fSieve.fXOffset+fSieve.fX[row]+Xsiv_tr;
    V3sieve_tr[1] = fSieve.fYOffset+fSieve.fY[col]+Ysiv_tr;
    V3sieve_tr[2] = fSieve.fZ;

    V3pd_tr[0] = V3sieve_tr[0]-Xreact_tr;
    V3pd_tr[1] = V3sieve_tr[1]-Yreact_tr;
    V3pd_tr[2] = V3sieve_tr[2]-Zreact_tr;

    double Thetareact_tr = atan(V3pd_tr[0]/V3pd_tr[2]);
    double Phireact_tr = atan(V3pd_tr[1]/V3pd_tr[2]);

    // Calculate delta based on angle
    double V3pd_lab[3];
    HRSTransTCSNHCS::X_TCS2HCS(V3pd_tr[0], V3pd_tr[1], V3pd_tr[2], fHRSAngle, V3pd_lab[0], V3pd_lab[1], V3pd_lab[2]);

    double cosscatangle = V3pd_lab[2]/(sqrt(V3pd_lab[0]*V3pd_lab[0]+V3pd_lab[1]*V3pd_lab[1]+V3pd_lab[2]*V3pd_lab[2]));

    double scatmom = (fTargetMass*fBeamEnergy)/(fTargetMass+fBeamEnergy-fBeamEnergy*cosscatangle);

    double Delta = scatmom/fHRSMomentum-1-fEnergyLoss/fHRSMomentum;

    V5react_tr[0] = Xreact_tr;
    V5react_tr[1] = Thetareact_tr;
    V5react_tr[2] = Yreact_tr;
    V5react_tr[3] = Phireact_tr;
    V5react_tr[4] = Delta;

#ifdef GUN_DEBUG
    printf("G2PGun: %e\t%e\t%e\t%e\t%e\n", V5react_tr[0], V5react_tr[1], V5react_tr[2], V5react_tr[3], V5react_tr[4]);
#endif

    return true;
}

bool G2PGun::ShootData(double* V5bpm_lab, double* V5react_tr, double* V5fp_tr)
{
    bool noerror = true;

    sData tempdata;

    if (fData.empty()) return false;
    tempdata = fData.back();
    V5bpm_lab[0] = tempdata.xb;
    V5bpm_lab[1] = tempdata.tb;
    V5bpm_lab[2] = tempdata.yb;
    V5bpm_lab[3] = tempdata.pb;
    V5bpm_lab[4] = tempdata.zb;
    V5fp_tr[0] = tempdata.xf;
    V5fp_tr[1] = tempdata.tf;
    V5fp_tr[2] = tempdata.yf;
    V5fp_tr[3] = tempdata.pf;
    V5fp_tr[4] = 0.0;

#ifdef GUN_DEBUG
    printf("G2PGun: %e\t%e\t%e\t%e\t%e\n", V5fp_tr[0], V5fp_tr[1], V5fp_tr[2], V5fp_tr[3], V5fp_tr[4]);
#endif

    fData.pop_back();
    return noerror;
}

bool G2PGun::ShootOpticsData(double* V5bpm_lab, double* V5react_tr, double* V5fp_tr)
{
    bool noerror = true;

    sData tempdata;

    if (fData.empty()) return false;
    tempdata = fData.back();
    int index = tempdata.ind;
    V5bpm_lab[0] = tempdata.xb;
    V5bpm_lab[1] = tempdata.tb;
    V5bpm_lab[2] = tempdata.yb;
    V5bpm_lab[3] = tempdata.pb;
    V5bpm_lab[4] = tempdata.zb;
    V5fp_tr[0] = tempdata.xf;
    V5fp_tr[1] = tempdata.tf;
    V5fp_tr[2] = tempdata.yf;
    V5fp_tr[3] = tempdata.pf;
    V5fp_tr[4] = 0.0;

    int col = index/(fSieve.nRow);
    int row = index%(fSieve.nRow);

    double Xreact_tr, Yreact_tr, Zreact_tr;
    HRSTransTCSNHCS::X_HCS2TCS(V5bpm_lab[0], V5bpm_lab[2], V5bpm_lab[4], fHRSAngle, Xreact_tr, Yreact_tr, Zreact_tr);

    double V3sieve_tr[3];
    double V3pd_tr[3];

    V3sieve_tr[0] = fSieve.fXOffset+fSieve.fX[row];
    V3sieve_tr[1] = fSieve.fYOffset+fSieve.fY[col];
    V3sieve_tr[2] = fSieve.fZ;

    V3pd_tr[0] = V3sieve_tr[0]-Xreact_tr;
    V3pd_tr[1] = V3sieve_tr[1]-Yreact_tr;
    V3pd_tr[2] = V3sieve_tr[2]-Zreact_tr;

    double Thetareact_tr = atan(V3pd_tr[0]/V3pd_tr[2]);
    double Phireact_tr = atan(V3pd_tr[1]/V3pd_tr[2]);

    // Calculate delta based on angle
    double V3pd_lab[3];
    HRSTransTCSNHCS::X_TCS2HCS(V3pd_tr[0], V3pd_tr[1], V3pd_tr[2], fHRSAngle, V3pd_lab[0], V3pd_lab[1], V3pd_lab[2]);

    double cosscatangle = V3pd_lab[2]/(sqrt(V3pd_lab[0]*V3pd_lab[0]+V3pd_lab[1]*V3pd_lab[1]+V3pd_lab[2]*V3pd_lab[2]));

    double scatmom = (fTargetMass*fBeamEnergy)/(fTargetMass+fBeamEnergy-fBeamEnergy*cosscatangle);

    double Delta = scatmom/fHRSMomentum-1-fEnergyLoss/fHRSMomentum;

    V5react_tr[0] = Xreact_tr;
    V5react_tr[1] = Thetareact_tr;
    V5react_tr[2] = Yreact_tr;
    V5react_tr[3] = Phireact_tr;
    V5react_tr[4] = Delta;

#ifdef GUN_DEBUG
    printf("G2PGun: %e\t%e\t%e\t%e\t%e\n", V5react_tr[0], V5react_tr[1], V5react_tr[2], V5react_tr[3], V5react_tr[4]);
#endif

    fData.pop_back();
    return noerror;
}

void G2PGun::GetReactPoint(double xb, double yb, double zb, double* V5)
{
    if (bUseField) {
        double x[3] = { xb, yb, 0.0 };
        double p[3] = { 0.0,
                        fBeamEnergy*sin(fBeamTiltAngle),
                        fBeamEnergy*cos(fBeamTiltAngle) };
        G2PDrift::Drift(x, p, zb, 10.0, x, p);
        V5[0] = x[0];
        V5[1] = acos(p[2]/fBeamEnergy);
        V5[2] = x[1];
        V5[3] = atan(p[1]/p[0]);
        if (p[1]/p[0]<0) V5[3]+= 180.0*kDEG;
        if (p[1]<0) V5[3]+= 180.0*kDEG;
        V5[4] = x[2];
    }
    else {
        V5[0] = xb;
        V5[1] = 0.0;
        V5[2] = yb;
        V5[3] = 0.0;
        V5[4] = zb;
    }
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
