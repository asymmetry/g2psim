// This file defines a class G2PGun.
// This class is used in G2PSim class as particle gun.
// The active gun is chosen during initializing.
// G2PSim class will call Shoot() to get kinematic variables. It is a virtual
//+function so you can rewrite it by inheriting this class.
//
// History:
//   Jan 2013, C. Gu, First public version.
//   Jan 2013, C. Gu, Add ShootSieve() method.
//   Feb 2013, C. Gu, Rewrite the sieve and data parts
//

#ifndef G2P_GUN_H
#define G2P_GUN_H

#include <cstdio>
#include <vector>

#include "TObject.h"

using namespace std;

class G2PGun : public TObject
{
public:
    G2PGun();
////////////////////////////////////////////////////////////////////////
// Definition of guns:
// gaus : gaussian distribution
// flat : flat distribution
// sieve : use sieve as mask
// fastsieve : directly calculation from survey, only work without field
// data : production data
// opticsdata : use cut index to indicate sieve holes
////////////////////////////////////////////////////////////////////////
    G2PGun(const char* dist);
    ~G2PGun();

    typedef bool (G2PGun::*pf_Gun)(double*, double*, double*);
    
    void SetHRSAngle(double value) { fHRSAngle = value; }
    void SetHRSMomentum(double value) { fHRSMomentum = value; }
    void SetBeamEnergy(double value) { fBeamEnergy = value; }
    
    void SetBeamX(double value) { fBeamX_lab = value; }
    void SetBeamY(double value) { fBeamY_lab = value; }
    void SetBeamR(double value) { fBeamR = value; }
    
    void SetReactZ(double value) { fReactZLow_lab = value; fReactZHigh_lab = value;}
    void SetReactZRange(double low, double high) { fReactZLow_lab = low; fReactZHigh_lab = high; }

    void SetTargetTh(double value) { fTargetThLow_tr = value; fTargetThHigh_tr = value; }
    void SetTargetThRange(double low, double high) { fTargetThLow_tr = low; fTargetThHigh_tr = high; }
    void SetTargetPh(double value) { fTargetPhLow_tr = value; fTargetPhHigh_tr = value; }
    void SetTargetPhRange(double low, double high) { fTargetPhLow_tr = low; fTargetPhHigh_tr = high; }

    void SetDelta(double value) { fDeltaLow = value; fDeltaHigh = value; }
    void SetDeltaRange(double low, double high) { fDeltaLow = low; fDeltaHigh = high; }
    
    void SetSigmaPosLab(double value) { fSigmaPos_lab = value; }
    void SetSigmaAngLab(double value) { fSigmaAng_lab = value; }
    void SetSigmaAngTr(double value) { fSigmaAng_tr = value; }
    void SetSigmaDelta(double value) { fSigmaDelta = value; }

    void SetDataFile(const char* name) { pFileName = name; }

    bool IsInit() { return bIsInit; }
    bool IsUsingData() { return bUseData; }

    int GetSetting() { return iSetting; }

    virtual void Init();
    virtual bool Shoot(double* V51, double* V52, double* V53 = NULL) { return (this->*pfGunSelector)(V51, V52, V53); }

private:
    void SetGun();
    void SetBeamTiltAngle();
    void SetSieve();

    bool ShootGaus(double* V5beam_lab, double* V5react_tr, double* reserved);
    bool ShootFlat(double* V5beam_lab, double* V5react_tr, double* reserved);
    bool ShootSieve(double* V5beam_lab, double* V5react_tr, double* reserved);
    bool ShootSieveFast(double* V5beam_lab, double* V5react_tr, double* reserved);
    bool ShootData(double* V5bpm_lab, double* V5react_tr, double* V5fp_tr);
    bool ShootOpticsData(double* V5bpm_lab, double* V5react_tr, double* V5fp_tr);

    void GetReactPoint(double x, double y, double z, double* V5);
    bool LoadData();

    bool bIsInit;
    
    int iSetting;
    bool bUseData;
    bool bUseField;

    double fHRSAngle;
    double fHRSMomentum;
    double fBeamEnergy;

    double fBeamX_lab, fBeamY_lab;
    double fBeamTiltAngle;
    double fBeamR;
    
    double fReactZLow_lab;
    double fReactZHigh_lab;

    double fTargetThLow_tr;
    double fTargetThHigh_tr;
    double fTargetPhLow_tr;
    double fTargetPhHigh_tr;

    double fDeltaLow; // in the unit of delta
    double fDeltaHigh;

    double fSigmaPos_lab;
    double fSigmaAng_lab;
    double fSigmaAng_tr;
    double fSigmaDelta;

    typedef struct {
        int nRow;
        int nCol;
        vector<double> fX;
        vector<double> fY;
        double fZ;
        double fXOffset;
        double fYOffset;
        int nLargerHole;
        vector<int> iLargerHole;
        vector<bool> bOpen;
        double fDHole;
        double fDLargerHole;
        double fThreshold;
    } sSieve;

    sSieve fSieve;
    double fTargetMass;
    double fEnergyLoss;

    typedef struct {
        int ind;
        double xb, tb, yb, pb, zb, xf, tf, yf, pf;
    } sData;
    
    vector<sData> fData;

    const char* pFileName;
    
    pf_Gun pfGunSelector;
    
    ClassDef(G2PGun,1);
};

#endif
