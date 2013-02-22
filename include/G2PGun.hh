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
    G2PGun(const char* dist);
    ~G2PGun();

    typedef bool (G2PGun::*pf_Gun)(double*, double*, double*);
    
    void SetHRSAngle(double value) { fHRSAngle = value; }
    void SetHRSMomentum(double value) { fHRSMomentum = value; }
    void SetBeamEnergy(double value) { fBeamEnergy = value; }
    
    void SetBeamX(double value) { fBeamX_lab = value; }
    void SetBeamY(double value) { fBeamY_lab = value; }
    void SetBeamTh(double value) { fBeamTh_lab = value; }
    void SetBeamPh(double value) { fBeamPh_lab = value; }
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

    void SetBPMZLab(double value) { fBPMZ_lab = value; }
    void SetBPMPosRes(double value) { fBPMPosRes = value; }
    void SetBPMAngRes(double value) { fBPMAngRes = value; }
   
    void SetDataFile(const char* name) { pFileName = name; }

    bool IsInit() { return bIsInit; }
    bool IsUsingData() { return bUseData; }

    int GetSetting() { return iSetting; }
    double GetBPMPosRes() { return fBPMPosRes; }
    double GetBPMAngRes() { return fBPMAngRes; }

    virtual void Init();
    virtual bool Shoot(double* V5beam_lab, double* V5bpm_lab, double* V5tg_tr) { return (this->*pfGunSelector)(V5beam_lab, V5bpm_lab, V5tg_tr); }
    virtual void GetFP(double* V5fp_tr);

private:
    void SetGun();   

    bool ShootDelta(double* V5beam_lab, double* V5bpm_lab, double* V5tg_tr);
    bool ShootGaus(double* V5beam_lab, double* V5bpm_lab, double* V5tg_tr);
    bool ShootFlat(double* V5beam_lab, double* V5bpm_lab, double* V5tg_tr);
    bool ShootTest(double* V5beam_lab, double* V5bpm_lab, double* V5tg_tr);
    bool ShootSieve(double* V5beam_lab, double* V5bpm_lab, double* V5tg_tr);
    bool ShootData(double* V5beam_lab, double* V5bpm_lab, double* V5tg_tr);

    void GetBPMValue(const double* V5beam_lab, double* V5bpm_lab);
    bool LoadData();

    bool bIsInit;
    
    int iSetting;
    bool bUseData;
    bool bUseField;

    double fHRSAngle;
    double fHRSMomentum;
    double fBeamEnergy;

    double fBeamX_lab, fBeamY_lab;
    double fBeamTh_lab, fBeamPh_lab;
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

    double fBPMZ_lab;
    double fBPMPosRes;
    double fBPMAngRes;

    typedef struct {
        int ind;
        double xb, tb, yb, pb, zb, xf, tf, yf, pf;
    } sData;
    
    vector<sData> fData;
    sData fDataAtIndex;
    int iIndex;

    const char* pFileName;
    
    pf_Gun pfGunSelector;
    
    ClassDef(G2PGun,1);
};

#endif
