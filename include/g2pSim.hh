#ifndef G2PSIM_H
#define G2PSIM_H

#include <cstdio>
#include <cstring>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

#include "HRSRecUseDB.hh"

// definition:
// iArm: Set Arm 0 means left arm, 1 means right arm
// pBPMRes: Set bpm resolution
// pHRSMomentum: Set up HRS Momentum
// iDirection: Set sim direction: 0 means forward, 1 means backward, 2 means both
// iSource: Set data source: 1,2,3 means use delta,flat,gaussian distribution for input position and angle, 0 means use real data in file input_tg.dat and input_fp.dat
// iSetting: Set experiment setting: 10 means normal 484816 septa, 11 means 484816 septa with shim, 12 means 403216 septa with shim, 13 means 400016 septa with shim

//void TestSNAKE(int iNEvent, int iArm, int iSetting, int iSource, int iDirection, double pHRSMomentum, double pBPMRes);

class g2pSim
{
public:
    g2pSim();
    ~g2pSim();

    void SetNEvent(int n) { iNEvent = n; }
    void SetArm(const char *label) { bIsLeftArm = (strcmp(label,"L")==0)?true:false; }
    void SetHRSAngle(double angle) { pHRSAngle = angle; }
    void SetHRSMomentum(double momentum) { pHRSMomentum = momentum; }
    void SetBPMRes(double res) { pBPMRes = res; }
    void SetHRSSetting(int setting) { iSetting = setting; }
    void SetDATASource(int setting) { iSource = setting; }
    void SetSimDirection(int setting) { iDirection = setting; }

    void Init();
    void Run();

private:
    void Clear();
    void ReadData(double *pV5);
    
    int    Index;
    int    iNEvent;
	bool   bIsLeftArm;
    int    iSource;
    int    iDirection;
    int    iSetting;
    
	double pHRSAngle;  // Set the same value to optics setting
    double pHRSMomentum;
    
	double pV3bpm_tr[3];
    double pV3bpm_lab[3];
    double pXtg_BPM_tr;
    double pBPMRes;

    double pV5fp_tr[5];
	double pV5tg_tr[5];
	double pV5rec_tr[5];
	double pV5tg_lab[5];
    double pV5rec_lab[5];
    double pV5dbrec_tr[5];

    HRSRecUseDB *pRecDB;

    TTree *pTree;
    TTree *pConfig;
    TFile *pFile;

    FILE *fp;
};

#endif
