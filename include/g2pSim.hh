#ifndef G2PSIM_H
#define G2PSIM_H

#include <cstdio>
#include <cstring>
#include <vector>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

#include "HRSGun.hh"
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

    typedef void (g2pSim::*run_ptr)();
    
    void SetNEvent(int n) { pNEvent = n; }
    void SetArm(const char *label) { pIsLeftArm = (strcmp(label,"L")==0)?true:false; }
    void SetHRSAngle(double angle) { pHRSAngel = angle; }
    void SetHRSMomentum(double momentum) { pHRSMomentum = momentum; }
    void SetHRSSetting(int setting) { pSetting = setting; }

    void SetRootFile(const char *name) { pFileName = name; }

    void SetGun(HRSGun &gun) { pGun = &gun; }

    bool IsInit() { return pIsInit; }

    void Init();
    void Run();
    void End();

private:
    void InitTree();

    void RunSim();
    void RunData();

    void Clear();

    bool pIsInit;

    int pIndex;
    int pNEvent;
	bool pIsLeftArm;
    int pSetting;
    int pGunSetting;

    double pHRSAngle;  // Set the same value to optics setting
    double pHRSMomentum;

    double pV3bpm_lab[3];

    double pV5fp_tr[5];
    double pV5fp_rot[5];
    double pV5fpdata_tr[5];
    double pV5fpdata_rot[5];
	double pV5tg_tr[5];
    double pV5tg_lab[5];
	double pV5rec_tr[5];
    double pV5rec_lab[5];
    double pV5recdb_tr[5];
    double pV5recdb_lab[5];

    double pCrossSection;

    HRSRecUseDB *pRecDB;
    
    HRSGun *pGun;

    TTree *pTree;
    TTree *pConfig;
    TFile *pFile;
    const char *pFileName;

    run_ptr pRunSelector;

    HRSRand *pRand;
};

#endif
