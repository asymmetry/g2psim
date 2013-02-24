// This file defines a class G2PSim.
// This class is the main class for this sim package.
//
// History:
//   Jan 2013, C. Gu, First public version.
//

#ifndef G2P_SIM_H
#define G2P_SIM_H

#include <cstdio>
#include <cstring>
#include <vector>

#include "TROOT.h"
#include "TObject.h"
#include "TFile.h"
#include "TTree.h"

#include "G2PGun.hh"
#include "HRSRecUseDB.hh"
#include "G2PTargetField.hh"

#include "G2PXS.hh"
#include "HRSTransport.hh"

class G2PSim : public TObject
{
public:
    G2PSim();
    ~G2PSim();

    typedef void (G2PSim::*pf_Run)();
    
    void SetNEvent(int n) { nEvent = n; }

    void SetArm(const char* label) { bIsLeftArm = (strcmp(label,"L")==0)?true:false; }
    void SetHRSAngle(double angle) { fHRSAngle = angle; }
    void SetHRSMomentum(double momentum) { fHRSMomentum = momentum; }
    void SetBeamEnergy(double value) { fBeamEnergy = value; }

    void SetEndPlaneZ(double value) { fEndPlaneZ = value; }

    void AddGun(G2PGun* gun) { pGunList.push_back(gun); }

    void SetHRSModel(HRSTransport* model) { pHRS = model; }
    void SetTargetField(G2PTargetField* model) { pField = model; bUseField = true; }
    void SetPhysModel(G2PXS* model) { pPhys=model; }

    void SetRootName(const char* name) { pFileName = name; }

    bool IsInit() { return bIsInit; }

    void Init();
    void Run();
    void Run(int n) { nEvent = n; Run(); }
    void End();

private:
    void Clear();
    void InitTree();

    void RunSim();
    void RunData();

    void VDCSmearing(double* V5_fp);
    double GetEffBPM(double xbpm_tr, double p);
    double DriftPath();

    bool bIsInit;

    bool bIsLeftArm;
    double fHRSAngle;  // Set the same value to optics setting
    double fHRSMomentum;
    double fBeamEnergy;
    bool bUseField;
    int nEvent;

    int nIndex;

    vector<G2PGun*> pGunList;
    G2PGun* pGun;
    int iGunSetting;
    double fV5beam_lab[5];
    double fV5bpm_lab[5];
    double fV5tg_tr[5];
    double fV5tg_lab[5];
    double fV5fpdata_tr[5];
    double fV5fpdata_rot[5];

    HRSTransport* pHRS;
    bool bIsGoodParticle;
    double fV5fp_tr[5];
    double fV5fp_rot[5];
	double fV5rec_tr[5];
    double fV5rec_lab[5];

    HRSRecUseDB* pRecUseDB;
    double fV5recdb_tr[5];
    double fV5recdb_lab[5];

    G2PTargetField *pField;
    double fEndPlaneZ;
    double fV5vb_tr[5];

    G2PXS* pPhys;
    double fXS;

    TTree* pTree;
    TTree* pConfig;

    TFile* pFile;
    const char* pFileName;

    pf_Run pfRunSelector;

    ClassDef(G2PSim,1);
};

#endif
