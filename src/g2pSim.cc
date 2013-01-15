#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

#include "HRSGun.hh"
#include "HRSRecUseDB.hh"
#include "HRSTransTCSNHCS.hh"
#include "HRSRand.hh"

#include "HRSTransport.hh"
#include "CrossSection.hh"

#include "g2pSim.hh"

const double deg = TMath::Pi()/180.0;

// void VDCSmearing(double* pV5_fp)
// {
//     double mWireChamberRes_x = 0.0013; //m;
//     double mWireChamberRes_y = 0.0013; //m;
//     double mWireChamberRes_theta = 0.0003; //rad;
//     double mWireChamberRes_phi = 0.0003; //rad;

//     pV5_fp[0] += fGausRand(0, mWireChamberRes_x);
//     pV5_fp[2] += fGausRand(0, mWireChamberRes_y);
//     pV5_fp[1] += fGausRand(0, mWireChamberRes_theta);
//     pV5_fp[3] += fGausRand(0, mWireChamberRes_phi);
// }

// definition:
// iArm: Set Arm 0 means left arm, 1 means right arm
// pBPMRes: Set bpm resolution
// pHRSMomentum: Set up HRS Momentum
// pDirection: Set sim direction: 0 means forward, 1 means backward, 2 means both
// iSource: Set data source: 1,2,3 means use delta, gaussian, flat distribution for input position and angle, 0 means use real data in file input_tg.dat and input_fp.dat
// iSetting: Set experiment setting: 10 means normal 484816 septa, 11 means
// 484816 septa with shim, 12 means 403216 septa with shim, 13 means 400016
// septa with shim

g2pSim::g2pSim()
    :pIsInit(false), pIndex(1), pNEvent(10000), pIsLeftArm(true),
     pSetting(11), pHRSAngle(5.767*deg), pHRSMomentum(2.251),
     pCrossSection(0), pRecDB(NULL), pGun(NULL), pRand(NULL)
{
    Clear();
}

g2pSim::~g2pSim()
{
    // Nothing to be done
}

void g2pSim::Init()
{
    pIsInit = true;
    
    pRand = new HRSRand();

    pGun->SetRand(pRand);
    pGun->SetHRSAngle(pHRSAngle);
    pGun->Init();

    if (pGun->IsInit()) {
        if (pIsLeftArm)
            pRecDB = new HRSRecUseDB("L","db_L.vdc.dat");
        else
            pRecDB = new HRSRecUseDB("R","db_R.vdc.dat");

        InitTree();

        if (pGun->IsUsingData())
            pRunSelector = &g2pSim::RunData;
        else
            pRunSelector = &g2pSim::RunSim;
    }
}

void g2pSim::Run()
{
    if (pGun==NULL) {
        //
    }
    else{
        Init();
        if (IsInit()) {
            (this->*pRunSelector)();
            End();
        }
    }
}

void g2pSim::RunSim()
{
    while (pIndex<=pNEvent) {
        pGun->Shoot(pV3bpm_lab, pV5tg_tr);

        pV5tg_tr[1] = tan(pV5tg_tr[1]);
        pV5tg_tr[3] = tan(pV5tg_tr[3]);

        pGoodParticle = SNAKEForward(pIsLeftArm, pSetting, pV5tg_tr, pV5fp_tr);
        pRecDB->TransTr2Rot(pV5fp_tr, pV5fp_rot);

        pV5fp_tr[4] = pV5tg_tr[0];
        pGoodParticle &= SNAKEBackward(pIsLeftArm, pSetting, pV5fp_tr, pV5rec_tr);

        pRecDB->CalcTargetCoords(pV5fp_rot, pV5recdb_tr);

        pTree->Fill();

        if ((pIndex%10000)==0) printf("%d\n", pIndex);
        pIndex++;
    }
}

void g2pSim::RunData()
{
    while (pIndex<=pNEvent) {
        pGun->Shoot(pV3bpm_lab, pV5fpdata_tr);
        if (!pGun->IsInit()) break;

        pRecDB->TransTr2Rot(pV5fpdata_tr, pV5fpdata_rot);
        pRecDB->CalcTargetCoords(pV5fpdata_rot, pV5tg_tr);

        pGoodParticle = SNAKEForward(pIsLeftArm, pSetting, pV5tg_tr, pV5fp_tr);
        pRecDB->TransTr2Rot(pV5fp_tr, pV5fp_rot);

        pV5fpdata_tr[4] = pV5tg_tr[0];
        
        SNAKEBackward(pIsLeftArm, pSetting, pV5fpdata_tr, pV5rec_tr);

        pRecDB->CalcTargetCoords(pV5fpdata_rot, pV5recdb_tr);

        pTree->Fill();

        if ((pIndex%10000)==0) printf("%d\n", pIndex);
        pIndex++;
    }
}

void g2pSim::End()
{
    pFile->Write();
    pFile->Close();
}

void g2pSim::InitTree()
{
    pFile = new TFile(pFileName, "recreate");
	pTree = new TTree("T", "sim result");
    pConfig = new TTree("config", "sim configure");

    pConfig->Branch("N", &pNEvent, "N/I");
    pConfig->Branch("HRSSetting", &pSetting, "HRSSetting/I");
 	pConfig->Branch("IsLeftArm", &pIsLeftArm, "IsLeftArm/O");   
    pConfig->Branch("HRSAngle", &pHRSAngle, "HRSAngle/D");
    pConfig->Branch("HRSMomentum", &pHRSMomentum, "HRSMomentum/D");

    pGunSetting = pGun->GetSetting();
    double posres = pGun->GetPosResolution();
    double angleres = pGun->GetAngleResolution();
    double deltares = pGun->GetDeltaResolution();
    pConfig->Branch("GunSetting", &pGunSetting, "GunSetting/I");
	pConfig->Branch("GunPosRes", &posres, "GunPosRes/D");
	pConfig->Branch("GunAngleRes", &angleres, "GunAngleRes/D");
	pConfig->Branch("GunDeltaRes", &deltares, "GunDeltaRes/D");

    pConfig->Fill();

    pTree->Branch("Index", &pIndex,"Index/I");
    pTree->Branch("IsGood", &pGoodParticle, "IsGood/O");
    
	pTree->Branch("Xfp_tr",&pV5fp_tr[0],"Xfp_tr/D");
	pTree->Branch("Thetafp_tr",&pV5fp_tr[1],"Thetafp_tr/D");
	pTree->Branch("Yfp_tr",&pV5fp_tr[2],"Yfp_tr/D");
	pTree->Branch("Phifp_tr",&pV5fp_tr[3],"Phifp_tr/D");

    pTree->Branch("Xfp_rot",&pV5fp_rot[0],"Xfp_rot/D");
	pTree->Branch("Thetafp_rot",&pV5fp_rot[1],"Thetafp_rot/D");
	pTree->Branch("Yfp_rot",&pV5fp_rot[2],"Yfp_rot/D");
	pTree->Branch("Phifp_rot",&pV5fp_rot[3],"Phifp_rot/D");

    pTree->Branch("Xfpdata_tr",&pV5fpdata_tr[0],"Xfpdata_tr/D");
	pTree->Branch("Thetafpdata_tr",&pV5fpdata_tr[1],"Thetafpdata_tr/D");
	pTree->Branch("Yfpdata_tr",&pV5fpdata_tr[2],"Yfpdata_tr/D");
	pTree->Branch("Phifpdata_tr",&pV5fpdata_tr[3],"Phifpdata_tr/D");

    pTree->Branch("Xfpdata_rot",&pV5fpdata_rot[0],"Xfpdata_rot/D");
	pTree->Branch("Thetafpdata_rot",&pV5fpdata_rot[1],"Thetafpdata_rot/D");
	pTree->Branch("Yfpdata_rot",&pV5fpdata_rot[2],"Yfpdata_rot/D");
	pTree->Branch("Phifpdata_rot",&pV5fpdata_rot[3],"Phifpdata_rot/D");

	pTree->Branch("Xtg_tr",&pV5tg_tr[0],"Xtg_tr/D");
	pTree->Branch("Thetatg_tr",&pV5tg_tr[1],"Thetatg_tr/D");
	pTree->Branch("Ytg_tr",&pV5tg_tr[2],"Ytg_tr/D");
	pTree->Branch("Phitg_tr",&pV5tg_tr[3],"Phitg_tr/D");
	pTree->Branch("Delta",&pV5tg_tr[4],"Delta/D");

	pTree->Branch("Xrec_tr",&pV5rec_tr[0],"Xrec_tr/D");
	pTree->Branch("Thetarec_tr",&pV5rec_tr[1],"Thetarec_tr/D");
	pTree->Branch("Yrec_tr",&pV5rec_tr[2],"Yrec_tr/D");
	pTree->Branch("Phirec_tr",&pV5rec_tr[3],"Phirec_tr/D");
	pTree->Branch("Deltarec",&pV5rec_tr[4],"Deltarec/D");

    pTree->Branch("Xrecdb_tr",&pV5recdb_tr[0],"Xrecdb_tr/D");
	pTree->Branch("Thetarecdb_tr",&pV5recdb_tr[1],"Thetarecdb_tr/D");
	pTree->Branch("Yrecdb_tr",&pV5recdb_tr[2],"Yrecdb_tr/D");
	pTree->Branch("Phirecdb_tr",&pV5recdb_tr[3],"Phirecdb_tr/D");
	pTree->Branch("Deltarecdb",&pV5recdb_tr[4],"Deltarecdb/D");

	pTree->Branch("Xtg_lab",&pV5tg_lab[0],"Xtg_lab/D");
	pTree->Branch("Thetatg_lab",&pV5tg_lab[1],"Thetatg_lab/D");
	pTree->Branch("Ytg_lab",&pV5tg_lab[2],"Ytg_lab/D");
	pTree->Branch("Phitg_lab",&pV5tg_lab[3],"Phitg_lab/D");
	pTree->Branch("Ztg_lab",&pV5tg_lab[4],"Ztg_lab/D");

	pTree->Branch("Xrec_lab",&pV5rec_lab[0],"Xrec_lab/D");
	pTree->Branch("Thetarec_lab",&pV5rec_lab[1],"Thetarec_lab/D");
	pTree->Branch("Yrec_lab",&pV5rec_lab[2],"Yrec_lab/D");
	pTree->Branch("Phirec_lab",&pV5rec_lab[3],"Phirec_lab/D");
	pTree->Branch("Zrec_lab",&pV5tg_lab[4],"Zrec_lab/D");

    pTree->Branch("Xrecdb_lab",&pV5recdb_lab[0],"Xrecdb_lab/D");
	pTree->Branch("Thetarecdb_lab",&pV5recdb_lab[1],"Thetarecdb_lab/D");
	pTree->Branch("Yrecdb_lab",&pV5recdb_lab[2],"Yrecdb_lab/D");
	pTree->Branch("Phirecdb_lab",&pV5recdb_lab[3],"Phirecdb_lab/D");
    pTree->Branch("Zrecdb_lab",&pV5recdb_lab[3],"Zrecdb_lab/D");

    pTree->Branch("XS", &pCrossSection, "XS/D");
}

void g2pSim::Clear()
{
    memset(pV3bpm_lab, 0, sizeof(pV3bpm_lab));
    memset(pV5fp_tr, 0, sizeof(pV5fp_tr));
    memset(pV5tg_tr, 0, sizeof(pV5tg_tr));
    memset(pV5rec_tr, 0, sizeof(pV5rec_tr));
    memset(pV5tg_lab, 0, sizeof(pV5tg_lab));
    memset(pV5rec_lab, 0, sizeof(pV5rec_lab));
    memset(pV5recdb_tr, 0, sizeof(pV5recdb_tr));
    pCrossSection = 0;
    pGoodParticle = false;
}
        
		//throw away this event if delta_rec>=1.0
		// {
		// 	Transform::X_TCS2HCS(pV5tg_tr[0],pV5tg_tr[2],0.0,pHRSAngle,pV5tg_lab[0],pV5tg_lab[2],pV5tg_lab[4]);
		// 	Transform::P_TCS2HCS(pV5tg_tr[1],pV5tg_tr[3],pHRSAngle,pV5tg_lab[1],pV5tg_lab[3]);

		// 	Transform::X_TCS2HCS(pV5rec_tr[0],pV5rec_tr[2],0.0,pHRSAngle,pV5rec_lab[0],pV5rec_lab[2],pV5rec_lab[4]);
		// 	Transform::P_TCS2HCS(pV5rec_tr[1],pV5rec_tr[3],pHRSAngle,pV5rec_lab[1],pV5rec_lab[3]);
            
		// 	//Reconstruct use optics database
		// 	pRecDB->CalcTargetCoords(pV5fp_tr,pV5rec_db_tr);

		// 	pTree->Fill();
		// 	pIndex++;
		// }

