#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

#include "HRSRecUseDB.hh"
#include "HRSTransport.hh"
#include "TransTCSNHCS.hh"
#include "CrossSection.hh"
#include "Rand.hh"

#include "g2pSim.hh"

const double deg = TMath::Pi()/180.0;

void VDCSmearing(double* pV5_fp)
{
    double mWireChamberRes_x = 0.0013; //m;
    double mWireChamberRes_y = 0.0013; //m;
    double mWireChamberRes_theta = 0.0003; //rad;
    double mWireChamberRes_phi = 0.0003; //rad;

    pV5_fp[0] += fGausRand(0, mWireChamberRes_x);
    pV5_fp[2] += fGausRand(0, mWireChamberRes_y);
    pV5_fp[1] += fGausRand(0, mWireChamberRes_theta);
    pV5_fp[3] += fGausRand(0, mWireChamberRes_phi);
}

// definition:
// iArm: Set Arm 0 means left arm, 1 means right arm
// pBPMRes: Set bpm resolution
// pHRSMomentum: Set up HRS Momentum
// iDirection: Set sim direction: 0 means forward, 1 means backward, 2 means both
// iSource: Set data source: 1,2,3 means use delta, gaussian, flat distribution for input position and angle, 0 means use real data in file input_tg.dat and input_fp.dat
// iSetting: Set experiment setting: 10 means normal 484816 septa, 11 means
// 484816 septa with shim, 12 means 403216 septa with shim, 13 means 400016
// septa with shim

g2pSim::g2pSim()
    :Index(0), bIsLeftArm(true), iSource(0), iDirection(0),
     iSetting(11), pHRSMomentum(2.251), pBPMRes(0.0), pRecDB(NULL)
{
    pHRSAngle = 5.767*deg;
    pXtg_BPM_tr = 0.0;

    Clear();
}

g2pSim::~g2pSim()
{
    pFile->Write();
    pFile->Close();
    
    fclose(fp);
}

void g2pSim::Init()
{
    if (bIsLeftArm) {
        pRecDB = new HRSRecUseDB("L","db_vdc.L.dat");
    }
    else{
        pRecDB = new HRSRecUseDB("R","db_vdc.R.dat");
    }

    pFile = new TFile(Form("result_E%02d_S%d.root", iSetting, iSource), "recreate");
	pTree = new TTree("T", "sim result");
    pConfig = new TTree("config", "sim configuration");

	pConfig->Branch("BPMRes",&pBPMRes,"BPMRes/D");
	pConfig->Branch("Left",&bIsLeftArm,"Left/B");
	pConfig->Branch("HRSAngle",&pHRSAngle,"HRSAngle/D");
	pConfig->Branch("P0",&pHRSMomentum,"P0/D");
	pConfig->Branch("Xtg_BPM_tr",&pXtg_BPM_tr,"Xtg_BPM_tr/D");
    pConfig->Branch("iSource",&iSource,"DATASource/D");
	pConfig->Branch("iDirection",&iDirection,"SimDirection/D");
	pConfig->Branch("iSetting",&iSetting,"HRSSetting/D");

    pConfig->Fill();

    pTree->Branch("Index",&Index,"Index/I");

	pTree->Branch("Xfp_tr",&pV5fp_tr[0],"Xfp_tr/D");
	pTree->Branch("Thetafp_tr",&pV5fp_tr[1],"Thetafp_tr/D");
	pTree->Branch("Yfp_tr",&pV5fp_tr[2],"Yfp_tr/D");
	pTree->Branch("Phifp_tr",&pV5fp_tr[3],"Phifp_tr/D");

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
    
	pTree->Branch("Xdbrec_tr",&pV5dbrec_tr[0],"Xdbrec_tr/D");
	pTree->Branch("Thetadbrec_tr",&pV5dbrec_tr[1],"Thetadbrec_tr/D");
	pTree->Branch("Ydbrec_tr",&pV5dbrec_tr[2],"Ydbrec_tr/D");
	pTree->Branch("Phidbrec_tr",&pV5dbrec_tr[3],"Phidbrec_tr/D");
	pTree->Branch("Deltadbrec",&pV5dbrec_tr[4],"Deltadbrec/D");

    if (iSource==0) {
        if ((iDirection==0)||(iDirection==2)) {
            fp=fopen("input_tg.dat","r");
        }
        else if (iDirection==1) {
            fp=fopen("input_fp.dat","r");
        }
    }
}

void g2pSim::Run()
{
    int NThrown=0;
    while (Index<iNEvent) {
        if (iSource==2)
        {
			pV5tg_tr[0] = 0.010 * 2 * (fRand()-0.5);
			pV5tg_tr[1] = 0.070 * 2 * (fRand()-0.5);
			pV5tg_tr[2] = 0.010 * 2 * (fRand()-0.5);
			pV5tg_tr[3] = 0.035 * 2 * (fRand()-0.5);
			pV5tg_tr[4] = 0.050 * 2 * (fRand()-0.5);
		}
		else if (iSource==3)
		{
			pV5tg_tr[0] = fGausRand(0,0.00010/2);
			pV5tg_tr[1] = fGausRand(0,0.010/2);
			pV5tg_tr[2] = fGausRand(0,0.00010/2);
			pV5tg_tr[3] = fGausRand(0,0.010/2);
			pV5tg_tr[4] = fGausRand(0,0.030/2);
		}
        else if (iSource==0)
        {
            ReadData(pV5tg_tr);
            pV5tg_tr[0] = 0.0;
        }

        bool bGoodParticle;

        bGoodParticle = SNAKEForward(bIsLeftArm, iSetting, pV5tg_tr, pV5fp_tr);

        if (iDirection==0) {
            pTree->Fill();
            continue;
        }

        if (iDirection==1) {
            ReadData(pV5fp_tr);
            pV5fp_tr[4] = 0.0;
        }
        
        bGoodParticle &= SNAKEBackward(bIsLeftArm, iSetting, pV5fp_tr, pV5tg_tr);
        pRecDB->CalcTargetCoords(pV5fp_tr, pV5dbrec_tr);

        if (bGoodParticle) {
            pTree->Fill();
            Index++;
        }

        NThrown++;
    }

    printf("%d/%d event saved into root file", Index, NThrown);
}

void g2pSim::ReadData(double *pV5)
{
    int temp;

    if ((iDirection==0)||(iDirection==2)) { // no x_tg
        fscanf(fp, "%d%lf%lf%lf%lf", &temp, &pV5[1], &pV5[2], &pV5[3], &pV5[4]);
    }
    else if (iDirection==1) { // no delta
        fscanf(fp, "%d%lf%lf%lf%lf", &temp, &pV5[0], &pV5[1], &pV5[2], &pV5[3]);
    }
}

void g2pSim::Clear()
{
    for(int i=0; i<3; i++) {
        pV3bpm_tr[i] = 0.0;
        pV3bpm_lab[i] = 0.0;
    }

    for(int i=0; i<5; i++) {
        pV5fp_tr[i] = 0.0;
        pV5tg_tr[i] = 0.0;
        pV5rec_tr[i] = 0.0;
        pV5tg_lab[i] = 0.0;
        pV5rec_lab[i] = 0.0;
        pV5dbrec_tr[i] = 0.0;
    }
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
		// 	Index++;
		// }

