#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

#include "HRSRecUseDB.hh"
#include "HRSTransport.hh"
#include "HRSTransform_TCSNHCS.hh"
#include "Rand.hh"

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
// iBPM: Set bpm resolution
// iDirection: Set sim direction: 0 means forward, 1 means backward, 2 means both
// iSource: Set data source: 1,2,3 means use delta,flat,gaussian distribution for input position and angle, 3 means use real data in file input_tg.dat and input_fp.dat
// iSetting: Set experiment setting: 10 means normal 484816 septa, 11 means 484816 septa with shim, 12 means 403216 septa with shim, 13 means 400016 septa with shim
void TestSNAKE(int iNEvent, int iArm, int iSetting, int iSource, int iDirection, double pHRSMomentum, double pBPMRes)
{
    const double deg = TMath::Pi()/180.0;

	int    Index=0;
	int    iLeftArm=1;
	double pHRSAngle=5.767*deg;  // Set the same value to optics setting
	double pXtg_BPM_tr=0;
	double pV5fp_tr[5]={0,0,0,0,0};
	double pV5tg_tr[5]={0,0,0,0,0};
	double pV5rec_tr[5]={0,0,0,0,0};
	double pV5tg_lab[5]={0,0,0,0,0};
    double pV5rec_lab[5]={0,0,0,0,0};

	TFile* pFile= new TFile(Form("snake_Exp%02d_X%d_Distr%d.root",
		iExperiment,iSmearX0,iSourceDistr),"recreate");
	TTree* pSnake=new TTree("s","snake reconstruction");
	pSnake->Branch("Index",&Index,"Index/I");

	pSnake->Branch("BPMRes",&pBPMRes,"BPMRes/D");
	pSnake->Branch("Left",&iLeftArm,"Left/I");
	pSnake->Branch("FieldRotation",&iFieldRotation,"FieldRotation/I");
	pSnake->Branch("HRSAngle",&pHRSAngle,"HRSAngle/D");
	pSnake->Branch("P0",&pHRSMomentum,"P0/D");
	pSnake->Branch("Xtg_BPM_tr",&pXtg_BPM_tr,"Xtg_BPM_tr/D");

	pSnake->Branch("Xfp_tr",&pV5fp_tr[0],"Xfp_tr/D");
	pSnake->Branch("Thetafp_tr",&pV5fp_tr[1],"Thetafp_tr/D");
	pSnake->Branch("Yfp_tr",&pV5fp_tr[2],"Yfp_tr/D");
	pSnake->Branch("Phifp_tr",&pV5fp_tr[3],"Phifp_tr/D");

	pSnake->Branch("Xtg_tr",&pV5tg_tr[0],"Xtg_tr/D");
	pSnake->Branch("Thetatg_tr",&pV5tg_tr[1],"Thetatg_tr/D");
	pSnake->Branch("Ytg_tr",&pV5tg_tr[2],"Ytg_tr/D");
	pSnake->Branch("Phitg_tr",&pV5tg_tr[3],"Phitg_tr/D");
	pSnake->Branch("Delta",&pV5tg_tr[4],"Delta/D");

	pSnake->Branch("Xrec_tr",&pV5rec_tr[0],"Xrec_tr/D");
	pSnake->Branch("Thetarec_tr",&pV5rec_tr[1],"Thetarec_tr/D");
	pSnake->Branch("Yrec_tr",&pV5rec_tr[2],"Yrec_tr/D");
	pSnake->Branch("Phirec_tr",&pV5rec_tr[3],"Phirec_tr/D");
	pSnake->Branch("Deltarec",&pV5rec_tr[4],"Deltarec/D");

	pSnake->Branch("Xtg_lab",&pV5tg_lab[0],"Xtg_lab/D");
	pSnake->Branch("Thetatg_lab",&pV5tg_lab[1],"Thetatg_lab/D");
	pSnake->Branch("Ytg_lab",&pV5tg_lab[2],"Ytg_lab/D");
	pSnake->Branch("Phitg_lab",&pV5tg_lab[3],"Phitg_lab/D");
	pSnake->Branch("Ztg_lab",&pV5tg_lab[4],"Ztg_lab/D");

	pSnake->Branch("Xrec_lab",&pV5rec_lab[0],"Xrec_lab/D");
	pSnake->Branch("Thetarec_lab",&pV5rec_lab[1],"Thetarec_lab/D");
	pSnake->Branch("Yrec_lab",&pV5rec_lab[2],"Yrec_lab/D");
	pSnake->Branch("Phirec_lab",&pV5rec_lab[3],"Phirec_lab/D");
	pSnake->Branch("Zrec_lab",&pV5tg_lab[4],"Zrec_lab/D");


	//rec by DB
	double pV5rec_db_tr[5];
	pSnake->Branch("Xrec_db_tr",&pV5rec_db_tr[0],"Xrec_db_tr/D");
	pSnake->Branch("Thetarec_db_tr",&pV5rec_db_tr[1],"Thetarec_db_tr/D");
	pSnake->Branch("Yrec_db_tr",&pV5rec_db_tr[2],"Yrec_db_tr/D");
	pSnake->Branch("Phirec_db_tr",&pV5rec_db_tr[3],"Phirec_db_tr/D");
	pSnake->Branch("Deltarec_db",&pV5rec_db_tr[4],"Deltarec_db/D");

	HRSRecUseDB *pRecDBL=new HRSRecUseDB("L","db_L.vdc.dat");
	HRSRecUseDB *pRecDBR=new HRSRecUseDB("R","db_R.vdc.dat");
	HRSRecUseDB *pRecDB=0;

	//do smear: 0: no smear; 1: flat; 2 gaus 
	int NThrown=0;
	Index=0;
	while(Index<iNEvent)
	{
		if(iSourceDistr==1)
		{
			pV5tg_tr[0] = 0.010 * 2 * (fRand()-0.5);
			pV5tg_tr[1] = 0.070 * 2 * (fRand()-0.5);
			pV5tg_tr[2] = 0.010 * 2 * (fRand()-0.5);
			pV5tg_tr[3] = 0.035 * 2 * (fRand()-0.5);
			pV5tg_tr[4] = 0.050 * 2 * (fRand()-0.5);
		}
		else if(iSourceDistr==2)
		{
			pV5tg_tr[0] = fRandGaus(0,0.00010/2);
			pV5tg_tr[1] = fRandGaus(0,0.010/2);
			pV5tg_tr[2] = fRandGaus(0,0.00010/2);
			pV5tg_tr[3] = fRandGaus(0,0.010/2);
			pV5tg_tr[4] = fRandGaus(0,0.030/2);
		}

		if(iSmearX0==0) pXtg_BPM_tr=pV5tg_tr[0];
		else if(iSmearX0==1) pXtg_BPM_tr=pV5tg_tr[0]+fRandGaus(0,pBPMRes);

		//iArm = 0,1 means left, right
		if(iArm==0) iLeftArm=1;
		else if(iArm==1) iLeftArm=0;

		pRecDB=(iLeftArm)?pRecDBL:pRecDBR;

		for(int j=0;j<5;j++) pV5rec_tr[j]=pV5tg_tr[j];
        bool bGoodParticle;
        if (iSourceDistr<=2){
            bGoodParticle=SNAKEThruHRS(iLeftArm, pHRSAngle, pXtg_BPM_tr, pHRSMomentum, iFieldRotation, iExperiment, pV5rec_tr, pV5fp_tr);
        }
        else{
            int temp;
            fin >> temp >> pV5fp_tr[0] >> pV5fp_tr[1] >> pV5fp_tr[2] >> pV5fp_tr[3];
            pXtg_BPM_tr=0;
            bGoodParticle=kTRUE;
            ReconstructHRS(iLeftArm, pHRSAngle, pXtg_BPM_tr, pHRSMomentum, iFieldRotation, iExperiment, pV5rec_tr, pV5fp_tr);
        }
        
		//throw away this event if delta_rec>=1.0
		if(bGoodParticle && pV5rec_tr[4]<1.0) 
		{
			Transform::X_TCS2HCS(pV5tg_tr[0],pV5tg_tr[2],0.0,pHRSAngle,pV5tg_lab[0],pV5tg_lab[2],pV5tg_lab[4]);
			Transform::P_TCS2HCS(pV5tg_tr[1],pV5tg_tr[3],pHRSAngle,pV5tg_lab[1],pV5tg_lab[3]);

			Transform::X_TCS2HCS(pV5rec_tr[0],pV5rec_tr[2],0.0,pHRSAngle,pV5rec_lab[0],pV5rec_lab[2],pV5rec_lab[4]);
			Transform::P_TCS2HCS(pV5rec_tr[1],pV5rec_tr[3],pHRSAngle,pV5rec_lab[1],pV5rec_lab[3]);
            
			//Reconstruct use optics database
			pRecDB->CalcTargetCoords(pV5fp_tr,pV5rec_db_tr);

			pSnake->Fill();
			Index++;
		}
		NThrown++;
	}
    
	pFile->Write();
	cout<<"\n"<<Index<<"/"<<NThrown<<" events saved into ntuple"<<endl;
	pFile->Delete();
    if (iSourceDistr==3) fin.close();
}
