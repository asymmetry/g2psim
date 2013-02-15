// This file defines a class G2PSim.
// This class is the main class for this sim package.
//
// History:
//   Jan 2013, C. Gu, First public version.
//

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cmath>

#include "TROOT.h"
#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

#include "G2PGun.hh"
#include "G2PRand.hh"
#include "HRSRecUseDB.hh"
#include "HRSTransTCSNHCS.hh"

#include "G2PXS.hh"
#include "HRSTransport.hh"

#include "G2PSim.hh"

//#define G2PSIM_DEBUG 1

using namespace std;

const double kDEG = TMath::Pi()/180.0;

ClassImp(G2PSim);

G2PSim::G2PSim()
    :bIsInit(false), pFile(NULL), pFileName(NULL), nIndex(1),
     nEvent(10000), bIsLeftArm(true), fHRSAngle(5.767*kDEG),
     fHRSMomentum(2.251), pGun(NULL),  pHRS(NULL), pRecUseDB(NULL),
     pPhys(NULL), pTree(NULL), pConfig(NULL), pfRunSelector(NULL)
{
    pGunList.clear();
    
    Clear();
}

G2PSim::~G2PSim()
{
    // Nothing to do
}

void G2PSim::Init()
{
    bool noerror = true;

    vector<G2PGun*>::iterator it = pGunList.begin();
    while (it!=pGunList.end()){
        (*it)->SetHRSAngle(fHRSAngle);
        (*it)->SetHRSMomentum(fHRSMomentum);
        (*it)->Init();
        if (!(*it)->IsInit()) noerror = false;
        it++;
    }

    pHRS->SetArm(bIsLeftArm);
    pHRS->SetHRSAngle(fHRSAngle);

    if (bIsLeftArm)
        pRecUseDB = new HRSRecUseDB("L","db_L.vdc.dat");
    else
        pRecUseDB = new HRSRecUseDB("R","db_R.vdc.dat");

    if (noerror) InitTree();

    bIsInit = noerror;
}

void G2PSim::Run()
{
    if (pGunList.empty()) { printf("No particle gun!\n"); exit(-1); }
    if (pHRS==NULL) { printf("No HRS Model!\n"); exit(-1); }

    Init();

    if (bIsInit) {
        vector<G2PGun*>::iterator it = pGunList.begin();
        while (it!=pGunList.end()){
            Clear();
            pGun = (*it);
            iGunSetting=pGun->GetSetting();
            if (pGun->IsUsingData())
                pfRunSelector = &G2PSim::RunData;
            else
                pfRunSelector = &G2PSim::RunSim;
            nIndex = 1;
            (this->*pfRunSelector)();
            it++;
        }
        End();
    }
}

void G2PSim::End()
{
    pFile->Write("", TObject::kOverwrite);
    pFile->Close();
    
    delete pRecUseDB;
}

void G2PSim::Clear()
{
    memset(fV3bpm_lab, 0, sizeof(fV3bpm_lab));
    memset(fV5tg_tr, 0, sizeof(fV5tg_tr));
    memset(fV5tg_lab, 0, sizeof(fV5tg_lab));
    memset(fV5fpdata_tr, 0, sizeof(fV5fpdata_tr));
    memset(fV5fpdata_rot, 0, sizeof(fV5fpdata_rot));

    bIsGoodParticle = false;
    memset(fV5fp_tr, 0, sizeof(fV5fp_tr));
    memset(fV5fp_rot, 0, sizeof(fV5fp_rot));
    memset(fV5rec_tr, 0, sizeof(fV5rec_tr));
    memset(fV5rec_lab, 0, sizeof(fV5rec_lab));
    
    memset(fV5recdb_tr, 0, sizeof(fV5recdb_tr));
    memset(fV5recdb_lab, 0, sizeof(fV5recdb_lab));

    fXS = 1;
}

void G2PSim::InitTree()
{
    pFile = new TFile(pFileName, "recreate");
	pTree = new TTree("T", "sim result");

    vector<G2PGun*>::iterator it = pGunList.begin();
    int i = 0;
    while (it!=pGunList.end()){
        pGun = (*it);
        pConfig = new TTree(Form("config%d", i), "sim configure");

        pConfig->Branch("NEvent", &nEvent, "NEvent/I");

        pConfig->Branch("IsLeftArm", &bIsLeftArm, "IsLeftArm/O");
        pConfig->Branch("HRSAngle", &fHRSAngle, "HRSAngle/D");
        pConfig->Branch("HRSMomentum", &fHRSMomentum, "HRSMomentum/D");

        int hrssetting = pHRS->GetModelIndex();
        pConfig->Branch("HRSModel", &hrssetting, "HRSModel/I");

        int gunsetting = pGun->GetSetting();
        double posres = pGun->GetPosResolution();
        double angleres = pGun->GetAngleResolution();
        double deltares = pGun->GetDeltaResolution();
        pConfig->Branch("GunSetting", &gunsetting, "GunSetting/I");
        pConfig->Branch("GunPosRes", &posres, "GunPosRes/D");
        pConfig->Branch("GunAngleRes", &angleres, "GunAngleRes/D");
        pConfig->Branch("GunDeltaRes", &deltares, "GunDeltaRes/D");

        pConfig->Fill();

        pFile->Write("", TObject::kOverwrite);

        delete pConfig;
        
        it++; i++;
    }
    
    pTree->Branch("Index", &nIndex,"Index/I");
    pTree->Branch("IsGood", &bIsGoodParticle, "IsGood/O");
    pTree->Branch("Gun", &iGunSetting, "Gun/I");

    pTree->Branch("Xbpm_lab", &fV3bpm_lab[0],"Xbpm_lab/D");
	pTree->Branch("Ybpm_lab", &fV3bpm_lab[1],"Ybpm_lab/D");
	pTree->Branch("Zbpm_lab", &fV3bpm_lab[2],"Zbpm_lab/D");

    pTree->Branch("Xtg_tr",&fV5tg_tr[0],"Xtg_tr/D");
	pTree->Branch("Thetatg_tr",&fV5tg_tr[1],"Thetatg_tr/D");
	pTree->Branch("Ytg_tr",&fV5tg_tr[2],"Ytg_tr/D");
	pTree->Branch("Phitg_tr",&fV5tg_tr[3],"Phitg_tr/D");
	pTree->Branch("Delta",&fV5tg_tr[4],"Delta/D");

    pTree->Branch("Xtg_lab",&fV5tg_lab[0],"Xtg_lab/D");
	pTree->Branch("Thetatg_lab",&fV5tg_lab[1],"Thetatg_lab/D");
	pTree->Branch("Ytg_lab",&fV5tg_lab[2],"Ytg_lab/D");
	pTree->Branch("Phitg_lab",&fV5tg_lab[3],"Phitg_lab/D");
	pTree->Branch("Ztg_lab",&fV5tg_lab[4],"Ztg_lab/D");

    pTree->Branch("Xfpdata_tr",&fV5fpdata_tr[0],"Xfpdata_tr/D");
	pTree->Branch("Thetafpdata_tr",&fV5fpdata_tr[1],"Thetafpdata_tr/D");
	pTree->Branch("Yfpdata_tr",&fV5fpdata_tr[2],"Yfpdata_tr/D");
	pTree->Branch("Phifpdata_tr",&fV5fpdata_tr[3],"Phifpdata_tr/D");

    pTree->Branch("Xfpdata_rot",&fV5fpdata_rot[0],"Xfpdata_rot/D");
	pTree->Branch("Thetafpdata_rot",&fV5fpdata_rot[1],"Thetafpdata_rot/D");
	pTree->Branch("Yfpdata_rot",&fV5fpdata_rot[2],"Yfpdata_rot/D");
	pTree->Branch("Phifpdata_rot",&fV5fpdata_rot[3],"Phifpdata_rot/D");
    
	pTree->Branch("Xfp_tr",&fV5fp_tr[0],"Xfp_tr/D");
	pTree->Branch("Thetafp_tr",&fV5fp_tr[1],"Thetafp_tr/D");
	pTree->Branch("Yfp_tr",&fV5fp_tr[2],"Yfp_tr/D");
	pTree->Branch("Phifp_tr",&fV5fp_tr[3],"Phifp_tr/D");

    pTree->Branch("Xfp_rot",&fV5fp_rot[0],"Xfp_rot/D");
	pTree->Branch("Thetafp_rot",&fV5fp_rot[1],"Thetafp_rot/D");
	pTree->Branch("Yfp_rot",&fV5fp_rot[2],"Yfp_rot/D");
	pTree->Branch("Phifp_rot",&fV5fp_rot[3],"Phifp_rot/D");

    pTree->Branch("Xrec_tr",&fV5rec_tr[0],"Xrec_tr/D");
	pTree->Branch("Thetarec_tr",&fV5rec_tr[1],"Thetarec_tr/D");
	pTree->Branch("Yrec_tr",&fV5rec_tr[2],"Yrec_tr/D");
	pTree->Branch("Phirec_tr",&fV5rec_tr[3],"Phirec_tr/D");
	pTree->Branch("Deltarec",&fV5rec_tr[4],"Deltarec/D");

    pTree->Branch("Xrec_lab",&fV5rec_lab[0],"Xrec_lab/D");
	pTree->Branch("Thetarec_lab",&fV5rec_lab[1],"Thetarec_lab/D");
	pTree->Branch("Yrec_lab",&fV5rec_lab[2],"Yrec_lab/D");
	pTree->Branch("Phirec_lab",&fV5rec_lab[3],"Phirec_lab/D");
	pTree->Branch("Zrec_lab",&fV5tg_lab[4],"Zrec_lab/D");

    pTree->Branch("Xrecdb_tr",&fV5recdb_tr[0],"Xrecdb_tr/D");
	pTree->Branch("Thetarecdb_tr",&fV5recdb_tr[1],"Thetarecdb_tr/D");
	pTree->Branch("Yrecdb_tr",&fV5recdb_tr[2],"Yrecdb_tr/D");
	pTree->Branch("Phirecdb_tr",&fV5recdb_tr[3],"Phirecdb_tr/D");
	pTree->Branch("Deltarecdb",&fV5recdb_tr[4],"Deltarecdb/D");

    pTree->Branch("Xrecdb_lab",&fV5recdb_lab[0],"Xrecdb_lab/D");
	pTree->Branch("Thetarecdb_lab",&fV5recdb_lab[1],"Thetarecdb_lab/D");
	pTree->Branch("Yrecdb_lab",&fV5recdb_lab[2],"Yrecdb_lab/D");
	pTree->Branch("Phirecdb_lab",&fV5recdb_lab[3],"Phirecdb_lab/D");
    pTree->Branch("Zrecdb_lab",&fV5recdb_lab[3],"Zrecdb_lab/D");

    pTree->Branch("XS", &fXS, "XS/D");
}

void G2PSim::RunSim()
{
    while (nIndex<=nEvent) {
        if (!pGun->Shoot(fV3bpm_lab, fV5tg_tr)) break;

        fV5tg_tr[1] = tan(fV5tg_tr[1]);
        fV5tg_tr[3] = tan(fV5tg_tr[3]);

        bIsGoodParticle = pHRS->Forward(fV5tg_tr, fV5fp_tr);
        pRecUseDB->TransTr2Rot(fV5fp_tr, fV5fp_rot);

#ifdef G2PSIM_DEBUG
        printf("G2PSim: %e\t%e\t%e\t%e\t%e\n", fV5fp_rot[0], fV5fp_rot[1], fV5fp_rot[2], fV5fp_rot[3], fV5fp_rot[4]);
#endif

        fV5fp_tr[4] = fV5tg_tr[0];
        bIsGoodParticle &= pHRS->Backward(fV5fp_tr, fV5rec_tr);
        pRecUseDB->CalcTargetCoords(fV5fp_rot, fV5recdb_tr);

#ifdef G2PSIM_DEBUG
        bIsGoodParticle = true;
        printf("G2PSim: %e\t%e\t%e\t%e\t%e\n", fV5rec_tr[0], fV5rec_tr[1], fV5rec_tr[2], fV5rec_tr[3], fV5rec_tr[4]);
        printf("G2PSim: %e\t%e\t%e\t%e\t%e\n", fV5recdb_tr[0], fV5recdb_tr[1], fV5recdb_tr[2], fV5recdb_tr[3], fV5recdb_tr[4]);
#endif

        pTree->Fill();

        if ((nIndex%10000)==0) printf("%d\n", nIndex);
        nIndex++;
    }
}

void G2PSim::RunData()
{
    while (nIndex<=nEvent) {
        if (!pGun->Shoot(fV3bpm_lab, fV5fpdata_tr)) break;

        double pV3[3];
        HRSTransTCSNHCS::X_HCS2TCS(fV3bpm_lab[0], fV3bpm_lab[1], fV3bpm_lab[2], fHRSAngle, pV3[0], pV3[1], pV3[2]);

        pRecUseDB->TransTr2Rot(fV5fpdata_tr, fV5fpdata_rot);
        pRecUseDB->CalcTargetCoords(fV5fpdata_rot, fV5tg_tr);

        HRSTransTCSNHCS::Project(pV3[0], pV3[1], pV3[2], -pV3[2], fV5tg_tr[1], fV5tg_tr[3]);
        fV5fpdata_tr[4] = pV3[0];
        fV5tg_tr[0] = pV3[0];

        bIsGoodParticle = pHRS->Forward(fV5tg_tr, fV5fp_tr);
        pRecUseDB->TransTr2Rot(fV5fp_tr, fV5fp_rot);

#ifdef G2PSIM_DEBUG
        printf("G2PSim: %e\t%e\t%e\t%e\t%e\n", fV5fp_rot[0], fV5fp_rot[1], fV5fp_rot[2], fV5fp_rot[3], fV5fp_rot[4]);
#endif

        bIsGoodParticle &= pHRS->Backward(fV5fpdata_tr, fV5rec_tr);
        pRecUseDB->CalcTargetCoords(fV5fpdata_rot, fV5recdb_tr);
        
#ifdef G2PSIM_DEBUG
        bIsGoodParticle = true;
        printf("G2PSim: %e\t%e\t%e\t%e\t%e\n", fV5rec_tr[0], fV5rec_tr[1], fV5rec_tr[2], fV5rec_tr[3], fV5rec_tr[4]);
        printf("G2PSim: %e\t%e\t%e\t%e\t%e\n", fV5recdb_tr[0], fV5recdb_tr[1], fV5recdb_tr[2], fV5recdb_tr[3], fV5recdb_tr[4]);
#endif

        pTree->Fill();

        if ((nIndex%10000)==0) printf("%d\n", nIndex);
        nIndex++;
    }
}

void G2PSim::VDCSmearing(double* V5_fp)
{
    double WireChamberResX = 0.0013; //m;
    double WireChamberResY = 0.0013; //m;
    double WireChamberResT = 0.0003; //rad;
    double WireChamberResP = 0.0003; //rad;

    V5_fp[0] += G2PRand::Gaus(0, WireChamberResX);
    V5_fp[2] += G2PRand::Gaus(0, WireChamberResY);
    V5_fp[1] += G2PRand::Gaus(0, WireChamberResT);
    V5_fp[3] += G2PRand::Gaus(0, WireChamberResP);
}

        //throw away this event if delta_rec>=1.0
		// {
		// 	HRSTransTCSNHCS::X_TCS2HCS(fV5tg_tr[0],fV5tg_tr[2],0.0,fHRSAngle,fV5tg_lab[0],fV5tg_lab[2],fV5tg_lab[4]);
		// 	HRSTransTCSNHCS::P_TCS2HCS(fV5tg_tr[1],fV5tg_tr[3],fHRSAngle,fV5tg_lab[1],fV5tg_lab[3]);

		// 	HRSTransTCSNHCS::X_TCS2HCS(fV5rec_tr[0],fV5rec_tr[2],0.0,fHRSAngle,fV5rec_lab[0],fV5rec_lab[2],fV5rec_lab[4]);
		// 	HRSTransTCSNHCS::P_TCS2HCS(fV5rec_tr[1],fV5rec_tr[3],fHRSAngle,fV5rec_lab[1],fV5rec_lab[3]);
            
		// 	//Reconstruct use optics database
		// 	pRecUseDB->CalcTargetCoords(fV5fp_tr,fV5rec_db_tr);

		// 	pTree->Fill();
		// 	nIndex++;
		// }
