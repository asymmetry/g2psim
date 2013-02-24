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

#include "G2PGun.hh"
#include "G2PRand.hh"
#include "HRSRecUseDB.hh"
#include "HRSTransTCSNHCS.hh"
#include "G2PDrift.hh"

#include "G2PXS.hh"
#include "HRSTransport.hh"

#include "G2PSim.hh"

//#define G2PSIM_DEBUG 1

using namespace std;

const double kDEG = 3.14159265358979323846/180.0;

ClassImp(G2PSim);

G2PSim::G2PSim()
    :bIsInit(false), bIsLeftArm(true), fHRSAngle(5.767*kDEG),
     fHRSMomentum(2.251), fBeamEnergy(2.254), bUseField(false),
     nEvent(10000), nIndex(1), pGun(NULL), iGunSetting(0), pHRS(NULL),
     bIsGoodParticle(false), pRecUseDB(NULL), pField(NULL),
     fEndPlaneZ(0), pPhys(NULL), fXS(0), pTree(NULL), pConfig(NULL),
     pFile(NULL), pFileName(NULL), pfRunSelector(NULL)
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

    if (bUseField) {
        pField->Init();
        if (!pField->IsInit()) noerror = false;
        else G2PDrift::SetField(pField);
    }

    vector<G2PGun*>::iterator it = pGunList.begin();
    while (it!=pGunList.end()){
        (*it)->SetHRSAngle(fHRSAngle);
        (*it)->SetHRSMomentum(fHRSMomentum);
        (*it)->SetBeamEnergy(fBeamEnergy);
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
    memset(fV5beam_lab, 0, sizeof(fV5beam_lab));
    memset(fV5bpm_lab, 0, sizeof(fV5bpm_lab));
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

    memset(fV5vb_tr, 0, sizeof(fV5recdb_lab));

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
        pConfig->Branch("BeamEnergy", &fBeamEnergy, "BeamEnergy/D");

        int hrssetting = pHRS->GetModelIndex();
        pConfig->Branch("HRSModel", &hrssetting, "HRSModel/I");

        int gunsetting = pGun->GetSetting();
        double posres = pGun->GetBPMPosRes();
        double angres = pGun->GetBPMAngRes();
        pConfig->Branch("GunSetting", &gunsetting, "GunSetting/I");
        pConfig->Branch("BPMPosRes", &posres, "BPMPosRes/D");
        pConfig->Branch("BPMAngRes", &angres, "BPMAngRes/D");

        pConfig->Fill();

        pFile->Write("", TObject::kOverwrite);

        delete pConfig;
        
        it++; i++;
    }
    
    pTree->Branch("Index", &nIndex,"Index/I");
    pTree->Branch("IsGood", &bIsGoodParticle, "IsGood/O");
    pTree->Branch("Gun", &iGunSetting, "Gun/I");

    pTree->Branch("Xbeam_lab", &fV5beam_lab[0],"Xbeam_lab/D");
    pTree->Branch("Thetabeam_lab", &fV5beam_lab[1],"Thetabeam_lab/D");
	pTree->Branch("Ybeam_lab", &fV5beam_lab[2],"Ybeam_lab/D");
    pTree->Branch("Phibeam_lab", &fV5beam_lab[3],"Phibeam_lab/D");
	pTree->Branch("Zbeam_lab", &fV5beam_lab[4],"Zbeam_lab/D");

    pTree->Branch("Xbpm_lab", &fV5bpm_lab[0],"Xbpm_lab/D");
    pTree->Branch("Thetabpm_lab", &fV5bpm_lab[1],"Thetabpm_lab/D");
	pTree->Branch("Ybpm_lab", &fV5bpm_lab[2],"Ybpm_lab/D");
    pTree->Branch("Phibpm_lab", &fV5bpm_lab[3],"Phibpm_lab/D");
	pTree->Branch("Zbpm_lab", &fV5bpm_lab[4],"Zbpm_lab/D");

    pTree->Branch("Xtg_tr",&fV5tg_tr[0],"Xtg_tr/D");
	pTree->Branch("Thetatg_tr",&fV5tg_tr[1],"Thetatg_tr/D");
	pTree->Branch("Ytg_tr",&fV5tg_tr[2],"Ytg_tr/D");
	pTree->Branch("Phitg_tr",&fV5tg_tr[3],"Phitg_tr/D");
	pTree->Branch("Delta",&fV5tg_tr[4],"Delta/D");

    pTree->Branch("Xtg_lab",&fV5tg_lab[0],"Xtg_lab/D");
	pTree->Branch("Thetatg_lab",&fV5tg_lab[1],"Thetatg_lab/D");
	pTree->Branch("Ytg_lab",&fV5tg_lab[2],"Ytg_lab/D");
	pTree->Branch("Phitg_lab",&fV5tg_lab[3],"Phitg_lab/D");

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

    pTree->Branch("Xrecdb_tr",&fV5recdb_tr[0],"Xrecdb_tr/D");
	pTree->Branch("Thetarecdb_tr",&fV5recdb_tr[1],"Thetarecdb_tr/D");
	pTree->Branch("Yrecdb_tr",&fV5recdb_tr[2],"Yrecdb_tr/D");
	pTree->Branch("Phirecdb_tr",&fV5recdb_tr[3],"Phirecdb_tr/D");
	pTree->Branch("Deltarecdb",&fV5recdb_tr[4],"Deltarecdb/D");

    pTree->Branch("Xrecdb_lab",&fV5recdb_lab[0],"Xrecdb_lab/D");
	pTree->Branch("Thetarecdb_lab",&fV5recdb_lab[1],"Thetarecdb_lab/D");
	pTree->Branch("Yrecdb_lab",&fV5recdb_lab[2],"Yrecdb_lab/D");
	pTree->Branch("Phirecdb_lab",&fV5recdb_lab[3],"Phirecdb_lab/D");

    pTree->Branch("XS", &fXS, "XS/D");
}

void G2PSim::RunSim()
{
    while (nIndex<=nEvent) {
        if (!pGun->Shoot(fV5beam_lab, fV5bpm_lab, fV5tg_tr)) break;

        double xtg_tr, ytg_tr, ztg_tr;
        HRSTransTCSNHCS::X_HCS2TCS(fV5beam_lab[0], fV5beam_lab[2], fV5beam_lab[4], fHRSAngle, xtg_tr, ytg_tr, ztg_tr);

        double xsave_tr;
        if (bUseField) {
            G2PDrift::Drift(fV5tg_tr, fHRSMomentum, fHRSAngle, ztg_tr, fEndPlaneZ, 10.0, fV5vb_tr);
            double virtualV5tg_tr[5];
            HRSTransTCSNHCS::Project(fV5vb_tr[0], fV5vb_tr[2], fEndPlaneZ, -fEndPlaneZ, fV5vb_tr[1], fV5vb_tr[3], virtualV5tg_tr[0], virtualV5tg_tr[2], virtualV5tg_tr[4]);
            virtualV5tg_tr[1] = fV5vb_tr[1];
            virtualV5tg_tr[3] = fV5vb_tr[3];
            virtualV5tg_tr[4] = fV5vb_tr[4];
            //printf("%e\t%e\n", fV5vb_tr[1], fV5vb_tr[3]);
            xsave_tr = virtualV5tg_tr[0];
            bIsGoodParticle = pHRS->Forward(virtualV5tg_tr, fV5fp_tr);
        }
        else {
            HRSTransTCSNHCS::Project(fV5tg_tr[0], fV5tg_tr[2], ztg_tr, -ztg_tr, fV5tg_tr[1], fV5tg_tr[3]);
            xsave_tr = fV5tg_tr[0];
            bIsGoodParticle = pHRS->Forward(fV5tg_tr, fV5fp_tr);
        }
        pRecUseDB->TransTr2Rot(fV5fp_tr, fV5fp_rot);

#ifdef G2PSIM_DEBUG
        printf("G2PSim: %e\t%e\t%e\t%e\t%e\n", fV5fp_rot[0], fV5fp_rot[1], fV5fp_rot[2], fV5fp_rot[3], fV5fp_rot[4]);
#endif
        VDCSmearing(fV5fp_tr);

        double xbpm_tr, ybpm_tr, zbpm_tr;
        HRSTransTCSNHCS::X_HCS2TCS(fV5bpm_lab[0], fV5bpm_lab[2], fV5bpm_lab[4], fHRSAngle, xbpm_tr, ybpm_tr, zbpm_tr);
        HRSTransTCSNHCS::Project(xbpm_tr, ybpm_tr, zbpm_tr, -zbpm_tr, fV5tg_tr[1], fV5tg_tr[3]);
        fV5fp_tr[4] = xbpm_tr;
        if (bUseField) {
            double xbpm_tr_eff = GetEffBPM(xbpm_tr, fHRSMomentum);
            fV5fp_tr[4] = xbpm_tr_eff;
            pHRS->Backward(fV5fp_tr, fV5rec_tr);
            double temprec=fHRSMomentum*(1.0+fV5rec_tr[4]);
            xbpm_tr_eff = GetEffBPM(xbpm_tr, temprec);
            fV5fp_tr[4] = xsave_tr;
            //printf("%e\t%e\n", xsave_tr, xbpm_tr_eff);
        }
        else {
            fV5fp_tr[4] = xsave_tr;
        }
        bIsGoodParticle &= pHRS->Backward(fV5fp_tr, fV5rec_tr);
        pRecUseDB->CalcTargetCoords(fV5fp_rot, fV5recdb_tr);

        if (bUseField) {
            HRSTransTCSNHCS::Project(fV5rec_tr[0], fV5rec_tr[2], 0.0, fEndPlaneZ, fV5rec_tr[1], fV5rec_tr[3], fV5vb_tr[0], fV5vb_tr[2], fV5vb_tr[4]);
            fV5vb_tr[1] = fV5rec_tr[1];
            fV5vb_tr[3] = fV5rec_tr[3];
            fV5vb_tr[4] = fV5rec_tr[4];
            G2PDrift::Drift(fV5vb_tr, fHRSMomentum, fHRSAngle, fEndPlaneZ, 0.0, 10.0, fV5rec_tr);
            double reactz = DriftPath();
            //double reactz = -13.6271e-3;
            //printf("%e\n", reactz);
            G2PDrift::Drift(fV5rec_tr, fHRSMomentum, fHRSAngle, 0.0, reactz, 10.0, fV5rec_tr);
        }

#ifdef G2PSIM_DEBUG
        bIsGoodParticle = true;
        printf("G2PSim: %e\t%e\t%e\t%e\t%e\n", fV5rec_tr[0], fV5rec_tr[1], fV5rec_tr[2], fV5rec_tr[3], fV5rec_tr[4]);
        printf("G2PSim: %e\t%e\t%e\t%e\t%e\n", fV5recdb_tr[0], fV5recdb_tr[1], fV5recdb_tr[2], fV5recdb_tr[3], fV5recdb_tr[4]);
#endif

        pTree->Fill();

        if ((nIndex%1000)==0) printf("%d\n", nIndex);
        nIndex++;
    }
}

void G2PSim::RunData()
{
    while (nIndex<=nEvent) {
        if (!pGun->Shoot(fV5beam_lab, fV5bpm_lab, fV5tg_tr)) break;
        pGun->GetFP(fV5fpdata_tr);

        fV5fpdata_tr[4] = fV5tg_tr[0];

        pRecUseDB->TransTr2Rot(fV5fpdata_tr, fV5fpdata_rot);
        pRecUseDB->CalcTargetCoords(fV5fpdata_rot, fV5recdb_tr);

        fV5recdb_tr[0] = fV5tg_tr[0];
        bIsGoodParticle = pHRS->Forward(fV5recdb_tr, fV5fp_tr);
        pRecUseDB->TransTr2Rot(fV5fp_tr, fV5fp_rot);

#ifdef G2PSIM_DEBUG
        printf("G2PSim: %e\t%e\t%e\t%e\t%e\n", fV5fp_rot[0], fV5fp_rot[1], fV5fp_rot[2], fV5fp_rot[3], fV5fp_rot[4]);
#endif

        bIsGoodParticle &= pHRS->Backward(fV5fpdata_tr, fV5rec_tr);
        
#ifdef G2PSIM_DEBUG
        bIsGoodParticle = true;
        printf("G2PSim: %e\t%e\t%e\t%e\t%e\n", fV5rec_tr[0], fV5rec_tr[1], fV5rec_tr[2], fV5rec_tr[3], fV5rec_tr[4]);
        printf("G2PSim: %e\t%e\t%e\t%e\t%e\n", fV5recdb_tr[0], fV5recdb_tr[1], fV5recdb_tr[2], fV5recdb_tr[3], fV5recdb_tr[4]);
#endif

        pTree->Fill();

        if ((nIndex%1000)==0) printf("%d\n", nIndex);
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

double G2PSim::GetEffBPM(double xbpm_tr, double p)
{
	//Apply the correction to get the effective X_BPM_tr at target plane
	//the following is the result of fitting "[0]/x+[1]" to 
	//(X_proj2tg_tr - Xtg_tr) vs Pvb @ 5.0T target field
	//need to fit for 2.5Ttarget field later
	// EXT PARAMETER                                   STEP         FIRST   
	//NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
	// 1  p0          -6.39438e+00   5.37354e-03   1.74320e-04  -2.30577e-05
	// 2  p1           8.06795e-04   4.92714e-03   1.02122e-05   6.92299e-07
	double xbpm_tr_eff=xbpm_tr;	
	if(bUseField)
    {
		xbpm_tr_eff += (-6.39438/p + 8.06795e-04)/1000 * pField->GetRatio();
	}
	return xbpm_tr_eff;
}


double G2PSim::DriftPath()
{   
    if (bUseField) {
        double x[3] = { fV5bpm_lab[0], fV5bpm_lab[2], fV5bpm_lab[4] };
        double p[3] = { fBeamEnergy*sin(fV5bpm_lab[1])*cos(fV5bpm_lab[3]),
                        fBeamEnergy*sin(fV5bpm_lab[1])*sin(fV5bpm_lab[3]),
                        fBeamEnergy*cos(fV5bpm_lab[1]) };

        double xrec[3];
        HRSTransTCSNHCS::X_TCS2HCS(fV5rec_tr[0], fV5rec_tr[2], 0.0, fHRSAngle, xrec[0], xrec[1], xrec[2]);
        double theta,phi;
        HRSTransTCSNHCS::P_TCS2HCS(fV5rec_tr[1], fV5rec_tr[3], fHRSAngle, theta, phi);
        double prec[3] = { fHRSMomentum*(1+fV5rec_tr[4])*sin(theta)*cos(phi),
                           fHRSMomentum*(1+fV5rec_tr[4])*sin(theta)*sin(phi),
                           fHRSMomentum*(1+fV5rec_tr[4])*cos(theta) };
        G2PDrift::Drift(xrec, prec, 0, 10.0, xrec, prec);

        double z = 0.0;
        double zmin = 0.0;
        double distmin = sqrt((x[0]-xrec[0])*(x[0]-xrec[0])+(x[1]-xrec[1])*(x[1]-xrec[1]));
        for (int i = 1; i<150; i++) {
            z = 0.1e-3*i;
            G2PDrift::Drift(x, p, z, 10.0, x, p);
            G2PDrift::Drift(xrec, prec, z, 10.0, xrec, prec);
            //printf("%e\t%e\t%e\t%e\t%e\n", z, x[0], x[1], xrec[0], xrec[1]);
            double distance = sqrt((x[0]-xrec[0])*(x[0]-xrec[0])+(x[1]-xrec[1])*(x[1]-xrec[1]));
            if (distance<distmin) {
                zmin = z;
                distmin = distance;
            }
        }

        // for (int i = 20; i<200; i++) {
        //     z = 1e-3*i;
        //     G2PDrift::Drift(xrec, prec, z, 10.0, xrec, prec);
        //     printf("%e\t%e\t%e\t%e\t%e\n", z, 1000.0, 1000.0, xrec[0], xrec[1]);
        // }

        for (int i = 1; i<150; i++) {
            z = -0.1e-3*i;
            G2PDrift::Drift(x, p, z, 10.0, x, p);
            G2PDrift::Drift(xrec, prec, z, 10.0, xrec, prec);
            //printf("%e\t%e\t%e\t%e\t%e\n", z, x[0], x[1], xrec[0], xrec[1]);
            double distance = sqrt((x[0]-xrec[0])*(x[0]-xrec[0])+(x[1]-xrec[1])*(x[1]-xrec[1]));
            if (distance<distmin) {
                zmin = z;
                distmin = distance;
            }
        }

        // for (int i = -20; i>-200; i--) {
        //     z = 1e-3*i;
        //     G2PDrift::Drift(x, p, z, 10.0, x, p);
        //     printf("%e\t%e\t%e\t%e\t%e\n", z, x[0], x[1], 1000.0, 1000.0);
        // }

        G2PDrift::Drift(x, p, zmin, 10.0, x, p);
        G2PDrift::Drift(xrec, prec, zmin, 10.0, xrec, prec);

        double x_tr,y_tr,z_tr;
        HRSTransTCSNHCS::X_HCS2TCS(xrec[0], xrec[1], xrec[2], fHRSAngle, x_tr, y_tr, z_tr);
        return (z_tr);
    }
    else return 0;
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
