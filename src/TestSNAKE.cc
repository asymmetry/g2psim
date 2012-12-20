#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#include <iostream>
#include <fstream>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TF1.h"
#include "TPad.h"
#include "TStyle.h"

#include "TText.h"
#include "TPaveText.h"

#include "HRSRecUseDB.hh"

#include "HRSTransport.hh"
#include "HRSTransform_TCSNHCS.hh"

using namespace std;

//flat random number generator between [0,1)
double fRand();
double fRandGaus(double m=0.0, double s=1.0);
//return random number in [low,High) following a*x+c prob density
double fLinearRand(double a=1.0,double c=0.0,double low=0.0,double high=1.0);

typedef double (*func)(double);

double fPol1(double x)
{
	double a=1.0,b=0.0;
	return a*x+b;
}

double fPol2(double x)
{
	double a=1.0,b=0.0,c=0.0;
	return a*x*x+b*x+c;
}

//generate a non-uniform x using the constrain(distribution) of f(x)
double RandomOfFunction(func f, double low,double high)
{
	//input:
	//	low,high	  lower,higher limit of x
	double x,y;
	double ylow=f(low),yhigh=f(high);
	if(ylow<yhigh) 
	{
		double tmp=ylow;
		ylow=yhigh;
		yhigh=tmp;
	}

	do{
		x=low+(high-low)*((double)rand())/(double(RAND_MAX));
		y=ylow+(yhigh-ylow)*((double)rand())/(double(RAND_MAX));

		if(f(x)>y) break;
	}while(true);
	return x;
}

bool SNAKEThruHRS(int pIsLeftArm, double pEndPlaneAngle, double pXtg_BPM_tr, 
				  double pHRSMomentum, int iFieldRotation, int iExperiment,
				  double* pV5tg_tr, double* pV5fp_tr);
void ReconstructHRS(int pIsLeftArm, double pEndPlaneAngle, double pX0_tg_m, 
				  double pHRSMomentum, int iFieldRotation, int iExperiment,
                    double* pV5_tg, double* pV5_fp);

//input:
//iExperiment=10: g2p normal, 11: g2p 484816_shim; 12 g2p 403216_shim; 13 g2p 400016_shim 
//iSmearX0=0,1,2,3 means Xtg_BPM_tr=0, Xtg_BPM_tr, 1mm_gaussian_smeared_XBMP_tr, 0-2.5mm_gaussian_smeared_XBMP_tr
//iSourceDistr=0,1,2 means pV5tg_tr[5] is in delta, flat, gaussian distribution
//iArm = 0,1,2 means random, left, right
void TestSNAKE(int iNEvent, int iExperiment, int iSmearX0, int iSourceDistr, int iArm)
{	
	srand(time(0));
    ifstream fin;
    
    if(iSourceDistr==3){
        fin.open("input.dat");
    }
    
	const double deg=acos(0.0)/90;

	int    Index=0;
	int    iLeftArm=1;
	double pHRSAngle=5.682*deg;  //6 deg
	double pXtg_BPM_tr=0;
	double pHRSMomentum=2.251;
	int    iFieldRotation=90;  //90 deg	
	double pV5fp_tr[5]={0,0,0,0,0};
	double pV5tg_tr[5]={0,0,0,0,0};
	double pV5rec_tr[5]={0,0,0,0,0};
	double pV5tg_lab[5]={0,0,0,0,0};
	double pV5rec_lab[5]={0,0,0,0,0};
	double pBPMRes=0.0005;  // 0.5 mm 

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

void ReconstructHRS(int pIsLeftArm, double pEndPlaneAngle, double pX0_tg_m, 
				  double pHRSMomentum, int iFieldRotation, int iExperiment,
				  double* pV5_tg, double* pV5_fp)
{
    const double deg = acos(0.0)/90.0;
    int iSnakeFlag=(iExperiment%10);
    double pV5_rec[5];
    for(int k=0;k<4;k++)pV5_rec[k]=pV5_fp[k];
    pV5_rec[4]=pX0_tg_m;

    if( iExperiment>=1 && iExperiment<20) {
        //use the transportation without target field, which is delta dependent, not momentum dependent
        if(pIsLeftArm) {
            if(cos(pEndPlaneAngle)>cos(12.0*deg)) {//with septa
                if(iSnakeFlag==0) ReconstructLeftHRS_g2_(pV5_rec);  //from JLL, No shim, Wrong Bx
                else if(iSnakeFlag==1) ReconstructLeftHRS_Shim_484816(pV5_rec);    //By John, 3cm raster,
                else if(iSnakeFlag==2) ReconstructLeftHRS_Shim_403216(pV5_rec);
                else if(iSnakeFlag==3) ReconstructLeftHRS_Shim_400016(pV5_rec);
                else if(iSnakeFlag==9) ReconstructLeftHRS_Shim_484816_WrongBx(pV5_rec);   //using Min's 5.69 version, 2 cm raster and Wrong Bx
            }
            else {//no septa
                //John does not provide 12.5 degree resonstruction for g2p
                if(iSnakeFlag==9) ReconstructLeftHRS_12p5_Min_(pV5_rec);  //using Min's tuned reconstruction, this package need to be updated
                else ReconstructLeftHRS(pV5_rec);
            }
        }
        else {
            if(cos(pEndPlaneAngle)>cos(12.0*deg)) {
                if(iSnakeFlag==0) ReconstructRightHRS_g2_(pV5_rec);  //from JLL, No shim, Wrong Bx
                else if(iSnakeFlag==1) ReconstructRightHRS_Shim_484816(pV5_rec);    //By John, 3cm raster,
                else if(iSnakeFlag==2) ReconstructRightHRS_Shim_403216(pV5_rec);
                else if(iSnakeFlag==3) ReconstructRightHRS_Shim_400016(pV5_rec);
                else if(iSnakeFlag==9) ReconstructRightHRS_Shim_484816_WrongBx(pV5_rec);   //using Min's 5.69 version, 2 cm raster and Wrong Bx
            }
            else { //John does not provide 12.5 degree resonstruction for g2p
                if(iSnakeFlag==9) ReconstructRightHRS_12p5_Min_(pV5_rec);  //using Min's tuned reconstruction, this package need to be updated
                else ReconstructRightHRS(pV5_rec);
            }
        }
    }
    else if( iExperiment == 20 ) {
        if(pIsLeftArm) {
            if(pEndPlaneAngle/deg<12.0)	ReconstructLeftHRS6_LargeX0(pV5_rec);
            else ReconstructLeftHRS(pV5_rec);
        }
        else {
            if(pEndPlaneAngle/deg<12.0)	ReconstructRightHRS6_LargeX0(pV5_rec);
            else ReconstructRightHRS(pV5_rec);
        }
    }
    else {
        if(pIsLeftArm) {
            if(pEndPlaneAngle/deg<12.0)	ReconstructLeftHRS6(pV5_rec);
            else ReconstructLeftHRS(pV5_rec);
        }
        else {
            if(pEndPlaneAngle/deg<12.0)	ReconstructRightHRS6(pV5_rec);
            else ReconstructRightHRS(pV5_rec);
        }
    }
    
    for(int k=0;k<5;k++)  pV5_tg[k]=pV5_rec[k];
}
