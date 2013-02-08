#include <cstdio>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TVirtualFitter.h"

TFile* f;
TTree* T;

Bool_t IsGood;
Int_t Gun;
Double_t Xfp_tr,Tfp_tr,Yfp_tr,Pfp_tr;
Double_t Xfpd_tr,Tfpd_tr,Yfpd_tr,Pfpd_tr;

Int_t NPara = 4;
Double_t startvalues[4] = { 0.0, 0.0, 0.0, 0.0 };
Double_t parerrors[4] = { 0.1, 0.1, 0.1, 0.1 };
Double_t parlimitslo[4] = { -0.1, -0.1, -0.1, -0.1 };
Double_t parlimitshi[4] = { 0.1, 0.1, 0.1, 0.1 };
Bool_t free[4] = { kTRUE, kTRUE, kTRUE, kTRUE };

Int_t iter = 1;

// Zero order calculation
void myfcn0(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    Int_t N = T->GetEntries();
    Double_t result[4] = { 0, 0, 0, 0 };
    Int_t NC = 1;

    printf("myfcn0: iter %d\n", iter++);

    printf("myfcn0: pars: %e\t%e\t%e\t%e\n", par[0], par[1], par[2], par[3]);

    for (Int_t i = 0; i<N; i++){
        T->GetEntry(i);
        if (IsGood&&Gun==6){
            result[0] += (Xfpd_tr-(Xfp_tr+par[0]))*(Xfpd_tr-(Xfp_tr+par[0]));
            result[1] += (Tfpd_tr-(Tfp_tr+par[1]))*(Tfpd_tr-(Tfp_tr+par[1]));
            result[2] += (Yfpd_tr-(Yfp_tr+par[2]))*(Yfpd_tr-(Yfp_tr+par[2]));
            result[3] += (Pfpd_tr-(Pfp_tr+par[3]))*(Pfpd_tr-(Pfp_tr+par[3]));
            NC++;
        }
    }
    Double_t f = (result[0]+result[1]+result[2]+result[3])/(NC-1);
    printf("myfcn0: NC = %d\t XCN = %e\t \n", NC-1, f);
}

void LoadTree(const Char_t* filename = "test.root"){
    // Will add TChain later
    f = new TFile(filename);
    T = (TTree *)f->Get("T");

    T->SetBranchAddress("IsGood", &IsGood);
    T->SetBranchAddress("Gun", &Gun);
    T->SetBranchAddress("Xfp_tr", &Xfp_tr);
    T->SetBranchAddress("Thetafp_tr", &Tfp_tr);
    T->SetBranchAddress("Yfp_tr", &Yfp_tr);
    T->SetBranchAddress("Phifp_tr", &Pfp_tr);
    T->SetBranchAddress("Xfpdata_tr", &Xfpd_tr);
    T->SetBranchAddress("Thetafpdata_tr", &Tfpd_tr);
    T->SetBranchAddress("Yfpdata_tr", &Yfpd_tr);
    T->SetBranchAddress("Phifpdata_tr", &Pfpd_tr);
}

void FitCor(const Char_t* filename = "test.root") 
{
    LoadTree(filename);
    
	TVirtualFitter::SetDefaultFitter("Minuit2");  //default is Minuit
	TVirtualFitter *fitter = TVirtualFitter::Fitter(NULL, NPara);
	fitter->SetFCN(myfcn0);

    for(Int_t i = 0; i<NPara; i++)
	{
        fitter->SetParameter(i, Form("Para%02d",i), startvalues[i], parerrors[i], parlimitslo[i], parlimitshi[i]);

		if (!free[i]) fitter->FixParameter(i);
	}

    fitter->Print();

    Double_t arglist[1] = {0};
    iter = 1;
	fitter->ExecuteCommand("MIGRAD", arglist, 0);
}
