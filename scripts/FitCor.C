#include <cstdio>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TVirtualFitter.h"

#define ORDER 1

TFile* f;
TTree* T;

Bool_t IsGood;
Int_t Gun;
Double_t Xfp_tr, Tfp_tr, Yfp_tr, Pfp_tr;
Double_t Xfpd_tr, Tfpd_tr, Yfpd_tr, Pfpd_tr;
Double_t Xtg_tr, Ttg_tr, Ytg_tr, Ptg_tr, Dtg_tr;

Double_t scale[4] = { 0.04, 0.04, 0.08, 0.04 };

#if ORDER == 0
Int_t NPara = 4;
Double_t startvalues[4] = { 0.0, 0.0, 0.0, 0.0 };
Double_t parerrors[4] = { 0.1, 0.1, 0.1, 0.1 };
Double_t parlimitslo[4] = { -0.1, -0.1, -0.1, -0.1 };
Double_t parlimitshi[4] = { 0.1, 0.1, 0.1, 0.1 };
Bool_t free[4] = { kTRUE, kTRUE, kTRUE, kTRUE };
#endif

#if ORDER == 1
Int_t NPara = 24;
// 484816R00
Double_t startvalues[24] = {  -5.91123e-3,  0.0,  0.0,  0.0,  0.0,  0.0,
                              -6.76146e-5,  0.0,  0.0,  0.0,  0.0,  0.0,
                              -1.89684e-2,  0.0,  0.0,  0.0,  0.0,  0.0,
                              -1.11524e-2,  0.0,  0.0,  0.0,  0.0,  0.0 };
// 484816
// Double_t startvalues[24] = {  1.99661e-3,  0.0,  0.0,  0.0,  0.0,  0.0,
//                               6.63373e-4,  0.0,  0.0,  0.0,  0.0,  0.0,
//                               -2.41468e-2,  0.0,  0.0,  0.0,  0.0,  0.0,
//                               -1.48819e-2,  0.0,  0.0,  0.0,  0.0,  0.0 };
// 484816 corr term
// Double_t startvalues[24] = {  1.99661e-3,  0.0,  0.0,  0.0,  0.0,  0.0,
//                               6.63373e-4,  0.0,  0.0,  0.0,  0.0,  0.0,
//                               -2.41468e-2,  0.0,  0.0,  0.0,  0.0,  0.0,
//                               -1.48819e-2,  0.0,  0.0,  0.0,  0.0,  0.0 };
Double_t parerrors[24] =   { 1.0e-4, 1.0e-4, 1.0e-4, 1.0e-4, 1.0e-4, 1.0e-4,
                             1.0e-4, 1.0e-4, 1.0e-4, 1.0e-4, 1.0e-4, 1.0e-4,
                             1.0e-4, 1.0e-4, 1.0e-4, 1.0e-4, 1.0e-4, 1.0e-4,
                             1.0e-4, 1.0e-4, 1.0e-4, 1.0e-4, 1.0e-4, 1.0e-4 };
Double_t parlimitslo[24] = { -0.1, -0.1, -0.1, -0.1, -0.1, -0.1,
                             -0.1, -0.1, -0.1, -0.1, -0.1, -0.1,
                             -0.1, -0.1, -0.1, -0.1, -0.1, -0.1,
                             -0.1, -0.1, -0.1, -0.1, -0.1, -0.1 };
Double_t parlimitshi[24] = {  0.1,  0.1,  0.1,  0.1,  0.1,  0.1,
                              0.1,  0.1,  0.1,  0.1,  0.1,  0.1,
                              0.1,  0.1,  0.1,  0.1,  0.1,  0.1,
                              0.1,  0.1,  0.1,  0.1,  0.1,  0.1 };
Double_t free[24] = { kFALSE,  kTRUE,  kTRUE, kFALSE, kFALSE,  kTRUE,
                      kFALSE,  kTRUE,  kTRUE, kFALSE, kFALSE,  kTRUE,
                      kFALSE, kFALSE, kFALSE,  kTRUE,  kTRUE,  kTRUE,
                      kFALSE, kFALSE, kFALSE,  kTRUE,  kTRUE,  kTRUE };
// Double_t free[24] = { kFALSE,  kTRUE,  kTRUE, kFALSE, kFALSE,  kTRUE,
//                       kFALSE,  kTRUE,  kTRUE, kFALSE, kFALSE,  kTRUE,
//                       kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
//                       kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE };
#endif

Int_t iter = 1;

// Zero order calculation
void myfcn0(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    Int_t N = T->GetEntries();
    Double_t result[4] = { 0, 0, 0, 0 };
    Double_t Xfpc_tr, Tfpc_tr, Yfpc_tr, Pfpc_tr;
    Int_t NC = 1;

    printf("myfcn0: iter %d\n", iter++);

    printf("myfcn0: pars: %e\t%e\t%e\t%e\n", par[0], par[1], par[2], par[3]);

    for (Int_t i = 0; i<N; i++){
        T->GetEntry(i);
        if (IsGood&&Gun==6){
            Xfpc_tr = Xfpd_tr-(Xfp_tr+par[0]);
            Tfpc_tr = Tfpd_tr-(Tfp_tr+par[1]);
            Yfpc_tr = Yfpd_tr-(Yfp_tr+par[2]);
            Pfpc_tr = Pfpd_tr-(Pfp_tr+par[3]);

            result[0] += Xfpc_tr*Xfpc_tr;
            result[1] += Tfpc_tr*Tfpc_tr;
            result[2] += Yfpc_tr*Yfpc_tr;
            result[3] += Pfpc_tr*Pfpc_tr;
            
            NC++;
        }
    }

    f = ( 0.0
        +result[0]/(scale[0]*scale[0])
        +result[1]/(scale[1]*scale[1])
        +result[2]/(scale[2]*scale[2])
        +result[3]/(scale[3]*scale[3])
        )/(NC-1);

    printf("myfcn0: NC = %d\t XCN = %e\t \n", NC-1, f);
}

// First order calculation
void myfcn1(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    Int_t N = T->GetEntries();
    Double_t result[4] = { 0, 0, 0, 0 };
    Double_t Xfpc_tr, Tfpc_tr, Yfpc_tr, Pfpc_tr;
    Int_t NC = 1;

    printf("myfcn0: iter %d\n", iter++);
    
    printf("myfcn0: pars:\n");

    for (Int_t i = 0; i<4; i++) {
        for (Int_t j = 0; j<6; j++) {
            printf("\t%e", par[i*6+j]);
        }
        printf("\n");
    }

    for (Int_t i = 0; i<N; i++){
        T->GetEntry(i);
        if (IsGood&&Gun==6){
            Xfpc_tr = Xfpd_tr-(Xfp_tr+par[ 0]+Xtg_tr*par[ 1]+Ttg_tr*par[ 2]+Ytg_tr*par[ 3]+Ptg_tr*par[ 4]+Dtg_tr*par[ 5]);
            Tfpc_tr = Tfpd_tr-(Tfp_tr+par[ 6]+Xtg_tr*par[ 7]+Ttg_tr*par[ 8]+Ytg_tr*par[ 9]+Ptg_tr*par[10]+Dtg_tr*par[11]);
            Yfpc_tr = Yfpd_tr-(Yfp_tr+par[12]+Xtg_tr*par[13]+Ttg_tr*par[14]+Ytg_tr*par[15]+Ptg_tr*par[16]+Dtg_tr*par[17]);
            Pfpc_tr = Pfpd_tr-(Pfp_tr+par[18]+Xtg_tr*par[19]+Ttg_tr*par[20]+Ytg_tr*par[21]+Ptg_tr*par[22]+Dtg_tr*par[23]);

            result[0] += Xfpc_tr*Xfpc_tr;
            result[1] += Tfpc_tr*Tfpc_tr;
            result[2] += Yfpc_tr*Yfpc_tr;
            result[3] += Pfpc_tr*Pfpc_tr;
            
            NC++;
        }
    }
    
    f = ( 0.0
        +result[0]/(scale[0]*scale[0])
        +result[1]/(scale[1]*scale[1])
        +result[2]/(scale[2]*scale[2])
        +result[3]/(scale[3]*scale[3])
        )/(NC-1);

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
    T->SetBranchAddress("Xtg_tr", &Xtg_tr);
    T->SetBranchAddress("Thetatg_tr", &Ttg_tr);
    T->SetBranchAddress("Ytg_tr", &Ytg_tr);
    T->SetBranchAddress("Phitg_tr", &Ptg_tr);
    T->SetBranchAddress("Delta", &Dtg_tr);
}

void FitCor(const Char_t* filename = "test.root") 
{
    LoadTree(filename);
    
	TVirtualFitter::SetDefaultFitter("Minuit2");  //default is Minuit
	TVirtualFitter *fitter = TVirtualFitter::Fitter(NULL, NPara);
#if ORDER == 0
	fitter->SetFCN(myfcn0);
#endif
#if ORDER == 1
    fitter->SetFCN(myfcn1);
#endif

    for(Int_t i = 0; i<NPara; i++)
	{
        fitter->SetParameter(i, Form("Para%02d",i), startvalues[i], parerrors[i], parlimitslo[i], parlimitshi[i]);

		if (!free[i]) fitter->FixParameter(i);
	}

    fitter->Print();

    Double_t arglist[1] = {0};
    iter = 1;
	fitter->ExecuteCommand("MIGRAD", arglist, 0);

#if ORDER == 0
    printf("V5fp[0]+=(%8.5e);\n", fitter->GetParameter(0));
    printf("V5fp[1]+=(%8.5e);\n", fitter->GetParameter(1));
    printf("V5fp[2]+=(%8.5e);\n", fitter->GetParameter(2));
    printf("V5fp[3]+=(%8.5e);\n", fitter->GetParameter(3));
#endif
#if ORDER == 1
    printf("V5fp[0]+=(%8.5e)+(%8.5e)*V5tg[0]+(%8.5e)*V5tg[1]+(%8.5e)*V5tg[2]+(%8.5e)*V5tg[3]+(%8.5e)*V5tg[4];\n", fitter->GetParameter(0), fitter->GetParameter(1), fitter->GetParameter(2), fitter->GetParameter(3), fitter->GetParameter(4), fitter->GetParameter(5));
    printf("V5fp[1]+=(%8.5e)+(%8.5e)*V5tg[0]+(%8.5e)*V5tg[1]+(%8.5e)*V5tg[2]+(%8.5e)*V5tg[3]+(%8.5e)*V5tg[4];\n", fitter->GetParameter(6), fitter->GetParameter(7), fitter->GetParameter(8), fitter->GetParameter(9), fitter->GetParameter(10), fitter->GetParameter(11));
    printf("V5fp[2]+=(%8.5e)+(%8.5e)*V5tg[0]+(%8.5e)*V5tg[1]+(%8.5e)*V5tg[2]+(%8.5e)*V5tg[3]+(%8.5e)*V5tg[4];\n", fitter->GetParameter(12), fitter->GetParameter(13), fitter->GetParameter(14), fitter->GetParameter(15), fitter->GetParameter(16), fitter->GetParameter(17));
    printf("V5fp[3]+=(%8.5e)+(%8.5e)*V5tg[0]+(%8.5e)*V5tg[1]+(%8.5e)*V5tg[2]+(%8.5e)*V5tg[3]+(%8.5e)*V5tg[4];\n", fitter->GetParameter(18), fitter->GetParameter(19), fitter->GetParameter(20), fitter->GetParameter(21), fitter->GetParameter(22), fitter->GetParameter(23));
#endif
}
