#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"

TFile* f;
TTree* T;

double fphlowlimit = -0.04, fphhilimit = 0.04;
double fthlowlimit = -0.03, fthhilimit = 0.03;
double fxlowlimit  = -0.06, fxhilimit  = 0.06;
double fylowlimit  = -0.06, fyhilimit  = 0.06;
double phlowlimit  = -0.03, phhilimit  = 0.03;
double thlowlimit  = -0.06, thhilimit  = 0.06;
double dplowlimit  = -0.01, dphilimit  = 0.01;

void PlotComFPThPh()
{
    TCanvas* c1 = new TCanvas("c1","FP Th vs Ph", 600, 600);

    gPad->SetGrid();
    TH2F* h11 = new TH2F("h11", "FP Th vs Ph", 200, fphlowlimit, fphhilimit, 200, fthlowlimit, fthhilimit);
    T->Project("h11", "Thetafpdata_rot:Phifpdata_rot","IsGood&&Gun==6");
        
    TH2F* h12 = new TH2F("h12", "FP Th vs Ph", 200, fphlowlimit, fphhilimit, 200, fthlowlimit, fthhilimit);
    h12->SetMarkerColor(kRed);
    T->Project("h12", "Thetafp_rot:Phifp_rot","IsGood&&Gun==6");

    h11->Draw();
    h12->Draw("same");
    c1->Update();
}

void PlotComFPThY()
{
    TCanvas* c2 = new TCanvas("c2","FP Th vs Y", 600, 600);

    gPad->SetGrid();
    TH2F* h21 = new TH2F("h21", "FP Th vs Y", 200, fylowlimit, fyhilimit, 200, fthlowlimit, fthhilimit);
    T->Project("h21", "Thetafpdata_rot:Yfpdata_rot","IsGood&&Gun==6");
    
    TH2F* h22 = new TH2F("h22", "FP Th vs Y", 200, fylowlimit, fyhilimit, 200, fthlowlimit, fthhilimit);
    h22->SetMarkerColor(kRed);
    T->Project("h22", "Thetafp_rot:Yfp_rot","IsGood&&Gun==6");

    h21->Draw();
    h22->Draw("same");    
    c2->Update();
}

void PlotComFPXY()
{
    TCanvas* c3 = new TCanvas("c3","FP X vs Y", 600, 600);

    gPad->SetGrid();
    TH2F* h31 = new TH2F("h31", "FP X vs Y", 200, fylowlimit, fyhilimit, 200, fxlowlimit, fxhilimit);
    T->Project("h31","Xfpdata_rot:Yfpdata_rot","IsGood&&Gun==6");
    c3->Update();
    
    TH2F* h32 = new TH2F("h32", "FP X vs Y", 200, fylowlimit, fyhilimit, 200, fxlowlimit, fxhilimit);
    h32->SetMarkerColor(kRed);
    T->Project("h32","Xfp_rot:Yfp_rot","IsGood&&Gun==6");

    h31->Draw();
    h32->Draw("same");
    c3->Update();
}


void LoadTree(const char* filename = "test.root"){
    // Will add TChain later
    f = new TFile(filename);
    T = (TTree *)f->Get("T");
}
