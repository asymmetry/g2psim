#include "TROOT.h"

TFile *f;
TTree *T;

double fphlowlimit = -0.04, fphhilimit = 0.04;
double fthlowlimit = -0.03, fthhilimit = 0.03;
double fylowlimit  = -60.0, fyhilimit  = 60.0;
double phlowlimit  = -0.03, phhilimit  = 0.03;
double thlowlimit  = -0.06, thhilimit  = 0.06;
double dplowlimit  = -0.01, dphilimit  = 0.01;

void PlotFPThPh(){
    TCanvas *c1 = new TCanvas("c1","FP Th vs Ph", 1200, 600);
    c1->Divide(2,1);
    c1->cd(1);
    TH2F* h11 = new TH2F("h11", "FP Th vs Ph (Data)", 200, fphlowlimit, fphhilimit, 200, fthlowlimit, fthhilimit);

    T->Draw("Thetafpdata_rot:Phifpdata_rot>>h11","IsGood","COLZ");
    c1->Update();

    c1->cd(2);
    TH2F* h12 = new TH2F("h12", "FP Th vs Ph (SNAKE)", 200, fphlowlimit, fphhilimit, 200, fthlowlimit, fthhilimit);

    T->Draw("Thetafp_rot:Phifp_rot>>h12","IsGood","COLZ");
    c1->Update();
}

void PlotFPThY(){
    TCanvas *c2 = new TCanvas("c2","FP Th vs Y", 600, 600);
    c2->cd();
    TH2F* h2 = new TH2F("h2", "FP Th vs Y", 200, fylowlimit, fyhilimit, 200, fthlowlimit, fthhilimit);

    T->Draw("Thetafp_tr:Yfp_tr>>h2","IsGood","COLZ");
    c2->Update();
}

void PlotTPThPh(){
    TCanvas *c3 = new TCanvas("c3","Target Th vs Ph", 1200, 600);
    c3->Divide(2,1);
    c3->cd(1);
    TH2F* h31 = new TH2F("h31", "Target Th vs Ph (Snake Model)", 200, phlowlimit, phhilimit, 200, thlowlimit, thhilimit);

    T->Draw("Thetarec_tr:Phirec_tr>>h31","IsGood","COLZ");
    c3->Update();

    c3->cd(2);
    TH2F* h32 = new TH2F("h32", "Target Th vs Ph (Database)", 200, phlowlimit, phhilimit, 200, thlowlimit, thhilimit);

    T->Draw("Thetarecdb_tr:Phirecdb_tr>>h32","IsGood","COLZ");
    c3->Update();
}

void PlotTPDp(){
    TCanvas *c4 = new TCanvas("c4","TP Dp", 600, 300);
    c4->cd();
    TH1F* h41 = new TH1F("h41", "TP Dp", 200, dplowlimit, dphilimit);
    
    h41->SetLineColor(2);
    
    T->Draw("Deltarecdb>>h41","IsGood");
    c4->Update();
    
    TH1F* h42 = new TH1F("h42", "TP Dp", 200, dplowlimit, dphilimit);

    h42->SetLineColor(1);
    
    T->Draw("Delta_rec>>h42","IsGood","same");
    c4->Update();
}

void PlotORThPh(){
    TCanvas *c0 = new TCanvas("c0","Origin Th vs Ph", 600, 600);
    c0->cd();
    TH2F* h0 = new TH2F("h0", "Origin Th vs Ph", 100, phlowlimit, phhilimit, 100, thlowlimit, thhilimit);

    T->Draw("Thetatg_tr:Phitg_tr>>h0","IsGood","COLZ");
    c0->Update();
}

void LoadTree(const char* filename = "test.root"){
    // Will add TChain later
    f = new TFile(filename);
    T = (TTree *)f->Get("T");
}
