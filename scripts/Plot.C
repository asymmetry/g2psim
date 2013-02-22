#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"

TFile *f;
TTree *T;

double fphlowlimit = -0.04, fphhilimit = 0.04;
double fthlowlimit = -0.03, fthhilimit = 0.03;
double fxlowlimit  = -0.06, fxhilimit  = 0.06;
double fylowlimit  = -0.06, fyhilimit  = 0.06;
double phlowlimit  = -0.03, phhilimit  = 0.03;
double thlowlimit  = -0.06, thhilimit  = 0.06;
double dplowlimit  = -0.01, dphilimit  = 0.01;
double xlowlimit  = -0.01, xhilimit  = 0.01;
double ylowlimit  = -0.01, yhilimit  = 0.01;

void PlotFPThPh(){
    TCanvas *c1 = new TCanvas("c1","FP Th vs Ph", 1200, 600);
    c1->Divide(2,1);
    c1->cd(1);
    gPad->SetGrid();
    TH2F* h11 = new TH2F("h11", "FP Th vs Ph (Data)", 200, fphlowlimit, fphhilimit, 200, fthlowlimit, fthhilimit);

    T->Draw("Thetafpdata_rot:Phifpdata_rot>>h11","IsGood","COLZ");
    //T->Draw("Thetafpdata_tr:Phifpdata_tr>>h11","IsGood","COLZ");
    c1->Update();

    c1->cd(2);
    gPad->SetGrid();
    TH2F* h12 = new TH2F("h12", "FP Th vs Ph (SNAKE)", 200, fphlowlimit, fphhilimit, 200, fthlowlimit, fthhilimit);

    T->Draw("Thetafp_rot:Phifp_rot>>h12","IsGood","COLZ");
    //T->Draw("Thetafp_tr:Phifp_tr>>h12","IsGood","COLZ");
    c1->Update();
}

void PlotFPThPhTR(){
    TCanvas *c10 = new TCanvas("c10","FP Th vs Ph", 1200, 600);
    c10->Divide(2,1);
    c10->cd(1);
    gPad->SetGrid();
    TH2F* h101 = new TH2F("h101", "FP Th vs Ph (Data)", 200, fphlowlimit, fphhilimit, 200, fthlowlimit, fthhilimit);

    //T->Draw("Thetafpdata_tr:Phifpdata_tr>>h101","IsGood&&Gun==6","COLZ");
    T->Draw("Thetafpdata_tr:Phifpdata_tr>>h101","IsGood","COLZ");
    c10->Update();

    c10->cd(2);
    gPad->SetGrid();
    TH2F* h102 = new TH2F("h102", "FP Th vs Ph (SNAKE)", 200, fphlowlimit, fphhilimit, 200, fthlowlimit, fthhilimit);

    //T->Draw("Thetafp_tr:Phifp_tr>>h102","IsGood&&Gun==6","COLZ");
    T->Draw("Thetafp_tr:Phifp_tr>>h102","IsGood","COLZ");
    c10->Update();
}

void PlotFPThY(){
    TCanvas *c2 = new TCanvas("c2","FP Th vs Y", 1200, 600);
    c2->Divide(2,1);
    c2->cd(1);
    gPad->SetGrid();
    TH2F* h21 = new TH2F("h21", "FP Th vs Y (Data)", 200, fylowlimit, fyhilimit, 200, fthlowlimit, fthhilimit);

    T->Draw("Thetafpdata_rot:Yfpdata_rot>>h21","IsGood","COLZ");
    c2->Update();

    c2->cd(2);
    gPad->SetGrid();
    TH2F* h22 = new TH2F("h22", "FP Th vs Y (SNAKE)", 200, fylowlimit, fyhilimit, 200, fthlowlimit, fthhilimit);

    T->Draw("Thetafp_rot:Yfp_rot>>h22","IsGood","COLZ");
    c2->Update();
}

void PlotFPXY(){
    TCanvas *c8 = new TCanvas("c8","FP X vs Y", 1200, 600);
    c8->Divide(2,1);
    c8->cd(1);
    gPad->SetGrid();
    TH2F* h81 = new TH2F("h81", "FP X vs Y (Data)", 200, fylowlimit, fyhilimit, 200, fxlowlimit, fxhilimit);

    T->Draw("Xfpdata_rot:Yfpdata_rot>>h81","IsGood&&Gun==6","COLZ");
    c8->Update();

    c8->cd(2);
    gPad->SetGrid();
    TH2F* h82 = new TH2F("h82", "FP X vs Y (SNAKE)", 200, fylowlimit, fyhilimit, 200, fxlowlimit, fxhilimit);

    T->Draw("Xfp_rot:Yfp_rot>>h82","IsGood&&Gun==6","COLZ");
    c8->Update();

    // c8->cd(3);
    // gPad->SetGrid();
    // TH2F* h83 = new TH2F("h83", "FP Th vs Y (SNAKE, Sim)", 200, fylowlimit, fyhilimit, 200, fxlowlimit, fxhilimit);

    // T->Draw("Xfp_rot:Yfp_rot>>h83","IsGood&&Gun==5","COLZ");
    // c8->Update();
}

void PlotTPThPh(){
    TCanvas *c3 = new TCanvas("c3","Target Th vs Ph", 1200, 600);
    c3->Divide(2,1);
    c3->cd(1);
    gPad->SetGrid();
    TH2F* h31 = new TH2F("h31", "Target Th vs Ph (SNAKE)", 200, phlowlimit, phhilimit, 200, thlowlimit, thhilimit);

    T->Draw("Thetarec_tr:Phirec_tr>>h31","IsGood","COLZ");
    c3->Update();

    c3->cd(2);
    gPad->SetGrid();
    TH2F* h32 = new TH2F("h32", "Target Th vs Ph (Database)", 200, phlowlimit, phhilimit, 200, thlowlimit, thhilimit);

    T->Draw("Thetarecdb_tr:Phirecdb_tr>>h32","IsGood","COLZ");
    c3->Update();
}

void PlotTPDp(){
    TCanvas *c4 = new TCanvas("c4","TP Dp", 600, 300);
    c4->cd();
    gPad->SetGrid();
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
    TCanvas *c0 = new TCanvas("c0","Origin Th vs Ph", 1200, 600);
    c0->Divide(2,1);
    c0->cd(1);
    gPad->SetGrid();
    TH2F* h01 = new TH2F("h01", "Origin Th vs Ph (CALCULATE)", 100, phlowlimit, phhilimit, 100, thlowlimit, thhilimit);

    T->Draw("Thetatg_tr:Phitg_tr>>h01","IsGood","COLZ");
    c0->Update();

    c0->cd(2);
    gPad->SetGrid();
    TH2F* h02 = new TH2F("h02", "Origin Th vs Ph (DATA, Database)", 100, phlowlimit, phhilimit, 100, thlowlimit, thhilimit);

    T->Draw("Thetarecdb_tr:Phirecdb_tr>>h02","IsGood","COLZ");
    c0->Update();
}

void PlotORYX(){
    TCanvas *c01 = new TCanvas("c0","Origin Y vs X", 1200, 600);
    c01->Divide(2,1);
    c01->cd(1);
    gPad->SetGrid();
    TH2F* h001 = new TH2F("h001", "Origin Y vs X", 100, xlowlimit, xhilimit, 100, ylowlimit, yhilimit);

    T->Draw("Xbeam_lab:Ybeam_lab>>h001","IsGood","COLZ");
    c01->Update();

    c01->cd(2);
    gPad->SetGrid();
    TH2F* h002 = new TH2F("h002", "Origin Y vs X (BPM)", 100, xlowlimit, xhilimit, 100, ylowlimit, yhilimit);

    T->Draw("Xbpm_lab:Ybpm_lab>>h002","IsGood","COLZ");
    c01->Update();
}

void LoadTree(const char* filename = "test.root"){
    // Will add TChain later
    f = new TFile(filename);
    T = (TTree *)f->Get("T");
}
