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
double bxlowlimit  = -0.01, bxhilimit  = 0.01;
double bylowlimit  = -0.01, byhilimit  = 0.01;

void PlotFPThPh(){
    TCanvas *c1 = new TCanvas("c1","FP Th vs Ph", 1200, 600);
    c1->Divide(2,1);
    c1->cd(1);
    gPad->SetGrid();
    TH2F* h11 = new TH2F("h11", "FP Th vs Ph (Data)", 200, fphlowlimit, fphhilimit, 200, fthlowlimit, fthhilimit);

    T->Draw("Tfpdata_rot:Pfpdata_rot>>h11","IsGood","COLZ");
    //T->Draw("Tfpdata_tr:Pfpdata_tr>>h11","IsGood","COLZ");
    c1->Update();

    c1->cd(2);
    gPad->SetGrid();
    TH2F* h12 = new TH2F("h12", "FP Th vs Ph (SNAKE)", 200, fphlowlimit, fphhilimit, 200, fthlowlimit, fthhilimit);

    T->Draw("Tfp_rot:Pfp_rot>>h12","IsGood","COLZ");
    //T->Draw("Tfp_tr:Pfp_tr>>h12","IsGood","COLZ");
    c1->Update();
}

void PlotFPThPhTR(){
    TCanvas *c10 = new TCanvas("c10","FP Th vs Ph", 1200, 600);
    c10->Divide(2,1);
    c10->cd(1);
    gPad->SetGrid();
    TH2F* h101 = new TH2F("h101", "FP Th vs Ph (Data)", 200, fphlowlimit, fphhilimit, 200, fthlowlimit, fthhilimit);

    //T->Draw("Tfpdata_tr:Pfpdata_tr>>h101","IsGood&&Gun==6","COLZ");
    T->Draw("Tfpdata_tr:Pfpdata_tr>>h101","IsGood","COLZ");
    c10->Update();

    c10->cd(2);
    gPad->SetGrid();
    TH2F* h102 = new TH2F("h102", "FP Th vs Ph (SNAKE)", 200, fphlowlimit, fphhilimit, 200, fthlowlimit, fthhilimit);

    //T->Draw("Tfp_tr:Pfp_tr>>h102","IsGood&&Gun==6","COLZ");
    T->Draw("Tfp_tr:Pfp_tr>>h102","IsGood","COLZ");
    c10->Update();
}

void PlotFPThY(){
    TCanvas *c2 = new TCanvas("c2","FP Th vs Y", 1200, 600);
    c2->Divide(2,1);
    c2->cd(1);
    gPad->SetGrid();
    TH2F* h21 = new TH2F("h21", "FP Th vs Y (Data)", 200, fylowlimit, fyhilimit, 200, fthlowlimit, fthhilimit);

    T->Draw("Tfpdata_rot:Yfpdata_rot>>h21","IsGood","COLZ");
    c2->Update();

    c2->cd(2);
    gPad->SetGrid();
    TH2F* h22 = new TH2F("h22", "FP Th vs Y (SNAKE)", 200, fylowlimit, fyhilimit, 200, fthlowlimit, fthhilimit);

    T->Draw("Tfp_rot:Yfp_rot>>h22","IsGood","COLZ");
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

    T->Draw("Trec_tr:Prec_tr>>h31","IsGood","COLZ");
    c3->Update();

    c3->cd(2);
    gPad->SetGrid();
    TH2F* h32 = new TH2F("h32", "Target Th vs Ph (Database)", 200, phlowlimit, phhilimit, 200, thlowlimit, thhilimit);

    T->Draw("Trecdb_tr:Precdb_tr>>h32","IsGood","COLZ");
    c3->Update();
}

void PlotTPThY(){
    TCanvas *c32 = new TCanvas("c32","Target Th vs Y", 1200, 600);
    c32->Divide(2,1);
    c32->cd(1);
    gPad->SetGrid();
    TH2F* h321 = new TH2F("h321", "Target Th vs Y (SNAKE)", 200, ylowlimit, yhilimit, 200, thlowlimit, thhilimit);

    T->Draw("Trec_tr:Yrec_tr>>h321","IsGood","COLZ");
    c32->Update();

    c32->cd(2);
    gPad->SetGrid();
    TH2F* h322 = new TH2F("h322", "Target Th vs Ph (Database)", 200, ylowlimit, yhilimit, 200, thlowlimit, thhilimit);

    T->Draw("Trecdb_tr:Yrecdb_tr>>h322","IsGood","COLZ");
    c32->Update();
}

void PlotTPDp(){
    TCanvas *c4 = new TCanvas("c4","TP Dp", 600, 300);
    c4->cd();
    gPad->SetGrid();
    TH1F* h41 = new TH1F("h41", "TP Dp", 200, dplowlimit, dphilimit);
    h41->SetLineColor(2);
    
    T->Draw("Delta>>h41","IsGood");
    c4->Update();
    
    TH1F* h42 = new TH1F("h42", "TP Dp", 200, dplowlimit, dphilimit);
    h42->SetLineColor(1);

    T->Draw("Deltarec>>h42","IsGood","same");
    c4->Update();
}

void PlotORThPh(){
    TCanvas *c0 = new TCanvas("c0","Origin Th vs Ph", 1200, 600);
    c0->Divide(2,1);
    c0->cd(1);
    gPad->SetGrid();
    TH2F* h01 = new TH2F("h01", "Origin Th vs Ph (CALCULATE)", 100, phlowlimit, phhilimit, 100, thlowlimit, thhilimit);

    T->Draw("Ttg_tr:Ptg_tr>>h01","IsGood","COLZ");
    c0->Update();

    c0->cd(2);
    gPad->SetGrid();
    TH2F* h02 = new TH2F("h02", "Origin Th vs Ph (DATA, Database)", 100, phlowlimit, phhilimit, 100, thlowlimit, thhilimit);

    T->Draw("Trecdb_tr:Precdb_tr>>h02","IsGood","COLZ");
    c0->Update();
}

void PlotORYX(){
    TCanvas *c01 = new TCanvas("c0","Origin Y vs X", 1200, 600);
    c01->Divide(2,1);
    c01->cd(1);
    gPad->SetGrid();
    TH2F* h001 = new TH2F("h001", "Origin Y vs X", 100, bxlowlimit, bxhilimit, 100, bylowlimit, byhilimit);

    T->Draw("Xbeam_lab:Ybeam_lab>>h001","IsGood","COLZ");
    c01->Update();

    c01->cd(2);
    gPad->SetGrid();
    TH2F* h002 = new TH2F("h002", "Origin Y vs X (BPM)", 100, bxlowlimit, bxhilimit, 100, bylowlimit, byhilimit);

    T->Draw("Xbpm_lab:Ybpm_lab>>h002","IsGood","COLZ");
    c01->Update();
}

void PlotRecError()
{
    TCanvas* c100 = new TCanvas("c100", "Reconstruction Uncertainty", 1200, 600);
    c100->Divide(2,2);
    c100->cd(1);
    TH1F* h1001 = new TH1F("h1001", "Delta (x 10e-4)", 100, dplowlimit/10.0*10000, dphilimit/10.0*10000);

    T->Draw("(Deltarec-Delta)*10000>>h1001","IsGood");

    c100->cd(2);
    TH1F* h1002 = new TH1F("h1002", "T (mrad)", 100, thlowlimit/5.0*1000, thhilimit/5.0*1000);

    T->Draw("(Trec_tr-Treact_tr)*1000>>h1002","IsGood");

    c100->cd(3);
    TH1F* h1003 = new TH1F("h1003", "Y (mm)", 100, ylowlimit/2.0*1000, yhilimit/2.0*1000);

    T->Draw("(Yrec_tr-Yreact_tr)*1000>>h1003","IsGood");

    c100->cd(4);
    TH1F* h1004 = new TH1F("h1004", "P (mrad)", 100, phlowlimit/10.0*1000, phhilimit/10.0*1000);

    T->Draw("(Prec_tr-Preact_tr)*1000>>h1004","IsGood");

    c100->Update();
}

void PlotXRecError()
{
    TCanvas* c297 = new TCanvas("c297", "Reconstruction Uncertainty", 600, 600);
    TH1F* h2971 = new TH1F("h2971", "Eff Xproj Uncertainty", 100, -2, 2);

    h2971->SetXTitle("(Xeff-Xreal)/mm");

    T->Draw("(Xeff_tr-Xreal_tr)*1000>>h2971","IsGood");
}

void PlotSNAKERecError()
{
    TCanvas* c101 = new TCanvas("c101", "Reconstruction Uncertainty", 1200, 600);
    c101->Divide(2,2);
    c101->cd(1);
    TH1F* h1011 = new TH1F("h1011", "Delta (x 10e-4)", 100, dplowlimit/10.0*10000, dphilimit/10.0*10000);

    T->Draw("(Deltarec-Delta)*10000>>h1011","IsGood");

    c101->cd(2);
    TH1F* h1012 = new TH1F("h1012", "T (mrad)", 100, thlowlimit/5.0*1000, thhilimit/5.0*1000);

    T->Draw("(Trectg_tr-Ttgproj_tr)*1000>>h1012","IsGood");

    c101->cd(3);
    TH1F* h1013 = new TH1F("h1013", "Y (mm)", 100, ylowlimit/2.0*1000, yhilimit/2.0*1000);

    T->Draw("(Yrectg_tr-Ytgproj_tr)*1000>>h1013","IsGood");

    c101->cd(4);
    TH1F* h1014 = new TH1F("h1014", "P (mrad)", 100, phlowlimit/10.0*1000, phhilimit/10.0*1000);

    T->Draw("(Prectg_tr-Ptgproj_tr)*1000>>h1014","IsGood");

    c101->Update();
}

void PlotXBPM()
{
    TCanvas* c110 = new TCanvas("c110", "BPM", 1200, 600);
    c110->Divide(2,2);
    c110->cd(1);
    TH2F* h1101 = new TH2F("h1101", "BPM Deviation vs Delta", 100, -0.02, 0.02, 100, 0.0, 0.005);

    T->Draw("(Xtg_tr-Xtgproj_tr):Delta>>h1101","");

    c110->cd(2);
    TH2F* h1102 = new TH2F("h1102", "BPM Deviation vs z", 100, -0.015, 0.015, 100, -0.001, 0.005);

    T->Draw("(Xbpm_tr-Xtg_tr)*(1+Delta):Zbeam_lab>>h1102","");

    // c110->cd(3);
    // TH1F* h1103 = new TH1F("h1103", "Y (mm)", 110, ylowlimit/2.0*1100, yhilimit/2.0*1100);

    // T->Draw("(Yrec_tr-Ytg_tr)*1100>>h1103","IsGood");

    // c110->cd(4);
    // TH1F* h1104 = new TH1F("h1104", "P (mrad)", 110, phlowlimit/10.0*1100, phhilimit/10.0*1100);

    // T->Draw("(Prec_tr-Ptg_tr)*1100>>h1104","IsGood");

    // c110->Update();
}

void LoadTree(const char* filename = "test.root"){
    // Will add TChain later
    f = new TFile(filename);
    T = (TTree *)f->Get("T");

    gStyle->SetStatH(0.3);
    gStyle->SetStatW(0.25);
    gStyle->SetOptStat("emr");
}
