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

void CalCorTPD()
{
    TCanvas* c1 = new TCanvas("c1","FP Th vs Ph", 1200, 600);
    c1->Divide(2,1);

    c1->cd(1);
    gPad->SetGrid();
    TH2F* h11 = new TH2F("h11", "FP Th vs Ph", 200, fphlowlimit, fphhilimit, 200, fthlowlimit, fthhilimit);

    T->Draw("Thetafpdata_tr:Phifpdata_tr>>h11","IsGood&&Gun==6","COLZ");
    c1->Update();

    TCutG* cutg = (TCutG*)(TVirtualPad::Pad()->WaitPrimitive("CUTG", "CutG"));
    c1->Update();

    cutg->SetName("FPCenter");
    cutg->SetVarX("Phifpdata_tr");
    cutg->SetVarY("Thetafpdata_tr");

    cutg->SetLineColor(kMagenta);
    cutg->SetLineWidth(2);
    cutg->Draw("PL");
    c1->Update();

    c1->cd(2);
    gPad->SetGrid();
    TH2F* h12 = new TH2F("h12", "FP Th vs Ph", 200, fphlowlimit, fphhilimit, 200, fthlowlimit, fthhilimit);

    T->Draw("Thetafpdata_tr:Phifpdata_tr>>h12","IsGood&&Gun==6&&FPCenter","COLZ");
    c1->Update();
}

void CalCorXYD()
{
    TCanvas* c2 = new TCanvas("c2","FP X vs Y", 1200, 600);
    c2->Divide(2,1);

    c2->cd(1);
    gPad->SetGrid();
    TH2F* h21 = new TH2F("h21", "FP X vs Y", 200, fylowlimit, fyhilimit, 200, fxlowlimit, fxhilimit);

    T->Draw("Xfpdata_tr:Yfpdata_tr>>h21","IsGood&&Gun==6","COLZ");
    c2->Update();

    TCutG* cutg = (TCutG*)(TVirtualPad::Pad()->WaitPrimitive("CUTG", "CutG"));
    c2->Update();

    cutg->SetName("FPCenter");
    cutg->SetVarX("Yfpdata_tr");
    cutg->SetVarY("Xfpdata_tr");

    cutg->SetLineColor(kMagenta);
    cutg->SetLineWidth(2);
    cutg->Draw("PL");
    c2->Update();

    c2->cd(2);
    gPad->SetGrid();
    TH2F* h22 = new TH2F("h22", "FP X vs Y", 200, fylowlimit, fyhilimit, 200, fxlowlimit, fxhilimit);

    T->Draw("Xfpdata_tr:Yfpdata_tr>>h22","IsGood&&Gun==6&&FPCenter","COLZ");
    c2->Update();
}

void CalCorTPS()
{
    TCanvas* c3 = new TCanvas("c3","FP Th vs Ph", 1200, 600);
    c3->Divide(2,1);

    c3->cd(1);
    gPad->SetGrid();
    TH2F* h31 = new TH2F("h31", "FP Th vs Ph", 200, fphlowlimit, fphhilimit, 200, fthlowlimit, fthhilimit);

    T->Draw("Thetafp_tr:Phifp_tr>>h31","IsGood&&Gun==6","COLZ");
    c3->Update();

    TCutG* cutg = (TCutG*)(TVirtualPad::Pad()->WaitPrimitive("CUTG", "CutG"));
    c3->Update();

    cutg->SetName("FPCenter");
    cutg->SetVarX("Phifp_tr");
    cutg->SetVarY("Thetafp_tr");

    cutg->SetLineColor(kMagenta);
    cutg->SetLineWidth(2);
    cutg->Draw("PL");
    c3->Update();

    c3->cd(2);
    gPad->SetGrid();
    TH2F* h32 = new TH2F("h32", "FP Th vs Ph", 200, fphlowlimit, fphhilimit, 200, fthlowlimit, fthhilimit);

    T->Draw("Thetafp_tr:Phifp_tr>>h32","IsGood&&Gun==6&&FPCenter","COLZ");
    c3->Update();
}

void CalCorXYS()
{
    TCanvas* c4 = new TCanvas("c4","FP X vs Y", 1200, 600);
    c4->Divide(2,1);

    c4->cd(1);
    gPad->SetGrid();
    TH2F* h41 = new TH2F("h41", "FP X vs Y", 200, fylowlimit, fyhilimit, 200, fxlowlimit, fxhilimit);

    T->Draw("Xfp_tr:Yfp_tr>>h41","IsGood&&Gun==6","COLZ");
    c4->Update();

    TCutG* cutg = (TCutG*)(TVirtualPad::Pad()->WaitPrimitive("CUTG", "CutG"));
    c4->Update();

    cutg->SetName("FPCenter");
    cutg->SetVarX("Yfp_tr");
    cutg->SetVarY("Xfp_tr");

    cutg->SetLineColor(kMagenta);
    cutg->SetLineWidth(2);
    cutg->Draw("PL");
    c4->Update();

    c4->cd(2);
    gPad->SetGrid();
    TH2F* h42 = new TH2F("h42", "FP X vs Y", 200, fylowlimit, fyhilimit, 200, fxlowlimit, fxhilimit);

    T->Draw("Xfp_tr:Yfp_tr>>h42","IsGood&&Gun==6&&FPCenter","COLZ");
    c4->Update();
}

void LoadTree(const char* filename = "test.root"){
    // Will add TChain later
    f = new TFile(filename);
    T = (TTree *)f->Get("T");
}
