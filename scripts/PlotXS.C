#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"

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

double IsGood;
double Delta, Deltarec;
double Thetainit, Thetarec;
double Xfp, Yfp, Tfp, Pfp;
double XSinit, XSrec;

void PlotDelta()
{
    TCanvas* c110 = new TCanvas("c110", "XS", 1200, 600);
    c110->Divide(2,1);

    c110->cd(1);

    TH1F* h1101 = new TH1F("h1101", "Init Delta (weight by XS)", 50, -0.05, 0.05);
    h1101->SetXTitle("Delta");

    for (int i = 0; i<T->GetEntries(); i++) {
        T->GetEntry(i);
        //if (IsGood&&TMath::Abs(Tfp)<0.05&&TMath::Abs(Pfp-1e-3)<0.01) {
        if (TMath::Abs(IsGood-1.0)<1e-8) {
            h1101->Fill(Delta, XSinit);
        }
    }

    h1101->Draw();

    c110->Update();

    c110->cd(2);

    TH1F* h1102 = new TH1F("h1102", "Rec Delta (weight by XS)", 50, -0.05, 0.05);
    h1102->SetXTitle("Delta_rec");

    for (int i = 0; i<T->GetEntries(); i++) {
        T->GetEntry(i);
        if (TMath::Abs(IsGood-1.0)<1e-8) {
            h1102->Fill(Deltarec, XSrec);
        }
    }

    h1102->Draw();

    c110->Update();
}

void PlotTheta()
{
    TCanvas* c111 = new TCanvas("c111", "XS", 1200, 600);
    c111->Divide(2,1);

    c111->cd(1);

    TH1F* h1111 = new TH1F("h1111", "Init Scat Angle (weight by XS)", 50, 0.08, 0.12);
    h1111->SetXTitle("Scat Angle");

    for (int i = 0; i<T->GetEntries(); i++) {
        T->GetEntry(i);
        //if (IsGood&&TMath::Abs(Tfp)<0.05&&TMath::Abs(Pfp-1e-3)<0.01) {
        if (TMath::Abs(IsGood-1.0)<1e-8) {
            h1111->Fill(Thetainit, XSinit);
        }
    }

    h1111->Draw();

    c111->Update();

    c111->cd(2);

    TH1F* h1112 = new TH1F("h1112", "Rec Scat Angle (weight by XS)", 50, 0.08, 0.12);
    h1112->SetXTitle("Scat Angle");

    for (int i = 0; i<T->GetEntries(); i++) {
        T->GetEntry(i);
        if (TMath::Abs(IsGood-1.0)<1e-8) {
            h1112->Fill(Thetarec, XSrec);
        }
    }

    h1112->Draw();

    c111->Update();
}

void LoadTree(const char* filename = "test.root"){
    // Will add TChain later
    f = new TFile(filename);
    T = (TTree *)f->Get("T");

    T->SetBranchAddress("fwd.isgood", &IsGood);
    T->SetBranchAddress("gun.react.d", &Delta);
    T->SetBranchAddress("bwd.rec.d", &Deltarec);
    T->SetBranchAddress("phy.react.xs", &XSinit);
    T->SetBranchAddress("phy.rec.xs", &XSrec);
    T->SetBranchAddress("fwd.focus.r_x", &Xfp);
    T->SetBranchAddress("fwd.focus.r_y", &Yfp);
    T->SetBranchAddress("fwd.focus.r_t", &Tfp);
    T->SetBranchAddress("fwd.focus.r_p", &Pfp);
    T->SetBranchAddress("phy.react.angle", &Thetainit);
    T->SetBranchAddress("phy.rec.angle", &Thetarec);

    gStyle->SetStatH(0.3);
    gStyle->SetStatW(0.25);
    gStyle->SetOptStat("emr");
}
