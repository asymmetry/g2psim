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
double Xinit, Yinit, Tinit, Pinit, Zinit;
double Xrec, Yrec, Trec, Prec;
double Xfp, Yfp, Tfp, Pfp;
double XSinit, XSrec;

void PlotDelta()
{
    TCanvas* c110 = new TCanvas("c110", "XS", 600, 600);

    c110->cd();

    TH1F* h1101 = new TH1F("h1101", "", 50, -0.06, 0.08);
    h1101->SetXTitle("Delta");
    h1101->SetLineColor(1);
    h1101->SetLineWidth(2);

    for (int i = 0; i<T->GetEntries(); i++) {
        T->GetEntry(i);
        if ((TMath::Abs(IsGood-1.0)<1e-8)&&(Zinit>5.0e-3)) {
            h1101->Fill(Delta, XSinit);
        }
    }

    h1101->Draw();

    TH1F* h1102 = new TH1F("h1102", "", 50, -0.06, -0.08);
    h1102->SetXTitle("Delta");
    h1102->SetLineColor(2);
    h1102->SetLineWidth(2);

    for (int i = 0; i<T->GetEntries(); i++) {
        T->GetEntry(i);
        if ((TMath::Abs(IsGood-1.0)<1e-8)&&(Zinit>5.0e-3)) {
            h1102->Fill(Deltarec, XSrec);
        }
    }

    h1102->Draw("same");

    double mean1 = h1101->GetMean();
    double mean2 = h1102->GetMean();

    printf("%e\t%e\t%e\n", mean1, mean2, fabs(mean2-mean1));

    c110->Update();
}

void PlotAngle()
{
    TCanvas* c111 = new TCanvas("c111", "XS", 600, 600);

    c111->cd();

    TH1F* h1111 = new TH1F("h1111", "", 50, 0.08, 0.12);
    h1111->SetXTitle("Scat Angle");

    for (int i = 0; i<T->GetEntries(); i++) {
        T->GetEntry(i);
        if (TMath::Abs(IsGood-1.0)<1e-8) {
            h1111->Fill(Thetainit, XSinit);
        }
    }

    h1111->Draw();

    TH1F* h1112 = new TH1F("h1112", "", 50, 0.08, 0.12);
    h1112->SetXTitle("Scat Angle");

    for (int i = 0; i<T->GetEntries(); i++) {
        T->GetEntry(i);
        if (TMath::Abs(IsGood-1.0)<1e-8) {
            h1112->Fill(Thetarec, XSrec);
        }
    }

    h1112->Draw("same");

    double mean1 = h1111->GetMean();
    double mean2 = h1112->GetMean();

    printf("%e\t%e\t%e\n", mean1, mean2, fabs(mean2-mean1));

    c111->Update();
}

void PlotTheta()
{
    TCanvas* c113 = new TCanvas("c113", "XS", 600, 600);

    c113->cd();

    TH1F* h1131 = new TH1F("h1131", "", 50, -0.07-0.015, -0.07+0.015);
    h1131->SetXTitle("Theta");
    h1131->SetLineColor(1);
    h1131->SetLineWidth(2);   

    for (int i = 0; i<T->GetEntries(); i++) {
        T->GetEntry(i);
        if ((TMath::Abs(IsGood-1.0)<1e-8)&&(fabs(Tinit+0.07)<0.003)) {
            h1131->Fill(Tinit, XSinit);
        }
    }

    h1131->Draw();

    TH1F* h1132 = new TH1F("h1132", "", 50, -0.07-0.015, -0.07+0.015);
    h1132->SetXTitle("Theta");
    h1132->SetLineColor(2);
    h1132->SetLineWidth(2);

    for (int i = 0; i<T->GetEntries(); i++) {
        T->GetEntry(i);
        if ((TMath::Abs(IsGood-1.0)<1e-8)&&(fabs(Tinit+0.07)<0.003)) {
            h1132->Fill(Trec, XSrec);
        }
    }

    h1132->Draw("same");

    double mean1 = h1131->GetMean();
    double mean2 = h1132->GetMean();

    printf("%e\t%e\t%e\n", mean1, mean2, fabs(mean2-mean1));

    c113->Update();
}

void PlotPhi()
{
    TCanvas* c114 = new TCanvas("c114", "XS", 600, 600);

    c114->cd();

    TH1F* h1141 = new TH1F("h1141", "", 50, -0.02-0.005, -0.02+0.005);
    h1141->SetXTitle("Phi");
    h1141->SetLineColor(1);
    h1141->SetLineWidth(2);

    for (int i = 0; i<T->GetEntries(); i++) {
        T->GetEntry(i);
        if ((TMath::Abs(IsGood-1.0)<1e-8)&&(fabs(Pinit+0.02)<0.002)) {
            h1141->Fill(Pinit, XSinit);
        }
    }

    h1141->Draw();

    TH1F* h1142 = new TH1F("h1142", "", 50, -0.02-0.005, -0.02+0.005);
    h1142->SetXTitle("Phi");
    h1142->SetLineColor(2);
    h1142->SetLineWidth(2);

    for (int i = 0; i<T->GetEntries(); i++) {
        T->GetEntry(i);
        if ((TMath::Abs(IsGood-1.0)<1e-8)&&(fabs(Pinit+0.02)<0.002))  {
            h1142->Fill(Prec, XSrec);
        }
    }

    h1142->Draw("same");

    double mean1 = h1141->GetMean();
    double mean2 = h1142->GetMean();

    printf("%e\t%e\t%e\n", mean1, mean2, fabs(mean2-mean1));

    c114->Update();
}

void PlotXS()
{
    TCanvas* c112 = new TCanvas("c112", "XS", 600, 600);

    c112->cd();

    TH1F* h1121 = new TH1F("h1121", "Cross section central value", 50, 0.00, 0.80);
    h1121->SetXTitle("XS(mb)");

    for (int i = 0; i<T->GetEntries(); i++) {
        T->GetEntry(i);
        //if (IsGood&&TMath::Abs(Tfp)<0.05&&TMath::Abs(Pfp-1e-3)<0.01) {
        if (TMath::Abs(IsGood-1.0)<1e-8) {
            h1121->Fill(XSinit/1000);
        }
    }

    h1121->Draw();

    TH1F* h1122 = new TH1F("h1122", "Cross section central value", 50, 0.00, 0.80);
    h1122->SetXTitle("XS(mb)");

    for (int i = 0; i<T->GetEntries(); i++) {
        T->GetEntry(i);
        if (TMath::Abs(IsGood-1.0)<1e-8) {
            h1122->Fill(XSrec/1000);
        }
    }

    h1122->Draw("same");

    double mean1 = h1121->GetMean();
    double mean2 = h1122->GetMean();

    printf("%e\t%e\t%e\n", mean1, mean2, mean2/mean1);

    c112->Update();
}

void LoadTree(const char* filename = "test.root"){
    // Will add TChain later
    f = new TFile(filename);
    T = (TTree *)f->Get("T");

    T->SetBranchAddress("isgood", &IsGood);
    
    T->SetBranchAddress("gun.react.d", &Delta);
    T->SetBranchAddress("gun.react.x", &Xinit);
    T->SetBranchAddress("gun.react.t", &Tinit);
    T->SetBranchAddress("gun.react.y", &Yinit);
    T->SetBranchAddress("gun.react.p", &Pinit);
    T->SetBranchAddress("gun.react.l_z",&Zinit);
    T->SetBranchAddress("bwd.rec.d", &Deltarec);
    T->SetBranchAddress("bwd.rec.x", &Xrec);
    T->SetBranchAddress("bwd.rec.t", &Trec);
    T->SetBranchAddress("bwd.rec.y", &Yrec);
    T->SetBranchAddress("bwd.rec.p", &Prec);
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
    //gStyle->SetOptStat("emr");
    gStyle->SetOptStat("");
}
