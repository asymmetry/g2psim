#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
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

void PlotDXvsP()
{
    TCanvas* c110 = new TCanvas("c110", "BPM", 600, 600);

    TH2F* h1101 = new TH2F("h1101", "(Xtg-Xproj) vs P", 200, 0.4, 2.4, 200, 0.0, 20.0);

    T->Project("h1101", "(Xtg_tr-Xtgproj_tr)*1000:(1+Delta)*2.249497","");

    TF1* f = new TF1("f1101", "[0]+[1]/x", 0.4, 2.4);

    h1101->Fit(f);

    h1101->SetXTitle("P/GeV");
    h1101->SetYTitle("(Xtg-Xproj)/mm");

    h1101->Draw();
}

void PlotDXPvsZ()
{
    TCanvas* c111 = new TCanvas("c111", "BPM", 600, 600);

    TH2F* h1111 = new TH2F("h1111", "(Xbpm-Xtg) vs reactz", 200, -0.015, 0.015, 200, -2.0, 2.0);

    //T->Project("h1111", "(Xbpm_tr-Xtg_tr)*1000*(1+Delta)*1.0:Zbeam_lab","IsGood");
    T->Project("h1111", "(Xbpm_tr-Xtg_tr)*1000:Zbeam_lab","IsGood");

    TF1* f = new TF1("f1111", "[0]+[1]*x", -0.15, 0.15);

    h1111->Fit(f);

    h1111->SetXTitle("reactz/m");
    h1111->SetYTitle("(Xbpm-Xtg)/mm");

    h1111->Draw();
}

void LoadTree(const char* filename = "test.root"){
    // Will add TChain later
    f = new TFile(filename);
    T = (TTree *)f->Get("T");

    //gStyle->SetStatH(0.3);
    //gStyle->SetStatW(0.25);
    gStyle->SetOptStat("emr");
}
