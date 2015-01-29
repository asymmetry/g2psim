#include "TROOT.h"
#include "TError.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TMath.h"
#include "TTree.h"

void FitBPM(const char *filename = "test_sim.root", TString type = "x")
{
    // Will add TChain later
    f = new TFile(filename);
    T = (TTree *)f->Get("T");

    gStyle->SetOptStat("emr");

    TH2F *h11[13];
    TF1 *f11[13];
    TGraph *g11[2];
    g11[0] = new TGraph(13);
    g11[1] = new TGraph(13);

    for (Int_t i = 0; i < 13; i++) {
        h11[i] = new TH2F(Form("h11%d", i), Form("h11%d", i), 400, 0, 3.0, 400, -10.0, 20.0);
        Double_t x = 2.0e-3 * i - 12.0e-3;

        if (type == "x")
            T->Project(Form("h11%d", i), "(gen.tp.x-fwd.tp.proj.x)*1000:(1+fwd.tp.proj.d)*2.0", Form("isgood&&TMath::Abs(gen.tp.y-%lf)<1.0e-3", x));
        else if (type == "y")
            T->Project(Form("h11%d", i), "(gen.tp.y-fwd.tp.proj.y)*1000:(1+fwd.tp.proj.d)*2.0", Form("isgood&&TMath::Abs(gen.tp.x-%lf)<1.0e-3", x));
        else return;

        f11[i] = new TF1(Form("f11%d", i), "[0]+[1]/x", 0.4, 3.0);
        h11[i]->Fit(f11[i], "", "", 0.4, 3.0);
        h11[i]->SetTitle(Form("diff = (%lf) + (%lf) / p", f11[i]->GetParameter(0), f11[i]->GetParameter(1)));
        g11[0]->SetPoint(i, x, f11[i]->GetParameter(0));
        g11[1]->SetPoint(i, x, f11[i]->GetParameter(1));
    }

    TCanvas *c11 = new TCanvas("c11", "BPM", 1200, 600);
    c11->Divide(2, 1);

    c11->cd(1);
    g11[0]->SetTitle("Par0");
    g11[0]->Draw("A*");
    TF1 *f110 = new TF1("f110", "pol0", -0.015, 0.015);
    g11[0]->Fit(f110);
    Double_t p0 = f110->GetParameter(0);
    c11->Update();

    c11->cd(2);
    g11[1]->SetTitle("Par1 & Par2");
    g11[1]->Draw("A*");
    TF1 *f111 = new TF1("f111", "pol1", -0.015, 0.015);
    g11[1]->Fit(f111);
    Double_t p1 = f111->GetParameter(0);
    Double_t p2 = f111->GetParameter(1);
    c11->Update();

    TCanvas *c01 = new TCanvas("c01", "BPM", 1200, 600);
    c01->Divide(2, 1);

    c01->cd(1);

    if (type == "x")
        T->Draw("(gen.tp.x-fwd.tp.proj.x)*1000:(1+fwd.tp.proj.d)*2.0:gen.tp.y", "isgood");

    if (type == "y")
        T->Draw("(gen.tp.y-fwd.tp.proj.y)*1000:(1+fwd.tp.proj.d)*2.0:gen.tp.x", "isgood");

    c01->Update();

    c01->cd(2);
    TF2 *f01 = new TF2("f01", Form("(%lf)+((%lf)+(%lf)*x)/y", p0, p1, p2), -0.015, 0.015, 0, 3.5);
    f01->Draw("surf");
    c01->Update();
}
