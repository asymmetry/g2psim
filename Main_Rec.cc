// -*- C++ -*-

/* RecMain.cc
 * Main function of the reconstruction program.
 * It parses command line parameters use getopt (GNU C Library).
 */

// History:
//   Jan 2013, C. Gu, First public version.
//   Mar 2013, C. Gu, Use Run.C script to do simulation
//   Sep 2013, C. Gu, Since G2PRun will take care of the configuration files, make it simpler here.
//   Feb 2014, C. Gu, Modified for reconstruction.
//   Jan 2015, C. Gu, Move this file back to the g2psim package.
//

#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cstring>
#include <ctype.h>
#include <getopt.h>
#include <time.h>
#include <unistd.h>
#include <vector>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TChain.h"
#include "THaEvent.h"

#include "G2PAppBase.hh"
#include "G2PAppList.hh"
#include "G2PGlobals.hh"
#include "G2PRec.hh"
#include "G2PRun.hh"
#include "G2PVarDef.hh"
#include "G2PVarList.hh"

//#define DEBUG

using namespace std;

int fRun = 0;
int fGEP = 0;
int fDebug = 1;
int fSeparate = 0;

const char *DBDIR = "./recdb";
const char *ROOTDIR = ".";

int fbpmavail = 0;
float fV5bpm_bpm[5];
float fV5bpmave_bpm[5];

double fV5tpmat_tr[5];
double fV5tpcorr_tr[5];

double frecz_lab = 0;
double fV5rec_tr[5];
double fV5rec_lab[5];

G2PRec *fRec = NULL;

int Begin();
int Process();
int Process_separate();
void Clear();
void Usage(int argc, char **argv);

int main(int argc, char **argv)
{
    int c;

    while (1) {
        static struct option long_options[] = {
            {"dbdir", required_argument, 0, 'd'},
            {"help", no_argument, 0, 'h'},
            {"rootdir", required_argument, 0, 'r'},
            {"separate", no_argument, 0, 's'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long(argc, argv, "d:hr:s", long_options, &option_index);

        if (c == -1) break;

        switch (c) {
        case 0:
            break;

        case 'd':
            DBDIR = Form("%s", optarg);
            break;

        case 'h':
            Usage(argc, argv);
            exit(0);

        case 'r':
            ROOTDIR = Form("%s", optarg);
            break;
            
        case 's':
            fSeparate = 1;
            break;

        case '?':
        default:
            Usage(argc, argv);
            exit(-1);
        }
    }

    if (optind < argc)
        fRun = atoi(argv[optind++]);
    else {
        Usage(argc, argv);
        exit(-1);
    }

    if (access(Form("%s/g2p_%d.root", ROOTDIR, fRun), F_OK) != 0) {
        Usage(argc, argv);
        exit(-1);
    }

    gG2PRun = new G2PRun();

    if (Begin() != 0) {
        Error("Main()", "Please check database!");
        exit(-1);
    }

    gG2PRun->Begin();

    ConfDef debug = {"run.debuglevel", "Debug Level", kINT, &fDebug};
    gG2PRun->GetConfig(&debug, "");
    ConfDef gep = {"run.gep", "GEP Flag", kINT, &fGEP};
    gG2PRun->GetConfig(&gep, "");

    fRec = new G2PRec();
    gG2PApps->Add(fRec);

    TIter next(gG2PApps);

    while (G2PAppBase *aobj = static_cast<G2PAppBase *>(next())) {
        if (aobj->IsZombie()) {
            gG2PApps->Remove(aobj);
            continue;
        }

        if (!aobj->IsInit()) {
            aobj->SetDebugLevel(fDebug);

            if (aobj->Begin() != 0)
                return -1;
        }
    }

    if (fDebug > 0)
        gG2PRun->Print();

    Clear();

    if (fDebug > 0)
        Info("Main()", "Ready to go!");

    if(fSeparate) Process_separate();
    else Process();

    next.Reset();

    while (G2PAppBase *aobj = static_cast<G2PAppBase *>(next()))
        aobj->End();

    if (fDebug > 0)
        Info("Main()", "Run finished!");

    return 0;
}

int Begin()
{
    static const char *const here = "Main::Begin()";

    TFile *f = new TFile(Form("%s/g2p_%d.root", ROOTDIR, fRun), "READ");
    TTree *t = (TTree *) f->Get("T");

    THaEvent *event = new THaEvent();
    t->SetBranchAddress("Event_Branch", &event);

    t->GetEntry(0);

    int label = int(event->GetHeader()->GetEvtTime() / 1.0e6);

    delete event;
    event = NULL;

    f->Close();

    void *dirp = gSystem->OpenDirectory(DBDIR);

    if (dirp == NULL) return -1;

    const char *result;
    string item;
    vector<int> time;

    while ((result = gSystem->GetDirEntry(dirp)) != NULL) {
        item = result;
        size_t pos;

        for (pos = 0; pos < item.length(); pos++) {
            if (!isdigit(item[pos])) break;
        }

        if (pos == item.length())
            time.push_back(atoi(item.c_str()));
    }

    gSystem->FreeDirectory(dirp);

    if (time.size() > 0) {
        sort(time.begin(), time.end());
        vector<int>::iterator it;

        for (it = time.begin(); it != time.end(); it++) {
            if (label < (*it)) {
                if (it == time.begin()) return -1;

                it--;
                break;
            }
        }

        if (it == time.end()) it--; // Assume the last directory is always valid

        // Only work for g2p
        if (fRun < 20000) {
            Info(here, "Use %s", Form("%s/%d/db_L.optics.cfg", DBDIR, (*it)));
            gG2PRun->SetConfigFile(Form("%s/%d/db_L.optics.cfg", DBDIR, (*it)));
        } else if (fRun < 40000) {
            Info(here, "Use %s", Form("%s/%d/db_R.optics.cfg", DBDIR, (*it)));
            gG2PRun->SetConfigFile(Form("%s/%d/db_R.optics.cfg", DBDIR, (*it)));
        } else return -1;
    }

    const char *rundbfile = Form("%s/db_rec.dat", DBDIR);

    if (access(rundbfile, F_OK) != 0)
        return -1;

    FILE *fp = fopen(rundbfile, "r");

    int id = 0;
    double e = 0.0, p = 0.0;
    bool found = false;

    while (!feof(fp)) {
        fscanf(fp, "%d%le%le", &id, &p, &e);

        if (id == fRun) {
            gG2PRun->SetBeamEnergy(e);
            gG2PRun->SetHRSMomentum(p);
            found = true;
            break;
        }
    }

    if (!found) return -1;

    fclose(fp);

    VarDef vars[] = {
        {"bpm.b_x", "BPM X", kFLOAT, &fV5bpm_bpm[0]},
        {"bpm.b_t", "BPM T", kFLOAT, &fV5bpm_bpm[1]},
        {"bpm.b_y", "BPM Y", kFLOAT, &fV5bpm_bpm[2]},
        {"bpm.b_p", "BPM P", kFLOAT, &fV5bpm_bpm[3]},
        {"bpm.b_z", "BPM Z", kFLOAT, &fV5bpm_bpm[4]},
        {"bpm.ave.b_x", "Average BPM X", kFLOAT, &fV5bpmave_bpm[0]},
        {"bpm.ave.b_y", "Average BPM Y", kFLOAT, &fV5bpmave_bpm[2]},
        {"tp.mat.x", "Matrix Rec to Target Plane X", kDOUBLE, &fV5tpmat_tr[0]},
        {"tp.mat.t", "Matrix Rec to Target Plane T", kDOUBLE, &fV5tpmat_tr[1]},
        {"tp.mat.y", "Matrix Rec to Target Plane Y", kDOUBLE, &fV5tpmat_tr[2]},
        {"tp.mat.p", "Matrix Rec to Target Plane P", kDOUBLE, &fV5tpmat_tr[3]},
        {"tp.mat.d", "Matrix Rec to Target Plane D", kDOUBLE, &fV5tpmat_tr[4]},
        {0}
    };

    if (!gG2PVars) {
        Error(here, "No global variable list.");
        return -1;
    } else
        gG2PVars->DefineVariables(vars, "");

    return 0;
}

int Process()
{
    static const char *const here = "Main::Process()";

    const char *arm = "R";

    if (fRun < 20000) arm = "L";
    else if (fRun < 40000) arm = "R";

    const char *filename = Form("%s/g2p_%d.root", ROOTDIR, fRun);
    int inc = 0;

    while (access(filename, F_OK) == 0) {
        TFile *f = new TFile(filename, "UPDATE");
        Info(here, "Opening existed rootfile %s ...", filename);

        TTree *t = (TTree *) f->Get("T");

        THaEvent *event = new THaEvent();
        t->SetBranchAddress("Event_Branch", &event);

        if (fGEP == 0) {
            if (t->FindBranch(Form("%srb.tgt_m13_x", arm))) {
                t->SetBranchAddress(Form("%srb.tgt_m13_x", arm), &fV5bpm_bpm[0]);
                t->SetBranchAddress(Form("%srb.tgt_m13_theta", arm), &fV5bpm_bpm[1]);
                t->SetBranchAddress(Form("%srb.tgt_m13_y", arm), &fV5bpm_bpm[2]);
                t->SetBranchAddress(Form("%srb.tgt_m13_phi", arm), &fV5bpm_bpm[3]);
                t->SetBranchAddress(Form("%srb.tgtave_m13_x", arm), &fV5bpmave_bpm[0]);
                t->SetBranchAddress(Form("%srb.tgtave_m13_y", arm), &fV5bpmave_bpm[2]);
                frecz_lab = -13.6271e-3;
            } else if (t->FindBranch(Form("%srb.tgt_m12_x", arm))) {
                t->SetBranchAddress(Form("%srb.tgt_m12_x", arm), &fV5bpm_bpm[0]);
                t->SetBranchAddress(Form("%srb.tgt_m12_theta", arm), &fV5bpm_bpm[1]);
                t->SetBranchAddress(Form("%srb.tgt_m12_y", arm), &fV5bpm_bpm[2]);
                t->SetBranchAddress(Form("%srb.tgt_m12_phi", arm), &fV5bpm_bpm[3]);
                t->SetBranchAddress(Form("%srb.tgtave_m12_x", arm), &fV5bpmave_bpm[0]);
                t->SetBranchAddress(Form("%srb.tgtave_m12_y", arm), &fV5bpmave_bpm[2]);
                frecz_lab = -12.5476e-3;
            } else if (t->FindBranch(Form("%srb.tgt_m10_x", arm))) {
                t->SetBranchAddress(Form("%srb.tgt_m10_x", arm), &fV5bpm_bpm[0]);
                t->SetBranchAddress(Form("%srb.tgt_m10_theta", arm), &fV5bpm_bpm[1]);
                t->SetBranchAddress(Form("%srb.tgt_m10_y", arm), &fV5bpm_bpm[2]);
                t->SetBranchAddress(Form("%srb.tgt_m10_phi", arm), &fV5bpm_bpm[3]);
                t->SetBranchAddress(Form("%srb.tgtave_m10_x", arm), &fV5bpmave_bpm[0]);
                t->SetBranchAddress(Form("%srb.tgtave_m10_y", arm), &fV5bpmave_bpm[2]);
                frecz_lab = -10.81e-3;
            } else {
                t->SetBranchAddress(Form("%srb.tgt_0_x", arm), &fV5bpm_bpm[0]);
                t->SetBranchAddress(Form("%srb.tgt_0_theta", arm), &fV5bpm_bpm[1]);
                t->SetBranchAddress(Form("%srb.tgt_0_y", arm), &fV5bpm_bpm[2]);
                t->SetBranchAddress(Form("%srb.tgt_0_phi", arm), &fV5bpm_bpm[3]);
                t->SetBranchAddress(Form("%srb.tgtave_0_x", arm), &fV5bpmave_bpm[0]);
                t->SetBranchAddress(Form("%srb.tgtave_0_y", arm), &fV5bpmave_bpm[2]);
                frecz_lab = 0.0;
            }

            t->SetBranchAddress(Form("%srb.bpmavail", arm), &fbpmavail);
        } else if (fGEP == 1)
            fbpmavail = 1;

        t->SetBranchAddress(Form("%s.gold.x", arm), &fV5tpmat_tr[0]);
        t->SetBranchAddress(Form("%s.gold.th", arm), &fV5tpmat_tr[1]);
        t->SetBranchAddress(Form("%s.gold.y", arm), &fV5tpmat_tr[2]);
        t->SetBranchAddress(Form("%s.gold.ph", arm), &fV5tpmat_tr[3]);
        t->SetBranchAddress(Form("%s.gold.dp", arm), &fV5tpmat_tr[4]);

        TList newBranch;

        newBranch.Add(t->Branch(Form("%s.rec.x", arm), &fV5rec_tr[0], Form("%s.rec.x/D", arm)));
        newBranch.Add(t->Branch(Form("%s.rec.th", arm), &fV5rec_tr[1], Form("%s.rec.th/D", arm)));
        newBranch.Add(t->Branch(Form("%s.rec.y", arm), &fV5rec_tr[2], Form("%s.rec.y/D", arm)));
        newBranch.Add(t->Branch(Form("%s.rec.ph", arm), &fV5rec_tr[3], Form("%s.rec.ph/D", arm)));
        newBranch.Add(t->Branch(Form("%s.rec.dp", arm), &fV5rec_tr[4], Form("%s.rec.dp/D", arm)));

        newBranch.Add(t->Branch(Form("%s.rec.l_x", arm), &fV5rec_lab[0], Form("%s.rec.l_x/D", arm)));
        newBranch.Add(t->Branch(Form("%s.rec.l_th", arm), &fV5rec_lab[1], Form("%s.rec.l_th/D", arm)));
        newBranch.Add(t->Branch(Form("%s.rec.l_y", arm), &fV5rec_lab[2], Form("%s.rec.l_y/D", arm)));
        newBranch.Add(t->Branch(Form("%s.rec.l_ph", arm), &fV5rec_lab[3], Form("%s.rec.l_ph/D", arm)));
        newBranch.Add(t->Branch(Form("%s.rec.l_z", arm), &fV5rec_lab[4], Form("%s.rec.l_z/D", arm)));

#ifdef DEBUG
        newBranch.Add(t->Branch(Form("%s.corr.x", arm), &fV5tpcorr_tr[0], Form("%s.corr.x/D", arm)));
        newBranch.Add(t->Branch(Form("%s.corr.th", arm), &fV5tpcorr_tr[1], Form("%s.corr.th/D", arm)));
        newBranch.Add(t->Branch(Form("%s.corr.y", arm), &fV5tpcorr_tr[2], Form("%s.corr.y/D", arm)));
        newBranch.Add(t->Branch(Form("%s.corr.ph", arm), &fV5tpcorr_tr[3], Form("%s.corr.ph/D", arm)));
        newBranch.Add(t->Branch(Form("%s.corr.dp", arm), &fV5tpcorr_tr[4], Form("%s.corr.dp/D", arm)));
#endif

        int N = t->GetEntries();
        int evnum;
        TIter next(&newBranch);

        for (int i = 0; i < N; i++) {
            t->GetEntry(i);
            evnum = Int_t(event->GetHeader()->GetEvtNum());

            if (fDebug > 1)
                Info(here, "Processing event %d ......", evnum);
            else if ((i % 10000 == 0) && (i != 0))
                Info(here, "%d events Processed ......", i);

            if ((fbpmavail < 0.5) || (fV5tpmat_tr[0] > 1.0e8) || (fV5tpmat_tr[1] > 1.0e8) || (fV5tpmat_tr[2] > 1.0e8) || (fV5tpmat_tr[3] > 1.0e8) || (fV5tpmat_tr[4] > 1.0e8)) {
                for (int i = 0; i < 5; i++) {
                    fV5rec_tr[i] = 1e+38;
                    fV5rec_lab[i] = 1e+38;
#ifdef DEBUG
                    fV5tpcorr_tr[i] = 1e+38;
#endif
                }
            } else {
                fV5bpm_bpm[0] /= 1000.0;
                fV5bpm_bpm[2] /= 1000.0;
                fV5bpm_bpm[4] = frecz_lab;
                fV5bpmave_bpm[0] /= 1000.0;
                fV5bpmave_bpm[2] /= 1000.0;

                fRec->Process();

                if (gG2PVars->FindSuffix("rec.x") && gG2PVars->FindSuffix("rec.l_x")
#ifdef DEBUG
                        && gG2PVars->FindSuffix("tp.corr.x")
#endif
                   ) {
                    fV5rec_tr[0] = gG2PVars->FindSuffix("rec.x")->GetValue();
                    fV5rec_tr[1] = tan(gG2PVars->FindSuffix("rec.t")->GetValue());
                    fV5rec_tr[2] = gG2PVars->FindSuffix("rec.y")->GetValue();
                    fV5rec_tr[3] = tan(gG2PVars->FindSuffix("rec.p")->GetValue());
                    fV5rec_tr[4] = gG2PVars->FindSuffix("rec.d")->GetValue();

                    fV5rec_lab[0] = gG2PVars->FindSuffix("rec.l_x")->GetValue();
                    fV5rec_lab[1] = gG2PVars->FindSuffix("rec.l_t")->GetValue();
                    fV5rec_lab[2] = gG2PVars->FindSuffix("rec.l_y")->GetValue();
                    fV5rec_lab[3] = gG2PVars->FindSuffix("rec.l_p")->GetValue();
                    fV5rec_lab[4] = gG2PVars->FindSuffix("rec.l_z")->GetValue();
#ifdef DEBUG
                    fV5tpcorr_tr[0] = gG2PVars->FindSuffix("tp.corr.x")->GetValue();
                    fV5tpcorr_tr[1] = tan(gG2PVars->FindSuffix("tp.corr.t")->GetValue());
                    fV5tpcorr_tr[2] = gG2PVars->FindSuffix("tp.corr.y")->GetValue();
                    fV5tpcorr_tr[3] = tan(gG2PVars->FindSuffix("tp.corr.p")->GetValue());
                    fV5tpcorr_tr[4] = gG2PVars->FindSuffix("tp,corr.d")->GetValue();
#endif
                } else
                    return -1;
            }

            next.Reset();

            while (TBranch *br = (TBranch *) next())
                br->Fill();
        }

        t->Write("", TObject::kOverwrite);
        f->Close();

        inc++;
        filename = Form("%s/g2p_%d_%d.root", ROOTDIR, fRun, inc);
    }

    Info(here, "No more rootfiles.");

    return 0;
}

int Process_separate()
{
    static const char *const here = "Main::Process_separate()";
    const char *arm = "R";
    if (fRun < 20000) arm = "L";
    else if (fRun < 40000) arm = "R";
    const char *optfilename=Form("%s/kin/optics_%d.root", ROOTDIR, fRun);
    
    TChain *t=new TChain("T");
    t->Add(Form("%s/g2p_%d*.root",ROOTDIR,fRun));
    t->AddFriend("T",Form("%s/kin/bpm_%d.root",ROOTDIR,fRun));
    
    TFile *newf=new TFile(optfilename, "RECREATE");
    TTree *newt=new TTree("T","optics");
    
    THaEvent *event = new THaEvent();
    t->SetBranchAddress("Event_Branch", &event);
    
    if (fGEP == 0) {
        if (t->FindBranch(Form("%srb.tgt_m13_x", arm))) {
            t->SetBranchAddress(Form("%srb.tgt_m13_x", arm), &fV5bpm_bpm[0]);
            t->SetBranchAddress(Form("%srb.tgt_m13_theta", arm), &fV5bpm_bpm[1]);
            t->SetBranchAddress(Form("%srb.tgt_m13_y", arm), &fV5bpm_bpm[2]);
            t->SetBranchAddress(Form("%srb.tgt_m13_phi", arm), &fV5bpm_bpm[3]);
            t->SetBranchAddress(Form("%srb.tgtave_m13_x", arm), &fV5bpmave_bpm[0]);
            t->SetBranchAddress(Form("%srb.tgtave_m13_y", arm), &fV5bpmave_bpm[2]);
            frecz_lab = -13.6271e-3;
        } else if (t->FindBranch(Form("%srb.tgt_m12_x", arm))) {
            t->SetBranchAddress(Form("%srb.tgt_m12_x", arm), &fV5bpm_bpm[0]);
            t->SetBranchAddress(Form("%srb.tgt_m12_theta", arm), &fV5bpm_bpm[1]);
            t->SetBranchAddress(Form("%srb.tgt_m12_y", arm), &fV5bpm_bpm[2]);
            t->SetBranchAddress(Form("%srb.tgt_m12_phi", arm), &fV5bpm_bpm[3]);
            t->SetBranchAddress(Form("%srb.tgtave_m12_x", arm), &fV5bpmave_bpm[0]);
            t->SetBranchAddress(Form("%srb.tgtave_m12_y", arm), &fV5bpmave_bpm[2]);
            frecz_lab = -12.5476e-3;
        } else if (t->FindBranch(Form("%srb.tgt_m10_x", arm))) {
            t->SetBranchAddress(Form("%srb.tgt_m10_x", arm), &fV5bpm_bpm[0]);
            t->SetBranchAddress(Form("%srb.tgt_m10_theta", arm), &fV5bpm_bpm[1]);
            t->SetBranchAddress(Form("%srb.tgt_m10_y", arm), &fV5bpm_bpm[2]);
            t->SetBranchAddress(Form("%srb.tgt_m10_phi", arm), &fV5bpm_bpm[3]);
            t->SetBranchAddress(Form("%srb.tgtave_m10_x", arm), &fV5bpmave_bpm[0]);
            t->SetBranchAddress(Form("%srb.tgtave_m10_y", arm), &fV5bpmave_bpm[2]);
            frecz_lab = -10.81e-3;
        } else {
            t->SetBranchAddress(Form("%srb.tgt_0_x", arm), &fV5bpm_bpm[0]);
            t->SetBranchAddress(Form("%srb.tgt_0_theta", arm), &fV5bpm_bpm[1]);
            t->SetBranchAddress(Form("%srb.tgt_0_y", arm), &fV5bpm_bpm[2]);
            t->SetBranchAddress(Form("%srb.tgt_0_phi", arm), &fV5bpm_bpm[3]);
            t->SetBranchAddress(Form("%srb.tgtave_0_x", arm), &fV5bpmave_bpm[0]);
            t->SetBranchAddress(Form("%srb.tgtave_0_y", arm), &fV5bpmave_bpm[2]);
            frecz_lab = 0.0;
        }
        
        t->SetBranchAddress(Form("%srb.bpmavail", arm), &fbpmavail);
    } else if (fGEP == 1)
        fbpmavail = 1;
        
    t->SetBranchAddress(Form("%s.gold.x", arm), &fV5tpmat_tr[0]);
    t->SetBranchAddress(Form("%s.gold.th", arm), &fV5tpmat_tr[1]);
    t->SetBranchAddress(Form("%s.gold.y", arm), &fV5tpmat_tr[2]);
    t->SetBranchAddress(Form("%s.gold.ph", arm), &fV5tpmat_tr[3]);
    t->SetBranchAddress(Form("%s.gold.dp", arm), &fV5tpmat_tr[4]);
    
    TList newBranch;

    newBranch.Add(newt->Branch(Form("%s.rec.x", arm), &fV5rec_tr[0], Form("%s.rec.x/D", arm)));
    newBranch.Add(newt->Branch(Form("%s.rec.th", arm), &fV5rec_tr[1], Form("%s.rec.th/D", arm)));
    newBranch.Add(newt->Branch(Form("%s.rec.y", arm), &fV5rec_tr[2], Form("%s.rec.y/D", arm)));
    newBranch.Add(newt->Branch(Form("%s.rec.ph", arm), &fV5rec_tr[3], Form("%s.rec.ph/D", arm)));
    newBranch.Add(newt->Branch(Form("%s.rec.dp", arm), &fV5rec_tr[4], Form("%s.rec.dp/D", arm)));
    
    newBranch.Add(newt->Branch(Form("%s.rec.l_x", arm), &fV5rec_lab[0], Form("%s.rec.l_x/D", arm)));
    newBranch.Add(newt->Branch(Form("%s.rec.l_th", arm), &fV5rec_lab[1], Form("%s.rec.l_th/D", arm)));
    newBranch.Add(newt->Branch(Form("%s.rec.l_y", arm), &fV5rec_lab[2], Form("%s.rec.l_y/D", arm)));
    newBranch.Add(newt->Branch(Form("%s.rec.l_ph", arm), &fV5rec_lab[3], Form("%s.rec.l_ph/D", arm)));
    newBranch.Add(newt->Branch(Form("%s.rec.l_z", arm), &fV5rec_lab[4], Form("%s.rec.l_z/D", arm)));

#ifdef DEBUG
    newBranch.Add(newt->Branch(Form("%s.corr.x", arm), &fV5tpcorr_tr[0], Form("%s.corr.x/D", arm)));
    newBranch.Add(newt->Branch(Form("%s.corr.th", arm), &fV5tpcorr_tr[1], Form("%s.corr.th/D", arm)));
    newBranch.Add(newt->Branch(Form("%s.corr.y", arm), &fV5tpcorr_tr[2], Form("%s.corr.y/D", arm)));
    newBranch.Add(newt->Branch(Form("%s.corr.ph", arm), &fV5tpcorr_tr[3], Form("%s.corr.ph/D", arm)));
    newBranch.Add(newt->Branch(Form("%s.corr.dp", arm), &fV5tpcorr_tr[4], Form("%s.corr.dp/D", arm)));
#endif    

    int N = t->GetEntries();
    int evnum;
    TIter next(&newBranch);

    for (int i = 0; i < N; i++) {
        t->GetEntry(i);
        evnum = Int_t(event->GetHeader()->GetEvtNum());
        
        if (fDebug > 1)
            Info(here, "Processing event %d ......", evnum);
        else if ((i % 10000 == 0) && (i != 0))
            Info(here, "%d events Processed ......", i);

        if ((fbpmavail < 0.5) || (fV5tpmat_tr[0] > 1.0e8) || (fV5tpmat_tr[1] > 1.0e8) || (fV5tpmat_tr[2] > 1.0e8) || (fV5tpmat_tr[3] > 1.0e8) || (fV5tpmat_tr[4] > 1.0e8)) {
            for (int i = 0; i < 5; i++) {
                fV5rec_tr[i] = 1e+38;
                fV5rec_lab[i] = 1e+38;
#ifdef DEBUG
                fV5tpcorr_tr[i] = 1e+38;
#endif
            }
        } else {
            fV5bpm_bpm[0] /= 1000.0;
            fV5bpm_bpm[2] /= 1000.0;
            fV5bpm_bpm[4] = frecz_lab;
            fV5bpmave_bpm[0] /= 1000.0;
            fV5bpmave_bpm[2] /= 1000.0;

            fRec->Process();

            if (gG2PVars->FindSuffix("rec.x") && gG2PVars->FindSuffix("rec.l_x")
#ifdef DEBUG
                    && gG2PVars->FindSuffix("tp.corr.x")
#endif             
               ) {
                fV5rec_tr[0] = gG2PVars->FindSuffix("rec.x")->GetValue();
                fV5rec_tr[1] = tan(gG2PVars->FindSuffix("rec.t")->GetValue());
                fV5rec_tr[2] = gG2PVars->FindSuffix("rec.y")->GetValue();
                fV5rec_tr[3] = tan(gG2PVars->FindSuffix("rec.p")->GetValue());
                fV5rec_tr[4] = gG2PVars->FindSuffix("rec.d")->GetValue();
                
                fV5rec_lab[0] = gG2PVars->FindSuffix("rec.l_x")->GetValue();
                fV5rec_lab[1] = gG2PVars->FindSuffix("rec.l_t")->GetValue();
                fV5rec_lab[2] = gG2PVars->FindSuffix("rec.l_y")->GetValue();
                fV5rec_lab[3] = gG2PVars->FindSuffix("rec.l_p")->GetValue();
                fV5rec_lab[4] = gG2PVars->FindSuffix("rec.l_z")->GetValue();
#ifdef DEBUG
                fV5tpcorr_tr[0] = gG2PVars->FindSuffix("tp.corr.x")->GetValue();
                fV5tpcorr_tr[1] = tan(gG2PVars->FindSuffix("tp.corr.t")->GetValue());
                fV5tpcorr_tr[2] = gG2PVars->FindSuffix("tp.corr.y")->GetValue();
                fV5tpcorr_tr[3] = tan(gG2PVars->FindSuffix("tp.corr.p")->GetValue());
                fV5tpcorr_tr[4] = gG2PVars->FindSuffix("tp,corr.d")->GetValue();
#endif
            } else
                return -1;
        }

        next.Reset();

        newt->Fill();
    }
    newt->Write("", TObject::kOverwrite);
    newf->Close();
    return 0;
}

void Clear()
{
    fbpmavail = 0;

    memset(fV5bpm_bpm, 0, sizeof(fV5bpm_bpm));
    memset(fV5bpmave_bpm, 0, sizeof(fV5bpmave_bpm));
    memset(fV5tpmat_tr, 0, sizeof(fV5tpmat_tr));
}

void Usage(int argc, char **argv)
{
    printf("usage: %s [options] RunNumber \n", argv[0]);
    printf("  -d, --dbdir=$PWD/recdb         Set db directory\n");
    printf("  -h, --help                     Print this small usage guide\n");
    printf("  -s, --separate                 Separate insert the optics information\n");
    printf("  -r, --rootdir=$PWD             Set db directory\n");
    
}
