// This file defines a class G2PSim.
// This class is the main class for this sim package.
//
// History:
//   Jan 2013, C. Gu, First public version.
//

#include <cstdlib>
#include <cmath>

#include "TROOT.h"
#include "TObject.h"
#include "TError.h"
#include "TFile.h"
#include "TList.h"

#include "G2PAppBase.hh"
#include "G2POutput.hh"
#include "G2PRunBase.hh"
#include "G2PVarDef.hh"
#include "G2PVarList.hh"

#include "G2PSim.hh"

TList* gG2PApps = new TList();
G2PRunBase* gG2PRun = NULL;
G2PVarList* gG2PVars = new G2PVarList();

G2PSim* G2PSim::pG2PSim = NULL;

G2PSim::G2PSim() :
    fDebug(1), fFile(NULL), pOutFile(NULL),
    nEvent(50000), nCounter(1), 
    fApps(NULL), pRun(NULL), pOutput(NULL)
{
    if (pG2PSim) {
        Error("G2PSim::G2PSim()", "Only one instance of G2PSim allowed.");
        MakeZombie();
        return;
    }
    pG2PSim = this;

    fApps = gG2PApps;
}

G2PSim::~G2PSim()
{
    if (pG2PSim==this) pG2PSim = NULL;

    TIter next(fApps);
    while (G2PAppBase* aobj = static_cast<G2PAppBase*>(next())) {
        fApps->Remove(aobj);
        aobj->Delete();
    }

    if (pOutput) delete pOutput;
}

void G2PSim::SetSeed(int n)
{
    G2PAppBase::SetSeed(n);
}

int G2PSim::Init()
{
    static const char* const here = "Init()";

    gG2PRun = pRun;
    gG2PVars->DefineByType("event", "Event number", &nCounter, kInt);

    TIter next(fApps);
    while (TObject* obj = next()) {
        if (obj->IsZombie()) gG2PApps->Remove(obj);
    }
    
    next.Reset();
    while (G2PAppBase* aobj = static_cast<G2PAppBase*>(next())) {
        aobj->SetDebug(fDebug);
    }
    pRun->SetDebug(fDebug);

    next.Reset();
    while (G2PAppBase* aobj = static_cast<G2PAppBase*>(next())) {
        if (aobj->Init()!=0) return -1;
    }

    if (pRun->Init()!=0) return -1;

    fFile = new TFile(pOutFile, "RECREATE");

    pOutput = new G2POutput();
    if (pOutput->Init()!=0) return -1;

    if (fDebug>0) Info(here, "Initialize done!");

    return 0;
}

int G2PSim::Begin()
{
    static const char* const here = "Begin()";

    if (fDebug>0) Info(here, "Beginning ......");
    
    TIter next(fApps);
    while (G2PAppBase* aobj = static_cast<G2PAppBase*>(next())) {
        if (aobj->Begin()!=0) return -1;
    }

    if (pRun->Begin()!=0) return -1;

    return 0;
}

int G2PSim::End()
{
    static const char* const here = "End()";

    if (fDebug>0) Info(here, "Cleaning ......");
    
    pOutput->End();
    fFile->Close();

    return 0;
}

void G2PSim::Run()
{
    static const char* const here = "Run()";

    if (Init()!=0) {
        Error(here, "Cannot initialize, program will stop.");
        return;
    }

    if (Begin()!=0) {
        Error(here, "Critical error, program will stop.");
        return;
    }

    fFile->cd();

    while (nCounter<=nEvent) {
        if (fDebug>1) Info(here, "Processing event %d ...", nCounter);
        pRun->Clear();
        if (pRun->Process()!=0) break;
        pOutput->Process();
        if ((nCounter%100==0)&&(fDebug>0)) Info(here, "%d events processed ...", nCounter);
        nCounter++;
    }

    End();
}

// void G2PSim::InitTree()
// {
//     pFile = new TFile(pFileName, "recreate");
// 	pTree = new TTree("T", "sim result");

//     vector<G2PGun*>::iterator it = pGunList.begin();
//     int i = 0;
//     while (it!=pGunList.end()){
//         pGun = (*it);
//         pConfig = new TTree(Form("config%d", i), "sim configure");

//         pConfig->Branch("NEvent", &nEvent, "NEvent/I");

//         pConfig->Branch("IsLeftArm", &bIsLeftArm, "IsLeftArm/O");
//         pConfig->Branch("HRSAngle", &fHRSAngle, "HRSAngle/D");
//         pConfig->Branch("HRSMomentum", &fHRSMomentum, "HRSMomentum/D");
//         pConfig->Branch("BeamEnergy", &fBeamEnergy, "BeamEnergy/D");

//         int hrssetting = pHRS->GetModelIndex();
//         pConfig->Branch("HRSModel", &hrssetting, "HRSModel/I");

//         int gunsetting = pGun->GetSetting();
//         double resa = pBPM->GetBPMARes();
//         double resb = pBPM->GetBPMBRes();
//         pConfig->Branch("GunSetting", &gunsetting, "GunSetting/I");
//         pConfig->Branch("BPMARes", &resa, "BPMARes/D");
//         pConfig->Branch("BPMBRes", &resb, "BPMBRes/D");

//         pConfig->Fill();

//         pFile->Write("", TObject::kOverwrite);

//         delete pConfig;
        
//         it++; i++;
//     }
// }

ClassImp(G2PSim)
