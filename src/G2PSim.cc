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
#include "TList.h"
#include "TFile.h"
#include "TTree.h"

#include "G2PAppsBase.hh"
#include "G2PDrift.hh"
#include "G2PRunBase.hh"

#include "G2PSim.hh"

TList* gG2PApps = new TList();
G2PRunBase* gG2PRun = NULL;

G2PSim* G2PSim::pG2PSim = NULL;

G2PSim::G2PSim() :
    nEvent(50000), nCounter(1), fDebug(1), pFile(NULL), pOutFileName(NULL),
    pTree(NULL), fApps(NULL), pRun(NULL)
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
    TIter next(fApps);

    while (G2PAppsBase* aobj = static_cast<G2PAppsBase*>(next())) {
        fApps->Remove(aobj);
        aobj->Delete();
    }

    if (pG2PSim==this) pG2PSim = NULL;
}

void G2PSim::SetSeed(int n)
{
    G2PAppsBase::SetSeed(n);
}

bool G2PSim::Init()
{
    static const char* const here = "Init()";

    pRun = G2PRunBase::GetInstance();
    gG2PRun = pRun;

    TIter next(fApps);
    while (TObject* obj = next()) {
        if (obj->IsZombie()) gG2PApps->Remove(obj);
    }
    next.Reset();
    while (G2PAppsBase* aobj = static_cast<G2PAppsBase*>(next())) {
        aobj->RegisterModel();
        aobj->SetDebug(fDebug);
    }
    next.Reset();
    while (G2PAppsBase* aobj = static_cast<G2PAppsBase*>(next())) {
        if (!aobj->IsInit()) {
            if (aobj->Init()) return false;
        }
    }

    if (fDebug>0) Info(here, "Initialize done!");

    return true;
}

void G2PSim::Begin()
{
    pFile = new TFile(pOutFileName, "recreate");
    pTree = new TTree("T", "Sim result");

    pTree->Branch("Index", &nCounter, "Index/I");

    TIter next(gG2PApps);
    while (G2PAppsBase* aobj = static_cast<G2PAppsBase*>(next())) {
        aobj->DefineVariables(pTree);
    }
}

void G2PSim::End()
{
    pFile->Write("", TObject::kOverwrite);
    pFile->Close();

    //delete pTree;
}

void G2PSim::Run()
{
    static const char* const here = "Run()";

    if (Init()) {
        Begin();
        while (nCounter<=nEvent) {
            if (fDebug>1) Info(here, "Processing event %d ...", nCounter);
            if (pRun->Run()) break;
            pTree->Fill();
            if ((nCounter%100==0)&&(fDebug>0)) Info(here, "%d events processed ...", nCounter);
            nCounter++;
        }
        End();
    }
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
