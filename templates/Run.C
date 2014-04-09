// -*- C++ -*-

/* script Run.C
 * Can directly run in ROOT as a script.
 * To do this, some include paths should be set by rootlogon.C.
 */

// History:
//   Jan 2013, C. Gu, First public version.
//

#include <cstdlib>
#include <cstdio>
#include <time.h>

#include "TROOT.h"

#include "G2PBPM.hh"
#include "G2PBwdProc.hh"
#include "G2PData.hh"
#include "G2PDBRec.hh"
#include "G2PField.hh"
#include "G2PFlatGun.hh"
#include "G2PGlobals.hh"
#include "G2PGun.hh"
#include "G2PHRS.hh"
#include "G2PPhys.hh"
#include "G2PProcBase.hh"
#include "G2PRun.hh"
#include "G2PSieveGun.hh"
#include "G2PSim.hh"

static const double kDEG = 3.14159265358979323846 / 180.0;
static const double e = 1.60217656535e-19;
static const double kU = 0.931494028;

int Run()
{
    ///////////////////////////////////////////////////////////////////////////
    // run parameters
    ///////////////////////////////////////////////////////////////////////////
    G2PRun* run = new G2PRun();
    run->SetConfigFile("sim.cfg");
    run->SetRunType("optics");
    //run->SetRunType("optics21"); // 125 mil carbon, no LHe
    //run->SetRunType("optics22"); // 40 mil carbon, with LHe
    //run->SetRunType("optics23"); // 125 mil carbon, with LHe
    //run->SetRunType("empty");
    run->SetHRSAngle(5.767 * kDEG);
    run->SetHRSMomentum(2.24949);
    //run->SetHRSMomentum(2.0); // eloss
    run->SetBeamEnergy(2.25327);
    run->SetTarget(6, 12);
    run->SetTargetMass(12.0107 * kU);
    //run->SetFieldRatio(0.0);
    run->SetFieldRatio(0.5);
    //run->SetFieldRatio(1.0);
    //run->SetSieve();

    ///////////////////////////////////////////////////////////////////////////
    // field
    ///////////////////////////////////////////////////////////////////////////
    G2PField* field = new G2PField();
    field->SetEulerAngle(90*kDEG,90*kDEG,-90*kDEG);
    //field->SetEulerAngle(0,0,0);
    gG2PApps->Add(field);

    ///////////////////////////////////////////////////////////////////////////
    // gun
    ///////////////////////////////////////////////////////////////////////////
    G2PFlatGun* gun = new G2PFlatGun();
    //G2PSieveGun* gun = new G2PSieveGun();
    //G2PData* gun = new G2PData("data.dat");

    //gun->SetBeamPos(-3.623e-3, 0.0); // 2.253GeV, 0.0T, 6deg
    //gun->SetBeamPos(0.5009e-3, 0.9205e-3); // 2.254GeV, 0.0T, 90deg
    //gun->SetBeamPos(-3.897e-3, 0.016e-3); // 2.254GeV, 2.5T, 90deg
    //gun->SetBeamPos(-4.249e-3, 6.619e-3); // 2.254GeV, 2.5T, 90deg
    //gun->SetBeamPos(-0.2357e-3, -3.855e-3); // 2.254GeV, 5.0T, 0deg
    //gun->SetBeamPos(2.594e-3, -6.961e-3); // 1.158GeV, 2.5T, 0deg
    gun->SetBeamPos(0, 0, 0);

    //gun->SetTiltAngle(0, 0);

    //gun->SetReactZ(0.0, 0.0);
    //gun->SetReactZ(-14.1350e-3, -13.1191e-3); // 40mil C12
    //gun->SetReactZ(-14.1350e-3, -10.9601e-3); // 125mil C12
    gun->SetReactZ(-14.1350e-3, 14.1350e-3); // production
    //gun->SetReactZ(-14.1350e-3, -14.1340e-3); // eloss

    //gun->SetRasterSize(14.0e-3); // production
    gun->SetRasterSize(0.0); // optics

    //gun->SetTargetTh(-100.0e-3, 10.0e-3); // 2.254GeV, 2.5T, 90deg
    //gun->SetTargetPh(-20.0e-3, 30.0e-3); // 2.254GeV, 2.5T, 90deg
    //gun->SetTargetTh(-160.0e-3, -50.0e-3); // 1.158GeV, 2.5T, 90deg
    //gun->SetTargetPh(-28.0e-3, 22.0e-3); // 1.158GeV, 2.5T, 90deg
    gun->SetTargetTh(-60.0e-3, 60.0e-3); // 2.254GeV, 5.0T, 0deg
    gun->SetTargetPh(-30.0e-3, 30.0e-3); // 2.254GeV, 5.0T, 0deg

    //gun->SetDelta(-0.04, 0.04);
    gun->SetDelta("elastic"); // optics
    gG2PApps->Add(gun);

    ///////////////////////////////////////////////////////////////////////////
    // BPM
    ///////////////////////////////////////////////////////////////////////////
    G2PBPM* bpm = new G2PBPM();
    bpm->SetBPMRes(0, 0);
    //bpm->SetBPMRes(0.2e-3, 0.4e-3);
    gG2PApps->Add(bpm);

    ///////////////////////////////////////////////////////////////////////////
    // HRS
    ///////////////////////////////////////////////////////////////////////////
    G2PHRS* hrs = new G2PHRS("484816R15");
    gG2PApps->Add(hrs);
    //G2PDBRec* rec = new G2PDBRec();
    //gG2PApps->Add(rec);

    ///////////////////////////////////////////////////////////////////////////
    // cross section
    ///////////////////////////////////////////////////////////////////////////
    //G2PPhys* model = new G2PPhys("qfs");
    //double par[1] = {2};
    //model->SetPars(par, 1);
    //gG2PApps->Add(model);

    ///////////////////////////////////////////////////////////////////////////
    // special commands
    ///////////////////////////////////////////////////////////////////////////
    //G2PProcBase* aobj;
    //if (aobj = static_cast<G2PProcBase*> (gG2PApps->Find("G2PBwdProc"))) gG2PApps->Remove(aobj);

    ///////////////////////////////////////////////////////////////////////////
    // simulation
    ///////////////////////////////////////////////////////////////////////////
    G2PSim *sim = new G2PSim();
    run->SetDebugLevel(1);
    run->SetSeed(1);
    //int N = 1;
    //int N = 6550;
    int N = 20000;
    //int N = 500000;
    sim->SetNEvent(N);
    //sim->SetOutFile("run_1051_data_22540090.root");
    //sim->SetOutFile("run_1060_optics_22545000.root");
    sim->SetOutFile("run_test.root");

    clock_t start = clock();
    sim->Run();
    clock_t end = clock();
#define CLOCKS_PER_SEC 1000000
    printf("Average calculation time for one event: %8.4f ms\n", (double) (end - start)*1000.0 / (double) CLOCKS_PER_SEC / N);

    return 0;
}
