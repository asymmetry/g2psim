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
#include "G2PData.hh"
#include "G2PDBRec.hh"
#include "G2PField.hh"
#include "G2PFlatGun.hh"
#include "G2PGlobals.hh"
#include "G2PGun.hh"
#include "G2PHRS.hh"
#include "G2PPhys.hh"
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
    //run->SetRunType("optics");
    run->SetRunType("optics21"); // 125 mil carbon, no LHe
    run->SetHRSAngle(5.767 * kDEG);
    run->SetHRSMomentum(2.24949);
    //run->SetHRSMomentum(1.5);
    run->SetBeamEnergy(2.253248);
    run->SetTarget(6, 12);
    run->SetTargetMass(12.0107 * kU);
    //run->SetEnergyLoss(1.009711e-3 + 0.501422e-3);
    //run->SetFieldRatio(0.0);
    //run->SetFieldRatio(0.5);
    run->SetFieldRatio(1.0);
    run->SetSieve();

    ///////////////////////////////////////////////////////////////////////////
    // field
    ///////////////////////////////////////////////////////////////////////////
    G2PField* field = new G2PField();
    //field->SetEulerAngle(90*kDEG,90*kDEG,-90*kDEG);
    field->SetEulerAngle(0,0,0);
    gG2PApps->Add(field);

    ///////////////////////////////////////////////////////////////////////////
    // gun
    ///////////////////////////////////////////////////////////////////////////
    G2PFlatGun* gun = new G2PFlatGun();
    //G2PSieveGun* gun = new G2PSieveGun();
    //G2PData* gun = new G2PData("data.dat");
    //gun->SetBeamPos(-3.623e-3, 0.0); // 0.0T, 6deg
    //gun->SetBeamPos(0.5009e-3, 0.9205e-3); // 0.0T, 90deg
    //gun->SetBeamPos(-3.813e-3, 0.083e-3); // 2.5T, 90deg
    gun->SetBeamPos(-0.2357e-3, -3.855e-3); // 5.0T, 0deg
    //gun->SetReactZ(0.0, 0.0);
    //gun->SetReactZ(-14.1350e-3, -13.1191e-3); // 40mil C12
    //gun->SetReactZ(-14.1350e-3, -10.9601e-3); // 125mil C12
    gun->SetReactZ(-14.1350e-3, -14.13e-3); // production
    //gun->SetReactZ(-14.1350e-3, 14.1350e-3); // production
    //gun->SetRasterSize(14.0e-3); // production
    gun->SetRasterSize(0.0); // optics
    //gun->SetTargetTh(-100.0e-3, 10.0e-3); // 2.5T, 90deg
    //gun->SetTargetPh(-20.0e-3, 30.0e-3); // 2.5T, 90deg
    //gun->SetTargetTh(-60.0e-3, 60.0e-3); // 5.0T, 0deg
    //gun->SetTargetPh(-30.0e-3, 30.0e-3); // 5.0T, 0deg
    gun->SetTargetTh(0.0e-3, 10.0e-3); // 5.0T, 0deg
    gun->SetTargetPh(-5.0e-3, 5.0e-3); // 5.0T, 0deg   
    //gun->SetDelta(-0.04, 0.04);
    gun->SetDelta("elastic");
    gG2PApps->Add(gun);

    ///////////////////////////////////////////////////////////////////////////
    // BPM
    ///////////////////////////////////////////////////////////////////////////
    G2PBPM* bpm = new G2PBPM();
    //bpm->SetBPMRes(0, 0);
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
    //G2PPhys* model = new G2PPhys("elastic");
    //double par[1] = {2};
    //model->SetPars(par, 1);
    //gG2PApps->Add(model);

    ///////////////////////////////////////////////////////////////////////////
    // simulation
    ///////////////////////////////////////////////////////////////////////////
    G2PSim *sim = new G2PSim();
    run->SetDebugLevel(1);
    run->SetSeed(1);
    //int N = 1;
    //int N = 6550;
    int N = 500000;
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
