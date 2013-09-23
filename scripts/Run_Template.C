// -*- C++ -*-

/* script Run.C
 * This file defines a function Run().
 * This file is able to run directly in ROOT as a script.
 * It is also included in Main.cc as a function so it could be compiled.
 */

// History:
//   Jan 2013, C. Gu, First public version.
//

#include <cstdlib>
#include <cstdio>
#include <time.h>
#include <set>

#include "TROOT.h"

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

int Run() {
    ///////////////////////////////////////////////////////////////////////////
    // run parameters
    ///////////////////////////////////////////////////////////////////////////
    G2PRun* run = new G2PRun();
    run->SetHRSAngle(5.65 * kDEG);
    run->SetHRSMomentum(2.249497);
    //run->SetHRSMomentum(1.5);
    run->SetParticleID(11);
    run->SetParticleMass(0.51099892811e-3);
    run->SetParticleCharge(-1 * e);
    run->SetBeamEnergy(2.253207);
    run->SetTarget(6, 12);
    run->SetTargetMass(12.0107 * kU);
    //run->SetEnergyLoss(1.009711e-3 + 0.501422e-3);
    run->SetFieldRatio(0.5);

    ///////////////////////////////////////////////////////////////////////////
    // field
    ///////////////////////////////////////////////////////////////////////////
    G2PField* field = new G2PField();
    field->SetEulerAngle(90, 90, -90); // transverse, g2p
    //field->SetEulerAngle(90, 6, -90);  // 6 deg, gep
    gG2PApps->Add(field);

    ///////////////////////////////////////////////////////////////////////////
    // gun
    ///////////////////////////////////////////////////////////////////////////
    G2PFlatGun* gun = new G2PFlatGun();
    //G2PSieveGun* gun = new G2PSieveGun();

    //G2PData* gun = new G2PData("input_fp_tr.dat");
    //gun->SetBeamPos(4.134e-3, 1.176e-3);     // 0T,6deg
    gun->SetBeamPos(0.0, 0.0);
    //gun->SetReactZ(0.0, 0.0);
    //gun->SetReactZ(-14.1350e-3, -13.1191e-3); //40mil
    //gun->SetReactZ(-14.1350e-3, -10.9600e-3); //125mil
    gun->SetReactZ(-14.135e-3, 14.135e-3); // full size
    gun->SetRasterSize(14.0e-3);
    gun->SetTargetTh(-75.0e-3, -55.0e-3);
    gun->SetTargetPh(-10.0e-3, 10.0e-3);
    gun->SetDelta(-0.04, 0.04);
    gG2PApps->Add(gun);

    ///////////////////////////////////////////////////////////////////////////
    // BPM
    ///////////////////////////////////////////////////////////////////////////
    G2PBPM* bpm = new G2PBPM();
    //bpm->SetBPMRes(0, 0);
    bpm->SetBPMRes(0.2e-3, 0.4e-3);
    gG2PApps->Add(bpm);

    ///////////////////////////////////////////////////////////////////////////
    // HRS
    ///////////////////////////////////////////////////////////////////////////
    G2PHRS* hrs = new G2PHRS("484816");
    gG2PApps->Add(hrs);
    //G2PDBRec* rec = new G2PDBRec();
    //gG2PApps->Add(rec);

    ///////////////////////////////////////////////////////////////////////////
    // cross section
    ///////////////////////////////////////////////////////////////////////////
    G2PPhys* model = new G2PPhys("qfs");
    gG2PApps->Add(model);

    ///////////////////////////////////////////////////////////////////////////
    // simulation
    ///////////////////////////////////////////////////////////////////////////
    G2PSim *sim = new G2PSim();
    run->SetDebugLevel(1);
    run->SetSeed(1);
    //int N = 1;
    int N = 50000;
    sim->SetNEvent(N);
    //sim->SetOutFile("run_1029_optics_484816.root");
    sim->SetOutFile("run_test.root");

    clock_t start = clock();
    sim->Run();
    clock_t end = clock();
#define CLOCKS_PER_SEC 1000000
    printf("Average calculation time for one event: %8.4f ms\n", (double) (end - start)*1000.0 / (double) CLOCKS_PER_SEC / N);

    return 0;
}
