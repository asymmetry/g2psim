// -*- C++ -*-

/* script Run.C
 * Can directly run in ROOT as a script.
 * To do this, some include paths should be set by rootlogon.C.
 */

// run list:
// 2000 2001 2002 2003 // LHRS
// 2010 2011 2012 2013 // RHRS

#include <cstdlib>
#include <cstdio>
#include <time.h>

#include "TROOT.h"

#include "G2PBPM.hh"
#include "G2PBwdDB.hh"
#include "G2PBwdHRS.hh"
#include "G2PData.hh"
#include "G2PEvGen.hh"
#include "G2PField.hh"
#include "G2PFwdHRS.hh"
#include "G2PFwdTarget.hh"
#include "G2PGlobals.hh"
#include "G2PPhys.hh"
#include "G2PRun.hh"
#include "G2PSim.hh"

static const double kDEG = 3.14159265358979323846 / 180.0;
static const double e = 1.60217656535e-19;
static const double kU = 0.931494028;

int Run_FitBPM()
{
    ///////////////////////////////////////////////////////////////////////////
    // run parameters
    ///////////////////////////////////////////////////////////////////////////
    G2PRun *run = new G2PRun();
    run->SetConfigFile("sim.cfg");

    run->SetBeamEnergy(3.5);
    run->SetHRSAngle(5.767 * kDEG); // LHRS
    //run->SetHRSAngle(-5.781 * kDEG); // RHRS
    run->SetHRSMomentum(2.0);

    run->SetTargetType("production");
    //run->SetTargetType("optics"); // 40 mil carbon, no LHe
    //run->SetTargetType("optics_C40He"); // 40 mil carbon, with LHe
    //run->SetTargetType("optics_C125"); // 125 mil carbon, no LHe
    //run->SetTargetType("optics_C125He"); // 125 mil carbon, with LHe
    //run->SetTargetType("dummy");
    //run->SetTargetType("dummy_nocap"); // without end cap
    //run->SetTargetType("empty");

    run->SetTarget(6, 12);
    run->SetTargetMass(12.0107 * kU);

    //run->SetFieldType("trans");
    run->SetFieldType("longitudinal");
    //run->SetFieldType("gep");

    //run->SetFieldRatio(0.0);
    //run->SetFieldRatio(0.5);
    run->SetFieldRatio(1.0);

    ///////////////////////////////////////////////////////////////////////////
    // field
    ///////////////////////////////////////////////////////////////////////////
    G2PField *field = new G2PField();
    //field->SetEulerAngle(90 * kDEG, 90 * kDEG, -90 * kDEG);
    //field->SetEulerAngle(90 * kDEG, 5.6 * kDEG, -90 * kDEG);
    field->SetEulerAngle(0, 0, 0);
    gG2PApps->Add(field);

    ///////////////////////////////////////////////////////////////////////////
    // gun
    ///////////////////////////////////////////////////////////////////////////
    G2PEvGen *gun = new G2PEvGen();
    gun->SetBeamPos(0, 0, 0);
    gun->SetReactZ(-14.1350e-3, 14.1350e-3); // production
    gun->SetRasterSize(15.0e-3); // production
    gun->SetTargetTh(-120.0e-3,  120.0e-3); // 5.0T, 0deg // 5.0T, 6deg
    //gun->SetTargetTh(-200.0e-3,   40.0e-3); // 2.5T, 90deg
    //gun->SetTargetTh(-440.0e-3,  -40.0e-3); // 5.0T, 90deg
    gun->SetTargetPh(-60.0e-3, 60.0e-3);
    gun->SetDelta(-0.8, 0.5); // effective bpm
    gG2PApps->Add(gun);

    ///////////////////////////////////////////////////////////////////////////
    // BPM
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    // special commands
    ///////////////////////////////////////////////////////////////////////////
    G2PFwdTarget *fwd = new G2PFwdTarget();
    gG2PApps->Add(fwd);

    ///////////////////////////////////////////////////////////////////////////
    // cross section
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    // simulation
    ///////////////////////////////////////////////////////////////////////////
    G2PSim *sim = new G2PSim();
    run->SetDebugLevel(1);
    run->SetSeed(1);
    int N = 100000;
    sim->SetNEvent(N);
    sim->SetOutFile("test_sim.root");

    clock_t start = clock();
    sim->Run();
    clock_t end = clock();
#define CLOCKS_PER_SEC 1000000
    printf("Average calculation time for one event: %8.4f ms\n", (double)(end - start) * 1000.0 / (double) CLOCKS_PER_SEC / N);

    return 0;
}
