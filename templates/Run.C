// -*- C++ -*-

/* script Run.C
 * Can directly run in ROOT as a script.
 * To do this, some include paths should be set by rootlogon.C.
 */

//

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

int Run()
{
    ///////////////////////////////////////////////////////////////////////////
    // run parameters
    ///////////////////////////////////////////////////////////////////////////
    G2PRun *run = new G2PRun();
    run->SetConfigFile("sim.cfg");

    run->SetBeamEnergy(2.2535);
    //run->SetBeamEnergy(1.7105);
    //run->SetBeamEnergy(1.1570);
    //run->SetBeamEnergy(1.1510); // 1.158GeV Optics

    run->SetHRSAngle(5.767 * kDEG);
    //run->SetHRSAngle(-5.781 * kDEG);

    //LHRS
    run->SetHRSMomentum(2.2495);
    //run->SetHRSMomentum(1.7044);
    //run->SetHRSMomentum(1.1573);
    //RHRS
    //run->SetHRSMomentum(2.2527);
    //run->SetHRSMomentum(1.7066);
    //run->SetHRSMomentum(1.1583);

    //run->SetTargetType("production");
    run->SetTargetType("optics"); // 40 mil carbon, no LHe
    //run->SetTargetType("optics_C40He"); // 40 mil carbon, with LHe
    //run->SetTargetType("optics_C125"); // 125 mil carbon, no LHe
    //run->SetTargetType("optics_C125He"); // 125 mil carbon, with LHe
    //run->SetTargetType("dummy");
    //run->SetTargetType("dummy_nocap"); // without end cap
    //run->SetTargetType("empty");

    run->SetTarget(6, 12);
    run->SetTargetMass(12.0107 * kU);

    run->SetFieldType("trans");
    //run->SetFieldType("longitudinal");
    //run->SetFieldType("gep");

    //run->SetFieldRatio(0.0);
    run->SetFieldRatio(0.5);
    //run->SetFieldRatio(1.0);

    ///////////////////////////////////////////////////////////////////////////
    // field
    ///////////////////////////////////////////////////////////////////////////
    G2PField *field = new G2PField();
    field->SetEulerAngle(90 * kDEG, 90 * kDEG, -90 * kDEG);
    //field->SetEulerAngle(90 * kDEG, 5.6 * kDEG, -90 * kDEG);
    //field->SetEulerAngle(0, 0, 0);
    gG2PApps->Add(field);

    ///////////////////////////////////////////////////////////////////////////
    // event generator
    ///////////////////////////////////////////////////////////////////////////
    G2PEvGen *gun = new G2PEvGen();
    //G2PData* gun = new G2PData("data.dat");

    //gun->SetBeamPos(0, 0, 0);
    //LHRS
    //gun->SetBeamPos( -0.8e-3,  2.3e-3, -13.6271e-3); // 2.253GeV, 0.0T, 6deg
    gun->SetBeamPos(-3.7e-3, -0.8e-3, -13.6271e-3);  // 2.253GeV, 2.5T, 90deg
    //gun->SetBeamPos( -6.2e-3,  7.3e-3, -13.6271e-3); // 2.253GeV, 2.5T, 90deg
    //gun->SetBeamPos(  1.9e-3, -0.8e-3, -13.6271e-3); // 1.706GeV, 2.5T, 90deg
    //gun->SetBeamPos(  2.6e-3, -7.4e-3, -13.6271e-3); // 1.158GeV, 2.5T, 90deg
    //gun->SetBeamPos(  0.6e-3, -2.3e-3, -12.5476e-3); // 2.253GeV, 5.0T, 0deg
    //RHRS
    //gun->SetBeamPos(  0.1e-3,  2.9e-3, -13.6271e-3); // 2.253GeV, 0.0T, 6deg
    //gun->SetBeamPos( -3.5e-3, -0.7e-3, -13.6271e-3); // 2.253GeV, 2.5T, 90deg
    //gun->SetBeamPos(-12.4e-3, -4.4e-3, -13.6271e-3); // 2.253GeV, 2.5T, 90deg
    //gun->SetBeamPos(  1.9e-3, -1.3e-3, -13.6271e-3); // 1.706GeV, 2.5T, 90deg
    //gun->SetBeamPos(  0.8e-3, -6.1e-3, -13.6271e-3); // 1.158GeV, 2.5T, 90deg

    //gun->SetTiltAngle(0, 0);

    //gun->SetReactZ( -1.0000e-8,   1.0000e-8);
    gun->SetReactZ(-14.1350e-3, -13.1191e-3); // 40mil C12
    //gun->SetReactZ(-14.1350e-3, -10.9601e-3); // 125mil C12
    //gun->SetReactZ(-14.1350e-3,  14.1350e-3); // production
    //gun->SetReactZ(-14.1350e-3, -14.1350e-3); // eloss

    //gun->SetRasterSize(15.0e-3); // production
    gun->SetRasterSize(0.5e-3); // optics

    //gun->SetTargetTh( -60.0e-3,   60.0e-3);
    gun->SetTargetTh(-100.0e-3,   20.0e-3); // 2.253GeV, 2.5T, 90deg
    //gun->SetTargetTh(-120.0e-3,    0.0e-3); // 1.706GeV, 2.5T, 90deg
    //gun->SetTargetTh(-150.0e-3,  -30.0e-3); // 1.158GeV, 2.5T, 90deg
    //gun->SetTargetTh( -60.0e-3,   60.0e-3); // 2.253GeV, 5.0T, 0deg
    //gun->SetTargetTh(-160.0e-3,  -50.0e-3); // 2.253GeV, 5.0T, 90deg
    gun->SetTargetPh(-30.0e-3,   30.0e-3);

    //gun->SetDelta(-0.04, 0.04); // normal
    gun->SetDelta("elastic"); // optics
    //gun->SetDelta(-0.8, 1.0); // effective bpm

    gG2PApps->Add(gun);

    ///////////////////////////////////////////////////////////////////////////
    // BPM
    ///////////////////////////////////////////////////////////////////////////
    //G2PBPM* bpm = new G2PBPM();
    //bpm->SetBPMRes(0, 0);
    //bpm->SetBPMRes(0.2e-3, 0.4e-3);
    //gG2PApps->Add(bpm);

    ///////////////////////////////////////////////////////////////////////////
    // major processes
    ///////////////////////////////////////////////////////////////////////////
    G2PFwdHRS *fwd = new G2PFwdHRS("484816");
    fwd->SetSieve("in");
    //G2PFwdTarget* fwd = new G2PFwdTarget();

    gG2PApps->Add(fwd);

    G2PBwdHRS *bwd = new G2PBwdHRS("484816");
    //G2PBwdDB* bwd = new G2PBwdDB();

    //LHRS
    //double parsx[3] = { 0.00079, -0.75038,  178.902}; // 5.0T, 0deg
    //double parsy[3] = { 0.15693, -0.39864, -176.049}; // 5.0T, 0deg
    double parsx[3] = { -0.00061,  3.13661, -4.87024}; // 2.5T, 90deg
    double parsy[3] = { -0.08772,  0.21862, -17.5404}; // 2.5T, 90deg
    //double parsx[3] = {-0.03853,  6.34323, -4.97020}; // 5.0T, 90deg
    //double parsy[3] = {-0.35369,  0.86343, -29.1603}; // 5.0T, 90deg
    //double parsx[3] = {-0.00179, -0.03300,  179.472}; // 5.0T, 6deg
    //double parsy[3] = { 0.00380, -0.03504, -176.711}; // 5.0T, 6deg
    //RHRS
    //double parsx[3] = {-0.00301,  0.76169,  177.993}; // 5.0T, 0deg
    //double parsy[3] = {-0.15435,  0.34665, -179.100}; // 5.0T, 0deg
    //double parsx[3] = { 0.00317,  3.13296,  4.63844}; // 2.5T, 90deg
    //double parsy[3] = { 0.08567, -0.21493,  16.5077}; // 2.5T, 90deg
    //double parsx[3] = {-0.05071,  6.36930,  5.63490}; // 5.0T, 90deg
    //double parsy[3] = { 0.34858, -0.85407,  30.4795}; // 5.0T, 90deg
    //double parsx[3] = { 0.02936,  1.41205,  175.656}; // 5.0T, 6deg
    //double parsy[3] = {-0.33310,  0.75685, -180.501}; // 5.0T, 6deg
    bwd->SetParsX(parsx);
    bwd->SetParsY(parsy);

    //bwd->SetRecZ(0.0); // production
    bwd->SetRecZ(-13.6271e-3); // 40mil C12
    //bwd->SetRecZ(-12.5476e-3); // 125mil C12

    gG2PApps->Add(bwd);

    ///////////////////////////////////////////////////////////////////////////
    // cross section
    ///////////////////////////////////////////////////////////////////////////
    G2PPhys *model = new G2PPhys("pbosted");
    //double par[2] = {0.0080,0.0080};
    //model->SetPars(par, 2);
    gG2PApps->Add(model);

    ///////////////////////////////////////////////////////////////////////////
    // simulation
    ///////////////////////////////////////////////////////////////////////////
    G2PSim *sim = new G2PSim();
    run->SetDebugLevel(1);
    run->SetSeed(1);
    //int N = 1;
    int N = 200000;
    sim->SetNEvent(N);
    sim->SetOutFile("test_sim.root");

    clock_t start = clock();
    sim->Run();
    clock_t end = clock();
#define CLOCKS_PER_SEC 1000000
    printf("Average calculation time for one event: %8.4f ms\n", (double)(end - start) * 1000.0 / (double) CLOCKS_PER_SEC / N);

    return 0;
}
