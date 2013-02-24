#include <cstdio>
#include <cstdlib>
#include <time.h>

#include "TROOT.h"

#include "G2PSim.hh"
#include "G2PGun.hh"
#include "G2PTargetField.hh"
#include "HRSTransport.hh"

int Run()
{
    clock_t start = clock();
    G2PSim *run = new G2PSim();
    
    //G2PGun *gun = new G2PGun("test");
    G2PGun *gun = new G2PGun("gaus");
    //G2PGun *gun = new G2PGun("sieve");
    //G2PGun *gun = new G2PGun("data");
    //gun->SetBeamX(4.134e-3);     // 0T,6deg
    //gun->SetBeamY(1.176e-3);     // 0T,6deg
    //gun->SetReactZ(-13.6271e-3); // 40mil
    //gun->SetReactZ(-12.5476e-3); // 125mil
    //gun->SetReactZRange(-14.1350e-3,-13.1191e-3); //40mil
    //gun->SetReactZRange(-14.1350e-3,-10.9600e-3); //125mil
    gun->SetBeamX(0.0e-3);
    gun->SetBeamY(0.0e-3);
    gun->SetBeamTh(0.0e-3);
    gun->SetBeamPh(0.0e-3);
    gun->SetReactZ(0.0e-3);
    gun->SetBeamR(0.0e-3);
    //gun->SetReactZRange(-14.135e-3,14.135e-3);
    //gun->SetTargetThRange(-10e-3,10e-3);
    gun->SetTargetTh(-0.005);
    //gun->SetTargetPhRange(-10e-3,10e-3);
    //gun->SetDeltaRange(-0.01,0.01);
    gun->SetSigmaPosLab(0.0e-3);
    gun->SetSigmaAngLab(0.0e-3);
    gun->SetSigmaAngTr(0.0e-3);
    gun->SetSigmaDelta(0.0e-3);
    gun->SetBPMPosRes(0.2e-3);
    gun->SetBPMAngRes(0.8e-3);
    //gun->SetDataFile("input_fp_tr.dat");
    run->AddGun(gun);

    HRSTransport *hrsmodel = new HRSTransport("484816");
    run->SetHRSModel(hrsmodel);

    G2PTargetField *field = new G2PTargetField("hallb");
    field->SetEulerAngle(90,90,-90); // transverse, g2p
    //field->SetEulerAngle(90,6,-90);  // 6 deg, gep
    field->SetRatio(0.5);
    run->SetTargetField(field);

    G2PXS *physmodel = new G2PXS("qfs");
    run->SetPhysModel(physmodel);

    run->SetArm("L");
    run->SetHRSAngle(5.65*3.14159265358979323846/180.0);
    run->SetBeamEnergy(2.253207);    // 0T,6deg
    run->SetHRSMomentum(2.249497);   // 0T,6deg

    //run->SetRootName("result_G4_484816R00.root");
    //run->SetRootName("result_G4_484816R00_noc.root"); // no correction
    //run->SetRootName("result_G4_484816R00_noc_200.root"); // fit

    run->SetRootName("result_G3_484816.root");
    //run->SetRootName("result_G4_484816_noc.root"); // no correction
    //run->SetRootName("result_G4_484816_noc_200.root"); // fit

    int N = 10000;
    run->Run(N);
    clock_t end = clock();
#define CLOCKS_PER_SEC 1000000
    printf("Average calcualtion time for one event: %8.4f ms\n", (double)(end-start)*1000.0/(double)CLOCKS_PER_SEC/N); 
    

//     // Test Drift
//     G2PTargetField *field = new G2PTargetField("hallb");
//     field->SetEulerAngle(90,90,-90);
//     field->SetRatio(0.5);
//     G2PDrift::SetField(field);
//     double x[3] = { 0.0, 0.0, 0.0 };
//     double p[3] = { 0.0, 0.0, 2.0 };
//     clock_t start = clock();
//     for (int i=0; i<1; i++) {
//         x[0] = 0.0; x[1] = 0.0; x[2] = 0.0;
//         p[0] = 0.0; p[1] = 0.0; p[2] = 2.0;
//         G2PDrift::Drift(x, p, 0.8, 10.0, x, p);
//         //printf("%e\t%e\t%e\t%e\t%e\t%e\n", x[0], x[1], x[2], p[0], p[1], p[2]);
//     }
//     clock_t end = clock();
// #define CLOCKS_PER_SEC 1000000
//     //printf("%f\t%f\t%f\n", (double)start, (double)end, (double)(end-start)*1000.0/(double)CLOCKS_PER_SEC/20000.0);
    
    return 0;
}
