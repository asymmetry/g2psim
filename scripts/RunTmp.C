#include <cstdio>
#include <cstdlib>
#include <time.h>

int Run()
{
    //G2PPointGun* gun = new G2PPointGun();
    G2PFlatGun* gun = new G2PFlatGun();
    //G2PSieveGun* gun = new G2PSieveGun();
    //gun->SetUseFast(true);
    //G2PDataGun* gun = new G2PDataGun("input_fp_tr.dat");
    //gun->SetOpticsData(true);
    //gun->SetBeamX(4.134e-3);     // 0T,6deg
    //gun->SetBeamY(1.176e-3);     // 0T,6deg
    //gun->SetBeamX(0.0);
    //gun->SetBeamY(0.0);
    gun->SetReactZ(0.0);
    //gun->SetReactZ(-13.6271e-3); // 40mil
    //gun->SetReactZ(-12.5476e-3); // 125mil
    //gun->SetReactZRange(-14.1350e-3,-13.1191e-3); //40mil
    //gun->SetReactZRange(-14.1350e-3,-10.9600e-3); //125mil
    gun->SetBeamR(15.0e-3);
    //gun->SetReactZRange(-14.135e-3, 14.135e-3);
    //gun->SetTargetTh(2.793e-02);
    gun->SetTargetThRange(-75.0e-3,-55.0e-3);
    //gun->SetTargetPh(-3.582e-03);
    gun->SetTargetPhRange(-10.0e-3,10.0e-3);
    //gun->SetDelta(-2.128e-03);
    gun->SetDeltaRange(-0.04, 0.04);
    //gun->SetSigmaPosLab(0.0e-3);
    //gun->SetSigmaAngLab(0.0e-3);
    //gun->SetSigmaAngTr(0.0e-3);
    //gun->SetSigmaDelta(0.0e-3);
    gG2PApps->Add(gun);

    G2PHRSTrans* hrs = new G2PHRSTrans("484816");
    gG2PApps->Add(hrs);
    G2PRecUseDB* rec = new G2PRecUseDB();
    gG2PApps->Add(rec);

    G2PHallBField* field = new G2PHallBField();
    field->SetEulerAngle(90, 90, -90); // transverse, g2p
    //field->SetEulerAngle(90, 6, -90);  // 6 deg, gep
    field->SetRatio(0.5);
    gG2PApps->Add(field);
    G2PDrift* drift = new G2PDrift();
    gG2PApps->Add(drift);

    G2PBPM* bpm = new G2PBPM();
    bpm->SetBPMRes(0.2e-3, 0.2e-3);
    gG2PApps->Add(bpm);

    G2PPhys* model = new G2PPhys("qfs");
    gG2PApps->Add(model);

    G2PRun* run = new G2PRun();
    run->SetHRSAngle(5.65*3.14159265358979323846/180.0);
    //run->SetHRSMomentum(2.249497);
    run->SetHRSMomentum(1.5);
    run->SetBeamEnergy(2.253207);
    run->SetTarget(6, 12);
    run->SetTargetMass(12.0107*0.931494028);
    run->SetEnergyLoss(1.009711e-3+0.501422e-3);
    run->SetParticlePID(11);
    run->SetParticleMass(0.51099892811e-3);
    run->SetParticleCharge(-1.60217656535e-19);
    //run->SetUseEffBPM(true);

    G2PSim *sim = new G2PSim();
    int N = 50000;
    sim->SetNEvent(N);
    sim->SetSeed(1);
    sim->SetDebug(1);
    //sim->SetOutFile("run1019_flat_484816.root");
    sim->SetOutFile("run_test.root");
    sim->SetRun(run);

    clock_t start = clock();
    sim->Run();
    clock_t end = clock();
#define CLOCKS_PER_SEC 1000000
    printf("Average calcualtion time for one event: %8.4f ms\n", (double)(end-start)*1000.0/(double)CLOCKS_PER_SEC/N);
    
    return 0;
}
