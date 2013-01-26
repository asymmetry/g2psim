#include <cstdio>
#include <cstdlib>

#include "TROOT.h"

#include "G2PSim.hh"
#include "G2PGun.hh"
#include "HRSTransport/HRSTransport.hh"

int Run()
{
    G2PSim *run = new G2PSim();
    //G2PGun *gun = new G2PGun("test");
    //G2PGun *gun = new G2PGun("data");
    G2PGun *gun = new G2PGun("sieve");
    HRSTransport *model = new HRSTransport(11);
    gun->SetDataFile("input_fp_tr.dat");
    gun->SetTargetX(0.0e-3);
    gun->SetTargetY(1.2e-3);
    gun->SetTargetZRange(-14.135e-3,-10.960e-3);
    gun->SetPositionRes(0.2e-3);
    gun->SetAngleRes(0.6e-3);
    gun->SetBeamEnergy(2.253207);
    
    run->SetGun(gun);
    run->SetHRSModel(model);
    run->SetRootName("result_G5_S11.root");
    //run->SetRootName("result_G6_S11.root"); // data
    //run->SetRootName("result_G4_S11.root");
    run->SetArm("L");
    run->SetHRSMomentum(2.24949710);
    
    run->Run(30000);

    return 0;
}
