#include <cstdio>
#include <cstdlib>

#include "TROOT.h"

#include "G2PSim.hh"
#include "G2PGun.hh"
#include "HRSTransport/HRSTransport.hh"
#include "G2PXSection/G2PXS.hh"

int Run()
{
    G2PSim *run = new G2PSim();
    
    //G2PGun *gun = new G2PGun("test");
    G2PGun *gun1 = new G2PGun("sieve");
    //gun1->SetTargetX(4.134e-3);
    gun1->SetTargetY(1.176e-3);
    gun1->SetTargetZRange(-14.135e-3,-10.960e-3);
    gun1->SetPositionRes(0.2e-3);
    gun1->SetAngleRes(0.6e-3);
    gun1->SetBeamEnergy(2.253207);    
    G2PGun *gun2 = new G2PGun("data");
    gun2->SetDataFile("input_fp_tr.dat");
    gun2->SetTargetZRange(-14.135e-3,-10.960e-3);
    run->AddGun(gun1);
    run->AddGun(gun2);

    HRSTransport *hrsmodel = new HRSTransport("484816R00");
    run->SetHRSModel(hrsmodel);

    G2PXS *physmodel = new G2PXS("qfs");
    run->SetPhysModel(physmodel);

    run->SetArm("L");
    run->SetHRSMomentum(2.24949710);

    //run->SetRootName("result_G4_484816R00.root"); //test
    run->SetRootName("result_G5G6_484816R00.root"); //data

    //run->SetRootName("result_G4_484816.root"); //test
    //run->SetRootName("result_G5G6_484816_noc.root"); //data

    run->Run(30000);

    return 0;
}
