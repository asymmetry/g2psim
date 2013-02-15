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
    //G2PGun *gun = new G2PGun("sieve");
    G2PGun *gun = new G2PGun("data");
    gun->SetTargetX(4.134e-3);
    gun->SetTargetY(1.176e-3);
    gun->SetTargetZ(-12.5475e-3);
    //gun->SetTargetZRange(-14.135e-3,-10.960e-3);
    gun->SetPositionRes(1.0e-3);
    gun->SetAngleRes(1.0e-3);
    gun->SetBeamEnergy(2.253207);
    gun->SetDataFile("input_fp_tr.dat");
    run->AddGun(gun);

    HRSTransport *hrsmodel = new HRSTransport("484816R00");
    run->SetHRSModel(hrsmodel);

    G2PXS *physmodel = new G2PXS("qfs");
    run->SetPhysModel(physmodel);

    run->SetArm("L");
    run->SetHRSMomentum(2.24949710);

    run->SetRootName("result_G6_484816R00.root");
    //run->SetRootName("result_G6_484816R00_noc.root"); // no correction
    //run->SetRootName("result_G6_484816R00_noc_200.root"); // fit

    //run->SetRootName("result_G6_484816.root");
    //run->SetRootName("result_G6_484816_noc.root"); // no correction
    //run->SetRootName("result_G6_484816_noc_200.root"); // fit

    run->Run(30000);

    return 0;
}
