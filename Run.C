#include <cstdio>
#include <cstdlib>

#include "TROOT.h"

#include "G2PSim.hh"
#include "HRSGun.hh"

int Run()
{
    G2PSim *run= new G2PSim();
    HRSGun *gun= new HRSGun("sieve");
    gun->SetDataFile("input_fp_tr.dat");
    gun->SetTargetX(0.0);
    gun->SetTargetY(0.0);
    gun->SetTargetZRange(-14.135e-3,-10.960e-3);
    gun->SetPositionRes(0.5e-3);

    run->SetGun(gun);
    run->SetRootName("result_test.root");
    run->SetArm("L");
    run->SetHRSMomentum(2.251);
    run->SetHRSSetting(11);
    
    run->Run(10000);
    
    
	return 0;
}
