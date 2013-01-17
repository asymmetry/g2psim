#include <cstdio>
#include <cstdlib>

#include "TROOT.h"

#include "G2PSim.hh"
#include "HRSGun.hh"

int Run()
{
    G2PSim *run= new G2PSim();
    HRSGun *gun= new HRSGun("data");
    gun->SetDataFile("input_fp_tr.dat");

    run->SetGun(gun);
    run->SetNEvent(10000);
    run->SetRootName("result_test.root");
    run->SetArm("L");
    run->SetHRSMomentum(2.251);
    run->SetHRSSetting(11);
    
    run->Run();
    
	return 0;
}
