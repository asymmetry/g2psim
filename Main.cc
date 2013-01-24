#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include "G2PSim.hh"
#include "G2PGun.hh"
#include "HRSTransport.hh"

void usage(int argc, char** argv);

int main(int argc, char** argv)
{
    int c;

    int nEvent = 50000;    // default 50000
    int iSetting = 11;     // default 484816 septa with shim
    int iArm = 0;          // default left arm

    while (1) {
        static struct option long_options[] = {
            {"help",        no_argument,       0, 'h'},
            {"setting",     required_argument, 0, 's'},
            {"arm",         required_argument, 0, 'a'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long(argc, argv, "a:hs:", long_options, &option_index);

        if (c==-1) break;

        switch (c) {
        case 'a':
            iArm = atoi(optarg);
            break;
        case 'h':
            usage(argc, argv);
            break;
        case 's':
            iSetting = atoi(optarg);
            break;
        case '?':
        default:
            usage(argc, argv);
        }
    }

    if (optind<argc) {
        nEvent = atoi(argv[optind++]);
    }
    else {
        usage(argc, argv);
        exit(-1);
    }

    G2PSim *run= new G2PSim();
    G2PGun *gun= new G2PGun("data");
    gun->SetDataFile("input_fp_tr.dat");
    HRSTransport *model = new HRSTransport(iSetting);

    run->SetGun(gun);
    run->SetHRSModel(model);
    run->SetNEvent(nEvent);
    run->SetRootName("result_test.root");
    if (iArm==0) run->SetArm("L");
    else if (iArm==1) run->SetArm("R");
    run->SetHRSMomentum(2.251);
    
    run->Run();
    
	return 0;
}

void usage(int argc, char** argv)
{
    printf("usage: %s [options] NEvent\n", argv[0]);
    printf("  -a, --arm=0           Set arm: 0 means left arm, 1 means right arm\n");
    printf("  -b, --bpm=0.0         Set bpm resolution in mm\n");
    printf("  -d, --direction=0     Set sim direction: 0 means forward, 1 means backward, 2 means both\n");
    printf("  -h, --help            This small usage guide\n");
    printf("  -i, --inputsource=0   Set data source: 1,2,3 means use delta,flat,gaussian distribution for input position and angle, 0 means use real data in file input_tg.dat and input_fp.dat\n");
    printf("  -m, --momentum=2.251  Set HRS momentum in GeV\n");
    printf("  -s, --setting=11      Set experiment setting: 10 means normal 484816 septa, 11 means 484816 septa with shim, 12 means 403216 septa with shim, 13 means 400016 septa with shim\n");
}
