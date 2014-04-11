// -*- C++ -*-

/* Main.cc
 * Main function of this simulation.
 * It parses command line parameters use getopt (GNU C Library).
 */

// History:
//   Jan 2013, C. Gu, First public version.
//   Mar 2013, C. Gu, Use Run.C script to do simulation
//   Sep 2013, C. Gu, Since G2PRun will take care of the configuration files, make it simpler here.
//

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <getopt.h>
#include <map>
#include <time.h>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"

#include "G2PAppList.hh"
#include "G2PBPM.hh"
#include "G2PDBBwd.hh"
#include "G2PFPData.hh"
#include "G2PField.hh"
#include "G2PFlatGun.hh"
#include "G2PGlobals.hh"
#include "G2PGun.hh"
#include "G2PHRSBwd.hh"
#include "G2PHRSFwd.hh"
#include "G2PPhys.hh"
#include "G2PRun.hh"
#include "G2PSieveGun.hh"
#include "G2PSim.hh"

// #define CLOCKS_PER_SEC 1000000

using namespace std;

void usage(int argc, char** argv);

int main(int argc, char** argv)
{
    int c;

    int fN = 50000; // default 50000
    int fDebug = 1;
    char fConfig[300] = "sim.cfg";
    char fOutput[300] = "run_test.root";
    int fBPM = 1;
    int fDB = 0;
    int fField = 1;
    char fGun[300] = "flat";
    map<string, int> fGunType;
    fGunType["flat"] = 1;
    fGunType["sieve"] = 2;
    char fHRS[300] = "484816";
    char fPhys[300] = "pbosted";

    while (1) {
        static struct option long_options[] = {
            {"nobpm", no_argument, &fBPM, 0},
            {"nofield", no_argument, &fField, 0},
            {"usedb", no_argument, &fDB, 1},
            {"config", required_argument, 0, 'c'},
            {"debug", required_argument, 0, 'd'},
            {"gun", required_argument, 0, 'g'},
            {"help", no_argument, 0, 'h'},
            {"hrsmodel", required_argument, 0, 'm'},
            {"output", required_argument, 0, 'o'},
            {"xsmodel", required_argument, 0, 'x'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long(argc, argv, "c:d:g:hm:o:x:", long_options, &option_index);

        if (c == -1) break;

        switch (c) {
        case 0:
            break;
        case 'c':
            strcpy(fConfig, optarg);
            break;
        case 'd':
            fDebug = atoi(optarg);
            break;
        case 'g':
            strcpy(fGun, optarg);
            break;
        case 'h':
            usage(argc, argv);
            exit(0);
        case 'm':
            strcpy(fHRS, optarg);
            break;
        case 'o':
            strcpy(fOutput, optarg);
            break;
        case 'x':
            strcpy(fPhys, optarg);
            break;
        case '?':
        default:
            usage(argc, argv);
            exit(-1);
        }
    }

    if (optind < argc) {
        fN = atoi(argv[optind++]);
    } else {
        usage(argc, argv);
        exit(-1);
    }

    G2PRun* run = new G2PRun();
    run->SetConfigFile(fConfig);
    if (fDebug != 1) run->SetDebugLevel(fDebug);

    G2PSim *sim = new G2PSim();
    sim->SetNEvent(fN);
    sim->SetOutFile(fOutput);

    if (fBPM) {
        G2PBPM* bpm = new G2PBPM();
        gG2PApps->Add(bpm);
    }

    G2PField* field = new G2PField();
    gG2PApps->Add(field);
    if (!fField) {
        run->SetFieldRatio(0.0);
    }

    if (fGunType.count(fGun) > 0) {
        switch (fGunType[fGun]) {
        case 1:
        {
            G2PGun* gun = new G2PFlatGun();
            gG2PApps->Add(gun);
            break;
        }
        case 2:
        {
            G2PGun* gun = new G2PSieveGun();
            gG2PApps->Add(gun);
            break;
        }
        }
    } else {
        G2PFPData* gun = new G2PFPData(fGun);
        gG2PApps->Add(gun);
    }

    if (strcmp(fHRS, "off") != 0) {
        G2PHRSFwd* fwd = new G2PHRSFwd(fHRS);
        gG2PApps->Add(fwd);
        if (fDB) {
            G2PDBBwd* db = new G2PDBBwd("db_L.vdc.dat");
            gG2PApps->Add(db);
        } else {
            G2PHRSBwd* bwd = new G2PHRSBwd(fHRS);
            gG2PApps->Add(bwd);
        }
    }

    if (strcmp(fPhys, "off") != 0) {
        G2PPhys* phys = new G2PPhys(fPhys);
        gG2PApps->Add(phys);
    }

    clock_t start = clock();
    sim->Run();
    clock_t end = clock();
    printf("Average calculation time for one event: %8.4f ms\n", (double) (end - start)*1000.0 / (double) CLOCKS_PER_SEC / fN);

    return 0;
}

void usage(int argc, char** argv)
{
    printf("usage: %s [options] NEvent\n", argv[0]);
    printf("  -c, --config=sim.cfg           Set configuration file\n");
    printf("  -d, --debug=1                  Set debug level\n");
    printf("  -g, --gun=flat                 Set gun type\n");
    printf("  -h, --help                     Print this small usage guide\n");
    printf("  -m, --hrsmodel=484816          Set HRS model\n");
    printf("  -o, --output=run_test.root     Set configuration file\n");
    printf("  -x, --xsmodel=pbosted          Set cross section model\n");
    printf("  --nobpm                        Turn off BPM\n");
    printf("  --nofield                      Turn off target field\n");
    printf("  --usedb                        Turn on database reconstruction\n");
}
