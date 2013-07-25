// -*- C++ -*-

/*
 * This file defines main function of this simulation.
 */

// History:
//   Jan 2013, C. Gu, First public version.
//   Mar 2013, C. Gu, Use Run.C script to do simulation
//

#include <cstdio>
#include <cstdlib>
#include <time.h>
#include <getopt.h>

#include "TROOT.h"

#include "G2PDBRec.hh"
#include "G2PGlobals.hh"
#include "G2PRun.hh"
#include "G2PRunBase.hh"
#include "G2PSim.hh"

#include "G2PField.hh"
#include "G2PHallBField.hh"
#include "G2PMapField.hh"
#include "G2PUniField.hh"

#include "G2PGun.hh"
#include "G2PDataGun.hh"
#include "G2PFlatGun.hh"
#include "G2PSieveGun.hh"

#include "G2PHRSTrans.hh"
#include "HRSTrans/G2PTrans400016/G2PTrans400016.hh"
#include "HRSTrans/G2PTrans484816/G2PTrans484816.hh"
#include "HRSTrans/HRSTransSTD/HRSTransSTD.hh"

#include "G2PPhys.hh"
#include "G2PPhys/G2PPhysEl/G2PPhysEl.hh"
#include "G2PPhys/G2PPhysPB/G2PPhysPB.hh"
#include "G2PPhys/G2PPhysQFS/G2PPhysQFS.hh"
#include "G2PPhys/G2PPhysWISER/G2PPhysWISER.hh"

void usage(int argc, char** argv);

int main(int argc, char** argv) {
    int c;

    int nEvent = 50000; // default 50000

    while (1) {
        static struct option long_options[] = {
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long(argc, argv, "h", long_options, &option_index);

        if (c == -1) break;

        switch (c) {
        case 'h':
            usage(argc, argv);
            break;
        case '?':
        default:
            usage(argc, argv);
        }
    }

    if (optind < argc) {
        nEvent = atoi(argv[optind++]);
    }
    else {
        usage(argc, argv);
        exit(-1);
    }

    return 0;
}

void usage(int argc, char** argv) {
    printf("usage: %s [options] NEvent\n", argv[0]);
    printf("  -h, --help            This small usage guide\n");
}
