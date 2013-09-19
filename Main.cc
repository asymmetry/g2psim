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
#include <getopt.h>

#include "TROOT.h"

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
