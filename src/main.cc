#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

void TestSNAKE(int iNEvent, int iArm, int iSetting, int iSource, int iDirection, double fBPM);
void usage(int argc, char** argv);

int main(int argc, char** argv)
{
    int c;

    int iNEvent = 50000;   // default 50000
    int iSetting = 11;  // default 484816 septa with shim
    double fBPM = 0.0;     // default no smear
    int iSource = 0;       // default delta distribution
    int iDirection = 0;    // default forward direction
    int iArm = 0;          // default left arm

    while (1) {
        static struct option long_options[] = {
            {"help",        no_argument,       0, 'h'},           
            {"inputsource", required_argument, 0, 'i'},
            {"direction",   required_argument, 0, 'd'},
            {"setting",     required_argument, 0, 's'},
            {"bpm",         required_argument, 0, 'b'},
            {"arm",         required_argument, 0, 'a'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long(argc, argv, "a:b:d:hs:", long_options, &option_index);

        if (c==-1) break;

        switch (c) {
        case 'a':
            iArm = atoi(optarg);
            break;
        case 'b':
            fBPM = atof(optarg);
            break;
        case 'd':
            iDirection = atoi(optarg);
            break;
        case 'h':
            usage(argc, argv);
            break;
        case 'i':
            iSource = atoi(optarg);
            break;
        case 's':
            iSetting = atoi(optarg);
            break;
        case '?':
            break;
        default:
            usage(argc, argv);
        }
    }

    if (optind<argc) {
        iNEvent = atoi(argv[optind++]);
    }
    else {
        usage(argc, argv);
        exit(-1);
    }

    TestSNAKE(iNEvent, iArm, iSetting, iSource, iDirection, fBPM);

	return 0;
}

void usage(int argc, char** argv)
{
    printf("usage: %s [options] NEvent\n", argv[0]);
    printf("  -a, --arm=0          Set arm: 0 means left arm, 1 means right arm\n");
    printf("  -b, --bpm=0.0        Set bpm resolution\n");
    printf("  -d, --direction=0    Set sim direction: 0 means forward, 1 means backward, 2 means both\n");
    printf("  -h, --help           This small usage guide\n");
    printf("  -i, --inputsource=0  Set data source: 0,1,2 means use delta,flat,gaussian distribution for input position and angle, 3 means use real data in file input_tg.dat and input_fp.dat\n");
    printf("  -s, --setting=11     Set experiment setting: 10 means normal 484816 septa, 11 means 484816 septa with shim, 12 means 403216 septa with shim, 13 means 400016 septa with shim\n");
}
