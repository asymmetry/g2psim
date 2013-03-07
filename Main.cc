#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <getopt.h>

#include "G2PSim.hh"

void usage(int argc, char** argv);

int main(int argc, char** argv)
{
    int c;

    int nEvent = 50000;    // default 50000
    char pSnake[300] = "484816";     // default 484816 septa with shim
    char pField[300] = "hallb";
    // char pExperiment[300] = "g2p";
    char pPhysics[300] = "qfs";
    char pGun[300] = "data";
    char pArm[300] = "L";          // default left arm

    while (1) {
        static struct option long_options[] = {
            {"help",        no_argument,       0, 'h'},
            {"snake",       required_argument, 0, 's'},
            {"arm",         required_argument, 0, 'a'},
            {"field",       required_argument, 0, 'f'},
            {"physics",     required_argument, 0, 'p'},
            {"gun",         required_argument, 0, 'g'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long(argc, argv, "a:f:g:hp:s:", long_options, &option_index);

        if (c==-1) break;

        switch (c) {
        case 'a':
            strcpy(pArm, optarg);
            break;
        case 'h':
            usage(argc, argv);
            break;
        case 's':
            strcpy(pSnake, optarg);
            break;
        case 'f':
            strcpy(pField, optarg);
            break;
        case 'p':
            strcpy(pPhysics, optarg);
            break;
        case 'g':
            strcpy(pGun, optarg);
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

    // clock_t start = clock();

    // G2PSim *run= new G2PSim();
    // G2PGun *gun= new G2PGun(pGun);
    // gun->SetDataFile("input_fp_tr.dat");
    // //G2PGun *gun = new G2PGun("sieve");
    // //G2PGun *gun = new G2PGun("data");
    // //gun->SetBeamX(4.134e-3);     // 0T,6deg
    // //gun->SetBeamY(1.176e-3);     // 0T,6deg
    // //gun->SetReactZ(-13.6271e-3); // 40mil
    // //gun->SetReactZ(-12.5476e-3); // 125mil
    // //gun->SetReactZRange(-14.1350e-3,-13.1191e-3); //40mil
    // //gun->SetReactZRange(-14.1350e-3,-10.9600e-3); //125mil
    // gun->SetBeamX(0.0e-3);
    // gun->SetBeamY(0.0e-3);
    // gun->SetReactZ(0.0e-3);
    // gun->SetSigmaPosLab(0.0e-3);
    // gun->SetSigmaAngLab(0.2e-3);
    // gun->SetSigmaAngTr(0.2e-3);
    // gun->SetSigmaDelta(0.0e-3);
    // //gun->SetDataFile("input_fp_tr.dat");
    // run->AddGun(gun);

    // HRSTransport *model = new HRSTransport(pSnake);
    // run->SetHRSModel(model);

    // G2PTargetField *field = new G2PTargetField(pField);
    // if (strcmp(pExperiment, "g2p")) {
    //     field->SetEulerAngle(90,90,-90); // transverse, g2p
    // }
    // else {
    //     field->SetEulerAngle(90,6,-90);  // 6 deg, gep
    // }
    // field->SetRatio(0.5);
    // run->SetTargetField(field);

    // G2PXS *physmodel = new G2PXS("qfs");
    // run->SetPhysModel(physmodel);

    // run->SetArm(pArm);
    // run->SetBeamEnergy(2.253207);    // 0T,6deg
    // run->SetHRSMomentum(2.249497);   // 0T,6deg

    // run->SetNEvent(nEvent);
    // run->SetRootName("result_test.root");
    
    // run->Run();

    // clock_t end = clock();

    // printf("Average calcualtion time for one event: %8.4f ms\n", (double)(end-start)*1000.0/(double)CLOCKS_PER_SEC/nEvent);

    // // Test Drift
    // G2PTargetField *field = new G2PTargetField("hallb");
    // field->SetEulerAngle(90,90,-90);
    // field->SetRatio(0.5);
    // field->Init();
    // G2PDrift::SetField(field);
    // double x[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
    // double p[3] = { 0.0, 0.0, 2.0 };
    // for (int i=0; i<1; i++) {
    //     x[0] = 0.01; x[1] = 0.01; x[2] = 0.0; x[3] = 0.0; x[4] = 0.0;
    //     p[0] = 0.2; p[1] = 0.0; p[2] = 2.0;
    //     G2PDrift::Drift(x, 2.251, 0.1, 0.0, 0.8, 10.0, x);
    //     printf("%e\t%e\t%e\t%e\t%e\n", x[0], x[1], x[2], x[3], x[4]);
    //     G2PDrift::Drift(x, 2.251, 0.1, 0.8, 0.0, 10.0, x);
    //     printf("%e\t%e\t%e\t%e\t%e\n", x[0], x[1], x[2], x[3], x[4]);
    // } 
    
    
	return 0;
}

void usage(int argc, char** argv)
{
    printf("usage: %s [options] NEvent\n", argv[0]);
    printf("  -a, --arm=L           Set arm: L or R\n");
    printf("  -f, --field=hallb     Set field: uniform or hallb \n");
    printf("  -g, --gun=data        Set gun: delta, gaus, flat, test, sieve or data\n");
    printf("  -h, --help            This small usage guide\n");
    printf("  -p, --physics=qfs     Set physics model: qfs\n");
    printf("  -s, --snake=484816    Set hrs model: STD, 484816, 483216, 400016, 484816R00, GDHSTD or GDHLargeX0\n");
}
