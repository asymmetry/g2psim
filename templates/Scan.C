#include <cstdio>
#include <cstdlib>

#include "TROOT.h"
#include "TError.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

void Scan(const char *filename = "test_sim.root", const char *defilename = "scan.def")
{
    // Will add TChain later
    TFile *f = new TFile(filename);
    TTree *T = (TTree *)f->Get("T");

    double isgood;
    T->SetBranchAddress("isgood", &isgood);

    FILE *defile = fopen(defilename, "r");

    double data[100];
    bool isint[100];
    bool usetan[100];
    char t1[300], t2[300], t3[300], temp[300];
    int NDATA = 0;

    while (!feof(defile)) {
        temp[0] = 0;
        fgets(temp, 300, defile);

        if (temp[0] == 0 || temp[0] == '#')
            continue;

        sscanf(temp, "%s%s%s", t1, t2, t3);

        if (t1[0] != 0) {
            T->SetBranchAddress(t1, &data[NDATA]);

            TString ts2 = t2;
            TString ts3 = t3;

            if (ts2 == "int")
                isint[NDATA] = true;
            else
                isint[NDATA] = false;

            if (ts3 == "tan")
                usetan[NDATA] = true;
            else
                usetan[NDATA] = false;

            NDATA++;
        }
    }

    fclose(defile);

    FILE *outfile = fopen("scan.dat", "w");

    int N = T->GetEntries();

    for (int i = 0; i < N; i++) {
        T->GetEntry(i);

        if (isgood > 0.5) {
            for (int j = 0; j < NDATA; j++) {
                if (isint[j])
                    fprintf(outfile, "%-5d  ", (int)data[j]);
                else if (usetan[j])
                    fprintf(outfile, "%13.6le  ", tan(data[j]));
                else
                    fprintf(outfile, "%13.6le  ", data[j]);
            }

            fprintf(outfile, "\n");
        }
    }

    fclose(outfile);
}
