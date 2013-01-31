#include <string>
#include <map>

#include "TROOT.h"
#include "TObject.h"

#include "G2PPModsQFS/G2PPModsQFS.hh"

#include "G2PXS.hh"

//#define PHYS_DEBUG 1

using namespace std;

ClassImp(G2PXS);

G2PXS::G2PXS()
    :iSetting(1), iZ(1), iA(1), fTb(0.0), fTa(0.0), fEPS(10.0),
     fEPSD(-10.0), fFP(220.0), pfModelSelector(NULL)
{
    // Nothing to do
}

G2PXS::G2PXS(const char *model)
    :iSetting(1), iZ(1), iA(1), fTb(0.0), fTa(0.0), fEPS(10.0),
     fEPSD(-10.0), fFP(220.0), pfModelSelector(NULL)
{
    map<string, int> model_map;
    model_map["qfs"] = 2;

    if (model_map[model]==0) {
        printf("Unknown phys model, set to elastic ...\n");
        iSetting = 1;
    }
    else {
        iSetting = model_map[model];
    }
}

G2PXS::~G2PXS()
{
    // Nothing to do
}

void G2PXS::SetModel()
{
    switch (iSetting) {
    case 2:
        pfModelSelector = &G2PXS::GetXSQFS;
        break;
    }
}

double G2PXS::GetXSQFS(double Eb, double Ef, double theta)
{
    double result;
    
    qfs(iZ, iA, Eb, Ef, theta, fEPS, fEPSD, fFP, fTb, fTa, &result);

    return result;
}
