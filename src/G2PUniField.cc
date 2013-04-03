#include <cstdlib>
#include <cmath>
#include <vector>

#include "TROOT.h"
#include "TObject.h"
#include "TError.h"

#include "G2PField.hh"
#include "G2PGlobals.hh"

#include "G2PUniField.hh"

using namespace std;

static const double kCM = 1.0e-2;

G2PUniField::G2PUniField()
{
    // Nothing to do
}

G2PUniField::~G2PUniField()
{
    // Nothing to do
}

int G2PUniField::Begin()
{
    static const char* const here = "Begin()";

    if (G2PField::Begin()!=0) return fStatus;

    fStatus = kERROR;
    if (CreateMap()==0) fStatus = kOK;
    else Error(here, "Cannot create field map.");

    return fStatus;
}

int G2PUniField::CreateMap()
{
    static const char* const here = "CreateMap()";

    if (!G2PField::CreateMap()) return -1;

    for (int i = 0; i<nR; i++) {
        for (int j = 0; j<nZ; j++) {
            fBField[i][j][0] = j*fZStep;
            fBField[i][j][1] = i*fRStep;
            fBField[i][j][2] = 1.0;
            fBField[i][j][3] = 0.0;
            fBField[i][j][4] = 1.0;

            if (fDebug>4) Info(here, "%10.3e %10.3e %10.3e %10.3e %10.3e",  fBField[i][j][0]/kCM, fBField[i][j][1]/kCM, fBField[i][j][2], fBField[i][j][3], fBField[i][j][4]);
        }
    }

    return 0;
}

ClassImp(G2PUniField)
