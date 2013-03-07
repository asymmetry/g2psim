#include <cstdlib>
#include <cmath>
#include <vector>

#include "TROOT.h"
#include "TObject.h"
#include "TError.h"

#include "G2PFieldBase.hh"
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

G2PAppsBase::EStatus G2PUniField::Init()
{
    static const char* const here = "Init()";

    if (G2PFieldBase::Init()) return fStatus;

    fStatus = kINITERROR;
    if (CreateMap()) fStatus = kOK;
    else Error(here, "Cannot initialize.");

    if (fDebug>3) SaveRootFile();

    return fStatus;
}

bool G2PUniField::CreateMap()
{
    static const char* const here = "CreateMap()";

    if (!G2PFieldBase::CreateMap()) return false;

    for (int i = 0; i<nNumR; i++) {
        for (int j = 0; j<nNumZ; j++) {
            fBField[i][j][0] = j*fStepZ;
            fBField[i][j][1] = i*fStepR;
            fBField[i][j][2] = 1.0;
            fBField[i][j][3] = 0.0;
            fBField[i][j][4] = 1.0;

            if (fDebug>3) Info(here, "%10.3e %10.3e %10.3e %10.3e %10.3e",  fBField[i][j][0]/kCM, fBField[i][j][1]/kCM, fBField[i][j][2], fBField[i][j][3], fBField[i][j][4]);
        }
    }

    return true;
}

ClassImp(G2PUniField)
