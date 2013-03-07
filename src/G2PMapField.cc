#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "TROOT.h"
#include "TObject.h"
#include "TError.h"

#include "G2PFieldBase.hh"
#include "G2PGlobals.hh"

#include "G2PMapField.hh"

using namespace std;

static const double kCM = 1.0e-2;

G2PMapField::G2PMapField(const char* name)
{
    pMapFileName = name;
}

G2PMapField::G2PMapField()
{
    // Nothing to do
}

G2PMapField::~G2PMapField()
{
    // Nothing to do
}

G2PAppsBase::EStatus G2PMapField::Init()
{
    static const char* const here = "Init()";

    if (G2PFieldBase::Init()) return fStatus;

    fStatus = kINITERROR;
    if (ReadMap()) fStatus = kOK;
    else Error(here, "Cannot initialize.");

    if (fDebug>3) SaveRootFile();

    return fStatus;
}

bool G2PMapField::ReadMap()
{
    static const char* const here = "ReadMap()";

    if (!G2PFieldBase::ReadMap()) return false;

    ifstream ifs;
    int count = 0;

    ifs.open(pMapFileName);
    if (ifs.fail()) return false;

    const int LEN = 300;
    char buff[LEN];

    ifs.getline(buff, LEN); // eat the first line

    while (ifs.getline(buff, LEN)!=NULL) {
        TString tmpline(buff);

        if (tmpline.EndsWith("\n")) tmpline.Chop();

        istringstream ist(tmpline.Data());
        string tmpstr;
        vector<string> line_spl;
        while (ist>>tmpstr) {
            line_spl.push_back(tmpstr);
        }

        if (line_spl.empty()) continue; // ignore empty lines

        double tempZ = atof(line_spl[0].c_str())*kCM;
        double tempR = atof(line_spl[1].c_str())*kCM;

        if ((tempZ>=fZMin)&&(tempZ<=fZMax)&&(tempR>=fRMin)&&(tempR<=fRMax)) {
            // store the value
            double tempBz = atof(line_spl[2].c_str());
            double tempBr = atof(line_spl[3].c_str());
            double tempB  = atof(line_spl[4].c_str());

            int indexZ = int((tempZ-fZMin)/fStepZ+1e-8);
            int indexR = int((tempR-fRMin)/fStepR+1e-8);

            fBField[indexR][indexZ][0] = tempZ;
            fBField[indexR][indexZ][1] = tempR;
            fBField[indexR][indexZ][2] = tempBz;
            fBField[indexR][indexZ][3] = tempBr;
            fBField[indexR][indexZ][4] = tempB;

            if (fDebug>3) Info(here, "%10.3e %10.3e %10.3e %10.3e %10.3e", tempZ, tempR, tempBz, tempBr, tempB);

            count++;
        }
    }

    ifs.close();

    if (count==0) return false;
    return true;
}

ClassImp(G2PMapField)
