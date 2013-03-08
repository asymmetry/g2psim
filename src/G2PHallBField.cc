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

#include "G2PHallBField.hh"

using namespace std;

extern "C"
{
    void sda_ptf_(float*, float*); // routine to calculate the field of HallB coil, which is not symmetric
}

static const double kCM = 1.0e-2;
static const double kKG = 1.0e-1;
static const double kSCALE = 4.9788476/5.0938709;

static void GetHallBField(const double* pos, double* field)
{
    float x[3], b[3];

    x[0] = pos[0]/kCM;
    x[1] = pos[1]/kCM;
    x[2] = pos[2]/kCM;

    sda_ptf_(x, b);

    // The routine will return fields in kG, maximum is 5.0938709T
    // Match it to the TOSCA map maximum 4.9788476T
    field[0] = b[0]*kKG*kSCALE;
    field[1] = b[1]*kKG*kSCALE;
    field[2] = b[2]*kKG*kSCALE;
}

G2PHallBField::G2PHallBField()
{
    pMapFileName = "hallbfield.map";
}

G2PHallBField::~G2PHallBField()
{
    // Nothing to do
}

G2PAppsBase::EStatus G2PHallBField::Init()
{
    static const char* const here = "Init()";

    if (G2PFieldBase::Init()) return fStatus;

    fStatus = kINITERROR;
    if (ReadMap()) fStatus = kOK;
    else if (CreateMap()) fStatus = kOK;
    else Error(here, "Cannot initialize.");

    if (fDebug>4) SaveRootFile();

    return fStatus;
}

bool G2PHallBField::ReadMap()
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

            if (fDebug>4) Info(here, "%10.3e %10.3e %10.3e %10.3e %10.3e", tempZ, tempR, tempBz, tempBr, tempB);

            count++;
        }
    }

    ifs.close();

    if (count==0) return false;
    return true;
}

bool G2PHallBField::CreateMap()
{
    static const char* const here = "CreateMap()";

    if (!G2PFieldBase::CreateMap()) return false;

    FILE* fp;

    if ((fp=fopen(pMapFileName,"w"))==NULL) return false;

    fprintf(fp, "   z        r       Bz              Br              Btot\n");
    double x[3], b[3];

    x[1] = 0;
    
    for (int i = 0; i<nNumR; i++) {
        x[0] = i*fStepR;
        for (int j = 0; j<nNumZ; j++) {
            x[2] = j*fStepZ;

            GetHallBField(x, b);

            fBField[i][j][0] = j*fStepZ;
            fBField[i][j][1] = i*fStepR;
            fBField[i][j][2] = b[2];
            fBField[i][j][3] = sqrt(b[0]*b[0]+b[1]*b[1]);
            fBField[i][j][4] = sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);
            
            fprintf(fp, "%8.3f %8.3f\t%e\t%e\t%e\n", x[2]/kCM, x[0]/kCM, fBField[i][j][2], fBField[i][j][3], fBField[i][j][4]);

            if (fDebug>4) Info(here, "%10.3e %10.3e %10.3e %10.3e %10.3e", x[2]/kCM, x[0]/kCM, fBField[i][j][2], fBField[i][j][3], fBField[i][j][4]);
        }
    }

    fclose(fp);

    return true;
}

ClassImp(G2PHallBField)
