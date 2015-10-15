// -*- C++ -*-

/* class G2PField
 * Generate field map for HallB magnets.
 * Calculate field strength of a particular point from the field map using Lagrange polynomial interpolation, default is 2nd order.
 * G2PProcBase classes will call GetField() to get field values.
 *
 * Field map may have an angle respect to the lab coordinates.
 * Use SetEulerAngle() to set this angle and the program will rotate the field map to correct direction.
 */

// History:
//   Mar 2013, C. Gu, First public version.
//   Sep 2013, C. Gu, Put HallB field map into G2PField.
//
// TODO: Rewrite it as a base class.
//

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"

#include "G2PAppBase.hh"
#include "G2PGlobals.hh"
#include "G2PVarDef.hh"

#include "G2PField.hh"

using namespace std;

static const double kCM = 1.0e-2;
static const double kKG = 1.0e-1;
static const double kSCALE = 4.9788476 / 5.0938709;

extern "C" {
    void sda_ptf_(float *, float *); // routine to calculate the field of HallB coil, which is not symmetric
}

static void GetHallBField(const double *pos, double *field)
{
    float x[3], b[3];

    x[0] = pos[0] / kCM;
    x[1] = pos[1] / kCM;
    x[2] = pos[2] / kCM;

    sda_ptf_(x, b);

    // The routine will return fields in kG, maximum is 5.0938709T
    // Match it to the TOSCA map maximum 4.9788476T
    field[0] = b[0] * kKG * kSCALE;
    field[1] = b[1] * kKG * kSCALE;
    field[2] = b[2] * kKG * kSCALE;
}

G2PField *G2PField::pG2PField = NULL;

G2PField::G2PField() : fMapFile(NULL), fZMin(0.0), fZMax(2.99), fRMin(0.0), fRMax(2.99), fZStep(0.01), fRStep(0.01), nZ(0), nR(0), fRotation(false), fRatio(0.0)
{
    if (pG2PField) {
        Error("G2PField()", "Only one instance of G2PField allowed.");
        MakeZombie();
        return;
    }

    pG2PField = this;

    memset(fOrigin, 0, sizeof(fOrigin));
    memset(fEulerAngle, 0, sizeof(fEulerAngle));
    memset(fRotationMatrix, 0, sizeof(fRotationMatrix));
    fBField.clear();

    fMapFile = "hallbfield.map";
}

G2PField::~G2PField()
{
    if (pG2PField == this)
        pG2PField = NULL;
}

int G2PField::Begin()
{
    static const char *const here = "Begin()";

    if (G2PAppBase::Begin() != 0)
        return (fStatus = kBEGINERROR);

    nZ = int((fZMax - fZMin) / fZStep + 1e-8) + 1;
    nR = int((fRMax - fRMin) / fRStep + 1e-8) + 1;

    fBField.resize(nR);

    for (int i = 0; i < nR; i++) {
        fBField[i].resize(nZ);

        for (int j = 0; j < nZ; j++)
            fBField[i][j].resize(5, 0.0);
    }

    EStatus status = kBEGINERROR;

    if (ReadMap() == 0)
        status = kOK;
    else if (CreateMap() == 0)
        status = kOK;
    else
        Error(here, "Cannot create field map.");

    SetRotationMatrix();

    return (fStatus = status);
}

int G2PField::End()
{
    typedef vector<vector<vector<double> > >::size_type size_t;

    for (size_t i = 0; i < fBField.size(); i++) {
        for (size_t j = 0; j < fBField.size(); j++)
            fBField[i][j].clear();

        fBField[i].clear();
    }

    fBField.clear();

    return (G2PAppBase::End());
}

void G2PField::GetField(const double *x, double *b)
{
    static const char *const here = "GetField()";

    double pos[3], field[3];

    TransLab2Field(x, pos);
    Interpolate(pos, field, 2);
    TransField2Lab(field, b);
    b[0] *= fRatio;
    b[1] *= fRatio;
    b[2] *= fRatio;

    if (fDebug > 9)
        Info(here, "%10.3e %10.3e %10.3e : %10.3e %10.3e %10.3e", pos[0], pos[1], pos[2], b[0], b[1], b[2]);
}

void G2PField::SetRotationMatrix()
{
    // The Euler angle is defined using Z-X'-Z" convention

    if ((fabs(fEulerAngle[0]) < 1e-5) && (fabs(fEulerAngle[1]) < 1e-5) && (fabs(fEulerAngle[2]) < 1e-5))
        fRotation = false;
    else {
        double s1 = sin(fEulerAngle[0]);
        double c1 = cos(fEulerAngle[0]);
        double s2 = sin(fEulerAngle[1]);
        double c2 = cos(fEulerAngle[1]);
        double s3 = sin(fEulerAngle[2]);
        double c3 = cos(fEulerAngle[2]);

        fRotationMatrix[0][0][0] = c1 * c3 - c2 * s1 * s3;
        fRotationMatrix[0][1][0] = -c1 * s3 - c2 * s1 * c3;
        fRotationMatrix[0][2][0] = s2 * s1;
        fRotationMatrix[0][0][1] = s1 * c3 + c2 * c1 * s3;
        fRotationMatrix[0][1][1] = -s1 * s3 + c2 * c1 * c3;
        fRotationMatrix[0][2][1] = -s2 * c1;
        fRotationMatrix[0][0][2] = s2 * s3;
        fRotationMatrix[0][1][2] = s2 * c3;
        fRotationMatrix[0][2][2] = c2;

        // inverse matrix
        fRotationMatrix[1][0][0] = c1 * c3 - c2 * s1 * s3;
        fRotationMatrix[1][0][1] = -c1 * s3 - c2 * c3 * s1;
        fRotationMatrix[1][0][2] = s1 * s2;
        fRotationMatrix[1][1][0] = c3 * s1 + c1 * c2 * s3;
        fRotationMatrix[1][1][1] = c1 * c2 * c3 - s1 * s3;
        fRotationMatrix[1][1][2] = -c1 * s2;
        fRotationMatrix[1][2][0] = s2 * s3;
        fRotationMatrix[1][2][1] = c3 * s2;
        fRotationMatrix[1][2][2] = c2;

        fRotation = true;
    }
}

int G2PField::ReadMap()
{
    static const char *const here = "ReadMap()";

    ifstream ifs;
    int count = 0;

    ifs.open(fMapFile);

    if (ifs.fail())
        return -1;

    const int LEN = 300;
    char buff[LEN];

    ifs.getline(buff, LEN); // eat the first line

    while (ifs.getline(buff, LEN) != 0) {
        TString tmpline(buff);

        if (tmpline.EndsWith("\n"))
            tmpline.Chop();

        istringstream ist(tmpline.Data());
        string tmpstr;
        vector<string> line_spl;

        while (ist >> tmpstr)
            line_spl.push_back(tmpstr);

        if (line_spl.empty()) {
            continue;    // ignore empty lines
        }

        double tempZ = atof(line_spl[0].c_str()) * kCM;
        double tempR = atof(line_spl[1].c_str()) * kCM;

        if ((tempZ >= fZMin) && (tempZ <= fZMax) && (tempR >= fRMin) && (tempR <= fRMax)) {
            // store the value
            double tempBz = atof(line_spl[2].c_str());
            double tempBr = atof(line_spl[3].c_str());
            double tempB = atof(line_spl[4].c_str());

            int indexZ = int((tempZ - fZMin) / fZStep + 1e-8);
            int indexR = int((tempR - fRMin) / fRStep + 1e-8);

            fBField[indexR][indexZ][0] = tempZ;
            fBField[indexR][indexZ][1] = tempR;
            fBField[indexR][indexZ][2] = tempBz;
            fBField[indexR][indexZ][3] = tempBr;
            fBField[indexR][indexZ][4] = tempB;

            if (fDebug > 9)
                Info(here, "%10.3e %10.3e %10.3e %10.3e %10.3e", tempZ, tempR, tempBz, tempBr, tempB);

            count++;
        }
    }

    ifs.close();

    if (count == 0)
        return -1;

    return 0;
}

int G2PField::CreateMap()
{
    static const char *const here = "CreateMap()";

    FILE *fp;

    if ((fp = fopen(fMapFile, "w")) == NULL)
        return -1;

    fprintf(fp, "   z        r       Bz              Br              Btot\n");

    double x[3], b[3];

    x[1] = 0;

    for (int i = 0; i < nR; i++) {
        x[0] = i * fRStep;

        for (int j = 0; j < nZ; j++) {
            x[2] = j * fZStep;

            GetHallBField(x, b);

            fBField[i][j][0] = j * fZStep;
            fBField[i][j][1] = i * fRStep;
            fBField[i][j][2] = b[2];
            fBField[i][j][3] = sqrt(b[0] * b[0] + b[1] * b[1]);
            fBField[i][j][4] = sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);

            fprintf(fp, "%8.3f %8.3f\t%e\t%e\t%e\n", x[2] / kCM, x[0] / kCM, fBField[i][j][2], fBField[i][j][3], fBField[i][j][4]);

            if (fDebug > 9)
                Info(here, "%10.3e %10.3e %10.3e %10.3e %10.3e", x[2] / kCM, x[0] / kCM, fBField[i][j][2], fBField[i][j][3], fBField[i][j][4]);
        }
    }

    fclose(fp);

    return 0;
}

int G2PField::Interpolate(const double *pos, double *b, int order)
{
    // Calculate the nth order Lagrange polynomial interpolation
    // 1) Find out (B[R0][Z0],...,B[Rn][Zn]), which is a (n+1)by(n+1) matrix
    // 2) Interpolate by Z first, get (B[R0][Z],...,B[Rn][Z]), which is a 1x(n+1)
    //+matrix
    // 3) Interpolate by R to get B

    // Only interpolate in the first octant
    double r = sqrt(pos[0] * pos[0] + pos[1] * pos[1]);
    double z = fabs(pos[2]);

    if ((r > fRMax) || (z > fZMax) || (r < fRMin) || (z < fZMin)) {
        b[0] = 0.0;
        b[1] = 0.0;
        b[2] = 0.0;
        return -1;
    }

    int indexZ0 = int((z - fZMin) / fZStep) - int((order - 1) / 2);
    int indexR0 = int((r - fRMin) / fRStep) - int((order - 1) / 2);

    if (indexZ0 < 0)
        indexZ0 = 0;

    if (indexR0 < 0)
        indexR0 = 0;

    if (indexZ0 + order > nZ - 1)
        indexZ0 = nZ - 1 - order;

    if (indexR0 + order > nR - 1)
        indexR0 = nR - 1 - order;

    // Lagrange polynomial interpolate on Z
    int indexZ, indexR, indexT;
    double tempBz[order + 1], tempBr[order + 1];

    for (int i = 0; i <= order; i++) {
        indexR = indexR0 + i;
        tempBz[i] = 0.0;
        tempBr[i] = 0.0;

        for (int j = 0; j <= order; j++) {
            double temp = 1.0;
            indexZ = indexZ0 + j;

            for (int k = 0; k <= order; k++) {
                indexT = indexZ0 + k;

                if (k != j)
                    temp *= (z - fBField[indexR][indexT][0]) / (fBField[indexR][indexZ][0] - fBField[indexR][indexT][0]);
            }

            tempBz[i] += temp * fBField[indexR][indexZ][2];
            tempBr[i] += temp * fBField[indexR][indexZ][3];
        }
    }

    // Lagrange polynomial interpolate on R
    double Bz = 0, Br = 0;
    indexZ = indexZ0;

    for (int i = 0; i <= order; i++) {
        indexR = indexR0 + i;
        double temp = 1.0;

        for (int j = 0; j <= order; j++) {
            indexT = indexR0 + j;

            if (j != i)
                temp *= (r - fBField[indexT][indexZ][1]) / (fBField[indexR][indexZ][1] - fBField[indexT][indexZ][1]);
        }

        Bz += temp * tempBz[i];
        Br += temp * tempBr[i];
    }

    // Sign has already been considered
    if (r == 0) {
        b[0] = 0.0;
        b[1] = 0.0;
        b[2] = Bz;
    } else if (z == 0) {
        b[0] = Br * (pos[0] / r);
        b[1] = Br * (pos[1] / r);
        b[2] = Bz;
    } else {
        b[0] = Br * (pos[0] / r) * (pos[2] / z);
        b[1] = Br * (pos[1] / r) * (pos[2] / z);
        b[2] = Bz;
    }

    return 0;
}

void G2PField::TransLab2Field(const double *x, double *xout)
{
    double temp[3] = {x[0], x[1], x[2]};

    for (int i = 0; i < 3; i++) temp[i] -= fOrigin[i];

    if (fRotation) {
        for (int i = 0; i < 3; i++) xout[i] = fRotationMatrix[0][i][0] * temp[0] + fRotationMatrix[0][i][1] * temp[1] + fRotationMatrix[0][i][2] * temp[2];
    } else {
        for (int i = 0; i < 3; i++) xout[i] = temp[i];
    }
}

void G2PField::TransField2Lab(const double *b, double *bout)
{
    if (fRotation) {
        for (int i = 0; i < 3; i++) bout[i] = fRotationMatrix[1][i][0] * b[0] + fRotationMatrix[1][i][1] * b[1] + fRotationMatrix[1][i][2] * b[2];
    } else {
        for (int i = 0; i < 3; i++) bout[i] = b[i];
    }
}

#ifdef DEBUGWITHROOT
void G2PField::SaveRootFile()
{
    TFile *file = new TFile("field.root", "RECREATE");
    TTree *field = new TTree("field", "field map");

    double x, y, z, r;
    double Bx, By, Bz, Br, Btot;

    field->Branch("x", &x, "x/D");
    field->Branch("y", &y, "y/D");
    field->Branch("r", &r, "r/D");
    field->Branch("z", &z, "z/D");
    field->Branch("Bx", &Bx, "Bx/D");
    field->Branch("By", &By, "By/D");
    field->Branch("Br", &Br, "Br/D");
    field->Branch("Bz", &Bz, "Bz/D");
    field->Branch("Btot", &Btot, "Btot/D");

    for (x = -1; x <= 1; x += 0.02)
        for (y = -1; y <= 1; y += 0.02)
            for (z = -1; z <= 1; z += 0.02) {
                double xx[3];
                xx[0] = x;
                xx[1] = y;
                xx[2] = z;
                double b[3];
                GetField(xx, b);
                Bx = b[0];
                By = b[1];
                Bz = b[2];
                Btot = sqrt(Bx * Bx + By * By + Bz * Bz);
                Br = sqrt(Bx * Bx + By * By);
                field->Fill();
            }

    field->Write();
    delete field;
    field = NULL;

    file->Close();
    delete file;
    file = NULL;
}
#endif

int G2PField::Configure(EMode mode)
{
    if ((mode == kREAD || mode == kTWOWAY) && fConfigured)
        return 0;

    if (G2PAppBase::Configure(mode) != 0)
        return -1;

    ConfDef confs[] = {
        {"ratio", "Field Ratio", kDOUBLE, &fRatio},
        {"origin.x", "Origin X", kDOUBLE, &fOrigin[0]},
        {"origin.y", "Origin Y", kDOUBLE, &fOrigin[1]},
        {"origin.z", "Origin Z", kDOUBLE, &fOrigin[2]},
        {"angle.alpha", "Euler Angle Alpha", kDOUBLE, &fEulerAngle[0]},
        {"angle.beta", "Euler Angle Beta", kDOUBLE, &fEulerAngle[1]},
        {"angle.gamma", "Euler Angle Gamma", kDOUBLE, &fEulerAngle[2]},
        {"r.min", "R Range", kDOUBLE, &fRMin},
        {"r.max", "R Range", kDOUBLE, &fRMax},
        {"r.step", "R Step", kDOUBLE, &fRStep},
        {"z.min", "Z Range", kDOUBLE, &fZMin},
        {"z.max", "Z Range", kDOUBLE, &fZMax},
        {"z.step", "Z Step", kDOUBLE, &fZStep},
        {0}
    };

    return ConfigureFromList(confs, mode);
}

void G2PField::MakePrefix()
{
    const char *base = "field";

    G2PAppBase::MakePrefix(base);
}

ClassImp(G2PField)
