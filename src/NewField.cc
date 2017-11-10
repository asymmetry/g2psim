// -*- C++ -*-

/* class NewField
 * Generate field map for New Oxford magnets for SoLID NH3 target.
 * Calculate field strength of a particular point from the field map using Lagrange polynomial interpolation, default is 2nd order.
 *
 * Field map may have an angle respect to the lab coordinates.
 * Use SetEulerAngle() to set this angle and the program will rotate the field map to correct direction.
 */

// History:
//   Sep 2013, C. Gu, Put New Oxford field map into NewField.
//

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"

#include "G2PField.hh"

#include "NewField.hh"

using namespace std;

static const double kCM = 1.0e-2;
static const double kKG = 1.0e-1;
static const double kSCALE = 1.0;

extern "C" {
    void sda_ptf_solid_nh3_(float *, float *); // routine to calculate the field of HallB coil, which is not symmetric
}

static void GetNewField(const double *pos, double *field)
{
    float x[3], b[3];

    x[0] = pos[0] / kCM;
    x[1] = pos[1] / kCM;
    x[2] = pos[2] / kCM;

    sda_ptf_solid_nh3_(x, b);

    // The routine will return fields in kG, maximum is 5.0938709T
    // Match it to the TOSCA map maximum 4.9788476T
    field[0] = b[0] * kKG * kSCALE;
    field[1] = b[1] * kKG * kSCALE;
    field[2] = b[2] * kKG * kSCALE;
}

NewField *NewField::pNewField = NULL;

NewField::NewField() : G2PField()
{
    if (pNewField) {
        Error("NewField()", "Only one instance of NewField allowed.");
        MakeZombie();
        return;
    }

    fMapFile = "newfield.map";
}

NewField::~NewField()
{
    if (pNewField == this)
        pNewField = NULL;
}

int NewField::CreateMap()
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

            GetNewField(x, b);

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

ClassImp(NewField)
