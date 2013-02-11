#include "G2PHallBField.hh"

extern "C"
{
    void sda_ptf_(float*, float*); // routine to calculate the field of HallB coil, which is not symmetric
}

const double kCM = 1.0e-2;
const double kKG = 1.0e-1;
const double kSCALE = 4.9788476/5.0938709;

G2PHallBField::G2PHallBField()
{
    // Nothing to do
}

G2PHallBField::~G2PHallBField()
{
    // Nothing to do
}

void G2PHallBField::operator() (const double* pos, double* field)
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
