#include <cmath>

#include "G2PTrans.hh"

G2PTrans::G2PTrans()
{
    // Nothing to do
}

G2PTrans::~G2PTrans()
{
    // Nothing to do
}

void G2PTrans::FPCorrection(double *v)
{
    // Nothing to do
}

void G2PTrans::CoordsCorrection(double angle, double* v)
{
    double cosangle = cos(angle);
    double sinangle = sin(angle);

    double x_p = v[0];
    double y_p = cosangle*v[2];
    double z_p = -sinangle*v[2];

    double p_z_1 = 1.0;
    double p_x_1 = v[1];
    double p_y_1 = v[3];

    double p_x_2 = p_x_1;
    double p_y_2 = p_y_1*cosangle+p_z_1*sinangle;
    double p_z_2 = -p_y_1*sinangle+p_z_1*cosangle;

    double t_p = p_x_2/p_z_2;
    double p_p = p_y_2/p_z_2;

    v[0] = x_p-t_p*z_p;
    v[1] = t_p;
    v[2] = y_p-p_p*z_p;
    v[3] = p_p;
}
