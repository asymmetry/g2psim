// -*- C++ -*-

/* class G2PTrans400016OLD
 * 400016 septum with shims, 5.65 central ray, no target field, 3cm raster
 * By M. Huang 11/12/2012
 */

// History:
//   Sep 2013, C. Gu, First public version.
//

#include <cmath>

#include "Fwd_l5p65_400016OLD.h"
#include "Bwd_l5p65_400016OLD.h"

#include "G2PTrans400016OLD.hh"

using namespace S400016OLD;

const float m2cm = 100.0;
const double kDEG = 3.14159265358979323846 / 180.0;

G2PTrans400016OLD::G2PTrans400016OLD()
{
    fModelAngle = 5.65 * kDEG;
}

G2PTrans400016OLD::~G2PTrans400016OLD()
{
    // Nothing to do
}

int G2PTrans400016OLD::TransLeftHRS(double* pV5, double* PlanePosX, double* PlanePosY)
{
    float vector_jjl[] = {pV5[0], pV5[1], pV5[2], pV5[3], pV5[4]};
    int ii = 5;

    float x_test, y_test;

    // Target to Septum ep5
    x_test = x_l5p65_ep5(vector_jjl, ii) * m2cm;
    y_test = y_l5p65_ep5(vector_jjl, ii) * m2cm;
    PlanePosX[5] = x_test;
    PlanePosY[5] = y_test;
    if (fabs(x_test) < 8.4 || fabs(x_test) > 38.8 || fabs(y_test) > 9.7)
        return 5;

    // Target to Septum ep7
    x_test = x_l5p65_ep7(vector_jjl, ii) * m2cm;
    y_test = y_l5p65_ep7(vector_jjl, ii) * m2cm;
    PlanePosX[7] = x_test;
    PlanePosY[7] = y_test;
    if (fabs(x_test) < 8.4 || fabs(x_test) > 38.8 || fabs(y_test) > 9.7)
        return 7;

    // Target to Q1 en ep10
    x_test = x_l5p65_ep10_q1en(vector_jjl, ii) * m2cm;
    y_test = y_l5p65_ep10_q1en(vector_jjl, ii) * m2cm;
    PlanePosX[10] = x_test;
    PlanePosY[10] = y_test;
    if (sqrt(x_test * x_test + y_test * y_test) > 14.92)
        return 10;

    // Target to Q1 ex ep13
    x_test = x_l5p65_ep13_q1ex(vector_jjl, ii) * m2cm;
    y_test = y_l5p65_ep13_q1ex(vector_jjl, ii) * m2cm;
    PlanePosX[13] = x_test;
    PlanePosY[13] = y_test;
    if (sqrt(x_test * x_test + y_test * y_test) > 30.)
        return 13;

    // Target to dipole exit, ep24
    // trapezoid, -46.19cm < x < 46.19cm, |y| < -0.0161 * x + 12.5
    x_test = x_l5p65_ep24_dex(vector_jjl, ii) * m2cm;
    y_test = y_l5p65_ep24_dex(vector_jjl, ii) * m2cm;
    PlanePosX[24] = x_test;
    PlanePosY[24] = y_test;
    if ((x_test < -46.19) || (x_test > 46.19) || fabs(y_test) > fabs(-0.0161 * x_test + 12.5))
        return 24;

    // Target to Q3 entrance
    // circle of radius 30.0 cm
    x_test = x_l5p65_ep26_q3en(vector_jjl, ii) * m2cm;
    y_test = y_l5p65_ep26_q3en(vector_jjl, ii) * m2cm;
    PlanePosX[26] = x_test;
    PlanePosY[26] = y_test;
    if (sqrt(x_test * x_test + y_test * y_test) > 30.0)
        return 26;

    // Target to Q3 exit
    // circle of radius 30.0 cm
    x_test = x_l5p65_ep29_q3ex(vector_jjl, ii) * m2cm;
    y_test = y_l5p65_ep29_q3ex(vector_jjl, ii) * m2cm;
    PlanePosX[29] = x_test;
    PlanePosY[29] = y_test;
    if (sqrt(x_test * x_test + y_test * y_test) > 30.0)
        return 29;

    // Successfully reach focus plane
    float x_fp = x_l5p65_fp(vector_jjl, ii);
    float theta_fp = t_l5p65_fp(vector_jjl, ii);
    float y_fp = y_l5p65_fp(vector_jjl, ii);
    float phi_fp = p_l5p65_fp(vector_jjl, ii);

    // Reset the vector and return it back to the caller
    pV5[0] = (double) x_fp;
    pV5[1] = (double) theta_fp;
    pV5[2] = (double) y_fp;
    pV5[3] = (double) phi_fp;
    //pV5[4] = (double)delta_fp; // delta is not change

    return 0;
}

int G2PTrans400016OLD::TransRightHRS(double* pV5, double* PlanePosX, double* PlanePosY)
{
    // Use left arm routines for right arm before right arm is ready

    pV5[2] *= -1.;
    pV5[3] *= -1.;
    int fGoodParticle = TransLeftHRS(pV5, PlanePosX, PlanePosY);
    pV5[2] *= -1.;
    pV5[3] *= -1.;

    return fGoodParticle;
}

void G2PTrans400016OLD::ReconLeftHRS(double* pV5)
{
    float vector_jjl[] = {pV5[0], pV5[1], pV5[2], pV5[3], pV5[4]};
    int ii = 5;

    vector_jjl[1] = vector_jjl[1] - txfit_l5p65(vector_jjl, ii);

    float x_or = vector_jjl[4];
    float delta_rec = delta_l5p65(vector_jjl, ii);
    float theta_rec = theta_l5p65(vector_jjl, ii);
    float phi_rec = phi_l5p65(vector_jjl, ii);
    float y_rec = y00_l5p65(vector_jjl, ii);

    // Reset the vector and return it back to the caller
    pV5[0] = (double) x_or;
    pV5[1] = (double) theta_rec;
    pV5[2] = (double) y_rec;
    pV5[3] = (double) phi_rec;
    pV5[4] = (double) delta_rec;
}

void G2PTrans400016OLD::ReconRightHRS(double* pV5)
{
    // In order to call left arm routines, need to flip y, phi
    pV5[2] *= -1;
    pV5[3] *= -1;
    ReconLeftHRS(pV5);
    pV5[2] *= -1;
    pV5[3] *= -1;
}
