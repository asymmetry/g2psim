// -*- C++ -*-

/* class G2PTrans484816R00
 * 484816 septum with shims, 5.767 central ray, no target field, no raster
 * By M. Huang 1/7/2013
 */

// History:
//   Sep 2013, C. Gu, First public version.
//

#include <cmath>

#include "Fwd_l5p77_484816R00.h"
#include "Bwd_l5p77_484816R00.h"

#include "G2PTrans484816R00.hh"

#define CORRECTION_ORDER -1

using namespace S484816R00;

const double kDEG = 3.14159265358979323846 / 180.0;

G2PTrans484816R00::G2PTrans484816R00()
{
    fModelAngle = 5.77 * kDEG;
}

G2PTrans484816R00::~G2PTrans484816R00()
{
    // Nothing to do
}

int G2PTrans484816R00::TransLeftHRS(double *pV5, double *PlanePosX, double *PlanePosY)
{
    float vector_jjl[] = {float(pV5[0]), float(pV5[1]), float(pV5[2]), float(pV5[3]), float(pV5[4])};
    int ii = 5;

    // Succesfully reach focus plane
    float x_fp = x_l5p77_fp(vector_jjl, ii);
    float theta_fp = t_l5p77_fp(vector_jjl, ii);
    float y_fp = y_l5p77_fp(vector_jjl, ii);
    float phi_fp = p_l5p77_fp(vector_jjl, ii);

    // Reset the vector and return it back to the caller
    pV5[0] = (double) x_fp;
    pV5[1] = (double) theta_fp;
    pV5[2] = (double) y_fp;
    pV5[3] = (double) phi_fp;
    //pV5[4] = (double)delta_fp;  // delta is not change

    return 0;
}

int G2PTrans484816R00::TransRightHRS(double *pV5, double *PlanePosX, double *PlanePosY)
{
    // Use right arm routines for left arm before left arm is ready

    pV5[2] *= -1.;
    pV5[3] *= -1.;
    int fGoodParticle = TransLeftHRS(pV5, PlanePosX, PlanePosY);
    pV5[2] *= -1.;
    pV5[3] *= -1.;

    return fGoodParticle;
}

void G2PTrans484816R00::ReconLeftHRS(double *pV5)
{
    float vector_jjl[] = {float(pV5[0]), float(pV5[1]), float(pV5[2]), float(pV5[3]), float(pV5[4])};
    int ii = 5;

    vector_jjl[1] = vector_jjl[1] - txfit_l5p77(vector_jjl, ii);

    float x_or = vector_jjl[4];
    float delta_rec = delta_l5p77(vector_jjl, ii);
    float theta_rec = theta_l5p77(vector_jjl, ii);
    float phi_rec = phi_l5p77(vector_jjl, ii);
    float y_rec = 0.0;

    // Reset the vector and return it back to the caller
    pV5[0] = (double) x_or;
    pV5[1] = (double) theta_rec;
    pV5[2] = (double) y_rec;
    pV5[3] = (double) phi_rec;
    pV5[4] = (double) delta_rec;
}

void G2PTrans484816R00::ReconRightHRS(double *pV5)
{
    // In order to call left arm routines, need to flip y, phi
    pV5[2] *= -1;
    pV5[3] *= -1;
    ReconLeftHRS(pV5);
    pV5[2] *= -1;
    pV5[3] *= -1;
}

void G2PTrans484816R00::FPCorrLeft(const double *V5tg, double *V5fp)
{
#if CORRECTION_ORDER == 0
    V5fp[0] += (-5.91123e-03);
    V5fp[1] += (-6.76146e-05);
    V5fp[2] += (-1.89684e-02);
    V5fp[3] += (-1.11524e-02);
#endif
#if CORRECTION_ORDER == 1
    //    V5fp[0] += (-5.91123e-03)+(-9.98542e-02) * V5tg[0]+(-2.96956e-02) * V5tg[1]+(0.00000e+00) * V5tg[2]+(0.00000e+00) * V5tg[3]+(-3.22433e-03) * V5tg[4];
    //    V5fp[1] += (-6.76146e-05)+(7.34029e-02) * V5tg[0]+(1.78179e-02) * V5tg[1]+(0.00000e+00) * V5tg[2]+(0.00000e+00) * V5tg[3]+(9.88179e-02) * V5tg[4];
    //    V5fp[2] += (-1.89684e-02)+(0.00000e+00) * V5tg[0]+(0.00000e+00) * V5tg[1]+(3.54528e-02) * V5tg[2]+(-9.47126e-04) * V5tg[3]+(-6.55867e-02) * V5tg[4];
    //    V5fp[3] += (-1.11524e-02)+(0.00000e+00) * V5tg[0]+(0.00000e+00) * V5tg[1]+(4.73617e-02) * V5tg[2]+(1.11894e-03) * V5tg[3]+(-6.39965e-02) * V5tg[4];
    //
    //    V5fp[0] += (-1.50923e-03)+(0.00000e+00) * V5tg[0]+(-2.98481e-02) * V5tg[1]+(9.98334e-03) * V5tg[2]+(-8.19924e-02) * V5tg[3]+(9.76081e+00) * V5tg[4];
    //    V5fp[1] += (1.35902e-03)+(0.00000e+00) * V5tg[0]+(2.29618e-02) * V5tg[1]+(-2.04595e-02) * V5tg[2]+(-2.20224e-03) * V5tg[3]+(5.30659e-01) * V5tg[4];
    //    V5fp[2] += (-1.78416e-02)+(0.00000e+00) * V5tg[0]+(9.90491e-03) * V5tg[1]+(-9.49504e-02) * V5tg[2]+(2.52013e-02) * V5tg[3]+(1.00987e+01) * V5tg[4];
    //    V5fp[3] += (-9.72573e-03)+(0.00000e+00) * V5tg[0]+(4.58426e-02) * V5tg[1]+(-9.99992e-02) * V5tg[2]+(1.55553e-02) * V5tg[3]+(1.00998e+01) * V5tg[4];

    V5fp[0] += (-1.50612e-03) + (0.00000e+00) * V5tg[0] + (-2.72655e-02) * V5tg[1] + (-9.48539e-03) * V5tg[2] + (-7.25842e-02) * V5tg[3] + (8.58806e+00) * V5tg[4];
    V5fp[1] += (1.35279e-03) + (0.00000e+00) * V5tg[0] + (2.23088e-02) * V5tg[1] + (1.21097e-02) * V5tg[2] + (-2.95080e-03) * V5tg[3] + (6.71229e-01) * V5tg[4];
    V5fp[2] += (-1.83197e-02) + (0.00000e+00) * V5tg[0] + (9.00558e-03) * V5tg[1] + (-5.02045e-02) * V5tg[2] + (-2.48741e-03) * V5tg[3] + (9.62701e+00) * V5tg[4];
    V5fp[3] += (-1.01323e-02) + (0.00000e+00) * V5tg[0] + (4.51187e-02) * V5tg[1] + (-9.76615e-02) * V5tg[2] + (-7.28207e-03) * V5tg[3] + (9.69853e+00) * V5tg[4];
#endif
}

void G2PTrans484816R00::FPCorrRight(const double *V5tg, double *V5fp)
{
    // Nothing to do
}
