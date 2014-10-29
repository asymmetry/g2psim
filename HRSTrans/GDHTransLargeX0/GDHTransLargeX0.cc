// -*- C++ -*-

/* class GDHTransLargeX0
 * 6 deg with septum, for GDH experiment (E97110), large X0
 * By J.J. LeRose 10/05/2012
 */

// History:
//   Sep 2013, J. Zhang, First public version.
//

#include <cmath>

#include "Bwd_R6_LargeX0_GDH.hh"
#include "Fwd_R6_LargeX0_GDH.hh"

#include "GDHTransLargeX0.hh"

const float m2cm = 100.0;
const double kDEG = 3.14159265358979323846 / 180.0;

GDHTransLargeX0::GDHTransLargeX0()
{
    fModelAngle = 6.0 * kDEG;
}

GDHTransLargeX0::~GDHTransLargeX0()
{
    // Nothing to do
}

int GDHTransLargeX0::TransLeftHRS(double* pV5)
{
    // Use right arm routines for left arm before left arm is ready

    pV5[2] *= -1.;
    pV5[3] *= -1.;
    int fGoodParticle = TransRightHRS(pV5);
    pV5[2] *= -1.;
    pV5[3] *= -1.;

    return fGoodParticle;
}

int GDHTransLargeX0::TransRightHRS(double* pV5)
{
    float vector_jjl[] = {pV5[0], pV5[1], pV5[2], pV5[3], pV5[4]};
    float x_test, y_test;
    int ii = 5;

    // the detail of input vector is
    // vector_jjl[0] = x_tr;
    // vector_jjl[1] = theta_tr;
    // vector_jjl[2] = y_tr;
    // vector_jjl[3] = phi_tr;
    // vector_jjl[4] = delta_tr;
    // the output is
    // vector_jjl[0] = x_fp;
    // vector_jjl[1] = theta_fp;
    // vector_jjl[2] = y_fp;
    // vector_jjl[3] = phi_fp;
    // vector_jjl[4] = delta_fp; // delta is not change

    // ep3: -0.1486233 < x < -0.08869672
    //          -0.110 < y < 0.110
    // ep4: -0.1792231 < x < -0.1089169
    //          -0.110 < y < 0.110
    // ep5: -0.2209211 < x < -0.1353789
    //          -0.110 < y < 0.110
    // ep6: -0.2763536 < x < -0.1697464
    //          -0.110 < y < 0.110
    // ep7: -0.3485396 < x < -0.2156404
    //          -0.110 < y < 0.110
    // q1ex is a circle of radius 0.1492 m
    // dent is a trapazoid:
    //                          -5.22008 < x < -4.98099
    // -(-0.192436784 * x - 0.192436784) < y < -0.192436784 * x - 0.192436784
    // dext is also a trapazoid:
    //                          -0.46188 < x < 0.46188
    //        -(-0.01610808 * x + 0.125) < y < -0.01610808 * x + 0.125
    // q3en is a circle of radius 0.3 m
    // q3ex is a circle of radius 0.3 m

    float y_min = -9.9;
    float y_max = 9.9;

    // Target to Septum entrance, -14.06cm < x < -8.87cm, -9.9cm < y < 9.9cm
    x_test = x_sr6_largex0_ep3_(vector_jjl, &ii) * m2cm;
    y_test = y_sr6_largex0_ep3_(vector_jjl, &ii) * m2cm;
    if ((x_test < -14.06) || (x_test > -8.87) || (y_test < y_min) || (y_test > y_max))
        return 5;

    // Target to 1/4 Septum, -17.12cm < x < -10.89cm, -9.9cm < y < 9.9cm
    x_test = x_sr6_largex0_ep4_(vector_jjl, &ii) * m2cm;
    y_test = y_sr6_largex0_ep4_(vector_jjl, &ii) * m2cm;
    if ((x_test < -17.12) || (x_test > -10.89) || (y_test < y_min) || (y_test > y_max))
        return 5;

    // Target to 1/2 Septum, -21.29cm < x < -13.54cm, -9.9cm < y < 9.9cm
    x_test = x_sr6_largex0_ep5_(vector_jjl, &ii) * m2cm;
    y_test = y_sr6_largex0_ep5_(vector_jjl, &ii) * m2cm;
    if ((x_test < -21.29) || (x_test > -13.54) || (y_test < y_min) || (y_test > y_max))
        return 5;

    // Target to 3/4 Septum, -26.84cm < x < -16.97cm, -9.9cm < y < 9.9cm
    x_test = x_sr6_largex0_ep6_(vector_jjl, &ii) * m2cm;
    y_test = y_sr6_largex0_ep6_(vector_jjl, &ii) * m2cm;
    if ((x_test < -26.84) || (x_test > -16.97) || (y_test < y_min) || (y_test > y_max))
        return 5;

    // Target to Septum exit, -34.05cm < x < -21.56cm, -9.9cm < y < 9.9cm
    x_test = x_sr6_largex0_ep7_(vector_jjl, &ii) * m2cm;
    y_test = y_sr6_largex0_ep7_(vector_jjl, &ii) * m2cm;
    if ((x_test < -34.05) || (x_test > -21.56) || (y_test < y_min) || (y_test > y_max))
        return 5;

    // Target to Q1 exit
    // circle of radius 14.92 cm
    x_test = x_sr6_largex0_q1ex_(vector_jjl, &ii) * m2cm;
    y_test = y_sr6_largex0_q1ex_(vector_jjl, &ii) * m2cm;
    x_test = x_test + 0.9;
    if ((x_test * x_test + y_test * y_test) > (14.92 * 14.92))
        return 5;

    // Target to dipole entrance
    // trapezoid, -522.0cm < x < -498.1cm, |y| < -0.1924 * x - 19.24
    x_test = x_sr6_largex0_dent_(vector_jjl, &ii) * m2cm;
    y_test = y_sr6_largex0_dent_(vector_jjl, &ii) * m2cm;
    if ((x_test < -522.0) || (x_test > -498.1) || fabs(y_test) > fabs(-0.1924 * x_test - 19.24))
        return 5;

    // Target to dipole exit
    // trapezoid, -46.19cm < x < 46.19cm, |y| < -0.0161 * x + 12.5
    x_test = x_sr6_largex0_dext_(vector_jjl, &ii) * m2cm;
    y_test = y_sr6_largex0_dext_(vector_jjl, &ii) * m2cm;
    if (fabs(x_test) > 46.19 || fabs(y_test) > fabs(-0.0161 * x_test + 12.5))
        return 5;

    // Target to Q3 entrance
    // circle of radius 30.0 cm
    x_test = x_sr6_largex0_q3en_(vector_jjl, &ii) * m2cm;
    y_test = y_sr6_largex0_q3en_(vector_jjl, &ii) * m2cm;
    if ((x_test * x_test + y_test * y_test) > (30.0 * 30.0))
        return 5;

    // Target to Q3 exit
    // circle of radius 30.0 cm -> 28.0cm
    x_test = x_sr6_largex0_q3ex_(vector_jjl, &ii) * m2cm;
    y_test = y_sr6_largex0_q3ex_(vector_jjl, &ii) * m2cm;
    x_test = (x_test - 1.0) / (28.0);
    y_test = y_test / (30.0);
    if ((x_test * x_test + y_test * y_test) > 1.0)
        return 5;

    // If we reach this point, it means the test was successful
    float x_fp = x_sr6_largex0_fp_(vector_jjl, &ii);
    float theta_fp = t_sr6_largex0_fp_(vector_jjl, &ii);
    float y_fp = y_sr6_largex0_fp_(vector_jjl, &ii);
    float phi_fp = p_sr6_largex0_fp_(vector_jjl, &ii);

    // Reset the vector and return it back to the caller
    pV5[0] = (double) x_fp;
    pV5[1] = (double) theta_fp;
    pV5[2] = (double) y_fp;
    pV5[3] = (double) phi_fp;
    //pV5[4] = (double)delta_fp;  // delta is not change

    return 0;
}

void GDHTransLargeX0::ReconLeftHRS(double *pV5)
{
    // I am using the right arm routines to do the reconstruction, JJL did not provide routines for the 2nd version of GDH experiment

    // In order to call right arm routines, need to flip y, phi
    float vector_jjl[] = {pV5[0], pV5[1], -pV5[2], -pV5[3], pV5[4]};
    int ii = 5, jj = 1;

    // Orthogonalize theta as JJL asks
    vector_jjl[1] = vector_jjl[1] - r6_largex0_txfit_(vector_jjl, &jj);
    float x_or = vector_jjl[4];
    float delta_rec = r6_largex0_delta_(vector_jjl, &ii);
    float theta_rec = r6_largex0_theta_(vector_jjl, &ii);
    float phi_rec = -r6_largex0_phi_(vector_jjl, &ii); // need to flip y, phi
    float y_rec = -r6_largex0_y00_(vector_jjl, &ii); // need to flip y, phi

    // Reset the vector and return it back to the caller
    pV5[0] = (double) x_or;
    pV5[1] = (double) theta_rec;
    pV5[2] = (double) y_rec;
    pV5[3] = (double) phi_rec;
    pV5[4] = (double) delta_rec;
}

void GDHTransLargeX0::ReconRightHRS(double *pV5)
{
    float vector_jjl[] = {pV5[0], pV5[1], pV5[2], pV5[3], pV5[4]};
    int ii = 5, jj = 1;

    // the detail of input vector in focus plane is
    // vector_jjl[0] = x_fp;
    // vector_jjl[1] = theta_fp;
    // vector_jjl[2] = y_fp;
    // vector_jjl[3] = phi_fp;
    // vector_jjl[4] = x_or;
    // the output is
    // vector_jjl[0] = x_or;
    // vector_jjl[1] = theta_rec;
    // vector_jjl[2] = y_rec;
    // vector_jjl[3] = phi_rec;
    // vector_jjl[4] = delta_rec;

    // Orthogonalize theta as JJL asks
    vector_jjl[1] = vector_jjl[1] - r6_largex0_txfit_(vector_jjl, &jj);
    float x_or = vector_jjl[4];
    float delta_rec = r6_largex0_delta_(vector_jjl, &ii);
    float theta_rec = r6_largex0_theta_(vector_jjl, &ii);
    float phi_rec = r6_largex0_phi_(vector_jjl, &ii);
    float y_rec = r6_largex0_y00_(vector_jjl, &ii);

    // Reset the vector and return it back to the caller
    pV5[0] = (double) x_or;
    pV5[1] = (double) theta_rec;
    pV5[2] = (double) y_rec;
    pV5[3] = (double) phi_rec;
    pV5[4] = (double) delta_rec;
}
