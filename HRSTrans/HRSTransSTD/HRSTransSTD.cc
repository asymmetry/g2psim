// -*- C++ -*-

/* class G2PTrans484816
 * Standard HRS transport functions
 * By J.J. LeRose
 */

// History:
//   Sep 2013, J. Zhang, First public version.
//

#include <cmath>

#include "Fwd_l12p5_STD.h"
#include "Fwd_r12p5_STD.h"
#include "Bwd_l12p5_STD.h"
#include "Bwd_r12p5_STD.h"

#include "HRSTransSTD.hh"

const float m2cm = 100.0;
const double kDEG = 3.14159265358979323846 / 180.0;

HRSTransSTD::HRSTransSTD()
{
    fModelAngle = 12.5 * kDEG;
}

HRSTransSTD::~HRSTransSTD()
{
    // Nothing to do
}

int HRSTransSTD::TransLeftHRS(double* pV5)
{
    float vector_jjl[] = {pV5[0], pV5[1], pV5[2], pV5[3], pV5[4]};
    int iii = 5;
    int *ii = &iii;

    float x_test, y_test;

    // Target to Q1 exit
    // circle of radius 14.92 cm
    x_test = x_l12p5_q1ex_(vector_jjl, ii) * m2cm;
    y_test = y_l12p5_q1ex_(vector_jjl, ii) * m2cm;
    if ((x_test * x_test + y_test * y_test) > (14.92 * 14.92))
        return 5;

    // Target to dipole entrance
    // trapezoid 522.0cm > x > 498.1cm, |y| < -0.1924 * x - 19.24
    //x_test = x_l12p5_dent_(vector_jjl, ii) * m2cm;
    //y_test = y_l12p5_dent_(vector_jjl, ii) * m2cm;
    //if ((x_test > 522.0) || (x_test < 498.1) || fabs(y_test) > fabs(-0.1924 * x_test - 19.24))
    //    return 5;

    // Target to dipole exit
    // trapezoid -46.19cm < x < 46.19cm, |y| < -0.0161 * x + 12.5
    x_test = x_l12p5_dext_(vector_jjl, ii) * m2cm;
    y_test = y_l12p5_dext_(vector_jjl, ii) * m2cm;
    if (fabs(x_test) > 46.19 || fabs(y_test) > fabs(-0.0161 * x_test + 12.5))
        return 5;

    // Target to Q3 entrance
    // circle of radius 30.0 cm
    x_test = x_l12p5_q3en_(vector_jjl, ii) * m2cm;
    y_test = y_l12p5_q3en_(vector_jjl, ii) * m2cm;
    if ((x_test * x_test + y_test * y_test) > (30.0 * 30.0))
        return 5;

    // Target to Q3 exit
    // circle of radius 30.0 cm
    x_test = x_l12p5_q3ex_(vector_jjl, ii) * m2cm;
    y_test = y_l12p5_q3ex_(vector_jjl, ii) * m2cm;
    if ((x_test * x_test + y_test * y_test) > (30.0 * 30.0))
        return 5;

    // If we reach this point, it means the test was successful
    float x_fp = x_l12p5_fp_(vector_jjl, ii);
    float theta_fp = t_l12p5_fp_(vector_jjl, ii);
    float y_fp = y_l12p5_fp_(vector_jjl, ii);
    float phi_fp = p_l12p5_fp_(vector_jjl, ii);

    // Reset the vector and return it back to the caller
    pV5[0] = (double) x_fp;
    pV5[1] = (double) theta_fp;
    pV5[2] = (double) y_fp;
    pV5[3] = (double) phi_fp;
    //pV5[4] = (double)delta_fp; // delta is not change

    return 0;
}

int HRSTransSTD::TransRightHRS(double* pV5)
{
    float vector_jjl[] = {pV5[0], pV5[1], pV5[2], pV5[3], pV5[4]};
    int iii = 5;
    int *ii = &iii;

    float x_test, y_test;

    // Target to Q1 exit
    // circle of radius 14.92 cm
    x_test = x_r12p5_q1ex_(vector_jjl, ii) * m2cm;
    y_test = y_r12p5_q1ex_(vector_jjl, ii) * m2cm;
    if ((x_test * x_test + y_test * y_test) > (14.92 * 14.92))
        return 5;

    // Target to dipole entrance
    // trapezoid 522.0cm > x > 498.1cm, |y| < -0.1924 * x - 19.24
    //x_test = x_r12p5_dent_(vector_jjl, ii) * m2cm;
    //y_test = y_r12p5_dent_(vector_jjl, ii) * m2cm;
    //if ((x_test < -522.0) || (x_test > -498.1) || fabs(y_test) > fabs(-0.1924 * x_test - 19.24))
    //    return 5;

    // Target to dipole exit
    // trapezoid -46.19cm < x < 46.19cm, |y| < -0.0161 * x + 12.5
    x_test = x_r12p5_dext_(vector_jjl, ii) * m2cm;
    y_test = y_r12p5_dext_(vector_jjl, ii) * m2cm;
    if (fabs(x_test) > 46.19 || fabs(y_test) > fabs(-0.0161 * x_test + 12.5))
        return 5;

    // Target to Q3 entrance
    // circle of radius 30.0 cm

    x_test = x_r12p5_q3en_(vector_jjl, ii) * m2cm;
    y_test = y_r12p5_q3en_(vector_jjl, ii) * m2cm;
    if ((x_test * x_test + y_test * y_test) > (30.0 * 30.0))
        return 5;

    // Target to Q3 exit
    // circle of radius 30.0 cm
    x_test = x_r12p5_q3ex_(vector_jjl, ii) * m2cm;
    y_test = y_r12p5_q3ex_(vector_jjl, ii) * m2cm;
    if ((x_test * x_test + y_test * y_test) > (30.0 * 30.0))
        return 5;

    // If we reach this point, it means the test was successful
    float x_fp = x_r12p5_fp_(vector_jjl, ii);
    float theta_fp = t_r12p5_fp_(vector_jjl, ii);
    float y_fp = y_r12p5_fp_(vector_jjl, ii);
    float phi_fp = p_r12p5_fp_(vector_jjl, ii);

    // Reset the vector and return it back to the caller
    pV5[0] = (double) x_fp;
    pV5[1] = (double) theta_fp;
    pV5[2] = (double) y_fp;
    pV5[3] = (double) phi_fp;
    //pV5[4] = (double)delta_fp; // delta is not change

    return 0;
}

void HRSTransSTD::ReconLeftHRS(double* pV5)
{
    float vector_jjl[] = {pV5[0], pV5[1], pV5[2], pV5[3], pV5[4]};
    int iii = 5;
    int *ii = &iii;
    int jjj = 1;
    int *jj = &jjj;

    vector_jjl[1] = vector_jjl[1] - txfit_l12p5_(vector_jjl, jj);
    float x_or = vector_jjl[4];
    float delta_rec = delta_l12p5_(vector_jjl, ii);
    float theta_rec = theta_l12p5_(vector_jjl, ii);
    float phi_rec = phi_l12p5_(vector_jjl, ii);
    float y_rec = y00_l12p5_(vector_jjl, ii);

    // Reset the vector and return it back to the caller
    pV5[0] = (double) x_or;
    pV5[1] = (double) theta_rec;
    pV5[2] = (double) y_rec;
    pV5[3] = (double) phi_rec;
    pV5[4] = (double) delta_rec;
}

void HRSTransSTD::ReconRightHRS(double* pV5)
{
    float vector_jjl[] = {pV5[0], pV5[1], pV5[2], pV5[3], pV5[4]};
    int iii = 5;
    int *ii = &iii;
    int jjj = 1;
    int *jj = &jjj;

    // Orthogonalize theta as JJL asks
    vector_jjl[1] = vector_jjl[1] - txfit_r12p5_(vector_jjl, jj);
    float x_or = vector_jjl[4];
    float delta_rec = delta_r12p5_(vector_jjl, ii);
    float theta_rec = theta_r12p5_(vector_jjl, ii);
    float phi_rec = phi_r12p5_(vector_jjl, ii);
    float y_rec = y00_r12p5_(vector_jjl, ii);

    // Reset the vector and return it back to the caller
    pV5[0] = (double) x_or;
    pV5[1] = (double) theta_rec;
    pV5[2] = (double) y_rec;
    pV5[3] = (double) phi_rec;
    pV5[4] = (double) delta_rec;
}
