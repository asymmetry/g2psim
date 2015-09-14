// -*- C++ -*-

/* class GDHTransSTD
 * 6 deg with septum, for GDH experiment (E97110), small X0
 * By J.J. LeRose 10/05/2012
 */

// History:
//   Sep 2013, J. Zhang, First public version.
//

#include <cmath>

#include "Bwd_L6_GDH.hh"
#include "Bwd_R6_GDH.hh"
#include "Fwd_L6_GDH.hh"
#include "Fwd_R6_GDH.hh"

#include "GDHTransSTD.hh"

const float m2cm = 100.0;
const double kDEG = 3.14159265358979323846 / 180.0;

GDHTransSTD::GDHTransSTD()
{
    fModelAngle = 6.0 * kDEG;
}

GDHTransSTD::~GDHTransSTD()
{
    // Nothing to do
}

int GDHTransSTD::TransLeftHRS(double *pV5, double *PlanePosX, double *PlanePosY)
{
    float vector_jjl[] = {float(pV5[0]), float(pV5[1]), float(pV5[2]), float(pV5[3]), float(pV5[4])};
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
    x_test = x_sl_ep3_(vector_jjl, &ii) * m2cm;
    y_test = y_sl_ep3_(vector_jjl, &ii) * m2cm;
    PlanePosX[5] = x_test;
    PlanePosY[5] = y_test;

    if ((x_test < -14.06) || (x_test > -8.87) || (y_test < y_min) || (y_test > y_max))
        return 5;

    // Target to 1/4 Septum, -17.12cm < x < -10.89cm, -9.9cm < y < 9.9cm
    x_test = x_sl_ep4_(vector_jjl, &ii) * m2cm;
    y_test = y_sl_ep4_(vector_jjl, &ii) * m2cm;
    PlanePosX[6] = x_test;
    PlanePosY[6] = y_test;

    if ((x_test < -17.12) || (x_test > -10.89) || (y_test < y_min) || (y_test > y_max))
        return 6;

    // Target to 1/2 Septum, -21.29cm < x < -13.54cm, -9.9cm < y < 9.9cm
    x_test = x_sl_ep5_(vector_jjl, &ii) * m2cm;
    y_test = y_sl_ep5_(vector_jjl, &ii) * m2cm;
    PlanePosX[7] = x_test;
    PlanePosY[7] = y_test;

    if ((x_test < -21.29) || (x_test > -13.54) || (y_test < y_min) || (y_test > y_max))
        return 7;

    // Target to 3/4 Septum, -26.84cm < x < -16.97cm, -9.9cm < y < 9.9cm
    x_test = x_sl_ep6_(vector_jjl, &ii) * m2cm;
    y_test = y_sl_ep6_(vector_jjl, &ii) * m2cm;
    PlanePosX[8] = x_test;
    PlanePosY[8] = y_test;

    if ((x_test < -26.84) || (x_test > -16.97) || (y_test < y_min) || (y_test > y_max))
        return 8;

    // Target to Septum exit, -34.05cm < x < -21.56cm, -9.9cm < y < 9.9cm
    x_test = x_sl_ep7_(vector_jjl, &ii) * m2cm;
    y_test = y_sl_ep7_(vector_jjl, &ii) * m2cm;
    PlanePosX[9] = x_test;
    PlanePosY[9] = y_test;

    if ((x_test < -34.05) || (x_test > -21.56) || (y_test < y_min) || (y_test > y_max))
        return 9;

    // Target to Q1 exit
    // circle of radius 14.92 cm
    x_test = x_sl_q1ex_(vector_jjl, &ii) * m2cm;
    y_test = y_sl_q1ex_(vector_jjl, &ii) * m2cm;
    PlanePosX[13] = x_test;
    PlanePosY[13] = y_test;
    x_test = x_test + 0.9;

    if ((x_test * x_test + y_test * y_test) > (14.92 * 14.92))
        return 13;

    // Target to dipole entrance
    // trapezoid, -522.0cm < x < -498.1cm, |y| < -0.1924 * x - 19.24
    x_test = x_sl_dent_(vector_jjl, &ii) * m2cm;
    y_test = y_sl_dent_(vector_jjl, &ii) * m2cm;
    PlanePosX[23] = x_test;
    PlanePosY[23] = y_test;

    if ((x_test < -522.0) || (x_test > -498.1) || fabs(y_test) > fabs(-0.1924 * x_test - 19.24))
        return 23;

    // Target to dipole exit
    // trapezoid, -46.19cm < x < 46.19cm, |y| < -0.0161 * x + 12.5
    x_test = x_sl_dext_(vector_jjl, &ii) * m2cm;
    y_test = y_sl_dext_(vector_jjl, &ii) * m2cm;
    PlanePosX[24] = x_test;
    PlanePosY[24] = y_test;

    if (fabs(x_test) > 46.19 || fabs(y_test) > fabs(-0.0161 * x_test + 12.5))
        return 24;

    // Target to Q3 entrance
    // circle of radius 30.0 cm
    x_test = x_sl_q3en_(vector_jjl, &ii) * m2cm;
    y_test = y_sl_q3en_(vector_jjl, &ii) * m2cm;
    PlanePosX[26] = x_test;
    PlanePosY[26] = y_test;

    if ((x_test * x_test + y_test * y_test) > (30.0 * 30.0))
        return 26;

    // Target to Q3 exit
    // circle of radius 30.0 cm -> 28.0cm
    x_test = x_sl_q3ex_(vector_jjl, &ii) * m2cm;
    y_test = y_sl_q3ex_(vector_jjl, &ii) * m2cm;
    PlanePosX[29] = x_test;
    PlanePosY[29] = y_test;
    x_test = (x_test - 1.0) / (28.0);
    y_test = y_test / (30.0);

    if ((x_test * x_test + y_test * y_test) > 1.0)
        return 29;

    // If we reach this point, it means the test was successful
    float x_fp = x_sl_fp_(vector_jjl, &ii);
    float theta_fp = t_sl_fp_(vector_jjl, &ii);
    float y_fp = y_sl_fp_(vector_jjl, &ii);
    float phi_fp = p_sl_fp_(vector_jjl, &ii);

    // Reset the vector and return it back to the caller
    pV5[0] = (double) x_fp;
    pV5[1] = (double) theta_fp;
    pV5[2] = (double) y_fp;
    pV5[3] = (double) phi_fp;
    //pV5[4] = (double)delta_fp; // delta is not change

    return 0;
}

int GDHTransSTD::TransRightHRS(double *pV5, double *PlanePosX, double *PlanePosY)
{
    float vector_jjl[] = {float(pV5[0]), float(pV5[1]), float(pV5[2]), float(pV5[3]), float(pV5[4])};
    float x_test, y_test;
    int ii = 5;

    float y_min = -9.9;
    float y_max = 9.9;

    // Target to Septum entrance, -14.06cm < x < -8.87cm, -9.9cm < y < 9.9cm
    x_test = x_sr_ep3_(vector_jjl, &ii) * m2cm;
    y_test = y_sr_ep3_(vector_jjl, &ii) * m2cm;
    PlanePosX[5] = x_test;
    PlanePosY[5] = y_test;

    if ((x_test < -14.06) || (x_test > -8.87) || (y_test < y_min) || (y_test > y_max))
        return 5;

    // Target to 1/4 Septum, -17.12cm < x < -10.89cm, -9.9cm < y < 9.9cm
    x_test = x_sr_ep4_(vector_jjl, &ii) * m2cm;
    y_test = y_sr_ep4_(vector_jjl, &ii) * m2cm;
    PlanePosX[6] = x_test;
    PlanePosY[6] = y_test;

    if ((x_test < -17.12) || (x_test > -10.89) || (y_test < y_min) || (y_test > y_max))
        return 6;

    // Target to 1/2 Septum, -21.29cm < x < -13.54cm, -9.9cm < y < 9.9cm
    x_test = x_sr_ep5_(vector_jjl, &ii) * m2cm;
    y_test = y_sr_ep5_(vector_jjl, &ii) * m2cm;
    PlanePosX[7] = x_test;
    PlanePosY[7] = y_test;

    if ((x_test < -21.29) || (x_test > -13.54) || (y_test < y_min) || (y_test > y_max))
        return 7;

    // Target to 3/4 Septum, -26.84cm < x < -16.97cm, -9.9cm < y < 9.9cm
    x_test = x_sr_ep6_(vector_jjl, &ii) * m2cm;
    y_test = y_sr_ep6_(vector_jjl, &ii) * m2cm;
    PlanePosX[8] = x_test;
    PlanePosY[8] = y_test;

    if ((x_test < -26.84) || (x_test > -16.97) || (y_test < y_min) || (y_test > y_max))
        return 8;

    // Target to Septum exit, -34.05cm < x < -21.56cm, -9.9cm < y < 9.9cm
    x_test = x_sr_ep7_(vector_jjl, &ii) * m2cm;
    y_test = y_sr_ep7_(vector_jjl, &ii) * m2cm;
    PlanePosX[9] = x_test;
    PlanePosY[9] = y_test;

    if ((x_test < -34.05) || (x_test > -21.56) || (y_test < y_min) || (y_test > y_max))
        return 9;

    // Target to Q1 exit
    // circle of radius 14.92 cm
    x_test = x_sr_q1ex_(vector_jjl, &ii) * m2cm;
    y_test = y_sr_q1ex_(vector_jjl, &ii) * m2cm;
    PlanePosX[13] = x_test;
    PlanePosY[13] = y_test;
    x_test = x_test + 0.9;

    if ((x_test * x_test + y_test * y_test) > (14.92 * 14.92))
        return 13;

    // Target to dipole entrance
    // trapezoid, -522.0cm < x < -498.1cm, |y| < -0.1924 * x - 19.24
    x_test = x_sr_dent_(vector_jjl, &ii) * m2cm;
    y_test = y_sr_dent_(vector_jjl, &ii) * m2cm;
    PlanePosX[23] = x_test;
    PlanePosY[23] = y_test;

    if ((x_test < -522.0) || (x_test > -498.1) || fabs(y_test) > fabs(-0.1924 * x_test - 19.24))
        return 23;

    // Target to dipole exit
    // trapezoid, -46.19cm < x < 46.19cm, |y| < -0.0161 * x + 12.5
    x_test = x_sr_dext_(vector_jjl, &ii) * m2cm;
    y_test = y_sr_dext_(vector_jjl, &ii) * m2cm;
    PlanePosX[24] = x_test;
    PlanePosY[24] = y_test;

    if (fabs(x_test) > 46.19 || fabs(y_test) > fabs(-0.0161 * x_test + 12.5))
        return 24;

    // Target to Q3 entrance
    // circle of radius 30.0 cm
    x_test = x_sr_q3en_(vector_jjl, &ii) * m2cm;
    y_test = y_sr_q3en_(vector_jjl, &ii) * m2cm;
    PlanePosX[26] = x_test;
    PlanePosY[26] = y_test;

    if ((x_test * x_test + y_test * y_test) > (30.0 * 30.0))
        return 26;

    // Target to Q3 exit
    // circle of radius 30.0 cm -> 28.0cm
    x_test = x_sr_q3ex_(vector_jjl, &ii) * m2cm;
    y_test = y_sr_q3ex_(vector_jjl, &ii) * m2cm;
    PlanePosX[29] = x_test;
    PlanePosY[29] = y_test;
    x_test = (x_test - 1.0) / (28.0);
    y_test = y_test / (30.0);

    if ((x_test * x_test + y_test * y_test) > 1.0)
        return 29;

    // If we reach this point, it means the test was successful
    float x_fp = x_sr_fp_(vector_jjl, &ii);
    float theta_fp = t_sr_fp_(vector_jjl, &ii);
    float y_fp = y_sr_fp_(vector_jjl, &ii);
    float phi_fp = p_sr_fp_(vector_jjl, &ii);

    // Reset the vector and return it back to the caller
    pV5[0] = (double) x_fp;
    pV5[1] = (double) theta_fp;
    pV5[2] = (double) y_fp;
    pV5[3] = (double) phi_fp;
    //pV5[4] = (double)delta_fp;  // delta is not change

    return 0;
}

void GDHTransSTD::ReconLeftHRS(double *pV5)
{
    // In order to call right arm routines, need to flip y, phi
    float vector_jjl[] = {float(pV5[0]), float(pV5[1]), float(-pV5[2]), float(-pV5[3]), float(pV5[4])};
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
    vector_jjl[1] = vector_jjl[1] - r6_txfit_(vector_jjl, &jj);
    float x_or_r = vector_jjl[4];
    float delta_rec_r = r6_delta_(vector_jjl, &ii);
    float theta_rec_r = r6_theta_(vector_jjl, &ii);
    float phi_rec_r = -r6_phi_(vector_jjl, &ii); // flip y
    float y_rec_r = -r6_y00_(vector_jjl, &ii); // flip phi

    // Reset the vector and return it back to the caller
    pV5[0] = (double) x_or_r;
    pV5[1] = (double) theta_rec_r;
    pV5[2] = (double) y_rec_r;
    pV5[3] = (double) phi_rec_r;
    pV5[4] = (double) delta_rec_r;

    // for (int i = 0; i < 5; i++) vector_jjl[i] = pV5[i];
    //
    // // Orthogonalize theta as JJL asks
    // // It turns out that sl6_txfit_ has bugs
    // vector_jjl[1] = vector_jjl[1] - sl6_txfit_(vector_jjl, &jj);
    // float x_or = vector_jjl[4];
    // float delta_rec = sl6_delta_(vector_jjl, &ii);
    // float theta_rec = sl6_theta_(vector_jjl, &ii);
    // float phi_rec = sl6_phi_(vector_jjl, &ii);
    // float y_rec = sl6_y00_(vector_jjl, &ii);
    //
    // // Reset the vector and return it back to the caller
    // pV5[0] = (double) x_or;
    // pV5[1] = (double) theta_rec;
    // pV5[2] = (double) y_rec;
    // pV5[3] = (double) phi_rec;
    // pV5[4] = (double) delta_rec;
}

void GDHTransSTD::ReconRightHRS(double *pV5)
{
    float vector_jjl[] = {float(pV5[0]), float(pV5[1]), float(pV5[2]), float(pV5[3]), float(pV5[4])};
    int ii = 5, jj = 1;

    // Orthogonalize theta as JJL asks
    vector_jjl[1] = vector_jjl[1] - r6_txfit_(vector_jjl, &jj);
    float x_or = vector_jjl[4];
    float delta_rec = r6_delta_(vector_jjl, &ii);
    float theta_rec = r6_theta_(vector_jjl, &ii);
    float phi_rec = r6_phi_(vector_jjl, &ii);
    float y_rec = r6_y00_(vector_jjl, &ii);

    // Reset the vector and return it back to the caller
    pV5[0] = (double) x_or;
    pV5[1] = (double) theta_rec;
    pV5[2] = (double) y_rec;
    pV5[3] = (double) phi_rec;
    pV5[4] = (double) delta_rec;
}

