// -*- C++ -*-

/* class GDHTransLargeX0
 * 6 deg with septum, for GDH experiment (E97110), large X0
 * By J.J. LeRose 10/05/2012
 */

// History:
//   Sep 2013, J. Zhang, First public version.
//

#include "GDHTransLargeX0.hh"

#include <iostream>
#include <iomanip>
#include <math.h>

using namespace std;

//#include "GlobalDebuger.hh"
//#define DEBUG_HRS_TRANSPORT 2

//include the header of fortran routines
#include "Bwd_R6_LargeX0_GDH.hh"
#include "Fwd_R6_LargeX0_GDH.hh"

const float m2cm = 100.0;

const double kDEG = 3.14159265358979323846 / 180.0;

GDHTransLargeX0::GDHTransLargeX0() {
    fModelAngle = 6.0 * kDEG;
}

GDHTransLargeX0::~GDHTransLargeX0() {
    // Nothing to do
}

//the detail of input vector is
//vector_jjl[0] = x_tr;
//vector_jjl[1] = theta_tr;
//vector_jjl[2] = y_tr;
//vector_jjl[3] = phi_tr;
//vector_jjl[4] = delta_tr;
//the output is
//vector_jjl[0] = x_fp;
//vector_jjl[1] = theta_fp;
//vector_jjl[2] = y_fp;
//vector_jjl[3] = phi_fp;
//vector_jjl[4] = delta_fp;  // delta is not change

//c     ep3: -0.1486233 < x < -0.08869672
//c          -0.110 < y < 0.110
//c     ep4: -0.1792231 < x < -0.1089169
//c          -0.110 < y < 0.110
//c     ep5: -0.2209211 < x < -0.1353789
//c          -0.110 < y < 0.110
//c     ep6: -0.2763536 < x < -0.1697464
//c          -0.110 < y < 0.110
//c     ep7: -0.3485396 < x < -0.2156404
//c          -0.110 < y < 0.110
//c     q1ex is a circle of radius 0.1492 m
//c     dent is a trapazoid:
//c                                   -5.22008 < x < -4.98099
//c             -(-0.192436784*x -0.192436784) < y < -0.192436784*x -0.192436784
//c     dext is also a trapazoid:
//c                                   -0.46188 < x < 0.46188
//c                   -(-0.01610808*x + 0.125) < y < -0.01610808*x + 0.125
//c     q3en is a circle of radius 0.3 m
//c     q3ex is a circle of radius 0.3 m
//c

//C Transport electron thru septum magnet, rectangle defined by (jjl)
//! includes reduced aperatures due to bore cooler,
//! assumes 0.8 cm thick side and 1.1 cm thick top and bottom
//Target to Septum entrance, -14.06cm>x>8.87cm, -9.9cm<y<9.9cm

/////////////////////////////////////////////////////////////////////////////
//2nd version of 6 degree seting, with larger X0, for E97110
/////////////////////////////////////////////////////////////////////////////

bool GDHTransLargeX0::TransRightHRS(double* pV5) {
    float vector_jjl[] = {pV5[0], pV5[1], pV5[2], pV5[3], pV5[4]};
    float x_test, y_test;
    int ii = 5;

    //By Jixie: since all collimators have the same y_min and y_max, I use
    //some variable to represent them for convenience
    float y_min = -9.9;
    float y_max = 9.9;

    x_test = x_sr6_largex0_ep3_(vector_jjl, &ii) * m2cm;
    y_test = y_sr6_largex0_ep3_(vector_jjl, &ii) * m2cm;

    if ((x_test<-14.06) || (x_test>-8.87) || (y_test < y_min) || (y_test > y_max))
        return false;

    //Target to 1/4 Septum, -17.12cm<x<-10.89cm, -9.9cm<y<9.9cm
    x_test = x_sr6_largex0_ep4_(vector_jjl, &ii) * m2cm;
    y_test = y_sr6_largex0_ep4_(vector_jjl, &ii) * m2cm;

    if ((x_test<-17.12) || (x_test>-10.89) || (y_test < y_min) || (y_test > y_max))
        return false;


    //Target to 1/2 Septum, -21.29cm<x<-13.54cm, -9.9cm<y<9.9cm
    x_test = x_sr6_largex0_ep5_(vector_jjl, &ii) * m2cm;
    y_test = y_sr6_largex0_ep5_(vector_jjl, &ii) * m2cm;
    if ((x_test<-21.29) || (x_test>-13.54) || (y_test < y_min) || (y_test > y_max))
        return false;

    //Target to 3/4 Septum, -26.84cm<x<-16.97cm, -9.9cm<y<9.9cm
    x_test = x_sr6_largex0_ep6_(vector_jjl, &ii) * m2cm;
    y_test = y_sr6_largex0_ep6_(vector_jjl, &ii) * m2cm;
    if ((x_test<-26.84) || (x_test>-16.97) || (y_test < y_min) || (y_test > y_max))
        return false;

    //Target to Septum exit, -34.05cm<x<-21.56cm, -9.9cm<y<9.9cm
    x_test = x_sr6_largex0_ep7_(vector_jjl, &ii) * m2cm;
    y_test = y_sr6_largex0_ep7_(vector_jjl, &ii) * m2cm;
    if ((x_test<-34.05) || (x_test>-21.56) || (y_test < y_min) || (y_test > y_max))
        return false;

    //Target to Q1 exit, circle of radius 14.92 cm
    x_test = x_sr6_largex0_q1ex_(vector_jjl, &ii) * m2cm;
    y_test = y_sr6_largex0_q1ex_(vector_jjl, &ii) * m2cm;
    x_test = x_test + 0.9;
    if ((x_test * x_test + y_test * y_test) > (14.92 * 14.92))
        return false;

    //Target to dipole entrance, trapezoid -522.0cm<x<-498.1cm  |y| < -0.1924*x-19.24
    x_test = x_sr6_largex0_dent_(vector_jjl, &ii) * m2cm;
    y_test = y_sr6_largex0_dent_(vector_jjl, &ii) * m2cm;
    if ((x_test<-522.0) || (x_test>-498.1) || fabs(y_test) > fabs(-0.1924 * x_test - 19.24))
        return false;

    //Target to dipole exit, trapezoid -46.19cm<x<46.19cm  |y| < -0.0161*x+12.5
    x_test = x_sr6_largex0_dext_(vector_jjl, &ii) * m2cm;
    y_test = y_sr6_largex0_dext_(vector_jjl, &ii) * m2cm;
    if ((x_test<-46.19) || (x_test > 46.19) || fabs(y_test) > fabs(-0.0161 * x_test + 12.5))
        return false;

    //Target to Q3 entrance, circle of radius 30.0 cm
    x_test = x_sr6_largex0_q3en_(vector_jjl, &ii) * m2cm;
    y_test = y_sr6_largex0_q3en_(vector_jjl, &ii) * m2cm;
    if ((x_test * x_test + y_test * y_test) > (30.0 * 30.0))
        return false;

    //Target to Q3 exit, circle of radius 30.0 cm  -> 28.0cm
    x_test = x_sr6_largex0_q3ex_(vector_jjl, &ii) * m2cm;
    y_test = y_sr6_largex0_q3ex_(vector_jjl, &ii) * m2cm;
    x_test = (x_test - 1.0) / (28.0);
    y_test = y_test / (30.0);
    if ((x_test * x_test + y_test * y_test) > 1.0)
        return false;

    /////////////////////////////////////////////////////////////
    /* If we reach this point, it means the test was succesful */
    float x_fp = x_sr6_largex0_fp_(vector_jjl, &ii);
    float theta_fp = t_sr6_largex0_fp_(vector_jjl, &ii);
    float y_fp = y_sr6_largex0_fp_(vector_jjl, &ii);
    float phi_fp = p_sr6_largex0_fp_(vector_jjl, &ii);

    //reset the vector and return it back to the caller
    pV5[0] = (double) x_fp;
    pV5[1] = (double) theta_fp;
    pV5[2] = (double) y_fp;
    pV5[3] = (double) phi_fp;
    //pV5[4] = (double)delta_fp;  // delta is not change

    return true;
}

bool GDHTransLargeX0::TransLeftHRS(double* pV5) {
    //in order to call right arm routines, need to flip y, phi
    float vector_jjl[] = {pV5[0], pV5[1], -pV5[2], -pV5[3], pV5[4]};
    bool bGoodParticle = TransRightHRS(pV5);
    if (bGoodParticle) {
        //need to flip y, phi again
        pV5[0] = (double) vector_jjl[0];
        pV5[1] = (double) vector_jjl[1];
        pV5[2] = (double) vector_jjl[2]*-1.0;
        pV5[3] = (double) vector_jjl[3]*-1.0;
        //pV5[4] = (double)vector_jjl[4];
    }
    return bGoodParticle;
}

/////////////////////////////////////////////////////////////////////////////
//6 degree seting, for E94110 and E97110
/////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
/* the detail of  JJLerose vector in focus plan is*/
//vector_jjl[0] = x_fp;
//vector_jjl[1] = theta_fp;
//vector_jjl[2] = y_fp;
//vector_jjl[3] = phi_fp;
//vector_jjl[4] = x_or;
//the output is
//vector_jjl[0] = x_or;
//vector_jjl[1] = theta_rec;
//vector_jjl[2] = y_rec;
//vector_jjl[3] = phi_rec;
//vector_jjl[4] = delta_rec;

/////////////////////////////////////////////////////////////////////////////////////

void GDHTransLargeX0::ReconRightHRS(double *pV5) {
    float vector_jjl[] = {pV5[0], pV5[1], pV5[2], pV5[3], pV5[4]};
    int ii = 5, jj = 1;

    /* Orthogonalize theta as JJL asks*/
    vector_jjl[1] = vector_jjl[1] - r6_largex0_txfit_(vector_jjl, &jj);
    float x_or = vector_jjl[4];
    float delta_rec = r6_largex0_delta_(vector_jjl, &ii);
    float theta_rec = r6_largex0_theta_(vector_jjl, &ii);
    float phi_rec = r6_largex0_phi_(vector_jjl, &ii);
    float y_rec = r6_largex0_y00_(vector_jjl, &ii);

    //reset the vector and return it back to the caller
    pV5[0] = (double) x_or;
    pV5[1] = (double) theta_rec;
    pV5[2] = (double) y_rec;
    pV5[3] = (double) phi_rec;
    pV5[4] = (double) delta_rec;
}

void GDHTransLargeX0::ReconLeftHRS(double *pV5) {
    //I am using the right arm routines to do the reconstruction, JJL did not provide
    //routines for the 2nd version of GDH experiment

    //in order to call right arm routines, need to flip y, phi
    float vector_jjl[] = {pV5[0], pV5[1], -pV5[2], -pV5[3], pV5[4]};
    int ii = 5, jj = 1;

    /* Orthogonalize theta as JJL asks*/
    vector_jjl[1] = vector_jjl[1] - r6_largex0_txfit_(vector_jjl, &jj);
    float x_or = vector_jjl[4];
    float delta_rec = r6_largex0_delta_(vector_jjl, &ii);
    float theta_rec = r6_largex0_theta_(vector_jjl, &ii);
    float phi_rec = -r6_largex0_phi_(vector_jjl, &ii); // need to flip y, phi
    float y_rec = -r6_largex0_y00_(vector_jjl, &ii); // need to flip y, phi

    //reset the vector and return it back to the caller
    pV5[0] = (double) x_or;
    pV5[1] = (double) theta_rec;
    pV5[2] = (double) y_rec;
    pV5[3] = (double) phi_rec;
    pV5[4] = (double) delta_rec;
}

/////////////////////////////////////////////////////////////////////////////////////

