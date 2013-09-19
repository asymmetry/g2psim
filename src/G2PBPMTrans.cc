// -*- C++ -*-

/* namespace Orbit1, Orbit4, Orbit5, Orbit7, Orbit9
 * G2PBPM uses these functions to project BPMA and BPMB readouts to target center.
 *
 * The definition of orbit number:
 * Orbit        Field/T         Eb/GeV
 * 1            2.5             1.159
 * 4            2.5             1.706
 * 5            2.5             2.254
 * 7            5.0             2.254
 * 9            5.0             3.355
 */

// History:
//   Apr 2013, C. Gu, First public version.
//

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "G2PBPMTrans.hh"

////////////////////////////////////////////////////////////////////////
// Orbit 1
////////////////////////////////////////////////////////////////////////
namespace Orbit1 {

float target_x(float *x, int m) {
    //int ncoeff= 11;
    float avdat = 0.1774270E+01;
    float xmin[10] = {-0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18025E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10] = {0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17977E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10] = {0};
    float coeff[12] = {0.11933310E-03, -0.39622177E+02, 0.65532028E+02, 0.83443993E+00, 0.25200529E+01, -0.47372710E-01, 0.78216873E-01, -0.15225040E+01, -0.13823655E+01, -0.18708317E-01, 0.26466766E-01, 0.};
    int ientry = 0;
    int i;
    if (ientry == 0) {
        ientry = 1;
        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i])
                continue;
            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

    // normalize variables between -1 and +1
    float x1 = 1. + (x[0] - xmax[0]) * scale[0];
    float x2 = 1. + (x[1] - xmax[1]) * scale[1];
    float x3 = 1. + (x[2] - xmax[2]) * scale[2];
    float x4 = 1. + (x[3] - xmax[3]) * scale[3];

    // set up monomial functions
    float x11 = x1;
    float x21 = x2;
    float x22 = x21 * x2;
    float x31 = x3;
    float x41 = x4;

    // function
    float v_target_x = avdat + coeff[0] + coeff[1] * x11 + coeff[2] * x31 + coeff[3] * x11 * x21 + coeff[4] * x31 * x41 + coeff[5] * x21 + coeff[6] * x41 + coeff[7] * x21 * x31;
    v_target_x = v_target_x + coeff[8] * x11 * x41 + coeff[9] * x11 * x22 + coeff[10] * x11 * x21 * x41;
    return v_target_x;
}

float target_y(float *x, int m) {
    //int ncoeff= 20;
    float avdat = -0.7529873E+01;
    float xmin[10] = {-0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18025E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10] = {0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17977E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10] = {0};
    float coeff[21] = {-0.33809584E-01, -0.38952614E+02, 0.64402359E+02, -0.10177314E+01, 0.90655327E+00, -0.72384363E+00, -0.82957273E-03, -0.30791414E+00, 0.28529131E+00, 0.97626418E+00, -0.38929533E-01, 0.30307530E-03, -0.57311816E-03, 0.30283552E-01, 0.36454204E-01, -0.20252923E-01, -0.13093712E-01, 0.78416966E-01, -0.10658498E+00, 0.17995641E-01, 0.};
    int ientry = 0;
    int i;
    if (ientry == 0) {
        ientry = 1;
        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i])
                continue;
            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

    // normalize variables between -1 and +1
    float x1 = 1. + (x[0] - xmax[0]) * scale[0];
    float x2 = 1. + (x[1] - xmax[1]) * scale[1];
    float x3 = 1. + (x[2] - xmax[2]) * scale[2];
    float x4 = 1. + (x[3] - xmax[3]) * scale[3];

    // set up monomial functions
    float x11 = x1;
    float x12 = x11 * x1;
    float x13 = x12 * x1;
    float x14 = x13 * x1;
    float x21 = x2;
    float x22 = x21 * x2;
    float x31 = x3;
    float x32 = x31 * x3;
    float x41 = x4;
    float x42 = x41 * x4;


    float v_target_y = avdat + coeff[0] + coeff[1] * x21 + coeff[2] * x41 + coeff[3] * x21 * x41 + coeff[4] * x42 + coeff[5] * x32 + coeff[6] * x14 + coeff[7] * x12;
    v_target_y = v_target_y + coeff[8] * x22 + coeff[9] * x11 * x31 + coeff[10] * x31 + coeff[11] * x13 + coeff[12] * x12 * x31 + coeff[13] * x11 + coeff[14] * x21 * x32 + coeff[15] * x12 * x21 + coeff[16] * x22 * x41;
    v_target_y = v_target_y + coeff[17] * x11 * x31 * x41 + coeff[18] * x32 * x41 + coeff[19] * x21 * x42;
    return v_target_y;
}

float target_theta(float *x, int m) {
    //int ncoeff= 20;
    float avdat = 0.1120523E+00;
    float xmin[10] = {-0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18025E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10] = {0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17977E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10] = {0};
    float coeff[21] = {0.93548515E-04, -0.50340626E-01, 0.57996720E-01, -0.47768317E-02, 0.15887059E-02, 0.15522169E-01, -0.12489236E-01, -0.63411877E-02, 0.41331380E-03, -0.65351982E-03, 0.60341577E-02, 0.43789251E-03, -0.62632468E-03, -0.10483676E-03, -0.22410172E-03, 0.31096616E-03, -0.86420892E-04, 0.23729268E-04, 0.17314035E-03, -0.33568846E-04, 0.};
    int ientry = 0;
    int i;
    if (ientry == 0) {
        ientry = 1;
        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i])
                continue;
            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

    // normalize variables between -1 and +1
    float x1 = 1. + (x[0] - xmax[0]) * scale[0];
    float x2 = 1. + (x[1] - xmax[1]) * scale[1];
    float x3 = 1. + (x[2] - xmax[2]) * scale[2];
    float x4 = 1. + (x[3] - xmax[3]) * scale[3];

    // set up monomial functions
    float x11 = x1;
    float x12 = x11 * x1;
    float x21 = x2;
    float x22 = x21 * x2;
    float x23 = x22 * x2;
    float x31 = x3;
    float x32 = x31 * x3;
    float x41 = x4;
    float x42 = x41 * x4;
    float x43 = x42 * x4;

    // function
    float v_target_theta = avdat + coeff[0] + coeff[1] * x21 + coeff[2] * x41 + coeff[3] * x12 + coeff[4] * x22 + coeff[5] * x11 * x31 + coeff[6] * x32 + coeff[7] * x21 * x41;
    v_target_theta = v_target_theta + coeff[8] * x11 + coeff[9] * x31 + coeff[10] * x42 + coeff[11] * x21 * x32 * x41 + coeff[12] * x32 * x42 + coeff[13] * x43 + coeff[14] * x12 * x22 + coeff[15] * x12 * x21 * x41 + coeff[16] * x23;
    v_target_theta = v_target_theta + coeff[17] * x22 * x31 + coeff[18] * x22 * x41 + coeff[19] * x21 * x31 * x41;
    return v_target_theta;
}

float target_phi(float *x, int m) {
    //int ncoeff= 20;
    float avdat = 0.7490737E-03;
    float xmin[10] = {-0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18025E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10] = {0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17977E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10] = {0};
    float coeff[21] = {0.12031157E-06, 0.52352659E-02, 0.55798171E-02, -0.10130769E-01, -0.93395617E-02, 0.16881166E-01, -0.11295979E+01, -0.49663335E-01, 0.38291793E+01, 0.55737391E-01, -0.44885817E+01, 0.18031590E+01, 0.22394061E-01, 0.30083301E-01, 0.18411198E+01, -0.24231635E-01, -0.34712013E-01, -0.63662648E+01, 0.74261703E+01, -0.29101140E+01, 0.};
    int ientry = 0;
    int i;
    if (ientry == 0) {
        ientry = 1;
        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i])
                continue;
            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

    // normalize variables between -1 and +1
    float x1 = 1. + (x[0] - xmax[0]) * scale[0];
    float x2 = 1. + (x[1] - xmax[1]) * scale[1];
    float x3 = 1. + (x[2] - xmax[2]) * scale[2];
    float x4 = 1. + (x[3] - xmax[3]) * scale[3];

    // set up monomial functions
    float x11 = x1;
    float x12 = x11 * x1;
    float x13 = x12 * x1;
    float x14 = x13 * x1;
    float x21 = x2;
    float x22 = x21 * x2;
    float x23 = x22 * x2;
    float x24 = x23 * x2;
    float x31 = x3;
    float x32 = x31 * x3;
    float x33 = x32 * x3;
    float x34 = x33 * x3;
    float x41 = x4;

    // function
    float v_target_phi = avdat + coeff[0] + coeff[1] * x31 + coeff[2] * x11 * x21 + coeff[3] * x21 * x31 + coeff[4] * x11 * x41 + coeff[5] * x31 * x41 + coeff[6] * x13 + coeff[7] * x11 * x22;
    v_target_phi = v_target_phi + coeff[8] * x12 * x31 + coeff[9] * x22 * x31 + coeff[10] * x11 * x32 + coeff[11] * x33 + coeff[12] * x13 * x22 + coeff[13] * x11 * x24 + coeff[14] * x14 * x31 + coeff[15] * x12 * x22 * x31 + coeff[16] * x24 * x31;
    v_target_phi = v_target_phi + coeff[17] * x13 * x32 + coeff[18] * x12 * x33 + coeff[19] * x11 * x34;
    return v_target_phi;
}
}
////////////////////////////////////////////////////////////////////////
// End of Orbit 1
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// Orbit 4
////////////////////////////////////////////////////////////////////////
namespace Orbit4 {

float target_x(float *x, int m) {
    //int ncoeff= 11;
    float avdat = 0.7018812E+00;
    float xmin[10] = {-0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18018E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10] = {0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17984E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10] = {0};
    float coeff[12] = {-0.68359968E-05, -0.38941566E+02, 0.64789429E+02, -0.10212842E+01, 0.17022043E+01, 0.55410087E+00, -0.92462826E+00, 0.20685967E-01, -0.12420867E-01, -0.96249459E-02, 0.13531869E-01, 0.};
    int ientry = 0;
    int i;
    if (ientry == 0) {
        ientry = 1;
        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i])
                continue;
            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

    // normalize variables between -1 and +1
    float x1 = 1. + (x[0] - xmax[0]) * scale[0];
    float x2 = 1. + (x[1] - xmax[1]) * scale[1];
    float x3 = 1. + (x[2] - xmax[2]) * scale[2];
    float x4 = 1. + (x[3] - xmax[3]) * scale[3];

    // set up monomial functions
    float x11 = x1;
    float x21 = x2;
    float x22 = x21 * x2;
    float x31 = x3;
    float x41 = x4;

    // function
    float v_target_x = avdat + coeff[0] + coeff[1] * x11 + coeff[2] * x31 + coeff[3] * x21 * x31 + coeff[4] * x31 * x41 + coeff[5] * x11 * x21 + coeff[6] * x11 * x41 + coeff[7] * x41;
    v_target_x = v_target_x + coeff[8] * x21 + coeff[9] * x11 * x22 + coeff[10] * x11 * x21 * x41;
    return v_target_x;
}

float target_y(float *x, int m) {
    //int ncoeff= 20;
    float avdat = 0.2109603E+01;
    float xmin[10] = {-0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18018E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10] = {0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17984E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10] = {0};
    float coeff[21] = {-0.21812087E-01, -0.38611103E+02, 0.64224586E+02, 0.63122535E+00, 0.19239204E+00, -0.58644545E+00, -0.69882125E+00, 0.45095468E-02, -0.12007308E-01, 0.72982730E-02, -0.27014941E+00, 0.81429839E+00, 0.48026238E-02, -0.68415175E-02, 0.20527579E-01, -0.11970920E-01, 0.43771654E-01, -0.58269650E-01, -0.96899252E-02, 0.12425808E-01, 0.};
    int ientry = 0;
    int i;
    if (ientry == 0) {
        ientry = 1;
        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i])
                continue;
            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

    // normalize variables between -1 and +1
    float x1 = 1. + (x[0] - xmax[0]) * scale[0];
    float x2 = 1. + (x[1] - xmax[1]) * scale[1];
    float x3 = 1. + (x[2] - xmax[2]) * scale[2];
    float x4 = 1. + (x[3] - xmax[3]) * scale[3];

    // set up monomial functions
    float x11 = x1;
    float x12 = x11 * x1;
    float x13 = x12 * x1;
    float x14 = x13 * x1;
    float x21 = x2;
    float x22 = x21 * x2;
    float x31 = x3;
    float x32 = x31 * x3;
    float x41 = x4;
    float x42 = x41 * x4;
    float x43 = x42 * x4;

    // function
    float v_target_y = avdat + coeff[0] + coeff[1] * x21 + coeff[2] * x41 + coeff[3] * x42 + coeff[4] * x22 + coeff[5] * x32 + coeff[6] * x21 * x41 + coeff[7] * x14;
    v_target_y = v_target_y + coeff[8] * x13 * x31 + coeff[9] * x12 * x32 + coeff[10] * x12 + coeff[11] * x11 * x31 + coeff[12] * x11 + coeff[13] * x31 + coeff[14] * x21 * x32 + coeff[15] * x12 * x21 + coeff[16] * x11 * x31 * x41;
    v_target_y = v_target_y + coeff[17] * x32 * x41 + coeff[18] * x21 * x42 + coeff[19] * x43;
    return v_target_y;
}

float target_theta(float *x, int m) {
    //int ncoeff= 20;
    float avdat = 0.8201031E-01;
    float xmin[10] = {-0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18018E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10] = {0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17984E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10] = {0};
    float coeff[21] = {0.62649066E-04, -0.54196086E-01, 0.64258933E-01, -0.32147225E-02, 0.96104975E-03, 0.10427103E-01, -0.40180697E-02, -0.84015857E-02, 0.39361473E-02, 0.10668020E-03, -0.17160157E-03, -0.57252757E-04, 0.24173827E-04, -0.36565578E-03, 0.11950999E-02, 0.50279271E-03, 0.18978960E-03, 0.13606023E-04, -0.61446748E-03, -0.98281982E-03, 0.};
    int ientry = 0;
    int i;
    if (ientry == 0) {
        ientry = 1;
        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i])
                continue;
            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

    // normalize variables between -1 and +1
    float x1 = 1. + (x[0] - xmax[0]) * scale[0];
    float x2 = 1. + (x[1] - xmax[1]) * scale[1];
    float x3 = 1. + (x[2] - xmax[2]) * scale[2];
    float x4 = 1. + (x[3] - xmax[3]) * scale[3];

    // set up monomial functions
    float x11 = x1;
    float x12 = x11 * x1;
    float x21 = x2;
    float x22 = x21 * x2;
    float x31 = x3;
    float x32 = x31 * x3;
    float x33 = x32 * x3;
    float x34 = x33 * x3;
    float x41 = x4;
    float x42 = x41 * x4;

    // function
    float v_target_theta = avdat + coeff[0] + coeff[1] * x21 + coeff[2] * x41 + coeff[3] * x12 + coeff[4] * x22 + coeff[5] * x11 * x31 + coeff[6] * x21 * x41 + coeff[7] * x32;
    v_target_theta = v_target_theta + coeff[8] * x42 + coeff[9] * x11 + coeff[10] * x31 + coeff[11] * x32 * x41 + coeff[12] * x12 * x21 + coeff[13] * x22 * x32 + coeff[14] * x21 * x32 * x41 + coeff[15] * x12 * x42 + coeff[16] * x12 * x22;
    v_target_theta = v_target_theta + coeff[17] * x34 + coeff[18] * x12 * x21 * x41 + coeff[19] * x32 * x42;
    return v_target_theta;
}

float target_phi(float *x, int m) {
    //int ncoeff= 20;
    float avdat = 0.3971374E-03;
    float xmin[10] = {-0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18018E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10] = {0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17984E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10] = {0};
    float coeff[21] = {-0.40782746E-08, -0.57241894E-01, 0.69043346E-01, 0.37130348E-02, -0.67628366E-02, -0.62625464E-02, 0.11349871E-01, -0.80460253E-04, 0.13431822E-03, 0.63897925E-04, -0.92961913E-04, 0.10851461E-03, 0.13341644E-03, -0.12398859E-03, 0.17970247E-03, -0.37678954E-04, 0.97782700E-04, -0.20281586E-03, 0.28518474E-03, -0.40884584E-03, 0.};
    int ientry = 0;
    int i;
    if (ientry == 0) {
        ientry = 1;
        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i])
                continue;
            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

    // normalize variables between -1 and +1
    float x1 = 1. + (x[0] - xmax[0]) * scale[0];
    float x2 = 1. + (x[1] - xmax[1]) * scale[1];
    float x3 = 1. + (x[2] - xmax[2]) * scale[2];
    float x4 = 1. + (x[3] - xmax[3]) * scale[3];

    // set up monomial functions
    float x11 = x1;
    float x12 = x11 * x1;
    float x13 = x12 * x1;
    float x21 = x2;
    float x22 = x21 * x2;
    float x23 = x22 * x2;
    float x31 = x3;
    float x32 = x31 * x3;
    float x33 = x32 * x3;
    float x41 = x4;
    float x42 = x41 * x4;
    float x43 = x42 * x4;

    // function
    float v_target_phi = avdat + coeff[0] + coeff[1] * x11 + coeff[2] * x31 + coeff[3] * x11 * x21 + coeff[4] * x21 * x31 + coeff[5] * x11 * x41 + coeff[6] * x31 * x41 + coeff[7] * x21;
    v_target_phi = v_target_phi + coeff[8] * x41 + coeff[9] * x11 * x21 * x41 + coeff[10] * x11 * x42 + coeff[11] * x21 * x33 + coeff[12] * x11 * x43 + coeff[13] * x21 * x31 * x41 + coeff[14] * x31 * x42 + coeff[15] * x13 * x21 + coeff[16] * x11 * x23;
    v_target_phi = v_target_phi + coeff[17] * x11 * x22 * x41 + coeff[18] * x11 * x32 * x41 + coeff[19] * x33 * x41;
    return v_target_phi;
}
}
////////////////////////////////////////////////////////////////////////
// End of Orbit 4
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// Orbit 5
////////////////////////////////////////////////////////////////////////
namespace Orbit5 {

float target_x(float *x, int m) {
    //int ncoeff=  9;
    float avdat = 0.1820886E+01;
    float xmin[10] = {-0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18014E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10] = {0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17987E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10] = {0};
    float coeff[10] = {0.14974939E-03, -0.38482136E+02, 0.64213394E+02, 0.41052529E+00, -0.68654931E+00, -0.75923431E+00, 0.12685835E+01, 0.39473522E-01, -0.23631420E-01, 0.};
    int ientry = 0;
    int i;
    if (ientry == 0) {
        ientry = 1;
        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i])
                continue;
            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

    // normalize variables between -1 and +1
    float x1 = 1. + (x[0] - xmax[0]) * scale[0];
    float x2 = 1. + (x[1] - xmax[1]) * scale[1];
    float x3 = 1. + (x[2] - xmax[2]) * scale[2];
    float x4 = 1. + (x[3] - xmax[3]) * scale[3];

    // set up monomials functions
    float x11 = x1;
    float x21 = x2;
    float x31 = x3;
    float x41 = x4;

    // function
    float v_target_x = avdat + coeff[0] + coeff[1] * x11 + coeff[2] * x31 + coeff[3] * x11 * x21 + coeff[4] * x11 * x41 + coeff[5] * x21 * x31 + coeff[6] * x31 * x41 + coeff[7] * x41;
    v_target_x = v_target_x + coeff[8] * x21;
    return v_target_x;
}

float target_y(float *x, int m) {
    //int ncoeff= 16;
    float avdat = 0.1612694E+01;
    float xmin[10] = {-0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18014E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10] = {0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17987E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10] = {0};
    float coeff[17] = {-0.16093630E-01, -0.38296303E+02, 0.63893883E+02, 0.89169078E-01, -0.39361799E+00, -0.38331112E+00, 0.39724809E+00, -0.47729790E-03, 0.14334489E-03, -0.23719139E-01, -0.16308597E+00, 0.51654381E+00, 0.17781040E-01, 0.44934195E-02, -0.62706377E-02, 0.14249448E-02, 0.};
    int ientry = 0;
    int i;
    if (ientry == 0) {
        ientry = 1;
        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i])
                continue;
            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

    // normalize variables between -1 and +1
    float x1 = 1. + (x[0] - xmax[0]) * scale[0];
    float x2 = 1. + (x[1] - xmax[1]) * scale[1];
    float x3 = 1. + (x[2] - xmax[2]) * scale[2];
    float x4 = 1. + (x[3] - xmax[3]) * scale[3];

    // set up monomials functions
    float x11 = x1;
    float x12 = x11 * x1;
    float x13 = x12 * x1;
    float x14 = x13 * x1;
    float x21 = x2;
    float x22 = x21 * x2;
    float x31 = x3;
    float x32 = x31 * x3;
    float x33 = x32 * x3;
    float x34 = x33 * x3;
    float x41 = x4;
    float x42 = x41 * x4;
    float x43 = x42 * x4;

    // function
    float v_target_y = avdat + coeff[0] + coeff[1] * x21 + coeff[2] * x41 + coeff[3] * x22 + coeff[4] * x21 * x41 + coeff[5] * x32 + coeff[6] * x42 + coeff[7] * x14;
    v_target_y = v_target_y + coeff[8] * x34 + coeff[9] * x31 + coeff[10] * x12 + coeff[11] * x11 * x31 + coeff[12] * x11 + coeff[13] * x11 * x21 * x31 + coeff[14] * x21 * x32 + coeff[15] * x43;
    return v_target_y;
}

float target_theta(float *x, int m) {
    //int ncoeff= 20;
    float avdat = 0.6180951E-01;
    float xmin[10] = {-0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18014E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10] = {0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17987E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10] = {0};
    float coeff[21] = {0.43998152E-04, -0.55196326E-01, 0.65793619E-01, 0.74639604E-02, -0.60790703E-02, -0.27675012E-02, 0.61537535E-03, 0.21662402E-03, -0.34620706E-03, -0.22640435E-02, 0.28185088E-02, 0.54932359E-04, 0.25132557E-03, -0.60061220E-03, -0.23329834E-04, 0.10197179E-03, -0.15672242E-03, -0.73724154E-05, -0.11046945E-03, 0.40850369E-03, 0.};
    int ientry = 0;
    int i;
    if (ientry == 0) {
        ientry = 1;
        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i])
                continue;
            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

    // normalize variables between -1 and +1
    float x1 = 1. + (x[0] - xmax[0]) * scale[0];
    float x2 = 1. + (x[1] - xmax[1]) * scale[1];
    float x3 = 1. + (x[2] - xmax[2]) * scale[2];
    float x4 = 1. + (x[3] - xmax[3]) * scale[3];

    // set up monomials   functions
    float x11 = x1;
    float x12 = x11 * x1;
    float x21 = x2;
    float x22 = x21 * x2;
    float x31 = x3;
    float x32 = x31 * x3;
    float x41 = x4;
    float x42 = x41 * x4;

    // function
    float v_target_theta = avdat + coeff[0] + coeff[1] * x21 + coeff[2] * x41 + coeff[3] * x11 * x31 + coeff[4] * x32 + coeff[5] * x21 * x41 + coeff[6] * x22 + coeff[7] * x11;
    v_target_theta = v_target_theta + coeff[8] * x31 + coeff[9] * x12 + coeff[10] * x42 + coeff[11] * x21 * x32 + coeff[12] * x21 * x32 * x41 + coeff[13] * x32 * x42 + coeff[14] * x12 * x21 + coeff[15] * x11 * x31 * x41 + coeff[16] * x32 * x41;
    v_target_theta = v_target_theta + coeff[17] * x12 * x22 + coeff[18] * x12 * x21 * x41 + coeff[19] * x11 * x31 * x42;
    return v_target_theta;
}

float target_phi(float *x, int m) {
    //int ncoeff= 20;
    float avdat = 0.1156772E-02;
    float xmin[10] = {-0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18014E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10] = {0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17987E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10] = {0};
    float coeff[21] = {0.16362627E-06, 0.69261701E-02, -0.30966168E-02, 0.19902112E-02, 0.23972974E-02, -0.11254207E+01, -0.56836046E-01, 0.37969577E+01, 0.55131987E-01, -0.44339662E+01, 0.17743859E+01, 0.32073297E-01, 0.36992654E-01, 0.18466847E+01, -0.23962980E-01, -0.34303959E-01, -0.63650475E+01, 0.73991542E+01, -0.28891239E+01, -0.93376543E-02, 0.};
    int ientry = 0;
    int i;
    if (ientry == 0) {
        ientry = 1;
        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i])
                continue;
            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

    // normalize variables between -1 and +1
    float x1 = 1. + (x[0] - xmax[0]) * scale[0];
    float x2 = 1. + (x[1] - xmax[1]) * scale[1];
    float x3 = 1. + (x[2] - xmax[2]) * scale[2];
    float x4 = 1. + (x[3] - xmax[3]) * scale[3];

    // set up monomial functions
    float x11 = x1;
    float x12 = x11 * x1;
    float x13 = x12 * x1;
    float x14 = x13 * x1;
    float x21 = x2;
    float x22 = x21 * x2;
    float x23 = x22 * x2;
    float x24 = x23 * x2;
    float x31 = x3;
    float x32 = x31 * x3;
    float x33 = x32 * x3;
    float x34 = x33 * x3;
    float x41 = x4;

    // function
    float v_target_phi = avdat + coeff[0] + coeff[1] * x31 + coeff[2] * x11 * x21 + coeff[3] * x21 * x31 + coeff[4] * x11 * x41 + coeff[5] * x13 + coeff[6] * x11 * x22 + coeff[7] * x12 * x31;
    v_target_phi = v_target_phi + coeff[8] * x22 * x31 + coeff[9] * x11 * x32 + coeff[10] * x33 + coeff[11] * x13 * x22 + coeff[12] * x11 * x24 + coeff[13] * x14 * x31 + coeff[14] * x12 * x22 * x31 + coeff[15] * x24 * x31 + coeff[16] * x13 * x32;
    v_target_phi = v_target_phi + coeff[17] * x12 * x33 + coeff[18] * x11 * x34 + coeff[19] * x13 * x24;
    return v_target_phi;
}
}
////////////////////////////////////////////////////////////////////////
// End of Orbit 5
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// Orbit 7
////////////////////////////////////////////////////////////////////////
namespace Orbit7 {

float target_x(float *x, int m) {
    //int ncoeff= 11;
    float avdat = 0.1225998E+01;
    float xmin[10] = {-0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18029E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10] = {0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17972E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10] = {0};
    float coeff[12] = {0.26783593E-04, -0.38496861E+02, 0.64328911E+02, 0.81722581E+00, -0.13696008E+01, -0.15153360E+01, 0.25363481E+01, -0.29490298E-01, 0.49459714E-01, -0.10544261E-01, 0.14767365E-01, 0.};
    int ientry = 0;
    int i;
    if (ientry == 0) {
        ientry = 1;
        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i])
                continue;
            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

    // normalize variables between -1 and +1
    float x1 = 1. + (x[0] - xmax[0]) * scale[0];
    float x2 = 1. + (x[1] - xmax[1]) * scale[1];
    float x3 = 1. + (x[2] - xmax[2]) * scale[2];
    float x4 = 1. + (x[3] - xmax[3]) * scale[3];

    // set up monomial functions
    float x11 = x1;
    float x21 = x2;
    float x22 = x21 * x2;
    float x31 = x3;
    float x41 = x4;

    // function
    float v_target_x = avdat + coeff[0] + coeff[1] * x11 + coeff[2] * x31 + coeff[3] * x11 * x21 + coeff[4] * x11 * x41 + coeff[5] * x21 * x31 + coeff[6] * x31 * x41 + coeff[7] * x21;
    v_target_x = v_target_x + coeff[8] * x41 + coeff[9] * x11 * x22 + coeff[10] * x11 * x21 * x41;
    return v_target_x;
}

float target_y(float *x, int m) {
    //int ncoeff= 19;
    float avdat = 0.3077399E+01;
    float xmin[10] = {-0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18029E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10] = {0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17972E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10] = {0};
    float coeff[20] = {-0.35191756E-01, -0.38178669E+02, 0.63770287E+02, 0.28721189E+00, -0.10513144E+01, -0.84744918E+00, 0.95576644E+00, -0.95743773E-03, 0.29999032E-03, -0.32043450E-01, -0.37909397E+00, 0.11631805E+01, 0.22509055E-01, -0.12034442E-01, 0.20753311E-01, -0.26044870E-01, 0.11202721E-01, -0.75282492E-02, 0.10060057E-01, 0.};
    int ientry = 0;
    int i;
    if (ientry == 0) {
        ientry = 1;
        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i])
                continue;
            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

    // normalize variables between -1 and +1
    float x1 = 1. + (x[0] - xmax[0]) * scale[0];
    float x2 = 1. + (x[1] - xmax[1]) * scale[1];
    float x3 = 1. + (x[2] - xmax[2]) * scale[2];
    float x4 = 1. + (x[3] - xmax[3]) * scale[3];

    // set up monomial functions
    float x11 = x1;
    float x12 = x11 * x1;
    float x13 = x12 * x1;
    float x14 = x13 * x1;
    float x21 = x2;
    float x22 = x21 * x2;
    float x31 = x3;
    float x32 = x31 * x3;
    float x33 = x32 * x3;
    float x34 = x33 * x3;
    float x41 = x4;
    float x42 = x41 * x4;

    // function
    float v_target_y = avdat + coeff[0] + coeff[1] * x21 + coeff[2] * x41 + coeff[3] * x22 + coeff[4] * x21 * x41 + coeff[5] * x32 + coeff[6] * x42 + coeff[7] * x14;
    v_target_y = v_target_y + coeff[8] * x34 + coeff[9] * x31 + coeff[10] * x12 + coeff[11] * x11 * x31 + coeff[12] * x11 + coeff[13] * x11 * x21 * x31 + coeff[14] * x21 * x32 + coeff[15] * x32 * x41 + coeff[16] * x12 * x41;
    v_target_y = v_target_y + coeff[17] * x22 * x41 + coeff[18] * x21 * x42;
    return v_target_y;
}

float target_theta(float *x, int m) {
    //int ncoeff= 20;
    float avdat = 0.8200859E-03;
    float xmin[10] = {-0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18029E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10] = {0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17972E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10] = {0};
    float coeff[21] = {0.82136430E-04, -0.55444140E-01, 0.66205032E-01, -0.46876529E-02, 0.17395908E-02, 0.15340483E-01, -0.12436306E-01, -0.67721987E-02, 0.28879344E-03, 0.63901697E-02, -0.46958990E-03, 0.24596680E-03, -0.11652384E-03, 0.48341413E-03, -0.70989854E-03, -0.61864918E-03, 0.88314648E-03, 0.88900409E-03, -0.12663672E-02, 0.15426091E-04, 0.};
    int ientry = 0;
    int i;
    if (ientry == 0) {
        ientry = 1;
        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i])
                continue;
            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

    // normalize variables between -1 and +1
    float x1 = 1. + (x[0] - xmax[0]) * scale[0];
    float x2 = 1. + (x[1] - xmax[1]) * scale[1];
    float x3 = 1. + (x[2] - xmax[2]) * scale[2];
    float x4 = 1. + (x[3] - xmax[3]) * scale[3];

    // set up monomial functions
    float x11 = x1;
    float x12 = x11 * x1;
    float x21 = x2;
    float x22 = x21 * x2;
    float x31 = x3;
    float x32 = x31 * x3;
    float x41 = x4;
    float x42 = x41 * x4;
    float x43 = x42 * x4;

    // function
    float v_target_theta = avdat + coeff[0] + coeff[1] * x21 + coeff[2] * x41 + coeff[3] * x12 + coeff[4] * x22 + coeff[5] * x11 * x31 + coeff[6] * x32 + coeff[7] * x21 * x41;
    v_target_theta = v_target_theta + coeff[8] * x11 + coeff[9] * x42 + coeff[10] * x31 + coeff[11] * x21 * x32 + coeff[12] * x12 * x21 + coeff[13] * x11 * x31 * x41 + coeff[14] * x32 * x41 + coeff[15] * x11 * x21 * x31 * x41 + coeff[16] * x21 * x32 * x41;
    v_target_theta = v_target_theta + coeff[17] * x11 * x31 * x42 + coeff[18] * x32 * x42 + coeff[19] * x43;
    return v_target_theta;
}

float target_phi(float *x, int m) {
    //int ncoeff= 20;
    float avdat = 0.1190358E-02;
    float xmin[10] = {-0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18029E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10] = {0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17972E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10] = {0};
    float coeff[21] = {0.21864603E-07, -0.57748999E-01, 0.70092417E-01, 0.54349815E-02, -0.99676084E-02, -0.92149135E-02, 0.16823130E-01, -0.19285118E-03, 0.32475023E-03, -0.30417266E-03, 0.42958659E-03, -0.22883539E-03, 0.32734202E-03, -0.14062611E-03, 0.39565176E-03, -0.57434203E-03, 0.18514736E-03, 0.34364985E-03, -0.65984967E-03, 0.25425808E-03, 0.};
    int ientry = 0;
    int i;
    if (ientry == 0) {
        ientry = 1;
        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i])
                continue;
            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

    // normalize variables between -1 and +1
    float x1 = 1. + (x[0] - xmax[0]) * scale[0];
    float x2 = 1. + (x[1] - xmax[1]) * scale[1];
    float x3 = 1. + (x[2] - xmax[2]) * scale[2];
    float x4 = 1. + (x[3] - xmax[3]) * scale[3];

    // set up monomial functions
    float x11 = x1;
    float x12 = x11 * x1;
    float x21 = x2;
    float x22 = x21 * x2;
    float x31 = x3;
    float x32 = x31 * x3;
    float x33 = x32 * x3;
    float x41 = x4;
    float x42 = x41 * x4;

    // function
    float v_target_phi = avdat + coeff[0] + coeff[1] * x11 + coeff[2] * x31 + coeff[3] * x11 * x21 + coeff[4] * x21 * x31 + coeff[5] * x11 * x41 + coeff[6] * x31 * x41 + coeff[7] * x21;
    v_target_phi = v_target_phi + coeff[8] * x41 + coeff[9] * x22 * x31 + coeff[10] * x21 * x31 * x41 + coeff[11] * x12 * x21 * x31 + coeff[12] * x11 * x21 * x32 + coeff[13] * x22 * x31 * x41 + coeff[14] * x11 * x32 * x41 + coeff[15] * x33 * x41 + coeff[16] * x21 * x31 * x42;
    v_target_phi = v_target_phi + coeff[17] * x11 * x22 + coeff[18] * x11 * x21 * x41 + coeff[19] * x11 * x42;
    return v_target_phi;
}
}
////////////////////////////////////////////////////////////////////////
// End of Orbit 7
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// Orbit 9
////////////////////////////////////////////////////////////////////////
namespace Orbit9 {

float target_x(float *x, int m) {
    //int ncoeff=  9;
    float avdat = 0.2198843E+01;
    float xmin[10] = {-0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18020E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10] = {0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17981E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10] = {0};
    float coeff[10] = {0.29159064E-03, -0.38260330E+02, 0.63978661E+02, -0.10097252E+01, 0.16923758E+01, 0.54349035E+00, -0.34861766E-01, -0.91231066E+00, 0.58452848E-01, 0.};
    int ientry = 0;
    int i;
    if (ientry == 0) {
        ientry = 1;
        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i])
                continue;
            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

    // normalize variables between -1 and +1
    float x1 = 1. + (x[0] - xmax[0]) * scale[0];
    float x2 = 1. + (x[1] - xmax[1]) * scale[1];
    float x3 = 1. + (x[2] - xmax[2]) * scale[2];
    float x4 = 1. + (x[3] - xmax[3]) * scale[3];

    // set up monomial functions
    float x11 = x1;
    float x21 = x2;
    float x31 = x3;
    float x41 = x4;

    // function
    float v_target_x = avdat + coeff[0] + coeff[1] * x11 + coeff[2] * x31 + coeff[3] * x21 * x31 + coeff[4] * x31 * x41 + coeff[5] * x11 * x21 + coeff[6] * x21 + coeff[7] * x11 * x41;
    v_target_x = v_target_x + coeff[8] * x41;
    return v_target_x;
}

float target_y(float *x, int m) {
    //int ncoeff= 14;
    float avdat = 0.2720693E+01;
    float xmin[10] = {-0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18020E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10] = {0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17981E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10] = {0};
    float coeff[15] = {-0.22295034E-01, -0.38116352E+02, 0.63725323E+02, 0.61193955E+00, 0.17184985E+00, 0.26848450E-01, -0.38815312E-01, -0.60057771E+00, -0.65558451E+00, 0.23434667E-02, -0.67835855E-02, 0.39609550E-02, -0.27852637E+00, 0.83647746E+00, 0.};
    int ientry = 0;
    int i;
    if (ientry == 0) {
        ientry = 1;
        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i])
                continue;
            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

    // normalize variables between -1 and +1
    float x1 = 1. + (x[0] - xmax[0]) * scale[0];
    float x2 = 1. + (x[1] - xmax[1]) * scale[1];
    float x3 = 1. + (x[2] - xmax[2]) * scale[2];
    float x4 = 1. + (x[3] - xmax[3]) * scale[3];

    // set up monomial functions
    float x11 = x1;
    float x12 = x11 * x1;
    float x13 = x12 * x1;
    float x14 = x13 * x1;
    float x21 = x2;
    float x22 = x21 * x2;
    float x31 = x3;
    float x32 = x31 * x3;
    float x41 = x4;
    float x42 = x41 * x4;

    // function
    float v_target_y = avdat + coeff[0] + coeff[1] * x21 + coeff[2] * x41 + coeff[3] * x42 + coeff[4] * x22 + coeff[5] * x11 + coeff[6] * x31 + coeff[7] * x32;
    v_target_y = v_target_y + coeff[8] * x21 * x41 + coeff[9] * x14 + coeff[10] * x13 * x31 + coeff[11] * x12 * x32 + coeff[12] * x12 + coeff[13] * x11 * x31;
    return v_target_y;
}

float target_theta(float *x, int m) {
    //int ncoeff= 20;
    float avdat = 0.1141953E-02;
    float xmin[10] = {-0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18020E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10] = {0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17981E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10] = {0};
    float coeff[21] = {0.54978864E-05, 0.43244953E-02, 0.26680531E-02, -0.17784934E-02, -0.36398428E-02, 0.25598458E-02, -0.93145981E-01, -0.11179094E+01, 0.10739040E+00, 0.38970244E+01, -0.46267352E+01, 0.18546287E+01, 0.30168299E-01, 0.14267412E+01, -0.34811120E-01, -0.50512481E+01, 0.45458049E-01, 0.59980774E+01, -0.44600341E-01, -0.23854923E+01, 0.};
    int ientry = 0;
    int i;
    if (ientry == 0) {
        ientry = 1;
        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i])
                continue;
            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

    // normalize variables between -1 and +1
    float x1 = 1. + (x[0] - xmax[0]) * scale[0];
    float x2 = 1. + (x[1] - xmax[1]) * scale[1];
    float x3 = 1. + (x[2] - xmax[2]) * scale[2];
    float x4 = 1. + (x[3] - xmax[3]) * scale[3];

    // set up monomial functions
    float x11 = x1;
    float x12 = x11 * x1;
    float x13 = x12 * x1;
    float x14 = x13 * x1;
    float x21 = x2;
    float x22 = x21 * x2;
    float x23 = x22 * x2;
    float x24 = x23 * x2;
    float x31 = x3;
    float x41 = x4;
    float x42 = x41 * x4;
    float x43 = x42 * x4;

    // function
    float v_target_theta = avdat + coeff[0] + coeff[1] * x41 + coeff[2] * x12 + coeff[3] * x22 + coeff[4] * x11 * x31 + coeff[5] * x21 * x41 + coeff[6] * x12 * x21 + coeff[7] * x23;
    v_target_theta = v_target_theta + coeff[8] * x12 * x41 + coeff[9] * x22 * x41 + coeff[10] * x21 * x42 + coeff[11] * x43 + coeff[12] * x14 * x21 + coeff[13] * x12 * x23 + coeff[14] * x14 * x41 + coeff[15] * x12 * x22 * x41 + coeff[16] * x24 * x41;
    v_target_theta = v_target_theta + coeff[17] * x12 * x21 * x42 + coeff[18] * x23 * x42 + coeff[19] * x12 * x43;
    return v_target_theta;
}

float target_phi(float *x, int m) {
    //int ncoeff= 20;
    float avdat = 0.2305556E-02;
    float xmin[10] = {-0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18020E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float xmax[10] = {0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17981E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00};
    float scale[10] = {0};
    float coeff[21] = {0.31433748E-06, -0.57134841E-01, 0.68960786E-01, 0.36454808E-02, -0.66899294E-02, -0.61481497E-02, 0.11246344E-01, 0.38622235E-03, -0.22769075E-03, -0.56440800E-04, 0.79597266E-04, 0.10659136E-03, 0.11166979E-03, -0.67178926E-05, -0.37444988E-04, 0.11103924E-03, -0.25095875E-03, 0.28138197E-03, -0.40265819E-03, 0.68492249E-04, 0.};
    int ientry = 0;
    int i;
    if (ientry == 0) {
        ientry = 1;
        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i])
                continue;
            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

    // normalize variables between -1 and +1
    float x1 = 1. + (x[0] - xmax[0]) * scale[0];
    float x2 = 1. + (x[1] - xmax[1]) * scale[1];
    float x3 = 1. + (x[2] - xmax[2]) * scale[2];
    float x4 = 1. + (x[3] - xmax[3]) * scale[3];

    // set up monomial functions
    float x11 = x1;
    float x12 = x11 * x1;
    float x13 = x12 * x1;
    float x21 = x2;
    float x22 = x21 * x2;
    float x23 = x22 * x2;
    float x31 = x3;
    float x32 = x31 * x3;
    float x33 = x32 * x3;
    float x41 = x4;
    float x42 = x41 * x4;
    float x43 = x42 * x4;

    // function
    float v_target_phi = avdat + coeff[0] + coeff[1] * x11 + coeff[2] * x31 + coeff[3] * x11 * x21 + coeff[4] * x21 * x31 + coeff[5] * x11 * x41 + coeff[6] * x31 * x41 + coeff[7] * x41;
    v_target_phi = v_target_phi + coeff[8] * x21 + coeff[9] * x22 * x31 + coeff[10] * x21 * x31 * x41 + coeff[11] * x21 * x33 + coeff[12] * x11 * x43 + coeff[13] * x32 * x41 + coeff[14] * x13 * x21 + coeff[15] * x11 * x23 + coeff[16] * x11 * x22 * x41;
    v_target_phi = v_target_phi + coeff[17] * x11 * x32 * x41 + coeff[18] * x33 * x41 + coeff[19] * x21 * x31 * x42;
    return v_target_phi;
}
}
////////////////////////////////////////////////////////////////////////
// End of Orbit 9
////////////////////////////////////////////////////////////////////////
