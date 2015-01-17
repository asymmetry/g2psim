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
namespace Orbit1
{

double target_x(double *x, int m)
{
    //int ncoeff= 17;
    double avdat =  0.1743476E+01;
    double xmin[10] = {
        -0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18023E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double xmax[10] = {
        0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17979E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double scale[10] = {0};
    double coeff[ 18] = {
        -0.17636981E-03, -0.40583023E+02, 0.66730476E+02, 0.88647437E+00,
        -0.14606659E+01, -0.16024895E+01, 0.26380181E+01, -0.47752146E-01,
        0.78487180E-01, 0.33229034E-01, -0.46933763E-01, -0.14068650E-01,
        -0.28934056E-01, 0.70947282E-01, -0.78048818E-02, 0.10551029E-01,
        0.11877671E-02,
        0.
    };
    int ientry = 0;
    int i;

    if (ientry == 0) {
        ientry = 1;

        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i]) continue;

            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

//  normalize variables between -1 and +1
    double x1 = 1. + (x[  0] - xmax[  0]) * scale[  0];
    double x2 = 1. + (x[  1] - xmax[  1]) * scale[  1];
    double x3 = 1. + (x[  2] - xmax[  2]) * scale[  2];
    double x4 = 1. + (x[  3] - xmax[  3]) * scale[  3];
//  set up monomials   functions
    double x11 = x1;
    double x12 = x11 * x1;
    double x21 = x2;
    double x22 = x21 * x2;
    double x31 = x3;
    double x32 = x31 * x3;
    double x41 = x4;
    double x42 = x41 * x4;

//                 function

    double v_target_x = avdat
                        + coeff[  0]
                        + coeff[  1] * x11
                        + coeff[  2]        * x31
                        + coeff[  3] * x11 * x21
                        + coeff[  4] * x11        * x41
                        + coeff[  5]    * x21 * x31
                        + coeff[  6]        * x31 * x41
                        + coeff[  7]    * x21
                        ;
    v_target_x = v_target_x
                 + coeff[  8]            * x41
                 + coeff[  9] * x11 * x22
                 + coeff[ 10] * x11 * x21    * x41
                 + coeff[ 11]    * x22 * x31
                 + coeff[ 12]    * x21 * x31 * x41
                 + coeff[ 13]        * x31 * x42
                 + coeff[ 14] * x12    * x31
                 + coeff[ 15] * x11    * x32
                 + coeff[ 16] * x11    * x31 * x42
                 ;
    return v_target_x;
}

double target_y(double *x, int m)
{
    //int ncoeff= 20;
    double avdat = -0.4025083E+01;
    double xmin[10] = {
        -0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18023E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double xmax[10] = {
        0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17979E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double scale[10] = {0};
    double coeff[ 21] = {
        -0.32111365E-01, -0.39886692E+02, 0.65563179E+02, -0.10011243E+01,
        0.90763301E+00, -0.80648327E+00, -0.76016056E-03, 0.11819611E-03,
        -0.36110654E+00, 0.27183130E+00, 0.11079873E+01, -0.26308103E-01,
        0.19471901E-01, 0.35428759E-01, -0.19442575E-01, -0.13427667E-01,
        0.75223878E-01, -0.10322651E+00, 0.18335944E-01, -0.55389776E-03,
        0.
    };
    int ientry = 0;
    int i;

    if (ientry == 0) {
        ientry = 1;

        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i]) continue;

            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

//  normalize variables between -1 and +1
    double x1 = 1. + (x[  0] - xmax[  0]) * scale[  0];
    double x2 = 1. + (x[  1] - xmax[  1]) * scale[  1];
    double x3 = 1. + (x[  2] - xmax[  2]) * scale[  2];
    double x4 = 1. + (x[  3] - xmax[  3]) * scale[  3];
//  set up monomials   functions
    double x11 = x1;
    double x12 = x11 * x1;
    double x13 = x12 * x1;
    double x14 = x13 * x1;
    double x21 = x2;
    double x22 = x21 * x2;
    double x31 = x3;
    double x32 = x31 * x3;
    double x33 = x32 * x3;
    double x34 = x33 * x3;
    double x41 = x4;
    double x42 = x41 * x4;

//                 function

    double v_target_y = avdat
                        + coeff[  0]
                        + coeff[  1]    * x21
                        + coeff[  2]            * x41
                        + coeff[  3]    * x21    * x41
                        + coeff[  4]            * x42
                        + coeff[  5]        * x32
                        + coeff[  6] * x14
                        + coeff[  7]        * x34
                        ;
    v_target_y = v_target_y
                 + coeff[  8] * x12
                 + coeff[  9]    * x22
                 + coeff[ 10] * x11    * x31
                 + coeff[ 11]        * x31
                 + coeff[ 12] * x11
                 + coeff[ 13]    * x21 * x32
                 + coeff[ 14] * x12 * x21
                 + coeff[ 15]    * x22    * x41
                 + coeff[ 16] * x11    * x31 * x41
                 ;
    v_target_y = v_target_y
                 + coeff[ 17]        * x32 * x41
                 + coeff[ 18]    * x21    * x42
                 + coeff[ 19] * x11 * x21
                 ;
    return v_target_y;
}

double target_theta(double *x, int m)
{
    //int ncoeff= 20;
    double avdat =  0.1117651E+00;
    double xmin[10] = {
        -0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18023E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double xmax[10] = {
        0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17979E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double scale[10] = {0};
    double coeff[ 21] = {
        0.99764933E-04, -0.50699815E-01, 0.58698010E-01, -0.50125537E-02,
        0.15764602E-02, 0.16160807E-01, -0.12924738E-01, -0.63843061E-02,
        0.39825638E-03, -0.63266262E-03, 0.61046197E-02, -0.63638022E-03,
        0.89668925E-03, 0.89368224E-03, -0.12596308E-02, 0.27908969E-04,
        -0.50039322E-04, 0.23908231E-04, -0.61413166E-05, -0.32098735E-04,
        0.
    };
    int ientry = 0;
    int i;

    if (ientry == 0) {
        ientry = 1;

        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i]) continue;

            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

//  normalize variables between -1 and +1
    double x1 = 1. + (x[  0] - xmax[  0]) * scale[  0];
    double x2 = 1. + (x[  1] - xmax[  1]) * scale[  1];
    double x3 = 1. + (x[  2] - xmax[  2]) * scale[  2];
    double x4 = 1. + (x[  3] - xmax[  3]) * scale[  3];
//  set up monomials   functions
    double x11 = x1;
    double x12 = x11 * x1;
    double x21 = x2;
    double x22 = x21 * x2;
    double x31 = x3;
    double x32 = x31 * x3;
    double x41 = x4;
    double x42 = x41 * x4;

//                 function

    double v_target_theta = avdat
                            + coeff[  0]
                            + coeff[  1]    * x21
                            + coeff[  2]            * x41
                            + coeff[  3] * x12
                            + coeff[  4]    * x22
                            + coeff[  5] * x11    * x31
                            + coeff[  6]        * x32
                            + coeff[  7]    * x21    * x41
                            ;
    v_target_theta = v_target_theta
                     + coeff[  8] * x11
                     + coeff[  9]        * x31
                     + coeff[ 10]            * x42
                     + coeff[ 11] * x11 * x21 * x31 * x41
                     + coeff[ 12]    * x21 * x32 * x41
                     + coeff[ 13] * x11    * x31 * x42
                     + coeff[ 14]        * x32 * x42
                     + coeff[ 15] * x11 * x21 * x31
                     + coeff[ 16]        * x32 * x41
                     ;
    v_target_theta = v_target_theta
                     + coeff[ 17] * x11 * x21    * x41
                     + coeff[ 18] * x11        * x42
                     + coeff[ 19]        * x31 * x42
                     ;
    return v_target_theta;
}

double target_phi(double *x, int m)
{
    //int ncoeff= 20;
    double avdat =  0.7711709E-03;
    double xmin[10] = {
        -0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18023E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double xmax[10] = {
        0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17979E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double scale[10] = {0};
    double coeff[ 21] = {
        0.35720810E-09, 0.52877953E-02, 0.58736377E-02, -0.10563425E-01,
        -0.97676162E-02, 0.17493853E-01, -0.11322616E+01, -0.50466456E-01,
        0.38551836E+01, 0.56646805E-01, -0.45394521E+01, 0.18317829E+01,
        0.22757124E-01, 0.30556185E-01, 0.18211429E+01, -0.24625771E-01,
        -0.35256881E-01, -0.63213725E+01, 0.74047604E+01, -0.29144576E+01,
        0.
    };
    int ientry = 0;
    int i;

    if (ientry == 0) {
        ientry = 1;

        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i]) continue;

            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

//  normalize variables between -1 and +1
    double x1 = 1. + (x[  0] - xmax[  0]) * scale[  0];
    double x2 = 1. + (x[  1] - xmax[  1]) * scale[  1];
    double x3 = 1. + (x[  2] - xmax[  2]) * scale[  2];
    double x4 = 1. + (x[  3] - xmax[  3]) * scale[  3];
//  set up monomials   functions
    double x11 = x1;
    double x12 = x11 * x1;
    double x13 = x12 * x1;
    double x14 = x13 * x1;
    double x21 = x2;
    double x22 = x21 * x2;
    double x23 = x22 * x2;
    double x24 = x23 * x2;
    double x31 = x3;
    double x32 = x31 * x3;
    double x33 = x32 * x3;
    double x34 = x33 * x3;
    double x41 = x4;

//                 function

    double v_target_phi  = avdat
                           + coeff[  0]
                           + coeff[  1]        * x31
                           + coeff[  2] * x11 * x21
                           + coeff[  3]    * x21 * x31
                           + coeff[  4] * x11        * x41
                           + coeff[  5]        * x31 * x41
                           + coeff[  6] * x13
                           + coeff[  7] * x11 * x22
                           ;
    v_target_phi  = v_target_phi
                    + coeff[  8] * x12    * x31
                    + coeff[  9]    * x22 * x31
                    + coeff[ 10] * x11    * x32
                    + coeff[ 11]        * x33
                    + coeff[ 12] * x13 * x22
                    + coeff[ 13] * x11 * x24
                    + coeff[ 14] * x14    * x31
                    + coeff[ 15] * x12 * x22 * x31
                    + coeff[ 16]    * x24 * x31
                    ;
    v_target_phi  = v_target_phi
                    + coeff[ 17] * x13    * x32
                    + coeff[ 18] * x12    * x33
                    + coeff[ 19] * x11    * x34
                    ;
    return v_target_phi  ;
}

}
////////////////////////////////////////////////////////////////////////
// End of Orbit 1
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// Orbit 4
////////////////////////////////////////////////////////////////////////
namespace Orbit4
{

double target_x(double *x, int m)
{
    //int ncoeff=  9;
    double avdat =  0.6676099E+00;
    double xmin[10] = {
        -0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18017E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double xmax[10] = {
        0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17985E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double scale[10] = {0};
    double coeff[ 10] = {
        -0.23387236E-04, -0.39861843E+02, 0.65917770E+02, -0.10740857E+01,
        0.17790698E+01, 0.58905059E+00, -0.97671878E+00, 0.20039359E-01,
        -0.12104545E-01,
        0.
    };
    int ientry = 0;
    int i;

    if (ientry == 0) {
        ientry = 1;

        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i]) continue;

            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

//  normalize variables between -1 and +1
    double x1 = 1. + (x[  0] - xmax[  0]) * scale[  0];
    double x2 = 1. + (x[  1] - xmax[  1]) * scale[  1];
    double x3 = 1. + (x[  2] - xmax[  2]) * scale[  2];
    double x4 = 1. + (x[  3] - xmax[  3]) * scale[  3];
//  set up monomials   functions
    double x11 = x1;
    double x21 = x2;
    double x31 = x3;
    double x41 = x4;

//                 function

    double v_target_x = avdat
                        + coeff[  0]
                        + coeff[  1] * x11
                        + coeff[  2]        * x31
                        + coeff[  3]    * x21 * x31
                        + coeff[  4]        * x31 * x41
                        + coeff[  5] * x11 * x21
                        + coeff[  6] * x11        * x41
                        + coeff[  7]            * x41
                        ;
    v_target_x = v_target_x
                 + coeff[  8]    * x21
                 ;
    return v_target_x;
}

double target_y(double *x, int m)
{
    //int ncoeff= 20;
    double avdat =  0.4571130E+01;
    double xmin[10] = {
        -0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18017E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double xmax[10] = {
        0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17985E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double scale[10] = {0};
    double coeff[ 21] = {
        -0.22821445E-01, -0.39518459E+02, 0.65333115E+02, 0.67580760E+00,
        0.21418124E+00, -0.64283490E+00, -0.76141554E+00, 0.25919394E-02,
        -0.75578876E-02, 0.46052765E-02, -0.30748791E+00, 0.90598822E+00,
        0.60408171E-02, -0.83330376E-02, 0.18341703E-01, -0.98282220E-02,
        0.37695020E-01, -0.52101631E-01, -0.89309392E-02, 0.11589214E-01,
        0.
    };
    int ientry = 0;
    int i;

    if (ientry == 0) {
        ientry = 1;

        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i]) continue;

            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

//  normalize variables between -1 and +1
    double x1 = 1. + (x[  0] - xmax[  0]) * scale[  0];
    double x2 = 1. + (x[  1] - xmax[  1]) * scale[  1];
    double x3 = 1. + (x[  2] - xmax[  2]) * scale[  2];
    double x4 = 1. + (x[  3] - xmax[  3]) * scale[  3];
//  set up monomials   functions
    double x11 = x1;
    double x12 = x11 * x1;
    double x13 = x12 * x1;
    double x14 = x13 * x1;
    double x21 = x2;
    double x22 = x21 * x2;
    double x31 = x3;
    double x32 = x31 * x3;
    double x41 = x4;
    double x42 = x41 * x4;
    double x43 = x42 * x4;

//                 function

    double v_target_y = avdat
                        + coeff[  0]
                        + coeff[  1]    * x21
                        + coeff[  2]            * x41
                        + coeff[  3]            * x42
                        + coeff[  4]    * x22
                        + coeff[  5]        * x32
                        + coeff[  6]    * x21    * x41
                        + coeff[  7] * x14
                        ;
    v_target_y = v_target_y
                 + coeff[  8] * x13    * x31
                 + coeff[  9] * x12    * x32
                 + coeff[ 10] * x12
                 + coeff[ 11] * x11    * x31
                 + coeff[ 12] * x11
                 + coeff[ 13]        * x31
                 + coeff[ 14]    * x21 * x32
                 + coeff[ 15] * x12 * x21
                 + coeff[ 16] * x11    * x31 * x41
                 ;
    v_target_y = v_target_y
                 + coeff[ 17]        * x32 * x41
                 + coeff[ 18]    * x21    * x42
                 + coeff[ 19]            * x43
                 ;
    return v_target_y;
}

double target_theta(double *x, int m)
{
    //int ncoeff= 20;
    double avdat =  0.8208974E-01;
    double xmin[10] = {
        -0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18017E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double xmax[10] = {
        0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17985E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double scale[10] = {0};
    double coeff[ 21] = {
        0.59646551E-04, -0.54393612E-01, 0.64627349E-01, -0.33002095E-02,
        0.11229175E-02, 0.10709785E-01, -0.44625797E-02, -0.86077312E-02,
        0.42447313E-02, 0.10586443E-03, -0.16888343E-03, -0.19640998E-03,
        0.28225780E-03, -0.46756095E-03, -0.95556796E-04, 0.55147818E-03,
        -0.60149119E-03, -0.16241729E-03, 0.54406968E-03, 0.66117777E-05,
        0.
    };
    int ientry = 0;
    int i;

    if (ientry == 0) {
        ientry = 1;

        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i]) continue;

            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

//  normalize variables between -1 and +1
    double x1 = 1. + (x[  0] - xmax[  0]) * scale[  0];
    double x2 = 1. + (x[  1] - xmax[  1]) * scale[  1];
    double x3 = 1. + (x[  2] - xmax[  2]) * scale[  2];
    double x4 = 1. + (x[  3] - xmax[  3]) * scale[  3];
//  set up monomials   functions
    double x11 = x1;
    double x12 = x11 * x1;
    double x21 = x2;
    double x22 = x21 * x2;
    double x31 = x3;
    double x32 = x31 * x3;
    double x41 = x4;
    double x42 = x41 * x4;

//                 function

    double v_target_theta = avdat
                            + coeff[  0]
                            + coeff[  1]    * x21
                            + coeff[  2]            * x41
                            + coeff[  3] * x12
                            + coeff[  4]    * x22
                            + coeff[  5] * x11    * x31
                            + coeff[  6]    * x21    * x41
                            + coeff[  7]        * x32
                            ;
    v_target_theta = v_target_theta
                     + coeff[  8]            * x42
                     + coeff[  9] * x11
                     + coeff[ 10]        * x31
                     + coeff[ 11] * x11 * x21 * x31
                     + coeff[ 12]    * x21 * x32
                     + coeff[ 13]        * x32 * x42
                     + coeff[ 14] * x12        * x41
                     + coeff[ 15] * x11    * x31 * x41
                     + coeff[ 16]        * x32 * x41
                     ;
    v_target_theta = v_target_theta
                     + coeff[ 17] * x12 * x22
                     + coeff[ 18] * x11 * x21 * x31 * x41
                     + coeff[ 19]    * x22    * x41
                     ;
    return v_target_theta;
}

double target_phi(double *x, int m)
{
    //int ncoeff= 20;
    double avdat =  0.3947516E-03;
    double xmin[10] = {
        -0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18017E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double xmax[10] = {
        0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17985E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double scale[10] = {0};
    double coeff[ 21] = {
        -0.22348891E-07, -0.57507142E-01, 0.69509491E-01, 0.39617335E-02,
        -0.71257828E-02, -0.65887882E-02, 0.11824809E-01, -0.78381629E-04,
        0.13004107E-03, 0.62197192E-04, -0.26274030E-04, -0.25839644E-03,
        0.34174198E-03, 0.11750825E-03, 0.11776925E-03, -0.26552597E-03,
        -0.11908164E-03, 0.69350947E-03, -0.71351830E-03, 0.72032897E-04,
        0.
    };
    int ientry = 0;
    int i;

    if (ientry == 0) {
        ientry = 1;

        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i]) continue;

            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

//  normalize variables between -1 and +1
    double x1 = 1. + (x[  0] - xmax[  0]) * scale[  0];
    double x2 = 1. + (x[  1] - xmax[  1]) * scale[  1];
    double x3 = 1. + (x[  2] - xmax[  2]) * scale[  2];
    double x4 = 1. + (x[  3] - xmax[  3]) * scale[  3];
//  set up monomials   functions
    double x11 = x1;
    double x12 = x11 * x1;
    double x21 = x2;
    double x22 = x21 * x2;
    double x23 = x22 * x2;
    double x31 = x3;
    double x32 = x31 * x3;
    double x33 = x32 * x3;
    double x41 = x4;
    double x42 = x41 * x4;
    double x43 = x42 * x4;

//                 function

    double v_target_phi  = avdat
                           + coeff[  0]
                           + coeff[  1] * x11
                           + coeff[  2]        * x31
                           + coeff[  3] * x11 * x21
                           + coeff[  4]    * x21 * x31
                           + coeff[  5] * x11        * x41
                           + coeff[  6]        * x31 * x41
                           + coeff[  7]    * x21
                           ;
    v_target_phi  = v_target_phi
                    + coeff[  8]            * x41
                    + coeff[  9]        * x31 * x42
                    + coeff[ 10] * x11 * x22
                    + coeff[ 11] * x11 * x21 * x32
                    + coeff[ 12]    * x21 * x33
                    + coeff[ 13] * x11        * x43
                    + coeff[ 14] * x11 * x23
                    + coeff[ 15] * x11 * x22    * x41
                    + coeff[ 16] * x12    * x31 * x41
                    ;
    v_target_phi  = v_target_phi
                    + coeff[ 17] * x11    * x32 * x41
                    + coeff[ 18]        * x33 * x41
                    + coeff[ 19]    * x21 * x31 * x42
                    ;
    return v_target_phi  ;
}

}
////////////////////////////////////////////////////////////////////////
// End of Orbit 4
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// Orbit 5
////////////////////////////////////////////////////////////////////////
namespace Orbit5
{

double target_x(double *x, int m)
{
    //int ncoeff=  9;
    double avdat =  0.1790319E+01;
    double xmin[10] = {
        -0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18013E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double xmax[10] = {
        0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17989E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double scale[10] = {0};
    double coeff[ 10] = {
        0.63137586E-05, -0.39382095E+02, 0.65306908E+02, 0.43441033E+00,
        -0.72221911E+00, -0.79527116E+00, 0.13208510E+01, 0.39553825E-01,
        -0.23821346E-01,
        0.
    };
    int ientry = 0;
    int i;

    if (ientry == 0) {
        ientry = 1;

        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i]) continue;

            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

//  normalize variables between -1 and +1
    double x1 = 1. + (x[  0] - xmax[  0]) * scale[  0];
    double x2 = 1. + (x[  1] - xmax[  1]) * scale[  1];
    double x3 = 1. + (x[  2] - xmax[  2]) * scale[  2];
    double x4 = 1. + (x[  3] - xmax[  3]) * scale[  3];
//  set up monomials   functions
    double x11 = x1;
    double x21 = x2;
    double x31 = x3;
    double x41 = x4;

//                 function

    double v_target_x = avdat
                        + coeff[  0]
                        + coeff[  1] * x11
                        + coeff[  2]        * x31
                        + coeff[  3] * x11 * x21
                        + coeff[  4] * x11        * x41
                        + coeff[  5]    * x21 * x31
                        + coeff[  6]        * x31 * x41
                        + coeff[  7]            * x41
                        ;
    v_target_x = v_target_x
                 + coeff[  8]    * x21
                 ;
    return v_target_x;
}

double target_y(double *x, int m)
{
    //int ncoeff= 16;
    double avdat =  0.3526662E+01;
    double xmin[10] = {
        -0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18013E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double xmax[10] = {
        0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17989E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double scale[10] = {0};
    double coeff[ 17] = {
        -0.16621113E-01, -0.39190487E+02, 0.64979134E+02, 0.15603992E+00,
        -0.56043983E+00, -0.46052605E+00, 0.50152624E+00, -0.43107473E-03,
        0.53009666E-04, -0.20796459E-01, -0.21615823E+00, 0.64443403E+00,
        0.15324904E-01, 0.45976071E-02, -0.63788076E-02, 0.14582445E-02,
        0.
    };
    int ientry = 0;
    int i;

    if (ientry == 0) {
        ientry = 1;

        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i]) continue;

            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

//  normalize variables between -1 and +1
    double x1 = 1. + (x[  0] - xmax[  0]) * scale[  0];
    double x2 = 1. + (x[  1] - xmax[  1]) * scale[  1];
    double x3 = 1. + (x[  2] - xmax[  2]) * scale[  2];
    double x4 = 1. + (x[  3] - xmax[  3]) * scale[  3];
//  set up monomials   functions
    double x11 = x1;
    double x12 = x11 * x1;
    double x13 = x12 * x1;
    double x14 = x13 * x1;
    double x21 = x2;
    double x22 = x21 * x2;
    double x31 = x3;
    double x32 = x31 * x3;
    double x33 = x32 * x3;
    double x34 = x33 * x3;
    double x41 = x4;
    double x42 = x41 * x4;
    double x43 = x42 * x4;

//                 function

    double v_target_y = avdat
                        + coeff[  0]
                        + coeff[  1]    * x21
                        + coeff[  2]            * x41
                        + coeff[  3]    * x22
                        + coeff[  4]    * x21    * x41
                        + coeff[  5]        * x32
                        + coeff[  6]            * x42
                        + coeff[  7] * x14
                        ;
    v_target_y = v_target_y
                 + coeff[  8]        * x34
                 + coeff[  9]        * x31
                 + coeff[ 10] * x12
                 + coeff[ 11] * x11    * x31
                 + coeff[ 12] * x11
                 + coeff[ 13] * x11 * x21 * x31
                 + coeff[ 14]    * x21 * x32
                 + coeff[ 15]            * x43
                 ;
    return v_target_y;
}

double target_theta(double *x, int m)
{
    //int ncoeff= 20;
    double avdat =  0.6202480E-01;
    double xmin[10] = {
        -0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18013E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double xmax[10] = {
        0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17989E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double scale[10] = {0};
    double coeff[ 21] = {
        0.45493704E-04, -0.55311322E-01, 0.66010721E-01, 0.80056107E-02,
        -0.64322841E-02, -0.31208822E-02, 0.74379134E-03, 0.21637636E-03,
        -0.34485047E-03, -0.24715587E-02, 0.30562545E-02, -0.64829459E-04,
        0.10978914E-03, -0.12935656E-03, -0.31915744E-03, 0.45420384E-03,
        0.45561741E-03, -0.64743758E-03, 0.54771735E-04, -0.16959303E-05,
        0.
    };
    int ientry = 0;
    int i;

    if (ientry == 0) {
        ientry = 1;

        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i]) continue;

            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

//  normalize variables between -1 and +1
    double x1 = 1. + (x[  0] - xmax[  0]) * scale[  0];
    double x2 = 1. + (x[  1] - xmax[  1]) * scale[  1];
    double x3 = 1. + (x[  2] - xmax[  2]) * scale[  2];
    double x4 = 1. + (x[  3] - xmax[  3]) * scale[  3];
//  set up monomials   functions
    double x11 = x1;
    double x12 = x11 * x1;
    double x21 = x2;
    double x22 = x21 * x2;
    double x31 = x3;
    double x32 = x31 * x3;
    double x41 = x4;
    double x42 = x41 * x4;

//                 function

    double v_target_theta = avdat
                            + coeff[  0]
                            + coeff[  1]    * x21
                            + coeff[  2]            * x41
                            + coeff[  3] * x11    * x31
                            + coeff[  4]        * x32
                            + coeff[  5]    * x21    * x41
                            + coeff[  6]    * x22
                            + coeff[  7] * x11
                            ;
    v_target_theta = v_target_theta
                     + coeff[  8]        * x31
                     + coeff[  9] * x12
                     + coeff[ 10]            * x42
                     + coeff[ 11] * x11 * x21 * x31
                     + coeff[ 12]    * x21 * x32
                     + coeff[ 13]        * x32 * x41
                     + coeff[ 14] * x11 * x21 * x31 * x41
                     + coeff[ 15]    * x21 * x32 * x41
                     + coeff[ 16] * x11    * x31 * x42
                     ;
    v_target_theta = v_target_theta
                     + coeff[ 17]        * x32 * x42
                     + coeff[ 18] * x12        * x41
                     + coeff[ 19] * x11 * x21
                     ;
    return v_target_theta;
}

double target_phi(double *x, int m)
{
    //int ncoeff= 20;
    double avdat =  0.1157775E-02;
    double xmin[10] = {
        -0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18013E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double xmax[10] = {
        0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17989E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double scale[10] = {0};
    double coeff[ 21] = {
        0.10655735E-07, -0.57073589E-01, 0.68776883E-01, 0.28866159E-02,
        -0.52297642E-02, -0.48527080E-02, 0.87516392E-02, -0.15369957E-03,
        0.25751023E-03, -0.24750456E-04, 0.63995205E-04, -0.12249475E-03,
        0.17248289E-03, 0.15560185E-03, 0.21245069E-03, -0.30277797E-03,
        -0.43908594E-03, 0.31329971E-03, -0.23957589E-05, -0.17305043E-04,
        0.
    };
    int ientry = 0;
    int i;

    if (ientry == 0) {
        ientry = 1;

        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i]) continue;

            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

//  normalize variables between -1 and +1
    double x1 = 1. + (x[  0] - xmax[  0]) * scale[  0];
    double x2 = 1. + (x[  1] - xmax[  1]) * scale[  1];
    double x3 = 1. + (x[  2] - xmax[  2]) * scale[  2];
    double x4 = 1. + (x[  3] - xmax[  3]) * scale[  3];
//  set up monomials   functions
    double x11 = x1;
    double x12 = x11 * x1;
    double x21 = x2;
    double x22 = x21 * x2;
    double x31 = x3;
    double x32 = x31 * x3;
    double x33 = x32 * x3;
    double x41 = x4;
    double x42 = x41 * x4;
    double x43 = x42 * x4;

//                 function

    double v_target_phi  = avdat
                           + coeff[  0]
                           + coeff[  1] * x11
                           + coeff[  2]        * x31
                           + coeff[  3] * x11 * x21
                           + coeff[  4]    * x21 * x31
                           + coeff[  5] * x11        * x41
                           + coeff[  6]        * x31 * x41
                           + coeff[  7]    * x21
                           ;
    v_target_phi  = v_target_phi
                    + coeff[  8]            * x41
                    + coeff[  9]    * x21 * x31 * x41
                    + coeff[ 10]        * x31 * x42
                    + coeff[ 11] * x12 * x21 * x31
                    + coeff[ 12] * x11 * x21 * x32
                    + coeff[ 13]    * x22 * x31 * x41
                    + coeff[ 14] * x11    * x32 * x41
                    + coeff[ 15]        * x33 * x41
                    + coeff[ 16]    * x21 * x31 * x42
                    ;
    v_target_phi  = v_target_phi
                    + coeff[ 17]        * x31 * x43
                    + coeff[ 18] * x11 * x21 * x31
                    + coeff[ 19] * x11 * x21    * x41
                    ;
    return v_target_phi  ;
}

}
////////////////////////////////////////////////////////////////////////
// End of Orbit 5
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// Orbit 7
////////////////////////////////////////////////////////////////////////
namespace Orbit7
{

double target_x(double *x, int m)
{
    //int ncoeff= 10;
    double avdat =  0.1191767E+01;
    double xmin[10] = {
        -0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18027E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double xmax[10] = {
        0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17974E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double scale[10] = {0};
    double coeff[ 11] = {
        0.33222139E-04, -0.39416534E+02, 0.65456841E+02, 0.86727107E+00,
        -0.14435229E+01, -0.15901225E+01, 0.26438158E+01, -0.29742345E-01,
        0.49378525E-01, 0.56314259E-02,
        0.
    };
    int ientry = 0;
    int i;

    if (ientry == 0) {
        ientry = 1;

        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i]) continue;

            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

//  normalize variables between -1 and +1
    double x1 = 1. + (x[  0] - xmax[  0]) * scale[  0];
    double x2 = 1. + (x[  1] - xmax[  1]) * scale[  1];
    double x3 = 1. + (x[  2] - xmax[  2]) * scale[  2];
    double x4 = 1. + (x[  3] - xmax[  3]) * scale[  3];
//  set up monomials   functions
    double x11 = x1;
    double x21 = x2;
    double x31 = x3;
    double x41 = x4;
    double x42 = x41 * x4;
    double x43 = x42 * x4;
    double x44 = x43 * x4;

//                 function

    double v_target_x = avdat
                        + coeff[  0]
                        + coeff[  1] * x11
                        + coeff[  2]        * x31
                        + coeff[  3] * x11 * x21
                        + coeff[  4] * x11        * x41
                        + coeff[  5]    * x21 * x31
                        + coeff[  6]        * x31 * x41
                        + coeff[  7]    * x21
                        ;
    v_target_x = v_target_x
                 + coeff[  8]            * x41
                 + coeff[  9]        * x31 * x44
                 ;
    return v_target_x;
}

double target_y(double *x, int m)
{
    //int ncoeff= 20;
    double avdat =  0.4782038E+01;
    double xmin[10] = {
        -0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18027E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double xmax[10] = {
        0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17974E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double scale[10] = {0};
    double coeff[ 21] = {
        -0.34589302E-01, -0.39092152E+02, 0.64887482E+02, 0.30870292E+00,
        -0.11154264E+01, -0.88635027E+00, 0.10022482E+01, -0.93715941E-03,
        0.25946356E-03, -0.34537494E-01, -0.40457165E+00, 0.12258588E+01,
        0.24603691E-01, -0.36185205E-01, 0.48397776E-01, -0.11011421E+00,
        -0.24892015E-01, -0.79524554E-02, 0.11463103E+00, 0.10580366E-01,
        0.
    };
    int ientry = 0;
    int i;

    if (ientry == 0) {
        ientry = 1;

        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i]) continue;

            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

//  normalize variables between -1 and +1
    double x1 = 1. + (x[  0] - xmax[  0]) * scale[  0];
    double x2 = 1. + (x[  1] - xmax[  1]) * scale[  1];
    double x3 = 1. + (x[  2] - xmax[  2]) * scale[  2];
    double x4 = 1. + (x[  3] - xmax[  3]) * scale[  3];
//  set up monomials   functions
    double x11 = x1;
    double x12 = x11 * x1;
    double x13 = x12 * x1;
    double x14 = x13 * x1;
    double x21 = x2;
    double x22 = x21 * x2;
    double x31 = x3;
    double x32 = x31 * x3;
    double x33 = x32 * x3;
    double x34 = x33 * x3;
    double x41 = x4;
    double x42 = x41 * x4;

//                 function

    double v_target_y = avdat
                        + coeff[  0]
                        + coeff[  1]    * x21
                        + coeff[  2]            * x41
                        + coeff[  3]    * x22
                        + coeff[  4]    * x21    * x41
                        + coeff[  5]        * x32
                        + coeff[  6]            * x42
                        + coeff[  7] * x14
                        ;
    v_target_y = v_target_y
                 + coeff[  8]        * x34
                 + coeff[  9]        * x31
                 + coeff[ 10] * x12
                 + coeff[ 11] * x11    * x31
                 + coeff[ 12] * x11
                 + coeff[ 13] * x11 * x21 * x31
                 + coeff[ 14]    * x21 * x32
                 + coeff[ 15]        * x32 * x41
                 + coeff[ 16] * x12        * x41
                 ;
    v_target_y = v_target_y
                 + coeff[ 17]    * x22    * x41
                 + coeff[ 18] * x11    * x31 * x41
                 + coeff[ 19]    * x21    * x42
                 ;
    return v_target_y;
}

double target_theta(double *x, int m)
{
    //int ncoeff= 20;
    double avdat =  0.1004010E-02;
    double xmin[10] = {
        -0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18027E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double xmax[10] = {
        0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17974E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double scale[10] = {0};
    double coeff[ 21] = {
        0.17238733E-03, 0.64394772E-02, -0.49919263E-02, -0.27214340E-02,
        0.16223796E-01, -0.13074216E-01, 0.39000700E-02, -0.75487830E-01,
        -0.75713021E+00, 0.84095903E-01, 0.25366268E+01, -0.29726195E+01,
        0.11964098E+01, 0.45413155E-01, 0.64331770E-01, -0.52402634E-01,
        -0.10982844E+00, 0.68441957E-01, 0.43932419E-01, -0.67155078E-01,
        0.
    };
    int ientry = 0;
    int i;

    if (ientry == 0) {
        ientry = 1;

        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i]) continue;

            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

//  normalize variables between -1 and +1
    double x1 = 1. + (x[  0] - xmax[  0]) * scale[  0];
    double x2 = 1. + (x[  1] - xmax[  1]) * scale[  1];
    double x3 = 1. + (x[  2] - xmax[  2]) * scale[  2];
    double x4 = 1. + (x[  3] - xmax[  3]) * scale[  3];
//  set up monomials   functions
    double x11 = x1;
    double x12 = x11 * x1;
    double x13 = x12 * x1;
    double x14 = x13 * x1;
    double x21 = x2;
    double x22 = x21 * x2;
    double x23 = x22 * x2;
    double x24 = x23 * x2;
    double x31 = x3;
    double x32 = x31 * x3;
    double x41 = x4;
    double x42 = x41 * x4;
    double x43 = x42 * x4;

//                 function

    double v_target_theta = avdat
                            + coeff[  0]
                            + coeff[  1]            * x41
                            + coeff[  2] * x12
                            + coeff[  3]    * x22
                            + coeff[  4] * x11    * x31
                            + coeff[  5]        * x32
                            + coeff[  6]    * x21    * x41
                            + coeff[  7] * x12 * x21
                            ;
    v_target_theta = v_target_theta
                     + coeff[  8]    * x23
                     + coeff[  9] * x12        * x41
                     + coeff[ 10]    * x22    * x41
                     + coeff[ 11]    * x21    * x42
                     + coeff[ 12]            * x43
                     + coeff[ 13] * x14 * x21
                     + coeff[ 14] * x12 * x23
                     + coeff[ 15] * x14        * x41
                     + coeff[ 16] * x12 * x22    * x41
                     ;
    v_target_theta = v_target_theta
                     + coeff[ 17]    * x24    * x41
                     + coeff[ 18] * x12 * x21    * x42
                     + coeff[ 19]    * x23    * x42
                     ;
    return v_target_theta;
}

double target_phi(double *x, int m)
{
    //int ncoeff= 20;
    double avdat =  0.1182557E-02;
    double xmin[10] = {
        -0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18027E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double xmax[10] = {
        0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17974E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double scale[10] = {0};
    double coeff[ 21] = {
        0.30169936E-07, -0.58040299E-01, 0.70596866E-01, 0.56541124E-02,
        -0.10317817E-01, -0.95822429E-02, 0.17370474E-01, -0.19255868E-03,
        0.32192009E-03, 0.12300612E-04, -0.94753799E-04, -0.57728703E-05,
        0.79150207E-03, 0.34391916E-04, -0.14829087E-03, 0.19474704E-03,
        0.15244866E-03, 0.24156134E-03, -0.76907064E-03, -0.33765318E-03,
        0.
    };
    int ientry = 0;
    int i;

    if (ientry == 0) {
        ientry = 1;

        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i]) continue;

            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

//  normalize variables between -1 and +1
    double x1 = 1. + (x[  0] - xmax[  0]) * scale[  0];
    double x2 = 1. + (x[  1] - xmax[  1]) * scale[  1];
    double x3 = 1. + (x[  2] - xmax[  2]) * scale[  2];
    double x4 = 1. + (x[  3] - xmax[  3]) * scale[  3];
//  set up monomials   functions
    double x11 = x1;
    double x12 = x11 * x1;
    double x13 = x12 * x1;
    double x21 = x2;
    double x22 = x21 * x2;
    double x23 = x22 * x2;
    double x31 = x3;
    double x32 = x31 * x3;
    double x33 = x32 * x3;
    double x41 = x4;
    double x42 = x41 * x4;
    double x43 = x42 * x4;
    double x44 = x43 * x4;

//                 function

    double v_target_phi  = avdat
                           + coeff[  0]
                           + coeff[  1] * x11
                           + coeff[  2]        * x31
                           + coeff[  3] * x11 * x21
                           + coeff[  4]    * x21 * x31
                           + coeff[  5] * x11        * x41
                           + coeff[  6]        * x31 * x41
                           + coeff[  7]    * x21
                           ;
    v_target_phi  = v_target_phi
                    + coeff[  8]            * x41
                    + coeff[  9]        * x31 * x44
                    + coeff[ 10]        * x33 * x41
                    + coeff[ 11] * x11 * x23    * x41
                    + coeff[ 12]        * x31 * x42
                    + coeff[ 13] * x13 * x21
                    + coeff[ 14]    * x22 * x31 * x41
                    + coeff[ 15]    * x21 * x31 * x42
                    + coeff[ 16]    * x22 * x31
                    ;
    v_target_phi  = v_target_phi
                    + coeff[ 17] * x11 * x21    * x41
                    + coeff[ 18]    * x21 * x31 * x41
                    + coeff[ 19] * x11        * x42
                    ;
    return v_target_phi  ;
}

}
////////////////////////////////////////////////////////////////////////
// End of Orbit 7
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// Orbit 9
////////////////////////////////////////////////////////////////////////
namespace Orbit9
{

double target_x(double *x, int m)
{
    //int ncoeff=  9;
    double avdat =  0.2168621E+01;
    double xmin[10] = {
        -0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18018E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double xmax[10] = {
        0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17983E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double scale[10] = {0};
    double coeff[ 10] = {
        0.68857335E-04, -0.39155697E+02, 0.65064911E+02, -0.10592775E+01,
        0.17628694E+01, 0.57712805E+00, -0.96148050E+00, -0.35512444E-01,
        0.58989104E-01,
        0.
    };
    int ientry = 0;
    int i;

    if (ientry == 0) {
        ientry = 1;

        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i]) continue;

            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

//  normalize variables between -1 and +1
    double x1 = 1. + (x[  0] - xmax[  0]) * scale[  0];
    double x2 = 1. + (x[  1] - xmax[  1]) * scale[  1];
    double x3 = 1. + (x[  2] - xmax[  2]) * scale[  2];
    double x4 = 1. + (x[  3] - xmax[  3]) * scale[  3];
//  set up monomials   functions
    double x11 = x1;
    double x21 = x2;
    double x31 = x3;
    double x41 = x4;

//                 function

    double v_target_x = avdat
                        + coeff[  0]
                        + coeff[  1] * x11
                        + coeff[  2]        * x31
                        + coeff[  3]    * x21 * x31
                        + coeff[  4]        * x31 * x41
                        + coeff[  5] * x11 * x21
                        + coeff[  6] * x11        * x41
                        + coeff[  7]    * x21
                        ;
    v_target_x = v_target_x
                 + coeff[  8]            * x41
                 ;
    return v_target_x;
}

double target_y(double *x, int m)
{
    //int ncoeff= 16;
    double avdat =  0.3957247E+01;
    double xmin[10] = {
        -0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18018E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double xmax[10] = {
        0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17983E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double scale[10] = {0};
    double coeff[ 17] = {
        -0.24037272E-01, -0.39005726E+02, 0.64802704E+02, 0.66658729E+00,
        0.20311919E+00, 0.27329957E-01, -0.39420936E-01, -0.58310479E+00,
        -0.73849553E+00, 0.14566478E-02, -0.45085372E-02, 0.25216362E-02,
        -0.26587474E+00, 0.80660069E+00, 0.46335957E-02, -0.64587002E-02,
        0.
    };
    int ientry = 0;
    int i;

    if (ientry == 0) {
        ientry = 1;

        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i]) continue;

            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

//  normalize variables between -1 and +1
    double x1 = 1. + (x[  0] - xmax[  0]) * scale[  0];
    double x2 = 1. + (x[  1] - xmax[  1]) * scale[  1];
    double x3 = 1. + (x[  2] - xmax[  2]) * scale[  2];
    double x4 = 1. + (x[  3] - xmax[  3]) * scale[  3];
//  set up monomials   functions
    double x11 = x1;
    double x12 = x11 * x1;
    double x13 = x12 * x1;
    double x14 = x13 * x1;
    double x21 = x2;
    double x22 = x21 * x2;
    double x31 = x3;
    double x32 = x31 * x3;
    double x41 = x4;
    double x42 = x41 * x4;

//                 function

    double v_target_y = avdat
                        + coeff[  0]
                        + coeff[  1]    * x21
                        + coeff[  2]            * x41
                        + coeff[  3]            * x42
                        + coeff[  4]    * x22
                        + coeff[  5] * x11
                        + coeff[  6]        * x31
                        + coeff[  7]        * x32
                        ;
    v_target_y = v_target_y
                 + coeff[  8]    * x21    * x41
                 + coeff[  9] * x14
                 + coeff[ 10] * x13    * x31
                 + coeff[ 11] * x12    * x32
                 + coeff[ 12] * x12
                 + coeff[ 13] * x11    * x31
                 + coeff[ 14] * x11 * x21 * x31
                 + coeff[ 15]    * x21 * x32
                 ;
    return v_target_y;
}

double target_theta(double *x, int m)
{
    //int ncoeff= 20;
    double avdat =  0.1411881E-02;
    double xmin[10] = {
        -0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18018E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double xmax[10] = {
        0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17983E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double scale[10] = {0};
    double coeff[ 21] = {
        0.47373142E-05, 0.43588234E-02, 0.27191939E-02, -0.18280419E-02,
        -0.37026971E-02, 0.26203997E-02, -0.93293555E-01, -0.11177913E+01,
        0.10765440E+00, 0.38986669E+01, -0.46316066E+01, 0.18578516E+01,
        0.30255878E-01, 0.14239165E+01, -0.34911972E-01, -0.50463147E+01,
        0.45635439E-01, 0.59986792E+01, -0.44780765E-01, -0.23884270E+01,
        0.
    };
    int ientry = 0;
    int i;

    if (ientry == 0) {
        ientry = 1;

        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i]) continue;

            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

//  normalize variables between -1 and +1
    double x1 = 1. + (x[  0] - xmax[  0]) * scale[  0];
    double x2 = 1. + (x[  1] - xmax[  1]) * scale[  1];
    double x3 = 1. + (x[  2] - xmax[  2]) * scale[  2];
    double x4 = 1. + (x[  3] - xmax[  3]) * scale[  3];
//  set up monomials   functions
    double x11 = x1;
    double x12 = x11 * x1;
    double x13 = x12 * x1;
    double x14 = x13 * x1;
    double x21 = x2;
    double x22 = x21 * x2;
    double x23 = x22 * x2;
    double x24 = x23 * x2;
    double x31 = x3;
    double x41 = x4;
    double x42 = x41 * x4;
    double x43 = x42 * x4;

//                 function

    double v_target_theta = avdat
                            + coeff[  0]
                            + coeff[  1]            * x41
                            + coeff[  2] * x12
                            + coeff[  3]    * x22
                            + coeff[  4] * x11    * x31
                            + coeff[  5]    * x21    * x41
                            + coeff[  6] * x12 * x21
                            + coeff[  7]    * x23
                            ;
    v_target_theta = v_target_theta
                     + coeff[  8] * x12        * x41
                     + coeff[  9]    * x22    * x41
                     + coeff[ 10]    * x21    * x42
                     + coeff[ 11]            * x43
                     + coeff[ 12] * x14 * x21
                     + coeff[ 13] * x12 * x23
                     + coeff[ 14] * x14        * x41
                     + coeff[ 15] * x12 * x22    * x41
                     + coeff[ 16]    * x24    * x41
                     ;
    v_target_theta = v_target_theta
                     + coeff[ 17] * x12 * x21    * x42
                     + coeff[ 18]    * x23    * x42
                     + coeff[ 19] * x12        * x43
                     ;
    return v_target_theta;
}

double target_phi(double *x, int m)
{
    //int ncoeff= 20;
    double avdat =  0.2306944E-02;
    double xmin[10] = {
        -0.15000E+02, -0.15000E+02, -0.18000E+02, -0.18018E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double xmax[10] = {
        0.15000E+02, 0.15000E+02, 0.18000E+02, 0.17983E+02, 0.00000E+00,
        0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00
    };
    double scale[10] = {0};
    double coeff[ 21] = {
        0.68881825E-07, -0.57274848E-01, 0.69202870E-01, 0.38690062E-02,
        -0.70175147E-02, -0.65140794E-02, 0.11760153E-01, 0.38778083E-03,
        -0.23018305E-03, 0.26180225E-04, 0.11273946E-03, 0.14055864E-03,
        -0.17463984E-04, -0.68420736E-05, 0.27533701E-04, -0.40552914E-04,
        0.10284865E-03, -0.21352051E-03, 0.30051859E-03, -0.42662484E-03,
        0.
    };
    int ientry = 0;
    int i;

    if (ientry == 0) {
        ientry = 1;

        for (i = 0; i < m; i++) {
            if (xmin[i] == xmax[i]) continue;

            scale[i] = 2. / (xmax[i] - xmin[i]);
        }
    }

//  normalize variables between -1 and +1
    double x1 = 1. + (x[  0] - xmax[  0]) * scale[  0];
    double x2 = 1. + (x[  1] - xmax[  1]) * scale[  1];
    double x3 = 1. + (x[  2] - xmax[  2]) * scale[  2];
    double x4 = 1. + (x[  3] - xmax[  3]) * scale[  3];
//  set up monomials   functions
    double x11 = x1;
    double x12 = x11 * x1;
    double x13 = x12 * x1;
    double x21 = x2;
    double x22 = x21 * x2;
    double x23 = x22 * x2;
    double x24 = x23 * x2;
    double x31 = x3;
    double x32 = x31 * x3;
    double x33 = x32 * x3;
    double x41 = x4;
    double x42 = x41 * x4;
    double x43 = x42 * x4;
    double x44 = x43 * x4;

//                 function

    double v_target_phi  = avdat
                           + coeff[  0]
                           + coeff[  1] * x11
                           + coeff[  2]        * x31
                           + coeff[  3] * x11 * x21
                           + coeff[  4]    * x21 * x31
                           + coeff[  5] * x11        * x41
                           + coeff[  6]        * x31 * x41
                           + coeff[  7]            * x41
                           ;
    v_target_phi  = v_target_phi
                    + coeff[  8]    * x21
                    + coeff[  9]        * x31 * x44
                    + coeff[ 10]    * x21 * x33
                    + coeff[ 11] * x11        * x43
                    + coeff[ 12] * x11 * x24
                    + coeff[ 13]        * x32 * x41
                    + coeff[ 14]        * x31 * x42
                    + coeff[ 15] * x13 * x21
                    + coeff[ 16] * x11 * x23
                    ;
    v_target_phi  = v_target_phi
                    + coeff[ 17] * x11 * x22    * x41
                    + coeff[ 18] * x11    * x32 * x41
                    + coeff[ 19]        * x33 * x41
                    ;
    return v_target_phi  ;
}

}
////////////////////////////////////////////////////////////////////////
// End of Orbit 9
////////////////////////////////////////////////////////////////////////
