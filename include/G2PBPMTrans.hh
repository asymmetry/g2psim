// -*- C++ -*-

/* namespace Orbit1, Orbit2, Orbit3, Orbit4, Orbit5, Orbit7, Orbit8, Orbit9
 * G2PBPM uses these functions to project BPMA and BPMB readouts to target center.
 *
 * The definition of orbit number:
 * Orbit  Field/T  Angle/deg  Eb/GeV
 * 1      2.5      90         1.157
 * 2      5.0      6          1.157
 * 3      2.5      6          1.706
 * 4      2.5      90         1.706
 * 5      2.5      90         2.253
 * 7      5.0      90         2.253
 * 8      5.0      6          2.253
 * 9      5.0      90         3.355
 */

// History:
//   Apr 2013, C. Gu, First public version.
//   Oct 2015, C. Gu, Add orbit 2, 3, 8
//

#ifndef G2P_BPMTRANS_H
#define G2P_BPMTRANS_H

namespace Orbit1
{
double target_x(double *x, int m);
double target_y(double *x, int m);
double target_theta(double *x, int m);
double target_phi(double *x, int m);
}

namespace Orbit2
{
double target_x(double *x, int m);
double target_y(double *x, int m);
double target_theta(double *x, int m);
double target_phi(double *x, int m);
}

namespace Orbit3
{
double target_x(double *x, int m);
double target_y(double *x, int m);
double target_theta(double *x, int m);
double target_phi(double *x, int m);
}

namespace Orbit4
{
double target_x(double *x, int m);
double target_y(double *x, int m);
double target_theta(double *x, int m);
double target_phi(double *x, int m);
}

namespace Orbit5
{
double target_x(double *x, int m);
double target_y(double *x, int m);
double target_theta(double *x, int m);
double target_phi(double *x, int m);
}

namespace Orbit7
{
double target_x(double *x, int m);
double target_y(double *x, int m);
double target_theta(double *x, int m);
double target_phi(double *x, int m);
}

namespace Orbit8
{
double target_x(double *x, int m);
double target_y(double *x, int m);
double target_theta(double *x, int m);
double target_phi(double *x, int m);
}

namespace Orbit9
{
double target_x(double *x, int m);
double target_y(double *x, int m);
double target_theta(double *x, int m);
double target_phi(double *x, int m);
}

#endif
