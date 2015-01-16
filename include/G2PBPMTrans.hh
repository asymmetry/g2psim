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

#ifndef G2P_BPMTRANS_H
#define G2P_BPMTRANS_H

namespace Orbit1
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

namespace Orbit9
{
double target_x(double *x, int m);
double target_y(double *x, int m);
double target_theta(double *x, int m);
double target_phi(double *x, int m);
}

#endif
