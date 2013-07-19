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

namespace Orbit1 {
float target_x(float *x, int m);
float target_y(float *x, int m);
float target_theta(float *x, int m);
float target_phi(float *x, int m);
}

namespace Orbit4 {
float target_x(float *x, int m);
float target_y(float *x, int m);
float target_theta(float *x, int m);
float target_phi(float *x, int m);
}

namespace Orbit5 {
float target_x(float *x, int m);
float target_y(float *x, int m);
float target_theta(float *x, int m);
float target_phi(float *x, int m);
}

namespace Orbit7 {
float target_x(float *x, int m);
float target_y(float *x, int m);
float target_theta(float *x, int m);
float target_phi(float *x, int m);
}

namespace Orbit9 {
float target_x(float *x, int m);
float target_y(float *x, int m);
float target_theta(float *x, int m);
float target_phi(float *x, int m);
}

#endif
