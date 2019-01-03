// -*- C++ -*-

/* class G2PTrans403216
 * 403216 septum with shims, 5.77 central ray, no target field, 15mm beam xy size
 * By M. Huang 4/15/2014
 */

// History:
//   Apr 2014, M. Huang, add in this module
//

#ifndef HRSTRANS_G2P_403216_H
#define HRSTRANS_G2P_403216_H

#include "HRSTransBase.hh"

class G2PTrans403216 : public HRSTransBase {
public:
    G2PTrans403216();
    ~G2PTrans403216();

    int TransLeftHRS(double* vector_jjl, double* PlanePosX, double* PlanePosY);
    int TransRightHRS(double* vector_jjl, double* PlanePosX, double* PlanePosY);
    void ReconLeftHRS(double* vector_jjl);
    void ReconRightHRS(double* vector_jjl);
};

#endif
