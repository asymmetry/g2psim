// -*- C++ -*-

/* class G2PTrans484816
 * 484816 septum with shims, 5.77 central ray, no target field, 15mm beam xy size
 * By M. Huang 4/14/2014
 */

// History:
//   Sep 2013, M. Huang, first add this module
//   Apr 2014, M. Huang, complete the module with forward transport functions to multiple end-planes along the trajectory

#ifndef HRSTRANS_G2P_484816_H
#define HRSTRANS_G2P_484816_H

#include "HRSTransBase.hh"

class G2PTrans484816 : public HRSTransBase {
public:
    G2PTrans484816();
    ~G2PTrans484816();

    int TransLeftHRS(double* vector_jjl);
    bool TransRightHRS(double* vector_jjl);
    void ReconLeftHRS(double* vector_jjl);
    void ReconRightHRS(double* vector_jjl);

    void FPCorrLeft(const double* V5tg, double* V5fp);
    void FPCorrRight(const double* V5tg, double* V5fp);
};

#endif
