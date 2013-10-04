// -*- C++ -*-

/* class G2PTrans484816R15
 * 484816 septum with shims, 5.785 central ray, no target field, 15mm beam xy size
 * By M. Huang 9/30/2013
 */

// History:
//   Inheritated from Sep 2013, C. Gu, First public version.
//
//   Sep 30, M. Huang, first add this module

#ifndef HRSTRANS_G2P_484816R15_H
#define HRSTRANS_G2P_484816R15_H

#include "HRSTransBase.hh"

class G2PTrans484816R15 : public HRSTransBase {
public:
    G2PTrans484816R15();
    ~G2PTrans484816R15();

    bool TransLeftHRS(double* vector_jjl);
    bool TransRightHRS(double* vector_jjl);
    void ReconLeftHRS(double* vector_jjl);
    void ReconRightHRS(double* vector_jjl);

    void FPCorrLeft(const double* V5tg, double* V5fp);
    void FPCorrRight(const double* V5tg, double* V5fp);
};

#endif
