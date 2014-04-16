// -*- C++ -*-

/* class G2PTrans400016R15
 * 400016 septum with shims, 5.77 central ray, no target field, 15mm beam xy size
 * By M. Huang 4/15/2014
 */

// History:
//   Inheritated from Sep 2013, C. Gu, First public version.
//
//   Apr 15, 2014 M. Huang, add in this module

#ifndef HRSTRANS_G2P_400016R15_H
#define HRSTRANS_G2P_400016R15_H

#include "HRSTransBase.hh"

class G2PTrans400016R15 : public HRSTransBase {
public:
    G2PTrans400016R15();
    ~G2PTrans400016R15();

    bool TransLeftHRS(double* vector_jjl);
    bool TransRightHRS(double* vector_jjl);
    void ReconLeftHRS(double* vector_jjl);
    void ReconRightHRS(double* vector_jjl);

    void FPCorrLeft(const double* V5tg, double* V5fp);
    void FPCorrRight(const double* V5tg, double* V5fp);
};

#endif
