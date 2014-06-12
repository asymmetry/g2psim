// -*- C++ -*-

/* class G2PTrans400016
 * 400016 septum with shims, 5.77 central ray, no target field, 15mm beam xy size
 * By M. Huang 4/15/2014
 */

// History:
//   Apr 2014, M. Huang, add in this module
//

#ifndef HRSTRANS_G2P_400016_H
#define HRSTRANS_G2P_400016_H

#include "HRSTransBase.hh"

class G2PTrans400016 : public HRSTransBase {
public:
    G2PTrans400016();
    ~G2PTrans400016();

    int TransLeftHRS(double* vector_jjl);
    bool TransRightHRS(double* vector_jjl);
    void ReconLeftHRS(double* vector_jjl);
    void ReconRightHRS(double* vector_jjl);

    void FPCorrLeft(const double* V5tg, double* V5fp);
    void FPCorrRight(const double* V5tg, double* V5fp);
};

#endif
