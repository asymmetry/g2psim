// -*- C++ -*-

/* class G2PTrans484816OLD
 * 484816 septum with shims, 5.65 central ray, no target field, 3cm raster
 * By J.J. LeRose 10/05/2012
 */

// History:
//   Sep 2013, C. Gu, First public version.
//

#ifndef HRSTRANS_G2P_484816OLD_H
#define HRSTRANS_G2P_484816OLD_H

#include "HRSTransBase.hh"

class G2PTrans484816OLD : public HRSTransBase {
public:
    G2PTrans484816OLD();
    ~G2PTrans484816OLD();

    int TransLeftHRS(double* vector_jjl);
    bool TransRightHRS(double* vector_jjl);
    void ReconLeftHRS(double* vector_jjl);
    void ReconRightHRS(double* vector_jjl);

    void FPCorrLeft(const double* V5tg, double* V5fp);
    void FPCorrRight(const double* V5tg, double* V5fp);
};

#endif
