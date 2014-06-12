// -*- C++ -*-

/* class G2PTrans484816R00
 * 484816 septum with shims, 5.767 central ray, no target field, no raster
 * By M. Huang 1/7/2013
 */

// History:
//   Sep 2013, C. Gu, First public version.
//

#ifndef HRSTRANS_G2P_484816R00_H
#define HRSTRANS_G2P_484816R00_H

#include "HRSTransBase.hh"

class G2PTrans484816R00 : public HRSTransBase {
public:
    G2PTrans484816R00();
    ~G2PTrans484816R00();

    int TransLeftHRS(double* vector_jjl);
    bool TransRightHRS(double* vector_jjl);
    void ReconLeftHRS(double* vector_jjl);
    void ReconRightHRS(double* vector_jjl);

    void FPCorrLeft(const double* V5tg, double* V5fp);
    void FPCorrRight(const double* V5tg, double* V5fp);
};

#endif
