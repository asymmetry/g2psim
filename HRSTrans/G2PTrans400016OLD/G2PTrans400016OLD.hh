// -*- C++ -*-

/* class G2PTrans400016OLD
 * 400016 septum with shims, 5.65 central ray, no target field, 3cm raster
 * By M. Huang 11/12/2012
 */

// History:
//   Sep 2013, C. Gu, First public version.
//

#ifndef HRSTRANS_G2P_400016OLD_H
#define HRSTRANS_G2P_400016OLD_H

#include "HRSTransBase.hh"

class G2PTrans400016OLD : public HRSTransBase {
public:
    G2PTrans400016OLD();
    ~G2PTrans400016OLD();

    int TransLeftHRS(double* vector_jjl);
    int TransRightHRS(double* vector_jjl);
    void ReconLeftHRS(double* vector_jjl);
    void ReconRightHRS(double* vector_jjl);
};

#endif
