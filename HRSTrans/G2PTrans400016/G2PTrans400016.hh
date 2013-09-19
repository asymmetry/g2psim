// -*- C++ -*-

/* class G2PTrans400016
 * 400016 septum with shims, 5.65 central ray, no target field, 3cm raster
 * By M. Huang 11/12/2012
 */

// History:
//   Sep 2013, C. Gu, First public version.
//

#ifndef HRSTRANS_G2P_400016_H
#define HRSTRANS_G2P_400016_H

#include "HRSTransBase.hh"

class G2PTrans400016 : public HRSTransBase {
public:
    G2PTrans400016();
    ~G2PTrans400016();

    bool TransLeftHRS(double* vector_jjl);
    bool TransRightHRS(double* vector_jjl);
    void ReconLeftHRS(double* vector_jjl);
    void ReconRightHRS(double* vector_jjl);
};

#endif
