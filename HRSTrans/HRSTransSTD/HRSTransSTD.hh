// -*- C++ -*-

/* class G2PTrans484816
 * Standard HRS transport functions
 * By J.J. LeRose
 */

// History:
//   Sep 2013, J. Zhang, First public version.
//

#ifndef HRSTRANS_STD_H
#define HRSTRANS_STD_H

#include "HRSTransBase.hh"

class HRSTransSTD : public HRSTransBase {
public:
    HRSTransSTD();
    ~HRSTransSTD();

    bool TransLeftHRS(double* vector_jjl);
    bool TransRightHRS(double* vector_jjl);
    void ReconLeftHRS(double* vector_jjl);
    void ReconRightHRS(double* vector_jjl);
};

#endif
