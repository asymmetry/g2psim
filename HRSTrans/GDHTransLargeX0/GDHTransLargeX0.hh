// -*- C++ -*-

/* class GDHTransLargeX0
 * 6 deg with septum, for GDH experiment (E97110), large X0
 * By J.J. LeRose 10/05/2012
 */

// History:
//   Sep 2013, J. Zhang, First public version.
//

#ifndef HRSTRANS_GDH_LARGEX0_H
#define HRSTRANS_GDH_LARGEX0_H

#include "HRSTransBase.hh"

class GDHTransLargeX0 : public HRSTransBase {
public:
    GDHTransLargeX0();
    ~GDHTransLargeX0();

    int TransLeftHRS(double* vector_jjl, double* PlanePosX, double* PlanePosY);
    int TransRightHRS(double* vector_jjl, double* PlanePosX, double* PlanePosY);
    void ReconLeftHRS(double* vector_jjl);
    void ReconRightHRS(double* vector_jjl);
};

#endif
