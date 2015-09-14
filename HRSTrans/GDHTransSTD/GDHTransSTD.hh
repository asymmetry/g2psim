// -*- C++ -*-

/* class GDHTransSTD
 * 6 deg with septum, for GDH experiment (E97110), small X0
 * By J.J. LeRose 10/05/2012
 */

// History:
//   Sep 2013, J. Zhang, First public version.
//

#ifndef HRSTRANS_GDH_STD_H
#define HRSTRANS_GDH_STD_H

#include "HRSTransBase.hh"

class GDHTransSTD : public HRSTransBase
{
public:
    GDHTransSTD();
    ~GDHTransSTD();

    int TransLeftHRS(double *vector_jjl, double *PlanePosX, double *PlanePosY);
    int TransRightHRS(double *vector_jjl, double *PlanePosX, double *PlanePosY);
    void ReconLeftHRS(double *vector_jjl);
    void ReconRightHRS(double *vector_jjl);
};

#endif
