#ifndef HRS_TRANSPORT_H
#define HRS_TRANSPORT_H

#include <map>

#include "TROOT.h"
#include "TObject.h"

#include "G2PTrans.hh"

using namespace std;

class HRSTransport : public TObject
{
public:
    HRSTransport();
    HRSTransport(const char *name);
    HRSTransport(int setting);
    ~HRSTransport();
    
    virtual void RegisterModel();

    void SetArm(bool isleftarm) { bIsLeftArm = isleftarm; }

    int GetModelIndex() { return iModelIndex; }
    
///////////////////////////////////////////////////////////////////////////
// Definition of variables
// forward:
// V5_tg = {x_tg, theta_tg, y_tg, phi_tg, delta@tg};
// V5_fp = {x_fg, theta_fp, y_fp, phi_fp, delta@tg};
// delta does not change
// backward:
// V5_fp = {x_fg, theta_fp, y_fp, phi_fp, x_tg};
// V5_tg = {x_tg, theta_tg, y_tg, phi_tg, delta@tg};
// x_tg does not change
//
// Transport particles through HRS using SNAKE model
// Use iModelIndex to identify which SNAKE model to be used
// 1: 484816 with shim, 5.65 deg, 3 cm raster, by JJL 
// 2: 403216 with shim, 5.65 deg, SNAKE Model not ready yet 
// 3: 400016 with shim, 5.65 deg, 3 cm raster, by Min
// Index > 10 means test
// 11: 484816 with shim, 5.65 deg, Wrong Bx, 2 cm raster, by Min
// 12: 484816 with shim, 5.76 deg, no raster, by Min
// 13: 484816 no shim, 5.65 deg, Wrong Bx, by JJL (g2p test run)
// May add more HRS packages later
///////////////////////////////////////////////////////////////////////////

    bool Forward(const double* V5_tg, double* V5_fp);
    bool Backward(const double* V5_fp, double* V5_tg);

    void DeltaCorrection();
    void XtgCorrection();
    void ThetatgCorrection();
    void YtgCorrection();
    void PhitgCorrection();

private:
    int iModelIndex;
    bool bIsLeftArm;

    map<int, G2PTrans *> mModel;
    map<const char *, int> mModelIndex;
    G2PTrans *pModel;

    ClassDef(HRSTransport,1);
};

#endif

///////////////////////////////////////////////////////////////////////////
// Original comments
///////////////////////////////////////////////////////////////////////////
//C++ interface to call JLL's Fortran transportation routines
//By Jixie Zhang
//Jixie has tested and verified that one can use the right arm routines to
//do the transportation and reconstruction. One only need to flip y_tr
//and phi_tr
//at the following pposition:
//1) flip y_tg and phi_tg then call transportation
//2) flip y_fp and phi_fp in the focus plane then call the reconstruction
//3) need to flip y_rec and phi_rec again
//C Forward transfer float s for left hrs with septum based on
//c Leftseptum_dir.dat
//c                     -JJL 4/4/04
//
//
//c typical call: answer = float (x,5)
//c INPUTS: x = 5 or more element array 
//c              x(1)=x0  (meters)
//c              x(2)=theta0 (really tan(theta0))
//c              x(3)=y0   (meters)
//c              x(4)=phi0 (really tan(phi0))
//c              x(5)=delta (fractional value NOT percent)
//c         M=5
//c
//c OUTPUT: units are the same as inputs
//c 
//c NOMENCLATURE: float  name = prefix + _sr_ +suffix
//c           prefixes:     x means xfinal
//c                         t means thetafinal
//c                         y means yfinal
//c                         p means phifinal
//c                         l means path length
//c     
//c           suffixes:     fp means target to focus
//c                         q1ex means target to Q1 exit
//c                         dent means target to Dipole entrance
//c                         dext means target to dipole exit
//c                         q3en means target to Q3 entrance
//c                         q3ex means target to Q3 exit
//c                         ep3  means septum entrance
//c                         ep4  means one quarter the way into the septum
//c                         ep5  means halfway through the septum
//c                         ep6  means three quarters the way through the septum
//c                         ep7  means septum exit
//c
//c          _sl_ is for left hrs with septum
//c
//c APERTURES:
//c    Coordinate systems are different in the spectrometer with septum model compared to the spectrometer
//c     without septum. So some apertures even in the body of the spectrometer appear to be different. Numbers
//c    given here supercede any numbers given in regard to transfer float s for the spectrometers alone.
//c
//c     ep3: -0.1486233 < x < -0.08869672
//c          -0.110 < y < 0.110
//c     ep4: -0.1792231 < x < -0.1089169
//c          -0.110 < y < 0.110
//c     ep5: -0.2209211 < x < -0.1353789
//c          -0.110 < y < 0.110
//c     ep6: -0.2763536 < x < -0.1697464
//c          -0.110 < y < 0.110
//c     ep7: -0.3485396 < x < -0.2156404
//c          -0.110 < y < 0.110
//c     q1ex is a circle of radius 0.1492 m
//c     dent is a trapazoid:
//c                                   -5.22008 < x < -4.98099
//c             -(-0.192436784*x -0.192436784) < y < -0.192436784*x -0.192436784
//c     dext is also a trapazoid: 
//c                                   -0.46188 < x < 0.46188
//c                   -(-0.01610808*x + 0.125) < y < -0.01610808*x + 0.125
//c     q3en is a circle of radius 0.3 m
//c     q3ex is a circle of radius 0.3 m
//c
//In each file there are 5 fortran functions:
//
//TXFIT(x,m)
//
//delta(x,m)
//
//theta(x,m)
//
//phi(x,m)
//
//y0(x,m)
//
//Calling convention is as follows. x is an array in which x(1->5) are 
//xf, thetaf, yf, phif, and x0 (the detector coordinates in the transport 
//coordinate system and the vertical position of the beam at the target). 
//For TXFIT m=1 for all others m=5.
//
//Procedure:
//
//Convert theta to a new orthogonalized theta' as follows:
//
//theta'=theta-TXFIT(x,1)
//
//Put the new theta' into your x array (or x(2)=x(2)-TXFIT(x,1) should work)
//
//Then I think it's self-explanitory:
//
//delta0 = delta(x,5)
//
//theta0= theta(x,5)
//
//phi0=phi(x,5)
//
//y0=y0(x,5)
//
///////////////////////////////////////////////////////////////////////////
