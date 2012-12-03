//C++ interface to call JLL's Fortran transportation routines
//By Jixie Zhang
//Jixie has tested and verified that one can use the right arm routines to
//do the transportation and reconstruction. One only need to flip y_tr and phi_tr
//at the following pposition:
//1) flip y_tg and phi_tg then call transportation
//2) flip y_fp and phi_fp in the focus plane then call the reconstruction
//3) need to flip y_rec and phi_rec again 

//C Forward transfer float s for left hrs with septum based on Leftseptum_dir.dat
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
//y00(x,m)
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
//y0=y00(x,5)
//
#ifndef HRS_TRANSPORT_H
#define HRS_TRANSPORT_H

#include "HRSTransport_STD.hh"
#include "HRSTransport_GDH.hh"
#include "HRSTransport_G2P.hh"

//Try to swing particle through the HRS, if it does not go through
//return false, otherwise return true and set the parameters of the 
//focus plan into the array
//
//For transportation, the input array is 
//the detail of input vector is
//vector_jjl[0] = x_tr;
//vector_jjl[1] = theta_tr;
//vector_jjl[2] = y_tr;
//vector_jjl[3] = phi_tr;
//vector_jjl[4] = delta_tr;
//the output is 
//vector_jjl[0] = x_fp;
//vector_jjl[1] = theta_fp;
//vector_jjl[2] = y_fp;
//vector_jjl[3] = phi_fp;
//vector_jjl[4] = delta_fp;  // delta is not change
//all length in unit of meter
//The name ended with underscore are using fortran code

//bool TransportLeftHRS_ID(double* vector_jjl, int experiment);	
//bool TransportRightHRS_ID(double* vector_jjl, int experiment);	

////////////////////////////////////////////////////////////////////////////////////////

// For reconstruction
//vector_jjl[0] = x_fp;
//vector_jjl[1] = theta_fp;
//vector_jjl[2] = y_fp;
//vector_jjl[3] = phi_fp;6
//vector_jjl[4] = x_or;
//the output is 
//vector_jjl[0] = x_or;
//vector_jjl[1] = theta_rec;
//vector_jjl[2] = y_rec;
//vector_jjl[3] = phi_rec;
//vector_jjl[4] = delta_rec;
//all length in unit of meter

//void ReconstructLeftHRS_ID(double* vector_jjl, int experiment);	
//void ReconstructRightHRS_ID(double* vector_jjl,int experiment);	

///////////////////////////////////////////////////////////////////////////

//input:
//pVector_tg[0] = x_tr;
//pVector_tg[1] = theta_tr;
//pVector_tg[2] = y_tr;
//pVector_tg[3] = phi_tr;
//pVector_tg[4] = delta_tr;
//output
//pVector_fp[0] = x_fp;
//pVector_fp[1] = theta_fp;
//pVector_fp[2] = y_fp;
//pVector_fp[3] = phi_fp;
//pVector_fp[4] = delta_fp;		// = delta_tr; not change
//pVector_rec_tg[0] = x_rec ;	// = x_tr; not change
//pVector_rec_tg[1] = theta_rec;
//pVector_rec_tg[2] = y_rec;
//pVector_rec_tg[3] = phi_rec;
//pVector_rec_tg[4] = delta_rec;
void CMP_HRS_TRANSPORTATION(float *pV5_tg, float *pV5_g2p);
void CMP_HRS_TRANSPORTATION(double *pV5_tg, double *pV5_g2p);

void DeltaCorrection(double &pDelta, double &pP_rec);
void XtgCorrection(double &pX,double pP_rec);
void ThetatgCorrection(double &pTheta,double pP_rec);
void YtgCorrection(double &pY,double pP_rec);
void PhitgCorrection(double &pPhi,double pP_rec);


//transport particles through HRS, use iExperiment to identify which HRS packages
//will be use, iExperiment==[10,20) is for g2p, 20 for E97110 GDH iExperiment
//will add more HRS packages later
bool SNAKEThruHRS(int pIsLeftArm, double pEndPlaneAngle, double pX0_tg_m, 
				  double pHRSMomentum, int iFieldRotation, int iExperiment,
				  double* pV5_tg, double* pV5_fp);

#endif
