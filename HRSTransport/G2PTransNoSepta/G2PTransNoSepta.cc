////////////////////////////////////////////////////////////////////////
// No septa, 12.5 central ray, no target field
// 3cm raster
// by M. Huang 11/12/2012
////////////////////////////////////////////////////////////////////////

#include <cmath>

#include "TROOT.h"
#include "TMath.h"

#include "Fwd_r12p5_NoSepta.h"
#include "Bwd_r12p5_NoSepta.h"

#include "G2PTransNoSepta.hh"

//using namespace SNoSepta;

const float m2cm = 100.0;

G2PTransNoSepta::G2PTransNoSepta()
    :cModelAngle(12.50*TMath::Pi()/180.0)
{
    // Nothing to do
}

G2PTransNoSepta::~G2PTransNoSepta()
{
    // Nothing to do
}

bool G2PTransNoSepta::TransLeftHRS(double* pV5)
{
    //use right arm routines for left arm before left arm is ready
    //return TransportLeftHRS(pV5);
    pV5[2]*=-1.;
    pV5[3]*=-1.;
    bool bGoodParticle=TransRightHRS(pV5);
    pV5[2]*=-1.;
    pV5[3]*=-1.;

    return bGoodParticle;
}

bool G2PTransNoSepta::TransRightHRS(double* pV5)
{
	float vector_jjl[]={pV5[0],pV5[1],pV5[2],pV5[3],pV5[4]};
    int *ii = new int; (*ii) = 5;

	float x_test, y_test;

	//Target to Q1 exit, circle of radius 14.92 cm
	//y,0., 0. ,none,150.,150.,0.,0.,0.,0.
	x_test = x_r12p5_q1ex_(vector_jjl, ii)*m2cm;
	y_test = y_r12p5_q1ex_(vector_jjl, ii)*m2cm;
	if( (x_test*x_test + y_test*y_test) > (15.*15.) )
		return false;

	//in TCS,             xmin,xmax,ymax1,ymax2,ymin1,ymin2,  
	//Target to dipole entrance, trapezoid -616.999cm<x<-596.293cm  |y| < 14.55cm
	//y,-26812.93, -75.,none,-6169.99,-5962.93,145.5,145.5,-145.5,-145.5
	x_test = x_r12p5_dent_(vector_jjl, ii)*m2cm;
	y_test = y_r12p5_dent_(vector_jjl, ii)*m2cm;
	if( (x_test<-616.999) || (x_test>-596.293) || fabs(y_test) > 145.5 )
		return false;

	//Target to dipole exit, trapezoid -46.19cm<x<46.19cm  |y| < -0.0161*x+12.5
	//y,0., 0.,none,-461.88,461.88,132.44,117.56,-132.44,-117.56
	x_test = x_r12p5_dext_(vector_jjl, ii)*m2cm;
	y_test = y_r12p5_dext_(vector_jjl, ii)*m2cm;
	//By Jixie: this cut have problems
	//if( (x_test<-46.19) || (x_test>46.19) || fabs(y_test) > fabs(-0.0161*x_test+12.5) )
	//	return false;
	
	//y,1154.7, -30.,none,153.59,846.41,150.,150.,-150.,-150.
	//if( (x_test<15.359) || (x_test>84.641) || fabs(y_test) > 15.0 )
	//	return false;

	//Target to Q3 entrance, circle of radius 30.0 cm
	x_test = x_r12p5_q3en_(vector_jjl, ii)*m2cm;
	y_test = y_r12p5_q3en_(vector_jjl, ii)*m2cm;
	if( (x_test*x_test + y_test*y_test) > (30.0*30.0) )
		return false;

	//Min forgot to generate routines for this end plane
	//Target to Q3 exit, circle of radius 30.0 cm
	//x_test = x_r12p5_q3ex_(vector_jjl, ii)*m2cm;
	//y_test = y_r12p5_q3ex_(vector_jjl, ii)*m2cm;
	//if( (x_test*x_test + y_test*y_test) > (30.0*30.0))
	//	return false;

	/////////////////////////////////////////////////////////////
	// succesfully reach focus plane
	float x_fp     = x_r12p5_fp_(vector_jjl, ii);
	float theta_fp = t_r12p5_fp_(vector_jjl, ii);
	float y_fp     = y_r12p5_fp_(vector_jjl, ii);
	float phi_fp   = p_r12p5_fp_(vector_jjl, ii);

	//reset the vector and return it back to the caller
	pV5[0] = (double)x_fp;
	pV5[1] = (double)theta_fp;
	pV5[2] = (double)y_fp;
	pV5[3] = (double)phi_fp;
	//pV5[4] = (double)delta_fp;  // delta is not change

    delete ii;

	return true;
}

void G2PTransNoSepta::ReconLeftHRS(double* pV5)
{   
    //in order to call right arm routines, need to flip y, phi 
    pV5[2]*=-1;
    pV5[3]*=-1;
    ReconRightHRS(pV5);
    pV5[2]*=-1;
    pV5[3]*=-1;
}

void G2PTransNoSepta::ReconRightHRS(double* pV5)
{
	float vector_jjl[]={pV5[0],pV5[1],pV5[2],pV5[3],pV5[4]};
	int *ii = new int; (*ii) = 5;

	float x_or      = vector_jjl[4];
	float delta_rec = delta_r12p5_(vector_jjl, ii);
	float theta_rec = theta_r12p5_(vector_jjl, ii);
	float phi_rec   = phi_r12p5_(vector_jjl, ii); 
	float y_rec     = y00_r12p5_(vector_jjl, ii); 

	//reset the vector and return it back to the caller
	pV5[0] = (double)x_or;
	pV5[1] = (double)theta_rec;
	pV5[2] = (double)y_rec;
	pV5[3] = (double)phi_rec;
	pV5[4] = (double)delta_rec;

    delete ii;
} 
