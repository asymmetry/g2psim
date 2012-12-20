#include "HRSTransport_G2P.hh"
#include "GlobalDebuger.hh"

#include <iostream>
#include <iomanip>
#include <math.h>

using namespace std;

//#define DEBUG_HRS_TRANSPORT 2


static const float m2cm = 100.0;

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

//C Transport electron thru septum magnet, rectangle defined by (jjl)
//! includes reduced aperatures due to bore cooler, 
//! assumes 0.8 cm thick side and 1.1 cm thick top and bottom
//Target to Septum entrance, -14.06cm>x>8.87cm, -9.9cm<y<9.9cm

/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////
//JLL provide the package for 5.69 degrees without the target field
//5.69 deg, septa, no target field setting
bool TransportRightHRS_g2_(double* pV5)
{
	float vector_jjl[]={pV5[0],pV5[1],pV5[2],pV5[3],pV5[4]};
	int ii=5;

	float x_test, y_test;

	//By Jixie @ 20110705: The septum tunnel is 30.4cm(W) X 24.4cm (H), center at (x=+/-23.6,y=0,z=68.646)
	
	//Target to Septum entrance, 
	//y, -478.5, -2.402,none, 94.2,409.7,120.,120.,-120.,-120.
	x_test = x_g2_sen_(vector_jjl, &ii)*m2cm;
	y_test = y_g2_sen_(vector_jjl, &ii)*m2cm;
	//By Jixie: revise according to the real septum tunnel size, not SNAKE directive file
	if( fabs(x_test)<8.4 || fabs(x_test)>38.6 || fabs(y_test)>12.0 ) 
		return false;

	//Target to Septum exit, 
	//y, 500., -2.402,none,0.,0.,0.,0.,0.,0.
	x_test = x_g2_sex_(vector_jjl, &ii)*m2cm;
	y_test = y_g2_sex_(vector_jjl, &ii)*m2cm;
	//By Jixie: revise according to the real septum tunnel size, not SNAKE directive file
	if( fabs(x_test)<8.4 || fabs(x_test)>38.6 || fabs(y_test)>12.0 ) 
		return false;

	//Target to Q1 entrance, circle of radius 14.92 cm
	//y,0., 0. ,none,150.,150.,0.,0.,0.,0.
	x_test = x_g2_q1en_(vector_jjl, &ii)*m2cm;
	y_test = y_g2_q1en_(vector_jjl, &ii)*m2cm;
	if( sqrt(x_test*x_test + y_test*y_test) > 14.92 )
		return false;

	//Target to Q1 exit, circle of radius 14.92 cm
	//y,0., 0. ,none,150.,150.,0.,0.,0.,0.
	x_test = x_g2_q1ex_(vector_jjl, &ii)*m2cm;
	y_test = y_g2_q1ex_(vector_jjl, &ii)*m2cm;
	if( sqrt(x_test*x_test + y_test*y_test) > 15.0 )
		return false;

	//in TCS,             xmin,xmax,ymax1,ymax2,ymin1,ymin2,  
	//Target to dipole entrance, trapezoid -616.999cm<x<-596.293cm  |y| < 14.55cm
	//y,-26812.93, -75.,none,-6169.99,-5962.93,145.5,145.5,-145.5,-145.5
	x_test = x_g2_dent_(vector_jjl, &ii)*m2cm;
	y_test = y_g2_dent_(vector_jjl, &ii)*m2cm;
	//By Jixie: I need to adjusted this aperture
	//if( (x_test<-616.999) || (x_test>-596.293) || fabs(y_test)>145.5 )
	//	return false;
	//if( (x_test<-520) || (x_test>-495) || fabs(y_test)>15.5 )
	//	return false;

	//Target to dipole exit, trapezoid -46.19cm<x<46.19cm  |y| < -0.0161*x+12.5
	//y,0., 0.,none,-461.88,461.88,132.44,117.56,-132.44,-117.56
	x_test = x_g2_dext_(vector_jjl, &ii)*m2cm;
	y_test = y_g2_dext_(vector_jjl, &ii)*m2cm;
	if( (x_test<-46.19) || (x_test>46.19) || fabs(y_test) > fabs(-0.0161*x_test+12.5) )
		return false;

	//Target to Q3 entrance, circle of radius 30.0 cm
	x_test = x_g2_q3en_(vector_jjl, &ii)*m2cm;
	y_test = y_g2_q3en_(vector_jjl, &ii)*m2cm;
	if( sqrt(x_test*x_test + y_test*y_test) > 30.0 )
		return false;

	//Target to Q3 exit, circle of radius 30.0 cm
	x_test = x_g2_q3ex_(vector_jjl, &ii)*m2cm;
	y_test = y_g2_q3ex_(vector_jjl, &ii)*m2cm;
	if( sqrt(x_test*x_test + y_test*y_test) > 30.0)
		return false;

	/////////////////////////////////////////////////////////////
	// succesfully reach focus plane
	float x_fp     = x_g2_fp_(vector_jjl,&ii);
	float theta_fp = t_g2_fp_(vector_jjl,&ii);
	float y_fp     = y_g2_fp_(vector_jjl,&ii);
	float phi_fp   = p_g2_fp_(vector_jjl,&ii);

	//reset the vector and return it back to the caller
	pV5[0] = (double)x_fp;
	pV5[1] = (double)theta_fp;
	pV5[2] = (double)y_fp;
	pV5[3] = (double)phi_fp;
	//pV5[4] = (double)delta_fp;  // delta is not change

	return true;
}

bool TransportLeftHRS_g2_(double* pV5)
{
	//use right arm routines for left arm before left arm is ready
	//return TransportLeftHRS(pV5);
	pV5[2]*=-1.;
	pV5[3]*=-1.;
	bool bGoodParticle=TransportRightHRS_g2_(pV5);
	pV5[2]*=-1.;
	pV5[3]*=-1.;
	return bGoodParticle;
}

//*************new septum with shims, 5.65 central ray, no target field*****M.Huang 2/23/2012***********

bool TransportRightHRS_Shim_484816_WrongBx(double* pV5)
{
	float vector_jjl[]={pV5[0],pV5[1],pV5[2],pV5[3],pV5[4]};
	int ii=5;

	float x_test, y_test;

	//Target to 1st end plane in Septum, 
	//y, -480., 0.,none,84.0,388.,97.,97.,-97.,-97. !!!NOTE that x, y switched here!!!
	x_test = x_r5p65_ep5(vector_jjl, ii)*m2cm;
	y_test = y_r5p65_ep5(vector_jjl, ii)*m2cm;
	if( fabs(x_test)<8.4 || fabs(x_test)>38.8 || fabs(y_test)>9.7 ) 
		return false;

	//Target to 3rd end plane in Septum, 
	//y, 480., 0.,none,84.0,388.,97.,97.,-97.,-97.  !!!NOTE that x, y switched here!!!
	x_test = x_r5p65_ep7(vector_jjl, ii)*m2cm;
	y_test = y_r5p65_ep7(vector_jjl, ii)*m2cm;
	if( fabs(x_test)<8.4 || fabs(x_test)>38.8 || fabs(y_test)>9.7 ) 
		return false;


	//Target to Q1 en 
	//y,-194., 0.,none,170.,170., 0.,0.,0.,0.
	x_test = x_r5p65_ep11_q1en(vector_jjl, ii)*m2cm;
	y_test = y_r5p65_ep11_q1en(vector_jjl, ii)*m2cm;
	if( sqrt(x_test*x_test + y_test*y_test) > 17.0 )
		return false;

	//Target to Q1 mid plane,
	//y,0., 0. ,none,150.,150.,0.,0.,0.,0.
	x_test = x_r5p65_ep13_q1(vector_jjl, ii)*m2cm;
	y_test = y_r5p65_ep13_q1(vector_jjl, ii)*m2cm;
	if( sqrt(x_test*x_test + y_test*y_test) > 15.)
		return false;

	//Target to Q1 ex,
	//y,610.4, 0. ,none,149.2,149.2,0.,0.,0.,0.
	x_test = x_r5p65_ep14_q1ex(vector_jjl, ii)*m2cm;
	y_test = y_r5p65_ep14_q1ex(vector_jjl, ii)*m2cm;
	if( sqrt(x_test*x_test + y_test*y_test) > 14.92 )
		return false;

	//in TCS,             xmin,xmax,ymax1,ymax2,ymin1,ymin2,  
	//Target to dipole entrance, trapezoid -616.999cm<x<-596.293cm  |y| < 14.55cm
	//y,-26812.93, -75.,none,-6169.99,-5962.93,145.5,145.5,-145.5,-145.5
	x_test = x_r5p65_ep23_den(vector_jjl, ii)*m2cm;
	y_test = y_r5p65_ep23_den(vector_jjl, ii)*m2cm;
	if( (x_test<-616.999) || (x_test>-596.293) || fabs(y_test)>145.5 )
		return false;

	//Target to dipole exit, trapezoid -46.19cm<x<46.19cm  |y| < -0.0161*x+12.5
	//y,0., 0.,none,-461.88,461.88,132.44,117.56,-132.44,-117.56
	x_test = x_r5p65_ep25_dex(vector_jjl, ii)*m2cm;
	y_test = y_r5p65_ep25_dex(vector_jjl, ii)*m2cm;
	if( (x_test<-46.19) || (x_test>46.19) || fabs(y_test) > fabs(-0.0161*x_test+12.5) )
		return false;

	//Target to Q3 entrance, circle of radius 30.0 cm
	x_test = x_r5p65_ep27_q3en(vector_jjl, ii)*m2cm;
	y_test = y_r5p65_ep27_q3en(vector_jjl, ii)*m2cm;
	if( sqrt(x_test*x_test + y_test*y_test) > 30.0 )
		return false;

	//Target to Q3 exit, circle of radius 30.0 cm
	x_test = x_r5p65_ep30_q3ex(vector_jjl, ii)*m2cm;
	y_test = y_r5p65_ep30_q3ex(vector_jjl, ii)*m2cm;
	if( sqrt(x_test*x_test + y_test*y_test) > 30.0)
		return false;

	/////////////////////////////////////////////////////////////
	// succesfully reach focus plane
	float x_fp     = x_r5p65_fp(vector_jjl,ii);
	float theta_fp = t_r5p65_fp(vector_jjl,ii);
	float y_fp     = y_r5p65_fp(vector_jjl,ii);
	float phi_fp   = p_r5p65_fp(vector_jjl,ii);


	//reset the vector and return it back to the caller
	pV5[0] = (double)x_fp;
	pV5[1] = (double)theta_fp;
	pV5[2] = (double)y_fp;
	pV5[3] = (double)phi_fp;
	//pV5[4] = (double)delta_fp;  // delta is not change

	return true;
}

bool TransportLeftHRS_Shim_484816_WrongBx(double* pV5)
{
	//use right arm routines for left arm before left arm is ready
	//return TransportLeftHRS(pV5);
	pV5[2]*=-1.;
	pV5[3]*=-1.;
	bool bGoodParticle=TransportRightHRS_Shim_484816_WrongBx(pV5);
	pV5[2]*=-1.;
	pV5[3]*=-1.;
	return bGoodParticle;
}

///////////////////////////////////////////////////////////////////////////////

bool TransportRightHRS_Shim_484816(double* pV5)
{
	float vector_jjl[]={pV5[0],pV5[1],pV5[2],pV5[3],pV5[4]};
	int nvar=5;
	int *ii=&nvar;

	float x_test, y_test;

	//Target to Septum exit, 
	//y, 480., 0.,none,84.0,388.,97.,97.,-97.,-97.
	x_test = x_r5p65_484816_sepex_(vector_jjl, ii)*m2cm;
	y_test = y_r5p65_484816_sepex_(vector_jjl, ii)*m2cm;
	if( fabs(x_test)<8.4 || fabs(x_test)>38.8 || fabs(y_test)>9.7 ) 
		return false;


	//Target to Q1 en 
	//y,-200., 0.,none,150.,150.,0.,0.,0.,0.
	x_test = x_r5p65_484816_q1ent_(vector_jjl, ii)*m2cm;
	y_test = y_r5p65_484816_q1ent_(vector_jjl, ii)*m2cm;
	if( sqrt(x_test*x_test + y_test*y_test) > 15.0 )
		return false;

	//Target to Q1 ex,
	//y,610.4, 0. ,none,149.2,149.2,0.,0.,0.,0.
	x_test = x_r5p65_484816_q1ext_(vector_jjl, ii)*m2cm;
	y_test = y_r5p65_484816_q1ext_(vector_jjl, ii)*m2cm;
	if( sqrt(x_test*x_test + y_test*y_test) > 14.92 )
		return false;

	//Target to Q2 ex,
	//y,-3040.40, 30.,none, 259.81, 300., 1316.53, 0.,0.,0.
	x_test = x_r5p65_484816_q2ext_(vector_jjl, ii)*m2cm;
	y_test = y_r5p65_484816_q2ext_(vector_jjl, ii)*m2cm;
	if( sqrt(x_test*x_test + y_test*y_test) > 30.0 )
		return false;

	//in TCS,             xmin,xmax,ymax1,ymax2,ymin1,ymin2,  
	//Target to dipole entrance, 
	//y,15121.67, -105.,none,-5220.08,-4980.99,132.44,117.56,-132.44,-117.56
	x_test = x_r5p65_484816_den_(vector_jjl, ii)*m2cm;
	y_test = y_r5p65_484816_den_(vector_jjl, ii)*m2cm;
	if( (x_test<-522.0) || (x_test>-498.1) || fabs(y_test)>13.244 )
		return false;

	//Target to dipole exit, trapezoid -46.19cm<x<46.19cm  |y| < -0.0161*x+12.5
	//y,0., 0.,none,-461.88,461.88,132.44,117.56,-132.44,-117.56
	x_test = x_r5p65_484816_dex_(vector_jjl, ii)*m2cm;
	y_test = y_r5p65_484816_dex_(vector_jjl, ii)*m2cm;
	if( (x_test<-46.19) || (x_test>46.19) || fabs(y_test) > fabs(-0.0161*x_test+12.5) )
		return false;

	//Target to Q3 entrance, circle of radius 30.0 cm
	x_test = x_r5p65_484816_q3ent_(vector_jjl, ii)*m2cm;
	y_test = y_r5p65_484816_q3ent_(vector_jjl, ii)*m2cm;
	if( sqrt(x_test*x_test + y_test*y_test) > 30.0 )
		return false;

	//Target to Q3 exit, circle of radius 30.0 cm
	x_test = x_r5p65_484816_q3ext_(vector_jjl, ii)*m2cm;
	y_test = y_r5p65_484816_q3ext_(vector_jjl, ii)*m2cm;
	if( sqrt(x_test*x_test + y_test*y_test) > 30.0)
		return false;

	/////////////////////////////////////////////////////////////
	// succesfully reach focus plane
	float x_fp     = x_r5p65_484816_fp_(vector_jjl,ii);
	float theta_fp = t_r5p65_484816_fp_(vector_jjl,ii);
	float y_fp     = y_r5p65_484816_fp_(vector_jjl,ii);
	float phi_fp   = p_r5p65_484816_fp_(vector_jjl,ii);


	//reset the vector and return it back to the caller
	pV5[0] = (double)x_fp;
	pV5[1] = (double)theta_fp;
	pV5[2] = (double)y_fp;
	pV5[3] = (double)phi_fp;
	//pV5[4] = (double)delta_fp;  // delta is not change

	return true;
}

bool TransportLeftHRS_Shim_484816(double* pV5)
{
	float vector_jjl[]={pV5[0],pV5[1],pV5[2],pV5[3],pV5[4]};
	int nvar=5;
	int ii=5;

	float x_test, y_test;

	//Target to Septum ep5
	//y, -480., 0.,none,84.0,388.,97.,97.,-97.,-97.  ;84.388.
	x_test = x_sl5p65_400016_ep5(vector_jjl, ii)*m2cm;
	y_test = y_sl5p65_400016_ep5(vector_jjl, ii)*m2cm;
	if( fabs(x_test)<8.4 || fabs(x_test)>38.8 || fabs(y_test)>9.7 ) 
		return false;

	//Target to Septum ep7
	//y, 480., 0.,none,84.0,388.,97.,97.,-97.,-97.
	x_test = x_sl5p65_400016_ep7(vector_jjl, ii)*m2cm;
	y_test = y_sl5p65_400016_ep7(vector_jjl, ii)*m2cm;
	if( fabs(x_test)<8.4 || fabs(x_test)>38.8 || fabs(y_test)>9.7 ) 
		return false;

	//Target to Q1 en ep10
	//y,-200., 0.,none,149.2,149.2,0.,0.,0.,0.
	x_test = x_sl5p65_400016_ep10_q1en(vector_jjl, ii)*m2cm;
	y_test = y_sl5p65_400016_ep10_q1en(vector_jjl, ii)*m2cm;
	if( sqrt(x_test*x_test + y_test*y_test) > 14.92 )
		return false;

	//Target to Q1 ex ep13
	//y,675., 0. ,none,300.,300.,0.,0.,0.,0.
	x_test = x_sl5p65_400016_ep13_q1ex(vector_jjl, ii)*m2cm;
	y_test = y_sl5p65_400016_ep13_q1ex(vector_jjl, ii)*m2cm;
	if( sqrt(x_test*x_test + y_test*y_test) > 30. )
		return false;



	//Target to dipole exit, trapezoid -46.19cm<x<46.19cm  |y| < -0.0161*x+12.5 ep24
	//y,0., 0.,none,-461.88,461.88,132.44,117.56,-132.44,-117.56 
	x_test = x_sl5p65_400016_ep24_dex(vector_jjl, ii)*m2cm;
	y_test = y_sl5p65_400016_ep24_dex(vector_jjl, ii)*m2cm;
	if( (x_test<-46.19) || (x_test>46.19) || fabs(y_test) > fabs(-0.0161*x_test+12.5) )
		return false;

	//Target to Q3 entrance, circle of radius 30.0 cm ep26
	//y,915.,0.,none,300.,300.,0.,0.,0.,0.
	x_test = x_sl5p65_400016_ep26_q3en(vector_jjl, ii)*m2cm;
	y_test = y_sl5p65_400016_ep26_q3en(vector_jjl, ii)*m2cm;
	if( sqrt(x_test*x_test + y_test*y_test) > 30.0 )
		return false;

	//Target to Q3 exit, circle of radius 30.0 cm ep29
	//y,575., 0.,none,300.,300.,0.,0.,0.,0.
	x_test = x_sl5p65_400016_ep29_q3ex(vector_jjl, ii)*m2cm;
	y_test = y_sl5p65_400016_ep29_q3ex(vector_jjl, ii)*m2cm;
	if( sqrt(x_test*x_test + y_test*y_test) > 30.0)
		return false;

	/////////////////////////////////////////////////////////////
	// succesfully reach focus plane
	float x_fp     = x_sl5p65_484816_unrastered_fp(vector_jjl,ii);
	float theta_fp = t_sl5p65_484816_unrastered_fp(vector_jjl,ii);
	float y_fp     = y_sl5p65_484816_unrastered_fp(vector_jjl,ii);
	float phi_fp   = p_sl5p65_484816_unrastered_fp(vector_jjl,ii);


	//reset the vector and return it back to the caller
	pV5[0] = (double)x_fp;
	pV5[1] = (double)theta_fp;
	pV5[2] = (double)y_fp;
	pV5[3] = (double)phi_fp;
	//pV5[4] = (double)delta_fp;  // delta is not change

	return true;



// 	//use right arm routines for left arm before left arm is ready
// 	//return TransportLeftHRS(pV5);
// 	pV5[2]*=-1.;
// 	pV5[3]*=-1.;
// 	bool bGoodParticle=TransportRightHRS_Shim_484816(pV5);
// 	pV5[2]*=-1.;
// 	pV5[3]*=-1.;
// 	return bGoodParticle;
}

///////////////////////////////////////////////////////////////////////////////
//empty prototype, need to be filled by Min

bool TransportRightHRS_Shim_403216(double* pV5)
{
	return true;
}

bool TransportLeftHRS_Shim_403216(double* pV5)
{
	return true;
}

void ReconstructRightHRS_Shim_403216(double* pV5)
{
}

void ReconstructLeftHRS_Shim_403216(double* pV5)
{
}


///////////////////////////////////////////////////////////////////////////////

bool TransportRightHRS_Shim_400016(double* pV5)
{
	return true;
}

///////////////////////////////////////////////////////////////////////////////
//*******400016 septum with shims, 5.65 central ray, no target field*********
//*****M.Huang 11/12/2012********
//the source code can be found in Fwd_sl5p65_400016.cpp
//1.5cm raster

bool TransportLeftHRS_Shim_400016(double* pV5)
{

	float vector_jjl[]={pV5[0],pV5[1],pV5[2],pV5[3],pV5[4]};
	int nvar=5;
	int ii=5;

	float x_test, y_test;

	//Target to Septum ep5
	//y, -480., 0.,none,84.0,388.,97.,97.,-97.,-97.  ;84.388.
	x_test = x_sl5p65_400016_ep5(vector_jjl, ii)*m2cm;
	y_test = y_sl5p65_400016_ep5(vector_jjl, ii)*m2cm;
	if( fabs(x_test)<8.4 || fabs(x_test)>38.8 || fabs(y_test)>9.7 ) 
		return false;

	//Target to Septum ep7
	//y, 480., 0.,none,84.0,388.,97.,97.,-97.,-97.
	x_test = x_sl5p65_400016_ep7(vector_jjl, ii)*m2cm;
	y_test = y_sl5p65_400016_ep7(vector_jjl, ii)*m2cm;
	if( fabs(x_test)<8.4 || fabs(x_test)>38.8 || fabs(y_test)>9.7 ) 
		return false;

	//Target to Q1 en ep10
	//y,-200., 0.,none,149.2,149.2,0.,0.,0.,0.
	x_test = x_sl5p65_400016_ep10_q1en(vector_jjl, ii)*m2cm;
	y_test = y_sl5p65_400016_ep10_q1en(vector_jjl, ii)*m2cm;
	if( sqrt(x_test*x_test + y_test*y_test) > 14.92 )
		return false;

	//Target to Q1 ex ep13
	//y,675., 0. ,none,300.,300.,0.,0.,0.,0.
	x_test = x_sl5p65_400016_ep13_q1ex(vector_jjl, ii)*m2cm;
	y_test = y_sl5p65_400016_ep13_q1ex(vector_jjl, ii)*m2cm;
	if( sqrt(x_test*x_test + y_test*y_test) > 30. )
		return false;



	//Target to dipole exit, trapezoid -46.19cm<x<46.19cm  |y| < -0.0161*x+12.5 ep24
	//y,0., 0.,none,-461.88,461.88,132.44,117.56,-132.44,-117.56 
	x_test = x_sl5p65_400016_ep24_dex(vector_jjl, ii)*m2cm;
	y_test = y_sl5p65_400016_ep24_dex(vector_jjl, ii)*m2cm;
	if( (x_test<-46.19) || (x_test>46.19) || fabs(y_test) > fabs(-0.0161*x_test+12.5) )
		return false;

	//Target to Q3 entrance, circle of radius 30.0 cm ep26
	//y,915.,0.,none,300.,300.,0.,0.,0.,0.
	x_test = x_sl5p65_400016_ep26_q3en(vector_jjl, ii)*m2cm;
	y_test = y_sl5p65_400016_ep26_q3en(vector_jjl, ii)*m2cm;
	if( sqrt(x_test*x_test + y_test*y_test) > 30.0 )
		return false;

	//Target to Q3 exit, circle of radius 30.0 cm ep29
	//y,575., 0.,none,300.,300.,0.,0.,0.,0.
	x_test = x_sl5p65_400016_ep29_q3ex(vector_jjl, ii)*m2cm;
	y_test = y_sl5p65_400016_ep29_q3ex(vector_jjl, ii)*m2cm;
	if( sqrt(x_test*x_test + y_test*y_test) > 30.0)
		return false;

	/////////////////////////////////////////////////////////////
	// succesfully reach focus plane
	float x_fp     = x_sl5p65_400016_fp(vector_jjl,ii);
	float theta_fp = t_sl5p65_400016_fp(vector_jjl,ii);
	float y_fp     = y_sl5p65_400016_fp(vector_jjl,ii);
	float phi_fp   = p_sl5p65_400016_fp(vector_jjl,ii);


	//reset the vector and return it back to the caller
	pV5[0] = (double)x_fp;
	pV5[1] = (double)theta_fp;
	pV5[2] = (double)y_fp;
	pV5[3] = (double)phi_fp;
	//pV5[4] = (double)delta_fp;  // delta is not change

	return true;
}

void ReconstructRightHRS_Shim_400016(double* pV5)
{
}

void ReconstructLeftHRS_Shim_400016(double* pV5)
{
	float vector_jjl[]={pV5[0],pV5[1],pV5[2],pV5[3],pV5[4]};
	int ii=5;

	vector_jjl[1]   = vector_jjl[1] - txfit_sl5p65_400016(vector_jjl,ii);

	float x_or      = vector_jjl[4];
	float delta_rec = delta_sl5p65_400016(vector_jjl,ii);
	float theta_rec = theta_sl5p65_400016(vector_jjl,ii);
	float phi_rec   = phi_sl5p65_400016(vector_jjl,ii); 
	float y_rec     = y00_sl5p65_400016(vector_jjl,ii);    

	//reset the vector and return it back to the caller
	pV5[0] = (double)x_or;
	pV5[1] = (double)theta_rec;
	pV5[2] = (double)y_rec;
	pV5[3] = (double)phi_rec;
	pV5[4] = (double)delta_rec;
}


//////////////////////////////////////////////////////////////////////////////////
//empty prototype, need to be filled by Min

bool TransportRightHRS_12p5_Min_(double* pV5)
{
	float vector_jjl[]={pV5[0],pV5[1],pV5[2],pV5[3],pV5[4]};
	int ii=5;

	float x_test, y_test;

	//Target to Q1 exit, circle of radius 14.92 cm
	//y,0., 0. ,none,150.,150.,0.,0.,0.,0.
	x_test = x_r12p5_q1ex_(vector_jjl, &ii)*m2cm;
	y_test = y_r12p5_q1ex_(vector_jjl, &ii)*m2cm;
	if( (x_test*x_test + y_test*y_test) > (15.*15.) )
		return false;

	//in TCS,             xmin,xmax,ymax1,ymax2,ymin1,ymin2,  
	//Target to dipole entrance, trapezoid -616.999cm<x<-596.293cm  |y| < 14.55cm
	//y,-26812.93, -75.,none,-6169.99,-5962.93,145.5,145.5,-145.5,-145.5
	x_test = x_r12p5_dent_(vector_jjl, &ii)*m2cm;
	y_test = y_r12p5_dent_(vector_jjl, &ii)*m2cm;
	if( (x_test<-616.999) || (x_test>-596.293) || fabs(y_test) > 145.5 )
		return false;

	//Target to dipole exit, trapezoid -46.19cm<x<46.19cm  |y| < -0.0161*x+12.5
	//y,0., 0.,none,-461.88,461.88,132.44,117.56,-132.44,-117.56
	x_test = x_r12p5_dext_(vector_jjl, &ii)*m2cm;
	y_test = y_r12p5_dext_(vector_jjl, &ii)*m2cm;
	//By Jx: this cut have problems
	//if( (x_test<-46.19) || (x_test>46.19) || fabs(y_test) > fabs(-0.0161*x_test+12.5) )
	//	return false;
	
	//y,1154.7, -30.,none,153.59,846.41,150.,150.,-150.,-150.
	//if( (x_test<15.359) || (x_test>84.641) || fabs(y_test) > 15.0 )
	//	return false;

	//Target to Q3 entrance, circle of radius 30.0 cm
	x_test = x_r12p5_q3en_(vector_jjl, &ii)*m2cm;
	y_test = y_r12p5_q3en_(vector_jjl, &ii)*m2cm;
	if( (x_test*x_test + y_test*y_test) > (30.0*30.0) )
		return false;

	//Min forgot to generate routines for this end plane
	//Target to Q3 exit, circle of radius 30.0 cm
	//x_test = x_r12p5_q3ex_(vector_jjl, &ii)*m2cm;
	//y_test = y_r12p5_q3ex_(vector_jjl, &ii)*m2cm;
	//if( (x_test*x_test + y_test*y_test) > (30.0*30.0))
	//	return false;

	/////////////////////////////////////////////////////////////
	// succesfully reach focus plane
	float x_fp     = x_r12p5_fp_(vector_jjl,&ii);
	float theta_fp = t_r12p5_fp_(vector_jjl,&ii);
	float y_fp     = y_r12p5_fp_(vector_jjl,&ii);
	float phi_fp   = p_r12p5_fp_(vector_jjl,&ii);

	//reset the vector and return it back to the caller
	pV5[0] = (double)x_fp;
	pV5[1] = (double)theta_fp;
	pV5[2] = (double)y_fp;
	pV5[3] = (double)phi_fp;
	//pV5[4] = (double)delta_fp;  // delta is not change

	return true;
}


bool TransportLeftHRS_12p5_Min_(double* pV5)
{
	//use right arm routines for left arm before left arm is ready
	//return TransportLeftHRS(pV5);
	pV5[2]*=-1.;
	pV5[3]*=-1.;
	bool bGoodParticle=TransportRightHRS_12p5_Min_(pV5);
	pV5[2]*=-1.;
	pV5[3]*=-1.;
	return bGoodParticle;
}


//////////////////////////////////////////////////////////////////////////////////////
/* the detail of  JJLerose vector in focus plane is*/
//vector_jjl[0] = x_fp;
//vector_jjl[1] = theta_fp;
//vector_jjl[2] = y_fp;
//vector_jjl[3] = phi_fp;
//vector_jjl[4] = x_or;
//the output is 
//vector_jjl[0] = x_or;
//vector_jjl[1] = theta_rec;
//vector_jjl[2] = y_rec;
//vector_jjl[3] = phi_rec;
//vector_jjl[4] = delta_rec;
/////////////////////////////
void ReconstructRightHRS_g2_(double *pV5)
{
	float vector_jjl[]={pV5[0],pV5[1],pV5[2],pV5[3],pV5[4]};
	int ii=5,jj=1;

	vector_jjl[1]   = vector_jjl[1] - g2_txfit_(vector_jjl,&jj);

	float x_or      = vector_jjl[4];
	float delta_rec = g2_delta_(vector_jjl,&ii);
	float theta_rec = g2_theta_(vector_jjl,&ii);
	float phi_rec   = g2_phi_(vector_jjl,&ii); 
	float y_rec     = g2_y00_(vector_jjl,&ii); 

	//reset the vector and return it back to the caller
	pV5[0] = (double)x_or;
	pV5[1] = (double)theta_rec;
	pV5[2] = (double)y_rec;
	pV5[3] = (double)phi_rec;
	pV5[4] = (double)delta_rec;
}

void ReconstructLeftHRS_g2_(double *pV5)
{
	//in order to call right arm routines, need to flip y, phi 
	pV5[2]*=-1;
	pV5[3]*=-1;
	ReconstructRightHRS_g2_(pV5);
	pV5[2]*=-1;
	pV5[3]*=-1;
}


//*************new septum with shims, 5.65 central ray, no target field*****M.Huang 2/23/2012***********
//2cm raster, With wrong Bx in septum field
void ReconstructRightHRS_Shim_484816_WrongBx(double *pV5)
{

	float vector_jjl[]={pV5[0],pV5[1],pV5[2],pV5[3],pV5[4]};
	int ii=5;

	vector_jjl[1]   = vector_jjl[1] - r5p65_txfit(vector_jjl,ii);

	float x_or      = vector_jjl[4];
	float delta_rec = r5p65_delta(vector_jjl,ii);
	float theta_rec = r5p65_theta(vector_jjl,ii);
	float phi_rec   = r5p65_phi(vector_jjl,ii); 
	float y_rec     = r5p65_y00(vector_jjl,ii);    

	//reset the vector and return it back to the caller
	pV5[0] = (double)x_or;
	pV5[1] = (double)theta_rec;
	pV5[2] = (double)y_rec;
	pV5[3] = (double)phi_rec;
	pV5[4] = (double)delta_rec;
}

void ReconstructLeftHRS_Shim_484816_WrongBx(double *pV5)
{
	//in order to call right arm routines, need to flip y, phi 
	pV5[2]*=-1;
	pV5[3]*=-1;
	ReconstructRightHRS_Shim_484816_WrongBx(pV5);
	pV5[2]*=-1;
	pV5[3]*=-1;
}

void ReconstructRightHRS_Shim_484816(double *pV5)
{

	float vector_jjl[]={pV5[0],pV5[1],pV5[2],pV5[3],pV5[4]};
	int n=5;
	int *ii=&n;

	vector_jjl[1]   = vector_jjl[1] - txfit_r5p65_484816_(vector_jjl,ii);

	float x_or      = vector_jjl[4];
	float delta_rec = delta_r5p65_484816_(vector_jjl,ii);
	float theta_rec = theta_r5p65_484816_(vector_jjl,ii);
	float phi_rec   = phi_r5p65_484816_(vector_jjl,ii); 
	float y_rec     = y00_r5p65_484816_(vector_jjl,ii);    

	//reset the vector and return it back to the caller
	pV5[0] = (double)x_or;
	pV5[1] = (double)theta_rec;
	pV5[2] = (double)y_rec;
	pV5[3] = (double)phi_rec;
	pV5[4] = (double)delta_rec;
}

void ReconstructLeftHRS_Shim_484816(double *pV5)
{
	//in order to call right arm routines, need to flip y, phi 
	pV5[2]*=-1;
	pV5[3]*=-1;
	ReconstructRightHRS_Shim_484816(pV5);
	pV5[2]*=-1;
	pV5[3]*=-1;
}


/////////////////////////////
void ReconstructRightHRS_12p5_Min_(double *pV5)
{
	float vector_jjl[]={pV5[0],pV5[1],pV5[2],pV5[3],pV5[4]};
	int ii=5;

	//these 2 routines are not available yet, need Min's input
	//vector_jjl[1]   = vector_jjl[1] - r12p5_txfit_(vector_jjl,&jj);
	//vector_jjl[2]   = vector_jjl[2] - r12p5_pyfit_(vector_jjl,&jj);

	float x_or      = vector_jjl[4];
	float delta_rec = r12p5_delta_(vector_jjl,&ii);
	float theta_rec = r12p5_theta_(vector_jjl,&ii);
	float phi_rec   = r12p5_phi_(vector_jjl,&ii);  
	float y_rec     = r12p5_y00_(vector_jjl,&ii); 

	//reset the vector and return it back to the caller
	pV5[0] = (double)x_or;
	pV5[1] = (double)theta_rec;
	pV5[2] = (double)y_rec;
	pV5[3] = (double)phi_rec;
	pV5[4] = (double)delta_rec;
}

void ReconstructLeftHRS_12p5_Min_(double *pV5)
{
	pV5[2]*=-1;
	pV5[3]*=-1;
	ReconstructRightHRS_12p5_Min_(pV5);
	pV5[2]*=-1;
	pV5[3]*=-1;
}

/////////////////////////////////////////////////////////////////////////////////////

//Sieve Slit to Target reconstruction
//all belows need to be updated
void ReconstructRight_S2T_S_H90R1(double *pV5)
{
	float vector_jjl[]={pV5[0],pV5[1],pV5[2],pV5[3],pV5[4]};
	int ii=5,jj=1;

	vector_jjl[1]   = vector_jjl[1] - s2t_sr5p69_2gev_txfit_(vector_jjl,&jj);
	vector_jjl[2]   = vector_jjl[2] - s2t_sr5p69_2gev_pyfit_(vector_jjl,&jj);

	//delta is not changed
	float x_rec = 0; //waiting for Min's routine
	float delta_rec = vector_jjl[4];
	float theta_rec = s2t_sr5p69_theta_(vector_jjl,&ii);
	float phi_rec   = s2t_sr5p69_phi_(vector_jjl,&ii); 
	float y_rec     = s2t_sr5p69_y00_(vector_jjl,&ii); 

	//reset the vector and return it back to the caller
	pV5[0] = (double)x_rec;
	pV5[1] = (double)theta_rec;
	pV5[2] = (double)y_rec;
	pV5[3] = (double)phi_rec;
	pV5[4] = (double)delta_rec;
}


void ReconstructLeft_S2T_S_H90R1(double *pV5)
{
	//in order to call right arm routines, need to flip y, phi 
	pV5[2]*=-1;
	pV5[3]*=-1;
	ReconstructRight_S2T_S_H90R1(pV5);
	pV5[2]*=-1;
	pV5[3]*=-1;
}


//GEP
void ReconstructLeft_S2T_S_H07R1(double *pV5)
{
	;
}

void ReconstructRight_S2T_S_H07R1(double *pV5)
{
	;
}

