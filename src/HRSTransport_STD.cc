#include "HRSTransport_STD.hh"
#include "GlobalDebuger.hh"

#include <iostream>
#include <iomanip>
#include <math.h>

using namespace std;

#define DEBUG_HRS_TRANSPORT 2


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
//Standard HRS seting, no septum field, no target field
/////////////////////////////////////////////////////////////////////////////

bool TransportRightHRS(double* pV5)
{
	float vector_jjl[]={pV5[0],pV5[1],pV5[2],pV5[3],pV5[4]};
	float x_test, y_test;
	int ii=5;

	//Target to Q1 exit, circle of radius 14.92 cm
	x_test = x_h_q1ex_(vector_jjl, &ii)*m2cm;
	y_test = y_h_q1ex_(vector_jjl, &ii)*m2cm;
	//x_test = x_test + 0.9;
	if( (x_test*x_test + y_test*y_test) > (14.92*14.92) )
		return false;

	//Target to dipole entrance, trapezoid -522.0cm<x<-498.1cm  |y| < -0.1924*x-19.24
	//x_test = x_h_dent_(vector_jjl, &ii)*m2cm;
	//y_test = y_h_dent_(vector_jjl, &ii)*m2cm;
	//if( (x_test<-522.0) || (x_test>-498.1) || fabs(y_test) > fabs(-0.1924*x_test-19.24) )
	//	return false;

	//Target to dipole exit, trapezoid -46.19cm<x<46.19cm  |y| < -0.0161*x+12.5
	x_test = x_h_dext_(vector_jjl, &ii)*m2cm;
	y_test = y_h_dext_(vector_jjl, &ii)*m2cm;
	if( (x_test<-46.19) || (x_test>46.19) || fabs(y_test) > fabs(-0.0161*x_test+12.5) )
		return false;

	//Target to Q3 entrance, circle of radius 30.0 cm
	x_test = x_h_q3en_(vector_jjl, &ii)*m2cm;
	y_test = y_h_q3en_(vector_jjl, &ii)*m2cm;
	if( (x_test*x_test + y_test*y_test) > (30.0*30.0) )
		return false;

	//Target to Q3 exit, circle of radius 30.0 cm  -> 28.0cm
	x_test = x_h_q3ex_(vector_jjl, &ii)*m2cm;
	y_test = y_h_q3ex_(vector_jjl, &ii)*m2cm;
	//x_test = (x_test - 1.0) / (28.0);
	//y_test = y_test / (30.0);
	if( (x_test*x_test + y_test*y_test) > (30.0*30.0) )
		return false;

	/////////////////////////////////////////////////////////////
	/* If we reach this point, it means the test was succesful */
	float x_fp     = x_h_fp_(vector_jjl,&ii);
	float theta_fp = t_h_fp_(vector_jjl,&ii);
	float y_fp     = y_h_fp_(vector_jjl,&ii);
	float phi_fp   = p_h_fp_(vector_jjl,&ii);

	//reset the vector and return it back to the caller
	pV5[0] = (double)x_fp;
	pV5[1] = (double)theta_fp;
	pV5[2] = (double)y_fp;
	pV5[3] = (double)phi_fp;
	//pV5[4] = (double)delta_fp;  // delta is not change

	return true;
}


bool TransportLeftHRS(double* pV5)
{
	float vector_jjl[]={pV5[0],pV5[1],pV5[2],pV5[3],pV5[4]};
	float x_test, y_test;
	int ii=5;

	//Target to Q1 exit, circle of radius 14.92 cm
	x_test = x_e_q1ex_(vector_jjl, &ii)*m2cm;
	y_test = y_e_q1ex_(vector_jjl, &ii)*m2cm;
	//x_test = x_test - 0.9;
	if( (x_test*x_test + y_test*y_test) > (14.92*14.92) )
		return false;

	//Target to dipole entrance, trapezoid 522.0cm>x>498.1cm  |y| < -0.1924*x-19.24
	//x_test = x_e_dent_(vector_jjl, &ii)*m2cm;
	//y_test = y_e_dent_(vector_jjl, &ii)*m2cm;
	//if( (x_test>522.0) || (x_test<498.1) || fabs(y_test) > fabs(-0.1924*x_test-19.24) )
	//	return false;

	//Target to dipole exit, trapezoid -46.19cm<x<46.19cm  |y| < -0.0161*x+12.5
	x_test = x_e_dext_(vector_jjl, &ii)*m2cm;
	y_test = y_e_dext_(vector_jjl, &ii)*m2cm;
	//cout<<"dipole_exit:(x,y)=\t"<<x_test<<"\t "<<y_test<<endl;
	if( fabs(x_test)>46.19 || fabs(y_test) > fabs(-0.0161*x_test+12.5) )
		return false;

	//Target to Q3 entrance, circle of radius 30.0 cm
	x_test = x_e_q3en_(vector_jjl, &ii)*m2cm;
	y_test = y_e_q3en_(vector_jjl, &ii)*m2cm;
	if( (x_test*x_test + y_test*y_test) > (30.0*30.0) )
		return false;

	//Target to Q3 exit, circle of radius 30.0 cm  -> 28.0cm
	x_test = x_e_q3ex_(vector_jjl, &ii)*m2cm;
	y_test = y_e_q3ex_(vector_jjl, &ii)*m2cm;
	//x_test = (x_test + 1.0) / (28.0);
	//y_test = y_test / (30.0);
	if( (x_test*x_test + y_test*y_test) > (30.0*30.0) )
		return false;

	/////////////////////////////////////////////////////////////
	/* If we reach this point, it means the test was succesful */
	float x_fp     = x_e_fp_(vector_jjl,&ii);
	float theta_fp = t_e_fp_(vector_jjl,&ii);
	float y_fp     = y_e_fp_(vector_jjl,&ii);
	float phi_fp   = p_e_fp_(vector_jjl,&ii);

	//reset the vector and return it back to the caller
	pV5[0] = (double)x_fp;
	pV5[1] = (double)theta_fp;
	pV5[2] = (double)y_fp;
	pV5[3] = (double)phi_fp;
	//pV5[4] = (double)delta_fp;  // delta is not change

	return true;
}

//////////////////////////////////////////////////////////////////////////////////////
/* the detail of  JJLerose vector in focus plan is*/
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

/////////////////////////////////////////////////////////////////////////////////////
//This routine is translated from the fortran, which is always return NAN
//it was used to check if there is anything wrong with the fortran routine
float L6_txfit(float *x, int *m)
{
	m+=0; //to avoid warning like "unused variable"
	float avdat = -0.3438139E-02;
	float xmin[] = { -0.56268E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00 };
	float xmax[] = { 0.49365E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00 };
	float scale[] = { 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00 };
	float coeff[] = {  -0.24954679E-02, 0.90844929E-01,-0.95373677E-03,  0.};
	
	for(int i=0;i<5;i++)
	{
		if (fabs(xmin[i]-xmax[i])<1.0E-8) break;
		scale[i]=2./(xmax[i]-xmin[i]);
	}

	//normalize variables between -1 and +1
	float x1 =1.+(x[0]-xmax[0])*scale[0];
	float x11 = x1;
	float x12 = x11*x1;

	float txfit = avdat + coeff[0] + coeff[1]*x11 + coeff[2]*x12;
	return txfit;

}

/////////////////////////////////////////////////////////////////////////////////////
void ReconstructRightHRS(double *pV5)
{
	float vector_jjl[]={pV5[0],pV5[1],pV5[2],pV5[3],pV5[4]};
	int ii=5,jj=1;
	/* Orthogonalize theta as JJL asks*/
	vector_jjl[1]   = vector_jjl[1] - r_txfit_(vector_jjl,&jj);
	float x_or      = vector_jjl[4];
	float delta_rec = r_delta_(vector_jjl,&ii);
	float theta_rec = r_theta_(vector_jjl,&ii);
	float phi_rec   = r_phi_(vector_jjl,&ii);
	float y_rec     = r_y00_(vector_jjl,&ii);

	//reset the vector and return it back to the caller
	pV5[0] = (double)x_or;
	pV5[1] = (double)theta_rec;
	pV5[2] = (double)y_rec;
	pV5[3] = (double)phi_rec;
	pV5[4] = (double)delta_rec;
}

void ReconstructLeftHRS(double *pV5)
{
	//using the right arm to reconstruct left arm, just for debugging
	//in order to call right arm routines, need to flip y, phi 
	float vector_jjl[]={pV5[0],pV5[1],-pV5[2],-pV5[3],pV5[4]};
	int ii=5,jj=1;

	pV5[1]   = pV5[1] - r_txfit_(vector_jjl,&jj);
	float x_or_r      = pV5[4];
	float delta_rec_r = r_delta_(vector_jjl,&ii);
	float theta_rec_r = r_theta_(vector_jjl,&ii);
	float phi_rec_r   = -r_phi_(vector_jjl,&ii);		//flip y
	float y_rec_r     = -r_y00_(vector_jjl,&ii);		//flip phi

	pV5[0] = (double)x_or_r;
	pV5[1] = (double)theta_rec_r;
	pV5[2] = (double)y_rec_r;
	pV5[3] = (double)phi_rec_r;
	pV5[4] = (double)delta_rec_r;

#ifndef DEBUG_HRS_TRANSPORT
	return;
#endif

	for(int i=0;i<5;i++) vector_jjl[i]=pV5[i];
	//By Jixie @ 20110303  the left arm routine has problems
	//once finish debuging, I will remove the above

	//I found that fortran routine l_txfit_ return something wrong
	vector_jjl[1] = vector_jjl[1] - L6_txfit(vector_jjl,&jj);
	//vector_jjl[1]   = vector_jjl[1] - l_txfit_(vector_jjl,&jj);
	float x_or      = vector_jjl[4];
	float delta_rec = l_delta_(vector_jjl,&ii);
	float theta_rec = l_theta_(vector_jjl,&ii);
	float phi_rec   = l_phi_(vector_jjl,&ii);
	float y_rec     = l_y00_(vector_jjl,&ii);

	//reset the vector and return it back to the caller
	pV5[0] = (double)x_or;
	pV5[1] = (double)theta_rec;
	pV5[2] = (double)y_rec;
	pV5[3] = (double)phi_rec;
	pV5[4] = (double)delta_rec;

#ifdef DEBUG_HRS_TRANSPORT
	if( Global_Debug_Level >= 2 )
	{
		cout<<"Std Right rec: delta="<<setw(10)<<delta_rec_r<<"  theta="<<setw(10)<<theta_rec_r
			<<"  y="<<setw(10)<<y_rec_r<<"  phi="<<setw(10)<<phi_rec_r<<endl;
		cout<<"Std Left  rec: delta="<<setw(10)<<delta_rec<<"  theta="<<setw(10)<<theta_rec
			<<"  y="<<setw(10)<<y_rec<<"  phi="<<setw(10)<<phi_rec<<endl;
		if( Global_Debug_Level >= 3 ) STOP4DEBUG;
	}
#endif

}

///////////////////////////////////////////////////////////////////

void ReconstructRightHRS(double *pV5, double &delta_rec, double &theta_rec, double &phi_rec, double &y_rec)
{
	float vector_jjl[]={pV5[0],pV5[1],pV5[2],pV5[3],pV5[4]};
	int ii=5,jj=1;

	vector_jjl[1] = vector_jjl[1] - r_txfit_(vector_jjl,&jj);
	delta_rec = r_delta_(vector_jjl,&ii);
	theta_rec = r_theta_(vector_jjl,&ii);
	phi_rec   = r_phi_(vector_jjl,&ii);
	y_rec     = r_y00_(vector_jjl,&ii);
}

void ReconstructLeftHRS(double *pV5, double &delta_rec, double &theta_rec, double &phi_rec, double &y_rec)
{
	float vector_jjl[]={pV5[0],pV5[1],pV5[2],pV5[3],pV5[4]};
	int ii=5,jj=1;

	vector_jjl[1] = vector_jjl[1] - l_txfit_(vector_jjl,&jj);
	delta_rec = l_delta_(vector_jjl,&ii);
	theta_rec = l_theta_(vector_jjl,&ii);
	phi_rec   = l_phi_(vector_jjl,&ii);
	y_rec     = l_y00_(vector_jjl,&ii);
}

