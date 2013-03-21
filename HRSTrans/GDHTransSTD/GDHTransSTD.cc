////////////////////////////////////////////////////////////////////////
// 6 deg with septum, for GDH experiment (E97110)
// small X0
// by J.J. LeRose 10/05/2012
////////////////////////////////////////////////////////////////////////

#include "GDHTransSTD.hh"

#include <iostream>
#include <iomanip>
#include <math.h>

using namespace std;

//#include "GlobalDebuger.hh"
//#define DEBUG_HRS_TRANSPORT 2

//include the header of fortran routines
#include "Bwd_L6_GDH.hh"
#include "Bwd_R6_GDH.hh"
#include "Fwd_L6_GDH.hh"
#include "Fwd_R6_GDH.hh"

const float m2cm = 100.0;
const double kDEG = 3.14159265358979323846/180.0;

GDHTransSTD::GDHTransSTD()
    :cModelAngle(6.0*kDEG)
{
    // Nothing to do
}

GDHTransSTD::~GDHTransSTD()
{
    // Nothing to do
}

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
//6 degree seting, for E94110 and E97110 
/////////////////////////////////////////////////////////////////////////////
bool GDHTransSTD::TransRightHRS(double* pV5)
{
	float vector_jjl[]={pV5[0],pV5[1],pV5[2],pV5[3],pV5[4]};
	float x_test, y_test;
	int ii=5;

	//By Jixie: since all collimators have the same y_min and y_max, I use
	//some variable to represent them for convenience 
	float y_min=-9.9;
	float y_max= 9.9;

	x_test = x_sr_ep3_(vector_jjl, &ii)*m2cm;
	y_test = y_sr_ep3_(vector_jjl, &ii)*m2cm;

	if( (x_test<-14.06) || (x_test>-8.87) || (y_test< y_min) || (y_test>y_max) )
		return false;

	//Target to 1/4 Septum, -17.12cm<x<-10.89cm, -9.9cm<y<9.9cm
	x_test = x_sr_ep4_(vector_jjl, &ii)*m2cm;
	y_test = y_sr_ep4_(vector_jjl, &ii)*m2cm;

	if( (x_test<-17.12) || (x_test>-10.89) || (y_test< y_min) || (y_test>y_max) )
		return false;


	//Target to 1/2 Septum, -21.29cm<x<-13.54cm, -9.9cm<y<9.9cm
	x_test = x_sr_ep5_(vector_jjl, &ii)*m2cm;
	y_test = y_sr_ep5_(vector_jjl, &ii)*m2cm;
	if( (x_test<-21.29) || (x_test>-13.54) || (y_test< y_min) || (y_test>y_max) )
		return false;

	//Target to 3/4 Septum, -26.84cm<x<-16.97cm, -9.9cm<y<9.9cm
	x_test = x_sr_ep6_(vector_jjl, &ii)*m2cm;
	y_test = y_sr_ep6_(vector_jjl, &ii)*m2cm;
	if( (x_test<-26.84) || (x_test>-16.97) || (y_test< y_min) || (y_test>y_max) )
		return false;

	//Target to Septum exit, -34.05cm<x<-21.56cm, -9.9cm<y<9.9cm
	x_test = x_sr_ep7_(vector_jjl, &ii)*m2cm;
	y_test = y_sr_ep7_(vector_jjl, &ii)*m2cm;
	if( (x_test<-34.05) || (x_test>-21.56) || (y_test< y_min) || (y_test>y_max) )
		return false;

	//Target to Q1 exit, circle of radius 14.92 cm
	x_test = x_sr_q1ex_(vector_jjl, &ii)*m2cm;
	y_test = y_sr_q1ex_(vector_jjl, &ii)*m2cm;
	x_test = x_test + 0.9;
	if( (x_test*x_test + y_test*y_test) > (14.92*14.92) )
		return false;

	//Target to dipole entrance, trapezoid -522.0cm<x<-498.1cm  |y| < -0.1924*x-19.24
	x_test = x_sr_dent_(vector_jjl, &ii)*m2cm;
	y_test = y_sr_dent_(vector_jjl, &ii)*m2cm;
	if( (x_test<-522.0) || (x_test>-498.1) || fabs(y_test) > fabs(-0.1924*x_test-19.24) )
		return false;

	//Target to dipole exit, trapezoid -46.19cm<x<46.19cm  |y| < -0.0161*x+12.5
	x_test = x_sr_dext_(vector_jjl, &ii)*m2cm;
	y_test = y_sr_dext_(vector_jjl, &ii)*m2cm;
	if( (x_test<-46.19) || (x_test>46.19) || fabs(y_test) > fabs(-0.0161*x_test+12.5) )
		return false;

	//Target to Q3 entrance, circle of radius 30.0 cm
	x_test = x_sr_q3en_(vector_jjl, &ii)*m2cm;
	y_test = y_sr_q3en_(vector_jjl, &ii)*m2cm;
	if( (x_test*x_test + y_test*y_test) > (30.0*30.0) )
		return false;

	//Target to Q3 exit, circle of radius 30.0 cm  -> 28.0cm
	x_test = x_sr_q3ex_(vector_jjl, &ii)*m2cm;
	y_test = y_sr_q3ex_(vector_jjl, &ii)*m2cm;
	x_test = (x_test - 1.0) / (28.0);
	y_test = y_test / (30.0);
	if( (x_test*x_test + y_test*y_test) > 1.0 )
		return false;

	/////////////////////////////////////////////////////////////
	/* If we reach this point, it means the test was succesful */
	float x_fp     = x_sr_fp_(vector_jjl,&ii);
	float theta_fp = t_sr_fp_(vector_jjl,&ii);
	float y_fp     = y_sr_fp_(vector_jjl,&ii);
	float phi_fp   = p_sr_fp_(vector_jjl,&ii);

	//reset the vector and return it back to the caller
	pV5[0] = (double)x_fp;
	pV5[1] = (double)theta_fp;
	pV5[2] = (double)y_fp;
	pV5[3] = (double)phi_fp;
	//pV5[4] = (double)delta_fp;  // delta is not change

	return true;
}


bool GDHTransSTD::TransLeftHRS(double* pV5)
{
	float vector_jjl[]={pV5[0],pV5[1],pV5[2],pV5[3],pV5[4]};
	float x_test, y_test;
	int ii=5;

	//By Jixie: since all collimators have the same y_min and y_max, I use
	//some variable to represent them for convenience 
	float y_min=-9.9;
	float y_max= 9.9;


	x_test = -x_sl_ep3_(vector_jjl, &ii)*m2cm;
	y_test = -y_sl_ep3_(vector_jjl, &ii)*m2cm;

#ifdef DEBUG_HRS_TRANSPORT
	//in order to call right arm routines, need to flip y, phi 
	float pV5_tg[5]={vector_jjl[0],vector_jjl[1],-vector_jjl[2],-vector_jjl[3],vector_jjl[4]};
	float xx_test,yy_test;
	if(Global_Debug_Level>=4)
	{
		xx_test = x_sr_ep3_(pV5_tg, &ii)*m2cm;
		yy_test = y_sr_ep3_(pV5_tg, &ii)*m2cm;
		cout<<"Left  6deg ep3: x="<<setw(10)<<x_test<<"  y="<<setw(10)<<y_test<<endl;
		cout<<"Right 6deg ep3: x="<<setw(10)<<xx_test<<"  y="<<setw(10)<<yy_test<<endl<<endl;
	}
#endif

	if( (x_test>14.06) || (x_test<8.87) || (y_test< y_min) || (y_test>y_max) )
		return false;

	//Target to 1/4 Septum, 17.12cm>x>10.89cm, -9.9cm<y<9.9cm
	x_test = -x_sl_ep4_(vector_jjl, &ii)*m2cm;
	y_test = -y_sl_ep4_(vector_jjl, &ii)*m2cm;

	if( (x_test>17.12) || (x_test<10.89) || (y_test< y_min) || (y_test>y_max) )
		return false;


	//Target to 1/2 Septum, 21.29cm>x>13.54cm, -9.9cm<y<9.9cm
	x_test = -x_sl_ep5_(vector_jjl, &ii)*m2cm;
	y_test = -y_sl_ep5_(vector_jjl, &ii)*m2cm;
	if( (x_test>21.29) || (x_test<13.54) || (y_test< y_min) || (y_test>y_max) )
		return false;

	//Target to 3/4 Septum, 26.84cm>x>16.97cm, -9.9cm<y<9.9cm
	x_test = -x_sl_ep6_(vector_jjl, &ii)*m2cm;
	y_test = -y_sl_ep6_(vector_jjl, &ii)*m2cm;
	if( (x_test>26.84) || (x_test<16.97) || (y_test< y_min) || (y_test>y_max) )
		return false;

	//Target to Septum exit, 34.05cm>x>21.56cm, -9.9cm<y<9.9cm
	x_test = -x_sl_ep7_(vector_jjl, &ii)*m2cm;
	y_test = -y_sl_ep7_(vector_jjl, &ii)*m2cm;
	if( (x_test>34.05) || (x_test<21.56) || (y_test< y_min) || (y_test>y_max) )
		return false;


	//Target to Q1 exit, circle of radius 14.92 cm
	x_test = -x_sl_q1ex_(vector_jjl, &ii)*m2cm;
	y_test = -y_sl_q1ex_(vector_jjl, &ii)*m2cm;
	x_test = x_test - 0.9;
	if( (x_test*x_test + y_test*y_test) > (14.92*14.92) )
		return false;

	//Target to dipole entrance, trapezoid 522.0cm>x>498.1cm  |y| < -0.1924*x-19.24
	x_test = -x_sl_dent_(vector_jjl, &ii)*m2cm;
	y_test = -y_sl_dent_(vector_jjl, &ii)*m2cm;

#ifdef DEBUG_HRS_TRANSPORT
	if(Global_Debug_Level>=4)
	{
		//in order to call right arm routines, need to flip y, phi 
		//float pV5_tg[5]={vector_jjl[0],vector_jjl[1],-vector_jjl[2],-vector_jjl[3],vector_jjl[4]};
		xx_test = x_sr_dent_(pV5_tg, &ii)*m2cm;
		yy_test = y_sr_dent_(pV5_tg, &ii)*m2cm;
		cout<<"Left  6deg dent: x="<<setw(10)<<x_test<<"  y="<<setw(10)<<y_test<<endl;
		cout<<"Right 6deg dent: x="<<setw(10)<<xx_test<<"  y="<<setw(10)<<yy_test<<endl<<endl;
	}
#endif

	if( (x_test>522.0) || (x_test<498.1) || fabs(y_test) > fabs(0.1924*x_test+19.24) )
		return false;

	//Target to dipole exit, trapezoid -46.19cm<x<46.19cm  |y| < -0.0161*x+12.5
	x_test = -x_sl_dext_(vector_jjl, &ii)*m2cm;
	y_test = -y_sl_dext_(vector_jjl, &ii)*m2cm;
	if( (x_test<-46.19) || (x_test>46.19) || fabs(y_test) > fabs(-0.0161*x_test+12.5) )
		return false;

	//Target to Q3 entrance, circle of radius 30.0 cm
	x_test = -x_sl_q3en_(vector_jjl, &ii)*m2cm;
	y_test = -y_sl_q3en_(vector_jjl, &ii)*m2cm;
	if( (x_test*x_test + y_test*y_test) > (30.0*30.0) )
		return false;

	//Target to Q3 exit, circle of radius 30.0 cm  -> 28.0cm
	x_test = -x_sl_q3ex_(vector_jjl, &ii)*m2cm;
	y_test = -y_sl_q3ex_(vector_jjl, &ii)*m2cm;
	x_test = (x_test + 1.0) / (28.0);
	y_test = y_test / (30.0);
	if( (x_test*x_test + y_test*y_test) > 1.0 )
		return false;

	/////////////////////////////////////////////////////////////
	/* If we reach this point, it means the test was succesful */
	float x_fp     = x_sl_fp_(vector_jjl,&ii);
	float theta_fp = t_sl_fp_(vector_jjl,&ii);
	float y_fp     = y_sl_fp_(vector_jjl,&ii);
	float phi_fp   = p_sl_fp_(vector_jjl,&ii);

#ifdef DEBUG_HRS_TRANSPORT
	if(Global_Debug_Level>=4)
	{
		//just debug if the parameters in focus plane for left and right are in the same convention
		//the conclusion is that after transfer, we need to flip y and phi again
		//in order to call right arm routines, need to flip y, phi 

		//float pV5_tg[5]={vector_jjl[0],vector_jjl[1],-vector_jjl[2],-vector_jjl[3],vector_jjl[4]};
		float x_fp_r     = x_sr_fp_(pV5_tg,&ii);
		float theta_fp_r = t_sr_fp_(pV5_tg,&ii);
		float y_fp_r     = -y_sr_fp_(pV5_tg,&ii);	// need to flip y, phi 
		float phi_fp_r   = -p_sr_fp_(pV5_tg,&ii);	// need to flip y, phi 
		cout<<"Left  6deg fp: x="<<setw(10)<<x_fp<<"  theta="<<setw(10)<<theta_fp
			<<"  y="<<setw(10)<<y_fp<<"  phi="<<setw(10)<<phi_fp<<endl;
		cout<<"Right 6deg fp: x="<<setw(10)<<x_fp_r<<"  theta="<<setw(10)<<theta_fp_r
			<<"  y="<<setw(10)<<y_fp_r<<"  phi="<<setw(10)<<phi_fp_r<<endl<<endl;
	}
#endif

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
void GDHTransSTD::ReconRightHRS(double *pV5)
{
	float vector_jjl[]={pV5[0],pV5[1],pV5[2],pV5[3],pV5[4]};
	int ii=5,jj=1;
	/* Orthogonalize theta as JJL asks*/
	vector_jjl[1]   = vector_jjl[1] - r6_txfit_(vector_jjl,&jj);
	float x_or      = vector_jjl[4];
	float delta_rec = r6_delta_(vector_jjl,&ii);
	float theta_rec = r6_theta_(vector_jjl,&ii);
	float phi_rec   = r6_phi_(vector_jjl,&ii);
	float y_rec     = r6_y00_(vector_jjl,&ii);

	//reset the vector and return it back to the caller
	pV5[0] = (double)x_or;
	pV5[1] = (double)theta_rec;
	pV5[2] = (double)y_rec;
	pV5[3] = (double)phi_rec;
	pV5[4] = (double)delta_rec;
}

void GDHTransSTD::ReconLeftHRS(double *pV5)
{
	//in order to call right arm routines, need to flip y, phi 
	float vector_jjl[]={pV5[0],pV5[1],-pV5[2],-pV5[3],pV5[4]};
	int ii=5,jj=1;

	/* Orthogonalize theta as JJL asks*/
	pV5[1]   = pV5[1] - r6_txfit_(vector_jjl,&jj);
	float x_or_r      = pV5[4];
	float delta_rec_r = r6_delta_(vector_jjl,&ii);
	float theta_rec_r = r6_theta_(vector_jjl,&ii);
	float phi_rec_r   = -r6_phi_(vector_jjl,&ii);		//flip y
	float y_rec_r     = -r6_y00_(vector_jjl,&ii);		//flip phi

	pV5[0] = (double)x_or_r;
	pV5[1] = (double)theta_rec_r;
	pV5[2] = (double)y_rec_r;
	pV5[3] = (double)phi_rec_r;
	pV5[4] = (double)delta_rec_r;

#ifndef DEBUG_HRS_TRANSPORT
	return;
#endif


	//By Jixie @ 20110303  the left arm 6 deg routine has problems
	//once finish debuging, I will move it up
	for(int i=0;i<5;i++) vector_jjl[i]=pV5[i];
	/* Orthogonalize theta as JJL asks*/
	//float tt2   = vector_jjl[1] - r6_txfit_(vector_jjl,&jj);
	//It turns out that sl6_txfit_ has bugs 
	vector_jjl[1]   = vector_jjl[1] - sl6_txfit_(vector_jjl,&jj);
	//vector_jjl[1]   = vector_jjl[1] - SL6_txfit(vector_jjl,&jj);
	float x_or      = vector_jjl[4];
	float delta_rec = sl6_delta_(vector_jjl,&ii);
	float theta_rec = sl6_theta_(vector_jjl,&ii);
	float phi_rec   = sl6_phi_(vector_jjl,&ii);
	float y_rec     = sl6_y00_(vector_jjl,&ii);

	//reset the vector and return it back to the caller
	pV5[0] = (double)x_or;
	pV5[1] = (double)theta_rec;
	pV5[2] = (double)y_rec;
	pV5[3] = (double)phi_rec;
	pV5[4] = (double)delta_rec;

#ifdef DEBUG_HRS_TRANSPORT
	if( Global_Debug_Level >= 2 )
	{
		cout<<"Right 6deg rec: delta="<<setw(10)<<delta_rec_r<<"  theta="<<setw(10)<<theta_rec_r
			<<"  y="<<setw(10)<<y_rec_r<<"  phi="<<setw(10)<<phi_rec_r<<endl;
		cout<<"Left  6deg rec: delta="<<setw(10)<<delta_rec<<"  theta="<<setw(10)<<theta_rec
			<<"  y="<<setw(10)<<y_rec<<"  phi="<<setw(10)<<phi_rec<<endl;
		STOP4DEBUG;
	}
#endif

}

/////////////////////////////////////////////////////////////////////////////////////

