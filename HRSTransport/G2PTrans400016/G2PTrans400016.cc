////////////////////////////////////////////////////////////////////////
// 400016 septum with shims, 5.65 central ray, no target field
// 3cm raster
// by M. Huang 11/12/2012
////////////////////////////////////////////////////////////////////////

#include <cmath>

#include "Fwd_l5p65_400016.h"
#include "Bwd_l5p65_400016.h"

#include "G2PTrans400016.hh"

using namespace S400016;

const float m2cm = 100.0;
const double kDEG = 3.14159265358979323846/180.0;

G2PTrans400016::G2PTrans400016()
    :cModelAngle(5.65*kDEG)
{
    // Nothing to do
}

G2PTrans400016::~G2PTrans400016()
{
    // Nothing to do
}

bool G2PTrans400016::TransLeftHRS(double* pV5)
{
    float vector_jjl[]={pV5[0],pV5[1],pV5[2],pV5[3],pV5[4]};
    int ii=5;

    float x_test, y_test;

    //Target to Septum ep5
    //y, -480., 0.,none,84.0,388.,97.,97.,-97.,-97.  ;84.388.
    x_test = x_l5p65_ep5(vector_jjl, ii)*m2cm;
    y_test = y_l5p65_ep5(vector_jjl, ii)*m2cm;
    if( fabs(x_test)<8.4 || fabs(x_test)>38.8 || fabs(y_test)>9.7 )
        return false;

    //Target to Septum ep7
    //y, 480., 0.,none,84.0,388.,97.,97.,-97.,-97.
    x_test = x_l5p65_ep7(vector_jjl, ii)*m2cm;
    y_test = y_l5p65_ep7(vector_jjl, ii)*m2cm;
    if( fabs(x_test)<8.4 || fabs(x_test)>38.8 || fabs(y_test)>9.7 )
        return false;

    //Target to Q1 en ep10
    //y,-200., 0.,none,149.2,149.2,0.,0.,0.,0.
    x_test = x_l5p65_ep10_q1en(vector_jjl, ii)*m2cm;
    y_test = y_l5p65_ep10_q1en(vector_jjl, ii)*m2cm;
    if( sqrt(x_test*x_test + y_test*y_test) > 14.92 )
        return false;
    
    //Target to Q1 ex ep13
    //y,675., 0. ,none,300.,300.,0.,0.,0.,0.
    x_test = x_l5p65_ep13_q1ex(vector_jjl, ii)*m2cm;
    y_test = y_l5p65_ep13_q1ex(vector_jjl, ii)*m2cm;
    if( sqrt(x_test*x_test + y_test*y_test) > 30. )
        return false;

    //Target to dipole exit, trapezoid -46.19cm<x<46.19cm  |y| < -0.0161*x+12.5 ep24
    //y,0., 0.,none,-461.88,461.88,132.44,117.56,-132.44,-117.56 
    x_test = x_l5p65_ep24_dex(vector_jjl, ii)*m2cm;
    y_test = y_l5p65_ep24_dex(vector_jjl, ii)*m2cm;
    if( (x_test<-46.19) || (x_test>46.19) || fabs(y_test) > fabs(-0.0161*x_test+12.5) )
        return false;
    
    //Target to Q3 entrance, circle of radius 30.0 cm ep26
    //y,915.,0.,none,300.,300.,0.,0.,0.,0.
    x_test = x_l5p65_ep26_q3en(vector_jjl, ii)*m2cm;
    y_test = y_l5p65_ep26_q3en(vector_jjl, ii)*m2cm;
    if( sqrt(x_test*x_test + y_test*y_test) > 30.0 )
        return false;
    
    //Target to Q3 exit, circle of radius 30.0 cm ep29
    //y,575., 0.,none,300.,300.,0.,0.,0.,0.
    x_test = x_l5p65_ep29_q3ex(vector_jjl, ii)*m2cm;
    y_test = y_l5p65_ep29_q3ex(vector_jjl, ii)*m2cm;
    if( sqrt(x_test*x_test + y_test*y_test) > 30.0)
        return false;

    // succesfully reach focus plane
    float x_fp     = x_l5p65_fp(vector_jjl,ii);
    float theta_fp = t_l5p65_fp(vector_jjl,ii);
    float y_fp     = y_l5p65_fp(vector_jjl,ii);
    float phi_fp   = p_l5p65_fp(vector_jjl,ii);

    //reset the vector and return it back to the caller
    pV5[0] = (double)x_fp;
    pV5[1] = (double)theta_fp;
    pV5[2] = (double)y_fp;
    pV5[3] = (double)phi_fp;
    //pV5[4] = (double)delta_fp;  // delta is not change
    
    return true;
}

bool G2PTrans400016::TransRightHRS(double* pV5)
{
    //use right arm routines for left arm before left arm is ready
    //return TransportLeftHRS(pV5);
    pV5[2]*=-1.;
    pV5[3]*=-1.;
    bool bGoodParticle=TransLeftHRS(pV5);
    pV5[2]*=-1.;
    pV5[3]*=-1.;
    return bGoodParticle;
}

void G2PTrans400016::ReconLeftHRS(double* pV5)
{   
    float vector_jjl[]={pV5[0],pV5[1],pV5[2],pV5[3],pV5[4]};
    int ii=5;
    
    vector_jjl[1]   = vector_jjl[1] - txfit_l5p65(vector_jjl,ii);
    
    float x_or      = vector_jjl[4];
    float delta_rec = delta_l5p65(vector_jjl,ii);
    float theta_rec = theta_l5p65(vector_jjl,ii);
    float phi_rec   = phi_l5p65(vector_jjl,ii); 
    float y_rec     = y00_l5p65(vector_jjl,ii);    
    
    //reset the vector and return it back to the caller
    pV5[0] = (double)x_or;
    pV5[1] = (double)theta_rec;
    pV5[2] = (double)y_rec;
    pV5[3] = (double)phi_rec;
    pV5[4] = (double)delta_rec;
}

void G2PTrans400016::ReconRightHRS(double* pV5)
{
    //in order to call right arm routines, need to flip y, phi 
    pV5[2]*=-1;
    pV5[3]*=-1;
    ReconLeftHRS(pV5);
    pV5[2]*=-1;
    pV5[3]*=-1;
} 
