//This file is supposed to declare the forward and backward functions for G2P exp
//including both c++ and fortran routines

#ifndef HRS_G2P_H
#define HRS_G2P_H
//By Jixie: how to read the snake directive file to get the aperture size?
//ep5 to ep8 is end planes inside the septum
//the shape and size is defined by the last 6 parameters in the following lines
//y, -478.5, -2.402,none, 104.2,409.7,120.,120.,-120.,-120.
//y, 0., 0.,none,0.,0.,0.,0.,0.,0.
//y, 478.5., -2.402,none,130.5,389.64,120.,120.,-120.,-120.
//y, 500., -2.402,none,0.,0.,0.,0.,0.,0.
//in TCS,             xmin,xmax,ymax1,ymax2,ymin1,ymin2,  
//1 is in bottom face and 2 is in up face


////////////////////////////////////////////////////////////////
//*******new septum with shims, 5.65 central ray, no target field*****M.Huang 2/23/2012********
//the source code can be found in Fwd_r5p65_notg.cpp
//2cm raster, With wrong Bx in septum field
float x_r5p65_ep11_q1en                       (float *x,int m);
float t_r5p65_ep11_q1en                       (float *x,int m);
float y_r5p65_ep11_q1en                       (float *x,int m);
float p_r5p65_ep11_q1en                       (float *x,int m);
float l_r6_h90s_ep11_q1en                     (float *x,int m);
float x_r5p65_ep13_q1                         (float *x,int m);
float t_r5p65_ep13_q1                         (float *x,int m);
float y_r5p65_ep13_q1                         (float *x,int m);
float p_r5p65_ep13_q1                         (float *x,int m);
float l_r6_h90s_ep13_q1                       (float *x,int m);
float x_r5p65_ep14_q1ex                       (float *x,int m);
float t_r5p65_ep14_q1ex                       (float *x,int m);
float y_r5p65_ep14_q1ex                       (float *x,int m);
float p_r5p65_ep14_q1ex                       (float *x,int m);
float l_r6_h90s_ep14_q1ex                     (float *x,int m);
float x_r5p65_ep23_den                        (float *x,int m);
float t_r5p65_ep23_den                        (float *x,int m);
float y_r5p65_ep23_den                        (float *x,int m);
float p_r5p65_ep23_den                        (float *x,int m);
float l_r6_h90s_ep23_den                      (float *x,int m);
float x_r5p65_ep25_dex                        (float *x,int m);
float t_r5p65_ep25_dex                        (float *x,int m);
float y_r5p65_ep25_dex                        (float *x,int m);
float p_r5p65_ep25_dex                        (float *x,int m);
float l_r6_h90s_ep25_dex                      (float *x,int m);
float x_r5p65_ep27_q3en                       (float *x,int m);
float t_r5p65_ep27_q3en                       (float *x,int m);
float y_r5p65_ep27_q3en                       (float *x,int m);
float p_r5p65_ep27_q3en                       (float *x,int m);
float l_r6_h90s_ep27_q3en                     (float *x,int m);
float x_r5p65_ep30_q3ex                       (float *x,int m);
float t_r5p65_ep30_q3ex                       (float *x,int m);
float y_r5p65_ep30_q3ex                       (float *x,int m);
float p_r5p65_ep30_q3ex                       (float *x,int m);
float l_r6_h90s_ep30_q3ex                     (float *x,int m);
float x_r5p65_ep5                             (float *x,int m);
float t_r5p65_ep5                             (float *x,int m);
float y_r5p65_ep5                             (float *x,int m);
float p_r5p65_ep5                             (float *x,int m);
float l_r6_h90s_ep5                           (float *x,int m);
float x_r5p65_ep7                             (float *x,int m);
float t_r5p65_ep7                             (float *x,int m);
float y_r5p65_ep7                             (float *x,int m);
float p_r5p65_ep7                             (float *x,int m);
float l_r6_h90s_ep7                           (float *x,int m);
float x_r5p65_fp                              (float *x,int m);
float t_r5p65_fp                              (float *x,int m);
float y_r5p65_fp                              (float *x,int m);
float p_r5p65_fp                              (float *x,int m);
float l_r6_h90s_fp                            (float *x,int m);

//the source code can be found in Rev_r5p65_notg.cpp
float r5p65_txfit(float *x,int m);
float r5p65_delta(float *x,int m);
float r5p65_theta(float *x,int m);
float r5p65_phi(float *x,int m);
float r5p65_y00(float *x,int m);


////////////////////////////////////////////////////////////////
//*******400016 septum with shims, 5.65 central ray, no target field*********
//*****M.Huang 11/12/2012********
//the source code can be found in Fwd_sl5p65_400016.cpp
//1.5cm raster
float x_sl5p65_400016_ep10_q1en               (float *x,int m);
float t_sl5p65_400016_ep10_q1en               (float *x,int m);
float y_sl5p65_400016_ep10_q1en               (float *x,int m);
float p_sl5p65_400016_ep10_q1en               (float *x,int m);
float l_sl5p65_400016_ep10_q1en               (float *x,int m);
float x_sl5p65_400016_ep13_q1ex               (float *x,int m);
float t_sl5p65_400016_ep13_q1ex               (float *x,int m);
float y_sl5p65_400016_ep13_q1ex               (float *x,int m);
float p_sl5p65_400016_ep13_q1ex               (float *x,int m);
float l_sl5p65_400016_ep13_q1ex               (float *x,int m);
float x_sl5p65_400016_ep20_q2ex               (float *x,int m);
float t_sl5p65_400016_ep20_q2ex               (float *x,int m);
float y_sl5p65_400016_ep20_q2ex               (float *x,int m);
float p_sl5p65_400016_ep20_q2ex               (float *x,int m);
float l_sl5p65_400016_ep20_q2ex               (float *x,int m);
float x_sl5p65_400016_ep23_den                (float *x,int m);
float t_sl5p65_400016_ep23_den                (float *x,int m);
float y_sl5p65_400016_ep23_den                (float *x,int m);
float p_sl5p65_400016_ep23_den                (float *x,int m);
float l_sl5p65_400016_ep23_den                (float *x,int m);
float x_sl5p65_400016_ep24_dex                (float *x,int m);
float t_sl5p65_400016_ep24_dex                (float *x,int m);
float y_sl5p65_400016_ep24_dex                (float *x,int m);
float p_sl5p65_400016_ep24_dex                (float *x,int m);
float l_sl5p65_400016_ep24_dex                (float *x,int m);
float x_sl5p65_400016_ep26_q3en               (float *x,int m);
float t_sl5p65_400016_ep26_q3en               (float *x,int m);
float y_sl5p65_400016_ep26_q3en               (float *x,int m);
float p_sl5p65_400016_ep26_q3en               (float *x,int m);
float l_sl5p65_400016_ep26_q3en               (float *x,int m);
float x_sl5p65_400016_ep29_q3ex               (float *x,int m);
float t_sl5p65_400016_ep29_q3ex               (float *x,int m);
float y_sl5p65_400016_ep29_q3ex               (float *x,int m);
float p_sl5p65_400016_ep29_q3ex               (float *x,int m);
float l_sl5p65_400016_ep29_q3ex               (float *x,int m);
float x_sl5p65_400016_ep5                     (float *x,int m);
float t_sl5p65_400016_ep5                     (float *x,int m);
float y_sl5p65_400016_ep5                     (float *x,int m);
float p_sl5p65_400016_ep5                     (float *x,int m);
float l_sl5p65_400016_ep5                     (float *x,int m);
float x_sl5p65_400016_ep7                     (float *x,int m);
float t_sl5p65_400016_ep7                     (float *x,int m);
float y_sl5p65_400016_ep7                     (float *x,int m);
float p_sl5p65_400016_ep7                     (float *x,int m);
float l_sl5p65_400016_ep7                     (float *x,int m);
float x_sl5p65_400016_fp                      (float *x,int m);
float t_sl5p65_400016_fp                      (float *x,int m);
float y_sl5p65_400016_fp                      (float *x,int m);
float p_sl5p65_400016_fp                      (float *x,int m);
float l_sl5p65_400016_fp                      (float *x,int m);

//the source code can be found in Rev_sl5p65_400016.cpp
float txfit_sl5p65_400016                     (float *x,int m);
float delta_sl5p65_400016                     (float *x,int m);
float theta_sl5p65_400016                     (float *x,int m);
float phi_sl5p65_400016                       (float *x,int m);
float y00_sl5p65_400016                       (float *x,int m);

////////////////////////////////////////////////////////////////

extern "C"
{
	////////////////////////////////////////////////////////////////
	//*******no septum, 5.65 484816 with shim , no target field, By JJL 20121005********
	//The following 45 functions are from file Fwd_r5p65_484816.f 
	float x_r5p65_484816_sepex_(float *value, int* i);
	float t_r5p65_484816_sepex_(float *value, int* i);
	float y_r5p65_484816_sepex_(float *value, int* i);
	float p_r5p65_484816_sepex_(float *value, int* i);
	float l_r5p65_484816_sepex_(float *value, int* i);
	float x_r5p65_484816_q1ent_(float *value, int* i);
	float t_r5p65_484816_q1ent_(float *value, int* i);
	float y_r5p65_484816_q1ent_(float *value, int* i);
	float p_r5p65_484816_q1ent_(float *value, int* i);
	float l_r5p65_484816_q1ent_(float *value, int* i);
	float x_r5p65_484816_q1ext_(float *value, int* i);
	float t_r5p65_484816_q1ext_(float *value, int* i);
	float y_r5p65_484816_q1ext_(float *value, int* i);
	float p_r5p65_484816_q1ext_(float *value, int* i);
	float l_r5p65_484816_q1ext_(float *value, int* i);
	float x_r5p65_484816_q2ext_(float *value, int* i);
	float t_r5p65_484816_q2ext_(float *value, int* i);
	float y_r5p65_484816_q2ext_(float *value, int* i);
	float p_r5p65_484816_q2ext_(float *value, int* i);
	float l_r5p65_484816_q2ext_(float *value, int* i);
	float x_r5p65_484816_den_(float *value, int* i);
	float t_r5p65_484816_den_(float *value, int* i);
	float y_r5p65_484816_den_(float *value, int* i);
	float p_r5p65_484816_den_(float *value, int* i);
	float l_r5p65_484816_den_(float *value, int* i);
	float x_r5p65_484816_dex_(float *value, int* i);
	float t_r5p65_484816_dex_(float *value, int* i);
	float y_r5p65_484816_dex_(float *value, int* i);
	float p_r5p65_484816_dex_(float *value, int* i);
	float l_r5p65_484816_dex_(float *value, int* i);
	float x_r5p65_484816_q3ent_(float *value, int* i);
	float t_r5p65_484816_q3ent_(float *value, int* i);
	float y_r5p65_484816_q3ent_(float *value, int* i);
	float p_r5p65_484816_q3ent_(float *value, int* i);
	float l_r5p65_484816_q3ent_(float *value, int* i);
	float x_r5p65_484816_q3ext_(float *value, int* i);
	float t_r5p65_484816_q3ext_(float *value, int* i);
	float y_r5p65_484816_q3ext_(float *value, int* i);
	float p_r5p65_484816_q3ext_(float *value, int* i);
	float l_r5p65_484816_q3ext_(float *value, int* i);
	float x_r5p65_484816_fp_(float *value, int* i);
	float t_r5p65_484816_fp_(float *value, int* i);
	float y_r5p65_484816_fp_(float *value, int* i);
	float p_r5p65_484816_fp_(float *value, int* i);
	float l_r5p65_484816_fp_(float *value, int* i);

	//The following 5 functions are from file Rev_r5p65_484816.f 
	float txfit_r5p65_484816_(float *value, int* i);
	float delta_r5p65_484816_(float *value, int* i);
	float theta_r5p65_484816_(float *value, int* i);
	float phi_r5p65_484816_(float *value, int* i);
	float y00_r5p65_484816_(float *value, int* i);

	////////////////////////////////////////////////////////////////
	//function for 12.5 degrees without target field
	//The following 25 functions are from file Fwd_r12p5_Min.f 
	float x_r12p5_dent_(float *value, int* i);
	float y_r12p5_dent_(float *value, int* i);
	float t_r12p5_dent_(float *value, int* i);
	float p_r12p5_dent_(float *value, int* i);
	float l_r12p5_dent_(float *value, int* i);
	float x_r12p5_dext_(float *value, int* i);
	float y_r12p5_dext_(float *value, int* i);
	float t_r12p5_dext_(float *value, int* i);
	float p_r12p5_dext_(float *value, int* i);
	float l_r12p5_dext_(float *value, int* i);
	float x_r12p5_fp_(float *value, int* i);
	float y_r12p5_fp_(float *value, int* i);
	float t_r12p5_fp_(float *value, int* i);
	float p_r12p5_fp_(float *value, int* i);
	float l_r12p5_fp_(float *value, int* i);
	float x_r12p5_q1ex_(float *value, int* i);
	float y_r12p5_q1ex_(float *value, int* i);
	float t_r12p5_q1ex_(float *value, int* i);
	float p_r12p5_q1ex_(float *value, int* i);
	float l_r12p5_q1ex_(float *value, int* i);
	float x_r12p5_q3en_(float *value, int* i);
	float y_r12p5_q3en_(float *value, int* i);
	float t_r12p5_q3en_(float *value, int* i);
	float p_r12p5_q3en_(float *value, int* i);
	float l_r12p5_q3en_(float *value, int* i);

	//The following 4 functions are from file Rev_r12p5_Min.f 
	float r12p5_delta_(float *value, int* i);
	float r12p5_theta_(float *value, int* i);
	float r12p5_phi_(float *value, int* i);
	float r12p5_y00_(float *value, int* i);


	////////////////////////////////////////////////////////////////
	//*******old septum without shims, 5.69 central ray, no target field, By John ********
	//function for 5.69 degree HRS, symmetric septa without shims, without target field
	//By JJL, 2 cm raster and wrong Bx in septum field 

	//The following 45 functions are from file Fwd_sr5p69_NoShim_John.f 
	float x_g2_fp_(float *value, int* i);
	float t_g2_fp_(float *value, int* i);
	float y_g2_fp_(float *value, int* i);
	float p_g2_fp_(float *value, int* i);
	float l_g2_fp_(float *value, int* i);
	float x_g2_sen_(float *value, int* i);
	float t_g2_sen_(float *value, int* i);
	float y_g2_sen_(float *value, int* i);
	float p_g2_sen_(float *value, int* i);
	float l_g2_sen_(float *value, int* i);
	float x_g2_sex_(float *value, int* i);
	float t_g2_sex_(float *value, int* i);
	float y_g2_sex_(float *value, int* i);
	float p_g2_sex_(float *value, int* i);
	float l_g2_sex_(float *value, int* i);
	float x_g2_q1en_(float *value, int* i);
	float t_g2_q1en_(float *value, int* i);
	float y_g2_q1en_(float *value, int* i);
	float p_g2_q1en_(float *value, int* i);
	float l_g2_q1en_(float *value, int* i);
	float x_g2_q1ex_(float *value, int* i);
	float t_g2_q1ex_(float *value, int* i);
	float y_g2_q1ex_(float *value, int* i);
	float p_g2_q1ex_(float *value, int* i);
	float l_g2_q1ex_(float *value, int* i);
	float x_g2_dent_(float *value, int* i);
	float t_g2_dent_(float *value, int* i);
	float y_g2_dent_(float *value, int* i);
	float p_g2_dent_(float *value, int* i);
	float l_g2_dent_(float *value, int* i);
	float x_g2_dext_(float *value, int* i);
	float t_g2_dext_(float *value, int* i);
	float y_g2_dext_(float *value, int* i);
	float p_g2_dext_(float *value, int* i);
	float l_g2_dext_(float *value, int* i);
	float x_g2_q3en_(float *value, int* i);
	float t_g2_q3en_(float *value, int* i);
	float y_g2_q3en_(float *value, int* i);
	float p_g2_q3en_(float *value, int* i);
	float l_g2_q3en_(float *value, int* i);
	float x_g2_q3ex_(float *value, int* i);
	float t_g2_q3ex_(float *value, int* i);
	float y_g2_q3ex_(float *value, int* i);
	float p_g2_q3ex_(float *value, int* i);
	float l_g2_q3ex_(float *value, int* i);

	//The following 5 functions are from file Rev_sr5p69_NoShim_John.f 
	float g2_delta_(float *value, int* i);
	float g2_theta_(float *value, int* i);
	float g2_phi_(float *value, int* i);
	float g2_y00_(float *value, int* i);
	float g2_txfit_(float *value, int* i);

	////////////////////////////////////////////////////////////////
	//fitting includes 1, 2, 4GeV electrons, rays from SNAKE. matrix from sl to tg. 5T
	//From file ReV_s2t_H90R1.f
	float s2t_sr5p69_2gev_txfit_(float *value, int* i);
	float s2t_sr5p69_2gev_pyfit_(float *value, int* i);
	float s2t_sr5p69_theta_(float *value, int* i);
	float s2t_sr5p69_phi_(float *value, int* i);
	float s2t_sr5p69_y00_(float *value, int* i);
}

////////////////////////////////////////////////////////////////
//usage:
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
////////////////////////////////////////////////////////////////

bool TransportRightHRS_g2_(double* vector_jjl);		//5.69 deg, septa, no target field setting. From John LeRose
bool TransportLeftHRS_g2_(double* vector_jjl);		//5.69 deg, septa, no target field setting. From John LeRose
void ReconstructRightHRS_g2_(double* vector_jjl);	//5.69 deg, septa, no target field setting. From John LeRose
void ReconstructLeftHRS_g2_(double* vector_jjl);	//5.69 deg, septa, no target field setting. From John LeRose


/////////////////////////////////////////////////////////////////
//TODO:  I hate to use NEW as a key since something keep changing
//NEW for today may not be NEW one month later, especially for g2p experiment
//I have to remove this part when finish debuging
bool TransportRightHRS_12p5_Min_(double* vector_jjl);	//no septa, no target field, simple HRS setting, for check with previous functions
bool TransportLeftHRS_12p5_Min_(double* vector_jjl);	//no septa, no target field, simple HRS setting, for check with previous functions
void ReconstructRightHRS_12p5_Min_(double* vector_jjl);	//no septa, no target field, simple HRS setting, for check with previous functions
void ReconstructLeftHRS_12p5_Min_(double* vector_jjl);	//no septa, no target field, simple HRS setting, for check with previous functions

/////////////////////////////////////////////////////////////////
//this one is from Min but with only 2.0 cm raster and wrong Bx field
bool TransportRightHRS_Shim_484816_WrongBx(double* vector_jjl);		//5.69 deg, 48-48-16 septa with shims, no target field, c++ tuned version
bool TransportLeftHRS_Shim_484816_WrongBx(double* vector_jjl);		//5.69 deg, 48-48-16 septa with shims, no target field, c++ tuned version
void ReconstructRightHRS_Shim_484816_WrongBx(double* vector_jjl);	//5.69 deg, 48-48-16 septa with shims, no target field, c++ tuned version
void ReconstructLeftHRS_Shim_484816_WrongBx(double* vector_jjl);	//5.69 deg, 48-48-16 septa with shims, no target field, c++ tuned version


/////////////////////////////////////////////////////////////////
bool TransportRightHRS_Shim_403216(double* vector_jjl);		//5.69 deg, 40-32-16 septa with shims, no target field, c++ tuned version
bool TransportLeftHRS_Shim_403216(double* vector_jjl);		//5.69 deg, 40-32-16 septa with shims, no target field, c++ tuned version
void ReconstructRightHRS_Shim_403216(double* vector_jjl);	//5.69 deg, 40-32-16 septa with shims, no target field, c++ tuned version
void ReconstructLeftHRS_Shim_403216(double* vector_jjl);	//5.69 deg, 40-32-16 septa with shims, no target field, c++ tuned version


/////////////////////////////////////////////////////////////////
bool TransportRightHRS_Shim_400016(double* vector_jjl);		//5.69 deg, 40-00-16 septa with shims, no target field, c++ tuned version
bool TransportLeftHRS_Shim_400016(double* vector_jjl);		//5.69 deg, 40-00-16 septa with shims, no target field, c++ tuned version
void ReconstructRightHRS_Shim_400016(double* vector_jjl);	//5.69 deg, 40-00-16 septa with shims, no target field, c++ tuned version
void ReconstructLeftHRS_Shim_400016(double* vector_jjl);	//5.69 deg, 40-00-16 septa with shims, no target field, c++ tuned version


/////////////////////////////////////////////////////////////////
//
bool TransportRightHRS_Shim_484816(double* vector_jjl);	    //5.69 deg, 48-48-16 septa with shims, no target field, from JJL
bool TransportLeftHRS_Shim_484816(double* vector_jjl);		//5.69 deg, 48-48-16 septa with shims, no target field, from JJL
void ReconstructRightHRS_Shim_484816(double* vector_jjl);	//5.69 deg, 48-48-16 septa with shims, no target field, from JJL
void ReconstructLeftHRS_Shim_484816(double* vector_jjl);	//5.69 deg, 48-48-16 septa with shims, no target field, from JJL


/////////////////////////////////////////////////////////////////
//Sieve Slit to Target reconstruction
//Need to be updated
//G2P
void ReconstructLeft_S2T_S_H90R1(double *pV5);
void ReconstructRight_S2T_S_H90R1(double *pV5);

//GEP
void ReconstructLeft_S2T_S_H07R1(double *pV5);
void ReconstructRight_S2T_S_H07R1(double *pV5);


/////////////////////////////////////////////////////////////////

#endif
