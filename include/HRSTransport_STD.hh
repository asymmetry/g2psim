//This file is supposed to declare the forward and backward functions for standard HRS
//the HRS was placed at >= 12.5 degrees and only Q1, Q2, Dipole, and Q3 were used
//including both c++ and fortran routines


#ifndef HRS_STD_H
#define HRS_STD_H
extern "C"
{
	float x_e_q1ex_(float *Value, int* i);
	float t_e_q1ex_(float *Value, int* i);
	float y_e_q1ex_(float *Value, int* i);
	float p_e_q1ex_(float *Value, int* i);
	float l_e_q1ex_(float *Value, int* i);

	float x_e_dent_(float *Value, int* i);
	float t_e_dent_(float *Value, int* i);
	float y_e_dent_(float *Value, int* i);
	float p_e_dent_(float *Value, int* i);
	float l_e_dent_(float *Value, int* i);

	float x_e_dext_(float *Value, int* i);
	float t_e_dext_(float *Value, int* i);
	float y_e_dext_(float *Value, int* i);
	float p_e_dext_(float *Value, int* i);
	float l_e_dext_(float *Value, int* i);

	float x_e_q3en_(float *Value, int* i);
	float t_e_q3en_(float *Value, int* i);
	float y_e_q3en_(float *Value, int* i);
	float p_e_q3en_(float *Value, int* i);
	float l_e_q3en_(float *Value, int* i);

	float x_e_q3ex_(float *Value, int* i);
	float t_e_q3ex_(float *Value, int* i);
	float y_e_q3ex_(float *Value, int* i);
	float p_e_q3ex_(float *Value, int* i);
	float l_e_q3ex_(float *Value, int* i);

	float x_e_fp_(float *Value, int* i);
	float t_e_fp_(float *Value, int* i);
	float y_e_fp_(float *Value, int* i);
	float p_e_fp_(float *Value, int* i);
	float l_e_fp_(float *Value, int* i);

	float x_h_q1ex_(float *Value, int* i);
	float t_h_q1ex_(float *Value, int* i);
	float y_h_q1ex_(float *Value, int* i);
	float p_h_q1ex_(float *Value, int* i);
	float l_h_q1ex_(float *Value, int* i);

	float x_h_dent_(float *Value, int* i);
	float t_h_dent_(float *Value, int* i);
	float y_h_dent_(float *Value, int* i);
	float p_h_dent_(float *Value, int* i);
	float l_h_dent_(float *Value, int* i);

	float x_h_dext_(float *Value, int* i);
	float t_h_dext_(float *Value, int* i);
	float y_h_dext_(float *Value, int* i);
	float p_h_dext_(float *Value, int* i);
	float l_h_dext_(float *Value, int* i);

	float x_h_q3en_(float *Value, int* i);
	float t_h_q3en_(float *Value, int* i);
	float y_h_q3en_(float *Value, int* i);
	float p_h_q3en_(float *Value, int* i);
	float l_h_q3en_(float *Value, int* i);

	float x_h_q3ex_(float *Value, int* i);
	float t_h_q3ex_(float *Value, int* i);
	float y_h_q3ex_(float *Value, int* i);
	float p_h_q3ex_(float *Value, int* i);
	float l_h_q3ex_(float *Value, int* i);

	float x_h_fp_(float *Value, int* i);
	float t_h_fp_(float *Value, int* i);
	float y_h_fp_(float *Value, int* i);
	float p_h_fp_(float *Value, int* i);
	float l_h_fp_(float *Value, int* i);


	//reverse 
	float r_txfit_(float *Value, int* i);
	float r_delta_(float *Value, int* i);
	float r_theta_(float *Value, int* i);
	float r_phi_(float *Value, int* i);
	float r_y00_(float *Value, int* i);

	float l_txfit_(float *Value, int* i);
	float l_delta_(float *Value, int* i);
	float l_theta_(float *Value, int* i);
	float l_phi_(float *Value, int* i);
	float l_y00_(float *Value, int* i);
}

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
bool TransportRightHRS(double* vector_jjl);		//For Stardard HRS
bool TransportLeftHRS(double* vector_jjl);		//For Stardard HRS

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
void ReconstructRightHRS(double* vector_jjl);		//For Stardard HRS
void ReconstructLeftHRS(double* vector_jjl);		//For Stardard HRS
void ReconstructRightHRS(double* vector_jjl, double &delta_rec, double &theta_rec, double &phi_rec, double &y_rec);
void ReconstructLeftHRS(double* vector_jjl, double &delta_rec, double &theta_rec, double &phi_rec, double &y_rec);

#endif


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
