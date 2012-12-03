//This file is supposed to declare the forward and backward functions for GDH exp (E97110)
//In this experiment,HRS was at 6 degrees and a septum was used
//including both c++ and fortran routines

#ifndef HRS_GDH_H
#define HRS_GDH_H
extern "C"
{
	float x_sl_ep3_(float *Value, int* i);
	float t_sl_ep3_(float *Value, int* i);
	float y_sl_ep3_(float *Value, int* i);
	float p_sl_ep3_(float *Value, int* i);
	float l_sl_ep3_(float *Value, int* i);

	float x_sl_ep4_(float *Value, int* i);
	float t_sl_ep4_(float *Value, int* i);
	float y_sl_ep4_(float *Value, int* i);
	float p_sl_ep4_(float *Value, int* i);
	float l_sl_ep4_(float *Value, int* i);

	float x_sl_ep5_(float *Value, int* i);
	float t_sl_ep5_(float *Value, int* i);
	float y_sl_ep5_(float *Value, int* i);
	float p_sl_ep5_(float *Value, int* i);
	float l_sl_ep5_(float *Value, int* i);

	float x_sl_ep6_(float *Value, int* i);
	float t_sl_ep6_(float *Value, int* i);
	float y_sl_ep6_(float *Value, int* i);
	float p_sl_ep6_(float *Value, int* i);
	float l_sl_ep6_(float *Value, int* i);

	float x_sl_ep7_(float *Value, int* i);
	float t_sl_ep7_(float *Value, int* i);
	float y_sl_ep7_(float *Value, int* i);
	float p_sl_ep7_(float *Value, int* i);
	float l_sl_ep7_(float *Value, int* i);

	float x_sl_q1ex_(float *Value, int* i);
	float t_sl_q1ex_(float *Value, int* i);
	float y_sl_q1ex_(float *Value, int* i);
	float p_sl_q1ex_(float *Value, int* i);
	float l_sl_q1ex_(float *Value, int* i);

	float x_sl_dent_(float *Value, int* i);
	float t_sl_dent_(float *Value, int* i);
	float y_sl_dent_(float *Value, int* i);
	float p_sl_dent_(float *Value, int* i);
	float l_sl_dent_(float *Value, int* i);

	float x_sl_dext_(float *Value, int* i);
	float t_sl_dext_(float *Value, int* i);
	float y_sl_dext_(float *Value, int* i);
	float p_sl_dext_(float *Value, int* i);
	float l_sl_dext_(float *Value, int* i);

	float x_sl_q3en_(float *Value, int* i);
	float t_sl_q3en_(float *Value, int* i);
	float y_sl_q3en_(float *Value, int* i);
	float p_sl_q3en_(float *Value, int* i);
	float l_sl_q3en_(float *Value, int* i);

	float x_sl_q3ex_(float *Value, int* i);
	float t_sl_q3ex_(float *Value, int* i);
	float y_sl_q3ex_(float *Value, int* i);
	float p_sl_q3ex_(float *Value, int* i);
	float l_sl_q3ex_(float *Value, int* i);

	float x_sl_fp_(float *Value, int* i);
	float t_sl_fp_(float *Value, int* i);
	float y_sl_fp_(float *Value, int* i);
	float p_sl_fp_(float *Value, int* i);
	float l_sl_fp_(float *Value, int* i);

	float l6_txfit_(float *Value, int* i);
	float l6_delta_(float *Value, int* i);
	float l6_theta_(float *Value, int* i);
	float l6_phi_(float *Value, int* i);
	float l6_y00_(float *Value, int* i);

	////////////////////////////////////////////////

	float x_sr_ep3_(float *Value, int* i);
	float t_sr_ep3_(float *Value, int* i);
	float y_sr_ep3_(float *Value, int* i);
	float p_sr_ep3_(float *Value, int* i);
	float l_sr_ep3_(float *Value, int* i);

	float x_sr_ep4_(float *Value, int* i);
	float t_sr_ep4_(float *Value, int* i);
	float y_sr_ep4_(float *Value, int* i);
	float p_sr_ep4_(float *Value, int* i);
	float l_sr_ep4_(float *Value, int* i);

	float x_sr_ep5_(float *Value, int* i);
	float t_sr_ep5_(float *Value, int* i);
	float y_sr_ep5_(float *Value, int* i);
	float p_sr_ep5_(float *Value, int* i);
	float l_sr_ep5_(float *Value, int* i);

	float x_sr_ep6_(float *Value, int* i);
	float t_sr_ep6_(float *Value, int* i);
	float y_sr_ep6_(float *Value, int* i);
	float p_sr_ep6_(float *Value, int* i);
	float l_sr_ep6_(float *Value, int* i);

	float x_sr_ep7_(float *Value, int* i);
	float t_sr_ep7_(float *Value, int* i);
	float y_sr_ep7_(float *Value, int* i);
	float p_sr_ep7_(float *Value, int* i);
	float l_sr_ep7_(float *Value, int* i);

	float x_sr_q1ex_(float *Value, int* i);
	float t_sr_q1ex_(float *Value, int* i);
	float y_sr_q1ex_(float *Value, int* i);
	float p_sr_q1ex_(float *Value, int* i);
	float l_sr_q1ex_(float *Value, int* i);

	float x_sr_dent_(float *Value, int* i);
	float t_sr_dent_(float *Value, int* i);
	float y_sr_dent_(float *Value, int* i);
	float p_sr_dent_(float *Value, int* i);
	float l_sr_dent_(float *Value, int* i);

	float x_sr_dext_(float *Value, int* i);
	float t_sr_dext_(float *Value, int* i);
	float y_sr_dext_(float *Value, int* i);
	float p_sr_dext_(float *Value, int* i);
	float l_sr_dext_(float *Value, int* i);

	float x_sr_q3en_(float *Value, int* i);
	float t_sr_q3en_(float *Value, int* i);
	float y_sr_q3en_(float *Value, int* i);
	float p_sr_q3en_(float *Value, int* i);
	float l_sr_q3en_(float *Value, int* i);

	float x_sr_q3ex_(float *Value, int* i);
	float t_sr_q3ex_(float *Value, int* i);
	float y_sr_q3ex_(float *Value, int* i);
	float p_sr_q3ex_(float *Value, int* i);
	float l_sr_q3ex_(float *Value, int* i);

	float x_sr_fp_(float *Value, int* i);
	float t_sr_fp_(float *Value, int* i);
	float y_sr_fp_(float *Value, int* i);
	float p_sr_fp_(float *Value, int* i);
	float l_sr_fp_(float *Value, int* i);

	//reverse float s
	float r6_txfit_(float *Value, int* i);
	float r6_delta_(float *Value, int* i);
	float r6_theta_(float *Value, int* i);
	float r6_phi_(float *Value, int* i);
	float r6_y00_(float *Value, int* i);

	float x_sr6_v2_ep3_(float *Value, int* i);
	float t_sr6_v2_ep3_(float *Value, int* i);
	float y_sr6_v2_ep3_(float *Value, int* i);
	float p_sr6_v2_ep3_(float *Value, int* i);
	float l_sr6_v2_ep3_(float *Value, int* i);

	float x_sr6_v2_ep4_(float *Value, int* i);
	float t_sr6_v2_ep4_(float *Value, int* i);
	float y_sr6_v2_ep4_(float *Value, int* i);
	float p_sr6_v2_ep4_(float *Value, int* i);
	float l_sr6_v2_ep4_(float *Value, int* i);

	float x_sr6_v2_ep5_(float *Value, int* i);
	float t_sr6_v2_ep5_(float *Value, int* i);
	float y_sr6_v2_ep5_(float *Value, int* i);
	float p_sr6_v2_ep5_(float *Value, int* i);
	float l_sr6_v2_ep5_(float *Value, int* i);

	float x_sr6_v2_ep6_(float *Value, int* i);
	float t_sr6_v2_ep6_(float *Value, int* i);
	float y_sr6_v2_ep6_(float *Value, int* i);
	float p_sr6_v2_ep6_(float *Value, int* i);
	float l_sr6_v2_ep6_(float *Value, int* i);

	float x_sr6_v2_ep7_(float *Value, int* i);
	float t_sr6_v2_ep7_(float *Value, int* i);
	float y_sr6_v2_ep7_(float *Value, int* i);
	float p_sr6_v2_ep7_(float *Value, int* i);
	float l_sr6_v2_ep7_(float *Value, int* i);

	float x_sr6_v2_q1ex_(float *Value, int* i);
	float t_sr6_v2_q1ex_(float *Value, int* i);
	float y_sr6_v2_q1ex_(float *Value, int* i);
	float p_sr6_v2_q1ex_(float *Value, int* i);
	float l_sr6_v2_q1ex_(float *Value, int* i);

	float x_sr6_v2_dent_(float *Value, int* i);
	float t_sr6_v2_dent_(float *Value, int* i);
	float y_sr6_v2_dent_(float *Value, int* i);
	float p_sr6_v2_dent_(float *Value, int* i);
	float l_sr6_v2_dent_(float *Value, int* i);

	float x_sr6_v2_dext_(float *Value, int* i);
	float t_sr6_v2_dext_(float *Value, int* i);
	float y_sr6_v2_dext_(float *Value, int* i);
	float p_sr6_v2_dext_(float *Value, int* i);
	float l_sr6_v2_dext_(float *Value, int* i);

	float x_sr6_v2_q3en_(float *Value, int* i);
	float t_sr6_v2_q3en_(float *Value, int* i);
	float y_sr6_v2_q3en_(float *Value, int* i);
	float p_sr6_v2_q3en_(float *Value, int* i);
	float l_sr6_v2_q3en_(float *Value, int* i);

	float x_sr6_v2_q3ex_(float *Value, int* i);
	float t_sr6_v2_q3ex_(float *Value, int* i);
	float y_sr6_v2_q3ex_(float *Value, int* i);
	float p_sr6_v2_q3ex_(float *Value, int* i);
	float l_sr6_v2_q3ex_(float *Value, int* i);

	float x_sr6_v2_fp_(float *Value, int* i);
	float t_sr6_v2_fp_(float *Value, int* i);
	float y_sr6_v2_fp_(float *Value, int* i);
	float p_sr6_v2_fp_(float *Value, int* i);
	float l_sr6_v2_fp_(float *Value, int* i);

	//reverse float s
	float r6_v2_txfit_(float *Value, int* i);
	float r6_v2_delta_(float *Value, int* i);
	float r6_v2_theta_(float *Value, int* i);
	float r6_v2_phi_(float *Value, int* i);
	float r6_v2_y00_(float *Value, int* i);

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

bool TransportRightHRS6_LargeX0(double* vector_jjl);	//For GDH exp, E97110, with larger X0
bool TransportLeftHRS6_LargeX0(double* vector_jjl);	//For GDH exp, E97110, with larger X0
bool TransportRightHRS6(double* vector_jjl);	//For GDH exp, E97110
bool TransportLeftHRS6(double* vector_jjl);		//For GDH exp, E97110


void ReconstructRightHRS6_LargeX0(double* vector_jjl);	//For GDH exp, E97110, with larger X0
void ReconstructLeftHRS6_LargeX0(double* vector_jjl);	//For GDH exp, E97110, with larger X0
void ReconstructRightHRS6(double* vector_jjl);		//For GDH exp, E97110
void ReconstructLeftHRS6(double* vector_jjl);		//For GDH exp, E97110

void ReconstructRightHRS6_LargeX0(double* vector_jjl, double &delta_rec, double &theta_rec, double &phi_rec, double &y_rec);
void ReconstructLeftHRS6_LargeX0(double* vector_jjl, double &delta_rec, double &theta_rec, double &phi_rec, double &y_rec);
void ReconstructRightHRS6(double* vector_jjl, double &delta_rec, double &theta_rec, double &phi_rec, double &y_rec);
void ReconstructLeftHRS6(double* vector_jjl, double &delta_rec, double &theta_rec, double &phi_rec, double &y_rec);


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
