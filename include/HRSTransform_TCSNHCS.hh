// ********************************************************************
//
// $Id: HRSTransform_TCSNHCS.hh,v 1.0, 2010/12/26  12:14:12 HRS Exp $
// GEANT4 tag $Name: geant4-09-04 $ 
//
//..............................................................................

#ifndef Transform_TCSNHCS_h
#define Transform_TCSNHCS_h 1

namespace Transform
{
	//transorm from Hall coordinate to tranportation coordinate 
	//assuming the collimator center is at hall center, if it is not, need to do translation
	//all angles are in rad
	void X_HCS2TCS(double x, double y, double z,double EndPlaneTheta_hall, 
		double &x_tr, double &y_tr, double &z_tr);
	//transorm from Hall coordinate to tranpotrtation coordinate 
	//all angles are in rad
	void P_HCS2TCS(double Theta_hall, double Phi_hall, double EndPlaneTheta_hall, 
		double &Theta_tr, double &Phi_tr);

	//transorm from tranpotrtation coordinate to hall coordinate 
	//all angles are in rad
	void X_TCS2HCS(double x_tr, double y_tr, double z_tr,double EndPlaneTheta_hall,
		double &x, double &y, double &z);

	//all angles are in rad
	void P_TCS2HCS(double Theta_tr, double Phi_tr, double EndPlaneTheta_hall, 
		double &Theta_hall, double &Phi_hall);


	//project the trajectory back or forth along z
	void Project(double &x,double &y,double &z,double z_drift,double theta,double phi);
	void Project(double x,double y,double z,double z_drift,double theta,double phi,
		double &x_out, double &y_out, double &z_out);
}

#endif


