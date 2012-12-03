// ********************************************************************
//
// $Id: HRSTransform_TCSNHCS.cc,v 1.0, 2010/12/26  12:14:12 HRS Exp $
// GEANT4 tag $Name: geant4-09-04 $ 
//
//..............................................................................
#include "HRSTransform_TCSNHCS.hh"

#include "math.h"
#include <iostream>
using namespace std;

//#define TRANSFORM_DEBUG 1

#if defined TRANSFORM_DEBUG && (TRANSFORM_DEBUG>=2)
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#endif

namespace Transform
{
	static const double deg=acos(0.0)/90;

	//transorm from Hall coordinate to tranportation coordinate 
	//assuming the collimator center is at hall center, if it is not, need to do translation
	//all angles are in rad
	void  X_HCS2TCS(double x, double y, double z,double EndPlaneTheta_hall, 
		double &x_tr, double &y_tr, double &z_tr)
	{
		//definition
		//x_tr = -y;
		//y_tr =  x * cosHRS - z * sinHRS;
		//z_tr =  x * sinHRS + z * cosHRS;  
		//Theta_tr=atan(dx_tr/dz_tr)=atan(-y/(x * sinHRS + z * cosHRS))
		//Phi_tr=atan(dy_tr/dz_tr)=atan((x * cosHRS - z * sinHRS)/(x * sinHRS + z * cosHRS))
		double cosHRS=cos(EndPlaneTheta_hall), sinHRS=sin(EndPlaneTheta_hall);
		x_tr = -y;
		y_tr =  x * cosHRS - z * sinHRS;
		z_tr =  x * sinHRS + z * cosHRS; 
	}

	//transorm from Hall coordinate to tranpotrtation coordinate 
	//all angles are in rad
	void  P_HCS2TCS(double Theta_hall, double Phi_hall, double EndPlaneTheta_hall, 
		double &Theta_tr, double &Phi_tr)
	{
		//The definitions:
		//TCS: EndPlane is the XY plane, Y to the horizontal left and X vertical down
		//The Z is the moving dirction, which perpendicular to EndPlane
		//
		//vertical angle, out of plane angle  Theta_tr = atan(dx/dz) 
		//horizontal angle, in plane angle Phi_tr = atan(dy/dz)


#if defined TRANSFORM_DEBUG && (TRANSFORM_DEBUG>=2)
		//Method 1, using the rotation matrix
		//This method works fine but not efficient
		G4RotationMatrix mRotHCS2TCS;
		mRotHCS2TCS.rotateY(-EndPlaneTheta_hall);
		mRotHCS2TCS.rotateZ(90.0*deg);

		G4ThreeVector pV3;
		pV3.setRThetaPhi(1,Theta_hall,Phi_hall); 
		pV3.transform(mRotHCS2TCS);	
		Theta_tr=atan(pV3.x()/pV3.z());
		Phi_tr=atan(pV3.y()/pV3.z());
		G4cout<<"P_HCS2TCS(): debug 1: By rotation: Theta_tr="<<Theta_tr/deg<<"  Phi_tr="<<Phi_tr/deg<<G4endl;

#endif
		//Method 2, by definition
		//x_tr = -y;
		//y_tr =  x * cosHRS - z * sinHRS;
		//z_tr =  x * sinHRS + z * cosHRS;  
		//Theta_tr=atan(dx_tr/dz_tr)=atan(-y/(x * sinHRS + z * cosHRS))
		//Phi_tr=atan(dy_tr/dz_tr)=atan((x * cosHRS - z * sinHRS)/(x * sinHRS + z * cosHRS))
		double cosHRS=cos(EndPlaneTheta_hall), sinHRS=sin(EndPlaneTheta_hall);
		double x=sin(Theta_hall)*cos(Phi_hall);
		double y=sin(Theta_hall)*sin(Phi_hall);
		double z=cos(Theta_hall);
		double x_tr=-y;
		double y_tr=x * cosHRS - z * sinHRS;
		double z_tr=x * sinHRS + z * cosHRS;
		Theta_tr=atan(x_tr/z_tr);
		Phi_tr=atan(y_tr/z_tr);

#ifdef TRANSFORM_DEBUG
		cout<<"P_HCS2TCS() debug 2 By definition: Theta_tr="<<Theta_tr/deg<<"  Phi_tr="<<Phi_tr/deg<<endl;
#endif
	}


	//transorm from tranpotrtation coordinate to hall coordinate 
	//all angles are in rad
	void  X_TCS2HCS(double x_tr, double y_tr, double z_tr,double EndPlaneTheta_hall,
		double &x, double &y, double &z)
	{
		//definition
		//x_tr = -y;
		//y_tr =  x * cosHRS - z * sinHRS;
		//z_tr =  x * sinHRS + z * cosHRS;  
		//Theta_tr=atan(dx_tr/dz_tr)=atan(-y/(x * sinHRS + z * cosHRS))
		//Phi_tr=atan(dy_tr/dz_tr)=atan((x * cosHRS - z * sinHRS)/(x * sinHRS + z * cosHRS))
		double cosHRS=cos(EndPlaneTheta_hall), sinHRS=sin(EndPlaneTheta_hall);
		x =  y_tr * cosHRS + z_tr * sinHRS; 
		y = -x_tr;
		z =  z_tr * cosHRS - y_tr * sinHRS;
	}

	//transorm from transportation coordinate to hall coordinate 
	//all angles are in rad
	void  P_TCS2HCS(double Theta_tr, double Phi_tr, double EndPlaneTheta_hall, 
		double &Theta_hall, double &Phi_hall)
	{		
		//The definitions:
		//TCS: EndPlane is the XY plane, Y to the horizontal left and X vertical down
		//The Z is the moving dirction, which perpendicular to EndPlane
		//
		//vertical angle, out of plane angle  Theta_tr = atan(dx/dz) 
		//horizontal angle, in plane angle Phi_tr = atan(dy/dz)

#if defined TRANSFORM_DEBUG && (TRANSFORM_DEBUG>=2)

		//Method 1: By definition, 
		//Rotation of HCS to TCS is: Rotate about Y axis by HRS angle
		//then switch axis as y_tr_axis=x_axis and x_tr_axis=-y_axis. Therefore
		//x_tr = -y;
		//y_tr =  x * cosHRS - z * sinHRS;
		//z_tr =  x * sinHRS + z * cosHRS;  
		//Theta_tr=atan(dx_tr/dz_tr)=atan(-y/(x * sinHRS + z * cosHRS))
		//Phi_tr=atan(dy_tr/dz_tr)=atan((x * cosHRS - z * sinHRS)/(x * sinHRS + z * cosHRS))

		//==> 	
		//x =  y_tr * cosHRS + z_tr * sinHRS; 
		//y = -x_tr;
		//z =  z_tr * cosHRS - y_tr * sinHRS;

		//==>
		//Theta=acos(z/sqrt(x*x+y*y+z*z))
		//Phi=atan2(y,x)

		double z_tr=1.0,x_tr=tan(Theta_tr),y_tr=tan(Phi_tr); 

		double cosHRS=cos(EndPlaneTheta_hall), sinHRS=sin(EndPlaneTheta_hall);
		double x =  y_tr * cosHRS + z_tr * sinHRS; 
		double y = -x_tr;
		double z =  z_tr * cosHRS - y_tr * sinHRS;

		Theta_hall=acos(z/sqrt(x*x+y*y+z*z));
		Phi_hall=atan2(y,x);

		//only good for left arm, need to +180 for right arm
		if(sin(EndPlaneTheta_hall)<0) Phi_hall=180*deg+Phi_hall;

		cout<<"P_TCS2HCS() Debug 1: Theta_hall="<<Theta_hall/deg<<"  Phi_hall="<<Phi_hall/deg<<endl;	

		//Method 2, by rotation matrix
		G4RotationMatrix mRotHCS2TCS;
		mRotHCS2TCS.rotateX(-EndPlaneTheta_hall);
		mRotHCS2TCS.rotateZ(-90.0*deg);
		G4ThreeVector pV3(tan(Theta_tr),tan(Phi_tr),1.0); 
		pV3.transform(mRotHCS2TCS);	
		Theta_hall=pV3.theta();
		Phi_hall=pV3.phi();
		cout<<"P_TCS2HCS() Debug 2: Theta_hall="<<Theta_hall/deg<<"  Phi_hall="<<Phi_hall/deg<<endl;
#endif
		//Method 3: by calculation or projection //Also in Jixie's G2p Note. Page 29
		//By definition, to make a HCS R vector into TCS, we have to do these: 
		//a)Rotate about Y axis by HRS+Phi_tr angle
		//b)Rotate about new X axis by Theta_tr angle
		//then we have these relation:
		// z=R*cos(Theta_tr)*cos(HRS+Phi_tr)
		// x=R*cos(Theta_tr)*sin(HRS+Phi_tr)
		// y=-R*sin(Theta_tr)
		// Theta_hall=acos(z/R)=...
		// Phi_hall=atan(y/x)=...
		Theta_hall=acos(cos(EndPlaneTheta_hall+Phi_tr)*cos(Theta_tr));
		Phi_hall=atan( -tan(Theta_tr) / sin(EndPlaneTheta_hall+Phi_tr) );
		//this line also work, but return in (-pi,pi]
		//Phi_hall=atan2(sin(-Theta_tr),cos(-Theta_tr)*sin(EndPlaneTheta_hall+Phi_tr));

		//The formula above only good for left arm, 
		//for the right arm, the caller need to change sign for both HRS angle and phi_tr 
		//and I have to add 180 deg to the phi_hall
		if(sin(EndPlaneTheta_hall)<0) Phi_hall=180.*deg-Phi_hall;

#ifdef TRANSFORM_DEBUG
		cout<<"P_TCS2HCS() Debug 3: Theta_hall="<<Theta_hall/deg<<"  Phi_hall="<<Phi_hall/deg<<endl;
#endif
	}


	//project the trajectory back or forth along z
	void Project(double &x,double &y,double &z,double z_drift,double theta,double phi)
	{
		//x = x + theta * z_drift;
		//y = y + phi * z_drift;
		x = x + tan(theta) * z_drift;
		y = y + tan(phi) * z_drift;
		z = z + z_drift;
	}

	//project the trajectory back or forth along z
	void Project(double x,double y,double z,double z_drift,double theta,double phi,
		double &x_out, double &y_out, double &z_out)
	{
		//x = x + theta * z_drift;
		//y = y + phi * z_drift;
		x_out = x + tan(theta) * z_drift;
		y_out = y + tan(phi) * z_drift;
		z_out = z + z_drift;
	}


} //end of namespace Transform

