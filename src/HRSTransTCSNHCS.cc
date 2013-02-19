// This file defines a namespace Transform.
// This namespace is used in G2PSim, G2PGun class.
// It has 4 standard transform and 2 projection between Hall Coordinate and
//+Transport Coordinate
// Most of the functions is developed from J.X. Zhang's G2P geant4 simulation
//+package.
// 
// History:
//   Oct 2010, J.X. Zhang, First version in G2P geant4 simulation.
//   Jan 2013, C. Gu, First public version.
//

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "TROOT.h"
#include "TMath.h"

#include "HRSTransTCSNHCS.hh"

//#define TRANSFORM_DEBUG 1

const double cDeg = TMath::Pi()/180.0;

namespace HRSTransTCSNHCS
{
    // Transform from Hall Coordinate to Transport Coordinate 
	// Assuming the center of Transport Coordinate is at hall center
	void X_HCS2TCS(double x, double y, double z, double EndPlaneTheta_hall, double &x_tr, double &y_tr, double &z_tr)
	{
		double cosHRS = cos(EndPlaneTheta_hall);
        double sinHRS = sin(EndPlaneTheta_hall);
		x_tr = -y;
		y_tr =  x*cosHRS-z*sinHRS;
		z_tr =  x*sinHRS+z*cosHRS;

#ifdef TRANSFORM_DEBUG
        printf("Transform: %e\t%e\t%e\n", x_tr, y_tr, z_tr);
#endif
    }
    
    void P_HCS2TCS(double Theta_hall, double Phi_hall, double EndPlaneTheta_hall, double &Theta_tr, double &Phi_tr)
    {       
		// Method 1, using the rotation matrix
		// This method works fine but not efficient
        
		// G4RotationMatrix mRotHCS2TCS;
		// mRotHCS2TCS.rotateY(-EndPlaneTheta_hall);
		// mRotHCS2TCS.rotateZ(90.0*cDeg);
        // G4ThreeVector pV3;
		// pV3.setRThetaPhi(1,Theta_hall,Phi_hall); 
		// pV3.transform(mRotHCS2TCS);	
		// Theta_tr=atan(pV3.x()/pV3.z());
		// Phi_tr=atan(pV3.y()/pV3.z());

// #ifdef TRANSFORM_DEBUG
//         printf("Transform: %e\t%e\n", Theta_tr/cDeg, Phi_tr/cDeg);
// #endif
        
		// Method 2, by definition
		// x_tr = -y;
		// y_tr =  x * cosHRS - z * sinHRS;
		// z_tr =  x * sinHRS + z * cosHRS;  
		// Theta_tr = atan(dx_tr/dz_tr)=atan(-y/(x * sinHRS + z * cosHRS))
		// Phi_tr = atan(dy_tr/dz_tr)=atan((x * cosHRS - z * sinHRS)/(x * sinHRS + z * cosHRS))
        
		double cosHRS = cos(EndPlaneTheta_hall);
        double sinHRS = sin(EndPlaneTheta_hall);
		double x = sin(Theta_hall)*cos(Phi_hall);
		double y = sin(Theta_hall)*sin(Phi_hall);
		double z = cos(Theta_hall);
		double x_tr = -y;
		double y_tr =  x*cosHRS-z*sinHRS;
		double z_tr =  x*sinHRS+z*cosHRS;
		Theta_tr = atan(x_tr/z_tr);
		Phi_tr = atan(y_tr/z_tr);

#ifdef TRANSFORM_DEBUG
        printf("Transform: %e\t%e\n", Theta_tr/cDeg, Phi_tr/cDeg);
#endif
    }
    
	// Transform from Transport Coordinate to Hall Coordinate
	void X_TCS2HCS(double x_tr, double y_tr, double z_tr, double EndPlaneTheta_hall, double &x, double &y, double &z)
	{
		double cosHRS = cos(EndPlaneTheta_hall);
        double sinHRS = sin(EndPlaneTheta_hall);
		x =  y_tr*cosHRS+z_tr*sinHRS; 
		y = -x_tr;
		z =  z_tr*cosHRS-y_tr*sinHRS;

#ifdef TRANSFORM_DEBUG
        printf("Transform: %e\t%e\t%e\n", x, y, z);
#endif
	}
    
	void P_TCS2HCS(double Theta_tr, double Phi_tr, double EndPlaneTheta_hall, double &Theta_hall, double &Phi_hall)
	{
		// Method 1: By definition, 
        // Just invert the HCS2TCS formula.	
		// x =  y_tr * cosHRS + z_tr * sinHRS; 
		// y = -x_tr;
		// z =  z_tr * cosHRS - y_tr * sinHRS;
		// Theta=acos(z/sqrt(x*x+y*y+z*z))
		// Phi=atan2(y,x)

		// double z_tr=1.0, x_tr=tan(Theta_tr), y_tr=tan(Phi_tr);
		// double cosHRS=cos(EndPlaneTheta_hall), sinHRS=sin(EndPlaneTheta_hall);
		// double x =  y_tr * cosHRS + z_tr * sinHRS; 
		// double y = -x_tr;
		// double z =  z_tr * cosHRS - y_tr * sinHRS;
        // Theta_hall=acos(z/sqrt(x*x+y*y+z*z));
		// Phi_hall=atan2(y,x);
        // //only good for left arm, need to +180 for right arm
		// if(sin(EndPlaneTheta_hall)<0) Phi_hall=180*cDeg+Phi_hall;

// #ifdef TRANSFORM_DEBUG
//         printf("Transform: %e\t%e\n", Theta_tr/cDeg, Phi_tr/cDeg);
// #endif

		// Method 2, by rotation matrix
    	// G4RotationMatrix mRotHCS2TCS;
		// mRotHCS2TCS.rotateX(-EndPlaneTheta_hall);
		// mRotHCS2TCS.rotateZ(-90.0*cDeg);
		// G4ThreeVector pV3(tan(Theta_tr),tan(Phi_tr),1.0); 
		// pV3.transform(mRotHCS2TCS);	
		// Theta_hall=pV3.theta();
		// Phi_hall=pV3.phi();

// #ifdef TRANSFORM_DEBUG
//         printf("Transform: %e\t%e\n", Theta_tr/cDeg, Phi_tr/cDeg);
// #endif
        
		// Method 3: by calculation or projection]
        // Also in Jixie's G2p Note. Page 29
		// By definition, to make a HCS R vector into TCS, we have to do these: 
		// a)Rotate about Y axis by HRS+Phi_tr angle
		// b)Rotate about new X axis by Theta_tr angle
		// then we have these relation:
		// z=R*cos(Theta_tr)*cos(HRS+Phi_tr)
		// x=R*cos(Theta_tr)*sin(HRS+Phi_tr)
		// y=-R*sin(Theta_tr)
		// Theta_hall=acos(z/R)=...
		// Phi_hall=atan(y/x)=...
        
		Theta_hall = acos(cos(EndPlaneTheta_hall+Phi_tr)*cos(Theta_tr));
		Phi_hall = atan(-tan(Theta_tr)/sin(EndPlaneTheta_hall+Phi_tr));
        // The formula above is only good for left arm. 
		// For right arm, the caller need to change sign for HRS angle
		// And add 180 deg to phi_hall
        // By Chao: it was 180.*cDeg-Phi_hall, need check
		if (sin(EndPlaneTheta_hall)<0)
            Phi_hall = 180.0*cDeg+Phi_hall; 

#ifdef TRANSFORM_DEBUG
        printf("Transform: %e\t%e\n", Theta_tr/cDeg, Phi_tr/cDeg);
#endif
    }

    // Project the trajectory along z
	void Project(double &x, double &y, double &z, double z_drift, double theta, double phi)
	{
		x = x+tan(theta)*z_drift;
		y = y+tan(phi)*z_drift;
		z = z+z_drift;
	}

	void Project(double x, double y, double z, double z_drift, double theta, double phi, double &x_out, double &y_out, double &z_out)
	{
        x_out = x;
        y_out = y;
        z_out = z;

        Project(x_out, y_out, z_out, z_drift, theta, phi);
	}
} // End of namespace Transform
