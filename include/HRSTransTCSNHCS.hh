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

#ifndef HRS_TRANSFORM_H
#define HRS_TRANSFORM_H

namespace Transform
{
	// Transform from Hall Coordinate to Transport Coordinate 
	// Assuming the center of Transport Coordinate is at hall center
	void X_HCS2TCS(double x, double y, double z, double EndPlaneTheta_hall, double &x_tr, double &y_tr, double &z_tr);
	void P_HCS2TCS(double Theta_hall, double Phi_hall, double EndPlaneTheta_hall, double &Theta_tr, double &Phi_tr);

	// Transform from Transport Coordinate to Hall Coordinate
	void X_TCS2HCS(double x_tr, double y_tr, double z_tr, double EndPlaneTheta_hall, double &x, double &y, double &z);
	void P_TCS2HCS(double Theta_tr, double Phi_tr, double EndPlaneTheta_hall, double &Theta_hall, double &Phi_hall);

    // Project the trajectory along z
	void Project(double &x, double &y, double &z, double z_drift, double theta, double phi);
	void Project(double x, double y, double z, double z_drift, double theta, double phi, double &x_out, double &y_out, double &z_out);
}

#endif
