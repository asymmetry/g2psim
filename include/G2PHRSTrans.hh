// This file defines a class G2PHRSTrans.
// This class is used in G2PSim class as HRS model.
// The definition of variables and the list of available models can be found in
//+the comments in the body
// The active model is chosen during constructing.
// 
// History:
//   Jan 2013, C. Gu, First public version.
//   Feb 2013, C. Gu, Add correction function.
//

#ifndef G2P_HRSTRANS_H
#define G2P_HRSTRANS_H

#include "G2PAppBase.hh"

class HRSTransBase;

class G2PHRSTrans : public G2PAppBase
{
public:
	G2PHRSTrans(const char* name);
    G2PHRSTrans(int setting);
	~G2PHRSTrans();

    ///////////////////////////////////////////////////////////////////////////
	// Transport particles through HRS using SNAKE model
	// Use iModelIndex to identify which SNAKE model to be used
	// 1: 484816 with shim, 5.65 deg, 3 cm raster, by JJL 
	// 2: 403216 with shim, 5.65 deg, SNAKE Model not ready yet 
	// 3: 400016 with shim, 5.65 deg, 3 cm raster, by Min
	// Index > 10 means test
	// 11: 484816 with shim, 5.76 deg, no raster, by Min
	// May add more HRS packages later
	///////////////////////////////////////////////////////////////////////////
	int Init();
    int Begin();

	///////////////////////////////////////////////////////////////////////////
	// Definition of variables
	// forward:
	// V5_tg = {x_tg, theta_tg, y_tg, phi_tg, delta@tg};
	// V5_fp = {x_fg, theta_fp, y_fp, phi_fp, delta@tg};
	// delta does not change
	///////////////////////////////////////////////////////////////////////////
    bool Forward(const double* V5_tg, double* V5_fp);

    ///////////////////////////////////////////////////////////////////////////
	// backward:
	// V5_fp = {x_fp, theta_fp, y_fp, phi_fp, x_tg};
	// V5_tg = {x_tg, theta_tg, y_tg, phi_tg, delta@tg};
	// x_tg does not change
	///////////////////////////////////////////////////////////////////////////
    bool Backward(const double* V5_fp, double* V5_tg);

    static G2PHRSTrans* GetInstance() { return pG2PHRSTrans; }

protected:
    G2PHRSTrans(); // Only for ROOT I/O

    int iSetting;

    double fHRSAngle;
	double fModelAngle;

    HRSTransBase* pModel;

private:
    static G2PHRSTrans* pG2PHRSTrans;

    ClassDef(G2PHRSTrans, 1)
};

#endif
