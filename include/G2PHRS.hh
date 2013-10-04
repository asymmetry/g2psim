// -*- C++ -*-

/* class G2PHRS
 * Interface class of HRSTrans package.
 * It provides HRS transport functions.
 * G2PProcBase classes will call Forward() to get focus plane kinematics to focus plane,
 * will call Backward() to get target plane kinematics.
 *
 * The definition of variables and the list of available models can be found in the comments in the body.
 */

// History:
//   Jan 2013, C. Gu, First public version.
//   Feb 2013, C. Gu, Add correction function.
//   Sep 2013, M. Huang, Add 484816R15 module in the comment block
//

#ifndef G2P_HRS_H
#define G2P_HRS_H

#include "G2PAppBase.hh"

class HRSTransBase;

class G2PHRS : public G2PAppBase {
public:
    G2PHRS(const char* name);
    G2PHRS(int setting);
    ~G2PHRS();

    ///////////////////////////////////////////////////////////////////////////
    // Transport particles through HRS using SNAKE model
    // Use iModelIndex to identify which SNAKE model to be used
    // 1: 484816 with shim, 5.65 deg, 3 cm raster, by JJL
    // 2: 403216 with shim, 5.65 deg, SNAKE Model not ready yet
    // 3: 400016 with shim, 5.65 deg, 3 cm raster, by Min
    // Index > 10 means test
    // 11: 484816 with shim, 5.76 deg, no raster, by Min
    // 12: 484816 with shim, 5.785 deg, 3cm raster, by Min
    // May add more HRS packages later
    ///////////////////////////////////////////////////////////////////////////
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

protected:
    G2PHRS(); // Only for ROOT I/O

    int Configure(EMode mode = kTWOWAY);
    void MakePrefix();

    int fSetting;

    double fHRSAngle;
    double fModelAngle;

    HRSTransBase* pModel;

private:
    static G2PHRS* pG2PHRS;

    ClassDef(G2PHRS, 1)
};

#endif
