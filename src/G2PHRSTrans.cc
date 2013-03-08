// This file defines a class HRSTransport.
// This class is used in G2PSim class as HRS model.
// The definition of variables and the list of available models can be found in
//+the comments in the body
// The active model is chosen during constructing.
// 
// History:
//   Jan 2013, C. Gu, First public version.
//   Feb 2013, C. Gu, Add correction function.
//   Feb 2013, J.X. Zhang, Add VC support
//

#include <cstdlib>
#include <cmath>
#include <cstring>
#include <map>

#include "TROOT.h"
#include "TObject.h"
#include "TError.h"

#include "HRSTransBase.hh"
#include "G2PTrans400016/G2PTrans400016.hh"
#include "G2PTrans484816/G2PTrans484816.hh"
#include "G2PTrans484816R00/G2PTrans484816R00.hh"

#include "G2PAppsBase.hh"
#include "G2PGlobals.hh"
#include "G2PRunBase.hh"

#include "G2PHRSTrans.hh"

using namespace std;

const double kDEG = 3.14159265358979323846/180.0;

G2PHRSTrans* G2PHRSTrans::pG2PHRSTrans = NULL;

G2PHRSTrans::G2PHRSTrans() :
    iSetting(0), fHRSAngle(5.767*kDEG), fModelAngle(5.767*kDEG),
    pModel(NULL)
{
    // Nothing to do
}

G2PHRSTrans::G2PHRSTrans(const char* name) :
    iSetting(0), fHRSAngle(5.767*kDEG), fModelAngle(5.767*kDEG),
    pModel(NULL)
{
    if (pG2PHRSTrans) {
        Error("G2PHRSTrans()", "Only one instance of G2PHRSTrans allowed.");
        MakeZombie();
        return;
    }
    pG2PHRSTrans = this;

    map<string, int> model_map;
    model_map["484816"] = 1;
    model_map["400016"] = 3;
    model_map["484816R00"] = 11;

    iSetting = model_map[name];
}

G2PHRSTrans::G2PHRSTrans(int setting) :
    iSetting(setting), fHRSAngle(5.767*kDEG), fModelAngle(5.767*kDEG),
    pModel(NULL)
{
    if (pG2PHRSTrans) {
        Error("G2PHRSTrans()", "Only one instance of G2PHRSTrans allowed.");
        MakeZombie();
        return;
    }
    pG2PHRSTrans = this;
}

G2PHRSTrans::~G2PHRSTrans()
{
    if (pG2PHRSTrans==this) pG2PHRSTrans = NULL;
}

////////////////////////////////////////////////////////////////////////////////
// Transport particles through HRS using SNAKE model
// Use iModelIndex to identify which SNAKE model to be used
// 1: 484816 with shim, 5.65 deg, 3 cm raster, by JJL
// 2: 403216 with shim, 5.65 deg, SNAKE Model not ready yet
// 3: 400016 with shim, 5.65 deg, 3 cm raster, by Min
// Index > 10 means test
// 11: 484816 with shim, 5.76 deg, no raster, by Min
// May add more HRS packages later
////////////////////////////////////////////////////////////////////////////////

G2PAppsBase::EStatus G2PHRSTrans::Init()
{
    static const char* const here = "Init()";

    if (G2PAppsBase::Init()) return fStatus;

    fHRSAngle = gG2PRun->GetHRSAngle();

    switch (iSetting) {
    case 1:
        pModel = new G2PTrans484816();
        break;
    case 3:
        pModel = new G2PTrans400016();
        break;
    case 11:
        pModel = new G2PTrans484816R00();
        break;
    default:
        Error(here, "Cannot initialize, invalid setting.");
        return (fStatus = kINITERROR);
        break;
    }

    fModelAngle = pModel->GetAngle();

    if (fDebug>0) Info(here, "Model angle is %10.3e." , fModelAngle/kDEG);

    return (fStatus = kOK);
}

bool G2PHRSTrans::Forward(const double* V5_tg, double* V5_fp)
{
    static const char* const here = "Forward()";

    // Definition of variables
    // V5_tg = {x_tg, theta_tg, y_tg, phi_tg, delta@tg};
    // V5_fp = {x_fg, theta_fp, y_fp, phi_fp, delta@tg};
    // delta does not change
    double V5[5];
    
    V5[0] = V5_tg[0];
    V5[1] = tan(V5_tg[1]);
    V5[2] = V5_tg[2];
    V5[3] = tan(V5_tg[3]);
    V5[4] = V5_tg[4];

    bool bGoodParticle=false;

    if (fHRSAngle>0) {
        //pModel->CoordsCorrection(fHRSAngle-fModelAngle, V5);
        bGoodParticle = pModel->TransLeftHRS(V5);
        //pModel->FPCorrLeft(V5_tg, V5);
    }
    else {
        //pModel->CoordsCorrection(fHRSAngle+fModelAngle, V5);
        bGoodParticle = pModel->TransRightHRS(V5);
        //pModel->FPCorrRight(V5_tg, V5);
    }

    V5_fp[0] = V5[0];
    V5_fp[1] = atan(V5[1]);
    V5_fp[2] = V5[2];
    V5_fp[3] = atan(V5[3]);
    V5_fp[4] = V5[4];

    if (fDebug>1) Info(here, "%10.3e %10.3e %10.3e %10.3e %10.3e", V5_fp[0], V5_fp[1], V5_fp[2], V5_fp[3], V5_fp[4]);

    return bGoodParticle;
}

bool G2PHRSTrans::Backward(const double* V5_fp, double* V5_tg)
{
    static const char* const here = "Backward()";

    // Definition of variables
    // V5_fp = {x_fp, theta_fp, y_fp, phi_fp, x_tg};
    // V5_tg = {x_tg, theta_tg, y_tg, phi_tg, delta@tg};
    // x_tg does not change
    
    double V5[5];
    
    V5[0] = V5_fp[0];
    V5[1] = tan(V5_fp[1]);
    V5[2] = V5_fp[2];
    V5[3] = tan(V5_fp[3]);
    V5[4] = V5_fp[4];
    
    if (fHRSAngle>0) {
        pModel->ReconLeftHRS(V5);
        //pModel->CoordsCorrection(fModelAngle-fHRSAngle, V5);
    }
    else {
        pModel->ReconRightHRS(V5);
        //pModel->CoordsCorrection(-fModelAngle-fHRSAngle, V5);
    }

    V5_tg[0] = V5[0];
    V5_tg[1] = atan(V5[1]);
    V5_tg[2] = V5[2];
    V5_tg[3] = atan(V5[3]);
    V5_tg[4] = V5[4];

    bool bGoodParticle = false;

    if (V5_tg[4]<1.0) bGoodParticle = true;

    if (fDebug>1) Info(here, "%10.3e %10.3e %10.3e %10.3e %10.3e", V5_tg[0], V5_tg[1], V5_tg[2], V5_tg[3], V5_tg[4]);
    
    return bGoodParticle;
}

ClassImp(G2PHRSTrans)
