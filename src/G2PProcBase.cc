// -*- C++ -*-

/* class G2PProcBase
 * Abstract base class for g2p simulation processes.
 * It provides fundamental functions like variable registration.
 * No instance allowed for this class.
 * Derived class must set its own internal variables and register them to the global variable list.
 */

// History:
//   Apr 2013, C. Gu, First public version.
//   Sep 2013, C. Gu, Move DefineVariables() function from G2PAppBase to here.
//   Dec 2014, C. Gu, Add drifting functions and projection functions.
//

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"

#include "G2PAppBase.hh"
#include "G2PAppList.hh"
#include "G2PDrift.hh"
#include "G2PGeoBase.hh"
#include "G2PGlobals.hh"
#include "G2PMaterial.hh"
#include "G2PRand.hh"
#include "G2PVarDef.hh"
#include "G2PVarList.hh"

#include "G2PProcBase.hh"

G2PProcBase::G2PProcBase() : fStage(kREADY), fDefined(false), fHRSMomentum(0.0), pDrift(NULL)
{
    // Nothing to do
}

G2PProcBase::~G2PProcBase()
{
    // Nothing to do
}

int G2PProcBase::Begin()
{
    //static const char* const here = "Begin()";

    if (G2PAppBase::Begin() != 0)
        return (fStatus = kBEGINERROR);

    pDrift = static_cast<G2PDrift *>(gG2PApps->Find("G2PDrift"));

    if (DefineVariables(kDEFINE) == 0)
        fDefined = true;
    else
        fDefined = false;

    return (fStatus = (fDefined ? kOK : kBEGINERROR));
}

G2PProcBase::EStage G2PProcBase::GetStage()
{
    return fStage;
}

void G2PProcBase::SetStage(EStage stage)
{
    fStage = stage;
}

double G2PProcBase::Drift(const char *dir, const double *x, const double *p, double zf, double *xout, double *pout)
{
    // Drift in lab coordinates

    static const char *const here = "Drift()";

    double xsave[3] = {x[0], x[1], x[2]};
    double psave[3] = {p[0], p[1], p[2]};

    G2PDriftCondition stop(x[2], zf);

    double result = pDrift->Drift(dir, x, p, stop, xout, pout);

    if (fDebug > 2)
        Info(here, "%10.3e %10.3e %10.3e %10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e", xsave[0], xsave[1], xsave[2], psave[0], psave[1], psave[2], xout[0], xout[1], xout[2], pout[0], pout[1], pout[2]);

    return result;
}

double G2PProcBase::Drift(const char *dir, const double *V5_tr, double z_tr, double zf_tr, double *V5out_tr)
{
    // Drift in transport coordinates

    static const char *const here = "Drift()";

    double V5save_tr[5] = {V5_tr[0], V5_tr[1], V5_tr[2], V5_tr[3], V5_tr[4]};

    double V5_lab[5];
    TCS2HCS(V5_tr, z_tr, V5_lab);
    double x[3] = {V5_lab[0], V5_lab[2], V5_lab[4]};
    double pp = (1 + V5_tr[4]) * fHRSMomentum;
    double p[3] = {pp * sin(V5_lab[1]) *cos(V5_lab[3]), pp * sin(V5_lab[1]) *sin(V5_lab[3]), pp * cos(V5_lab[1])};

    G2PDriftCondition stop(z_tr, zf_tr, fHRSAngle);

    double xout[3] = {0.0, 0.0, 0.0};
    double pout[3] = {0.0, 0.0, 0.0};
    double result = pDrift->Drift(dir, x, p, stop, xout, pout);

    double temp;
    HCS2TCS(acos(pout[2] / pp), atan2(pout[1], pout[0]), V5out_tr[1], V5out_tr[3]);
    HCS2TCS(xout[0], xout[1], xout[2], V5out_tr[0], V5out_tr[2], temp);
    V5out_tr[4] = V5save_tr[4];

    if (fDebug > 2)
        Info(here, "%10.3e %10.3e %10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e %10.3e %10.3e", V5save_tr[0], V5save_tr[1], V5save_tr[2], V5save_tr[3], z_tr, V5out_tr[0], V5out_tr[1], V5out_tr[2], V5out_tr[3], zf_tr);

    return result;
}

double G2PProcBase::Drift(const char *dir, const double *V5_tr, double z_tr, double zf_lab, double *V5out_tr, double &zout_tr)
{
    // Drift in transport coordinates

    static const char *const here = "Drift()";

    double V5save_tr[5] = {V5_tr[0], V5_tr[1], V5_tr[2], V5_tr[3], V5_tr[4]};

    double V5_lab[5];
    TCS2HCS(V5_tr, z_tr, V5_lab);
    double x[3] = {V5_lab[0], V5_lab[2], V5_lab[4]};
    double pp = (1 + V5_tr[4]) * fHRSMomentum;
    double p[3] = {pp * sin(V5_lab[1]) *cos(V5_lab[3]), pp * sin(V5_lab[1]) *sin(V5_lab[3]), pp * cos(V5_lab[1])};

    G2PDriftCondition stop(x[2], zf_lab);

    double xout[3] = {0.0, 0.0, 0.0};
    double pout[3] = {0.0, 0.0, 0.0};
    double result = pDrift->Drift(dir, x, p, stop, xout, pout);

    HCS2TCS(acos(pout[2] / pp), atan2(pout[1], pout[0]), V5out_tr[1], V5out_tr[3]);
    HCS2TCS(xout[0], xout[1], xout[2], V5out_tr[0], V5out_tr[2], zout_tr);
    V5out_tr[4] = V5save_tr[4];

    if (fDebug > 2)
        Info(here, "%10.3e %10.3e %10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e %10.3e %10.3e", V5save_tr[0], V5save_tr[1], V5save_tr[2], V5save_tr[3], z_tr, V5out_tr[0], V5out_tr[1], V5out_tr[2], V5out_tr[3], zout_tr);

    return result;
}

double G2PProcBase::Drift(const char *dir, const double *x, const double *p, G2PGeoBase *geo, double *xout, double *pout)
{
    // Drift in lab coordinates

    static const char *const here = "Drift()";

    double xsave[3] = {x[0], x[1], x[2]};
    double psave[3] = {p[0], p[1], p[2]};

    G2PDriftCondition stop(geo);

    double result = pDrift->Drift(dir, x, p, stop, xout, pout);

    if (fDebug > 2) {
        Info(here, "%s", geo->GetMaterial() ? geo->GetMaterial()->GetName() : "vacuum");
        Info(here, "%10.3e %10.3e %10.3e %10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e", xsave[0], xsave[1], xsave[2], psave[0], psave[1], psave[2], xout[0], xout[1], xout[2], pout[0], pout[1], pout[2]);
    }

    return result;
}

double G2PProcBase::Drift(const char *dir, const double *V5_tr, double z_tr, G2PGeoBase *geo, double *V5out_tr, double &zout_tr)
{
    // Drift in transport coordinates

    static const char *const here = "Drift()";

    double V5save_tr[5] = {V5_tr[0], V5_tr[1], V5_tr[2], V5_tr[3], V5_tr[4]};

    double V5_lab[5];
    TCS2HCS(V5_tr, z_tr, V5_lab);
    double x[3] = {V5_lab[0], V5_lab[2], V5_lab[4]};
    double pp = (1 + V5_tr[4]) * fHRSMomentum;
    double p[3] = {pp * sin(V5_lab[1]) *cos(V5_lab[3]), pp * sin(V5_lab[1]) *sin(V5_lab[3]), pp * cos(V5_lab[1])};

    G2PDriftCondition stop(geo);

    double xout[3] = {0.0, 0.0, 0.0};
    double pout[3] = {0.0, 0.0, 0.0};
    double result = pDrift->Drift(dir, x, p, stop, xout, pout);

    HCS2TCS(acos(pout[2] / pp), atan2(pout[1], pout[0]), V5out_tr[1], V5out_tr[3]);
    HCS2TCS(xout[0], xout[1], xout[2], V5out_tr[0], V5out_tr[2], zout_tr);
    V5out_tr[4] = V5save_tr[4];

    if (fDebug > 2) {
        Info(here, "%s", geo->GetMaterial() ? geo->GetMaterial()->GetName() : "vacuum");
        Info(here, "%10.3e %10.3e %10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e %10.3e %10.3e", V5save_tr[0], V5save_tr[1], V5save_tr[2], V5save_tr[3], z_tr, V5out_tr[0], V5out_tr[1], V5out_tr[2], V5out_tr[3], zout_tr);
    }

    return result;
}

double G2PProcBase::Project(const double *x, const double *p, double zf, double *xout, double *pout)
{
    // Project along z direction

    static const char *const here = "Project()";

    double xsave[3] = {x[0], x[1], x[2]};

    xout[0] = xsave[0] + (zf - xsave[2]) * p[0] / p[2];
    xout[1] = xsave[1] + (zf - xsave[2]) * p[1] / p[2];
    xout[2] = zf;

    pout[0] = p[0];
    pout[1] = p[1];
    pout[2] = p[2];

    if (fDebug > 2)
        Info(here, "%10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e", xsave[0], xsave[1], xsave[2], xout[0], xout[1], xout[2]);

    double dx[3] = {xout[0] - xsave[0], xout[1] - xsave[1], xout[2] - xsave[2]};

    return sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
}

double G2PProcBase::Project(const double *V5_tr, double z_tr, double zf_tr, double *V5out_tr)
{
    // Project along z direction

    static const char *const here = "Project()";

    V5out_tr[1] = V5_tr[1];
    V5out_tr[3] = V5_tr[3];
    V5out_tr[4] = V5_tr[4];

    double xsave = V5_tr[0];
    double ysave = V5_tr[2];

    V5out_tr[0] = xsave + (zf_tr - z_tr) * tan(V5out_tr[1]);
    V5out_tr[2] = ysave + (zf_tr - z_tr) * tan(V5out_tr[3]);

    if (fDebug > 2)
        Info(here, "%10.3e %10.3e %10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e %10.3e %10.3e", xsave, V5_tr[1], ysave, V5_tr[3], z_tr, V5out_tr[0], V5out_tr[1], V5out_tr[2], V5out_tr[3], zf_tr);

    double dx = V5out_tr[0] - xsave;
    double dy = V5out_tr[2] - ysave;
    double dz = zf_tr - z_tr;

    return sqrt(dx * dx + dy * dy + dz * dz);
}

double G2PProcBase::Project(double x, double y, double z, double zout, double t, double p, double &xout, double &yout)
{
    // Project along z direction

    static const char *const here = "Project()";

    double xsave = x;
    double ysave = y;

    xout = xsave + (zout - z) * tan(t);
    yout = ysave + (zout - z) * tan(p);

    if (fDebug > 2)
        Info(here, "%10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e", xsave, ysave, z, xout, yout, zout);

    double dx = xout - xsave;
    double dy = yout - ysave;
    double dz = zout - z;

    return sqrt(dx * dx + dy * dy + dz * dz);
}

double G2PProcBase:: InterBremsstrahlung(double E, double angle)
{
    // Take from gener_cone MC, modified by R.M
    // Modified by Jie Liu, Nov 25 2014

    // L.Van Hoorebeke, University of Gent, e-mail: Luc.VanHoorebeke@UGent.be
    // This function generates internal radiation the distribution used
    // contains multiple emission effects
    // k: electron momentum (GeV)
    // nu: equivalent radiator length in units radiation length
    //
    //     this is actually half the equivalent radiator length
    //     because 1/2 placed before and 1/2 placed after proton

    static double kMe = 0.510998918e-3; // GeV

    double qsq = 2 * E * E * (1 - cos(angle));
    double alpha = (1. / 137.);
    double bval = 4. / 3.;
    double msq = 2.6112e-7; // mass electron squared (GeV^2)
    // This is the equivalent radiator used for internal bremsstrahlung.
    double nu = (alpha / (bval * 3.141592627)) * (log(qsq / msq) - 1);
    double cut, Ekin, prob, prob_sample, sample;

    // Initialization of lower limit of bremsstrahlung (1 keV)
    cut = 1e-6;
    Ekin = E - kMe;

    // Calculation of probability to have internal radiation effect above 1 keV. *
    prob = 1. - pow(cut / Ekin, nu);
    prob_sample = pRand->Uniform();  // Random sampling

    if (prob_sample > prob)
        return 0.;

    // bremsstrahlung has taken place! Generate photon energy
    sample = pRand->Uniform();

    double result = Ekin * pow(sample * prob + pow(cut / Ekin, nu), 1. / nu);

    if (result > (E - 2 * kMe))
        result = E - 2 * kMe;

    if ((result < 0) || (E < 2 * kMe))
        result = 0;

    return result;
}

int G2PProcBase::Configure(EMode mode)
{
    if ((mode == kREAD || mode == kTWOWAY) && fConfigured)
        return 0;

    if (G2PAppBase::Configure(mode) != 0)
        return -1;

    ConfDef confs[] = {
        {"run.hrs.p0", "HRS Momentum", kDOUBLE, &fHRSMomentum},
        {0}
    };

    return ConfigureFromList(confs, mode);
}

int G2PProcBase::DefineVariables(EMode mode)
{
    if ((mode == kDEFINE) && fDefined)
        return 0;

    return 0;
}

int G2PProcBase::DefineVarsFromList(const VarDef *list, EMode mode) const
{
    return DefineVarsFromList(fPrefix, list, mode);
}

int G2PProcBase::DefineVarsFromList(const char *prefix, const VarDef *list, EMode mode) const
{
    // Add or delete global variables in "list" to the global list

    static const char *const here = "DefineVarsFromList()";

    if (!gG2PVars) {
        Error(here, "No global variable list.");
        return -1;
    }

    if (mode == kDEFINE)
        gG2PVars->DefineVariables(list, prefix);
    else if (mode == kDELETE) {
        const VarDef *item;

        while ((item = list++) && item->name)
            gG2PVars->RemoveName(Form("%s%s", prefix, item->name));
    } else
        return -1;

    return 0;
}

int G2PProcBase::RemoveVariables()
{
    // Default method for removing global variables

    return DefineVariables(kDELETE);
}

ClassImp(G2PProcBase)
