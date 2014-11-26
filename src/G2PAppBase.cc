// -*- C++ -*-

/* class G2PAppBase
 * Abstract base class for g2p simulation tools.
 * It provides fundamental functions like coordinates transport.
 * No instance allowed for this class.
 * Many functions are modified from THaAnalysisObject. Thanks to O. Hansen.
 */

// History:
//   Mar 2013, C. Gu, First public version.
//   Sep 2013, C. Gu, Add configure functions.
//   Nov 2014, C. Gu, Set random seed in G2PRun class.
//

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <map>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"

#include "G2PAppList.hh"
#include "G2PGlobals.hh"
#include "G2PRand.hh"
#include "G2PRun.hh"
#include "G2PVarDef.hh"
#include "G2PVarList.hh"

#include "G2PAppBase.hh"

using namespace std;

static const double kLT[4] = { -1.001135e+00, -3.313373e-01, -4.290819e-02, +4.470852e-03};
static const double kLY[4] = { -8.060915e-03, +1.071977e-03, +9.019102e-04, -3.239615e-04};
static const double kLP[4] = { -2.861912e-03, -2.469069e-03, +8.427172e-03, +2.274635e-03};
static const double kRT[4] = { -1.004600e+00, -3.349200e-01, -4.078700e-02, +0.000000e+00};
static const double kRY[4] = { -5.157400e-03, +2.642400e-04, +2.234600e-03, +0.000000e+00};
static const double kRP[4] = { -7.435500e-04, -2.130200e-03, +1.195000e-03, +0.000000e+00};

G2PRand *G2PAppBase::pRand = G2PRand::GetInstance();

G2PAppBase::G2PAppBase() : fPrefix(NULL), fStatus(kNOTINIT), fIsInit(false), fIsSetup(false), fDebug(0), fPriority(0)
{
    fConfigIsSet.clear();
}

G2PAppBase::~G2PAppBase()
{
    delete[] fPrefix;
    fPrefix = NULL;

    fConfigIsSet.clear();
}

int G2PAppBase::Init()
{
    static const char *const here = "Init()";

    if (fDebug > 1)
        Info(here, "Initializing ...");

    if (IsZombie())
        return (fStatus = kNOTINIT);

    EStatus status = kOK;

    MakePrefix();

    return (fStatus = status);
}

int G2PAppBase::Begin()
{
    static const char *const here = "Begin()";

    if (fDebug > 1)
        Info(here, "Beginning ...");

    EStatus status = kOK;

    if (Configure(kTWOWAY))
        status = kINITERROR;

    return (fStatus = status);
}

int G2PAppBase::End()
{
    //static const char* const here = "End()";

    return (Configure(kWRITE));
}

void G2PAppBase::Clear(Option_t *option)
{
    // Default does nothing

    return;
}

G2PAppBase::EStatus G2PAppBase::Status() const
{
    return fStatus;
}

bool G2PAppBase::IsInit() const
{
    return (fStatus == kOK);
}

int G2PAppBase::GetDebugLevel() const
{
    return fDebug;
}

int G2PAppBase::GetPriority() const
{
    return fPriority;
}

void G2PAppBase::SetDebugLevel(int level)
{
    fDebug = level;
}

void G2PAppBase::TCS2HCS(double x_tr, double y_tr, double z_tr, double angle, double &x_lab, double &y_lab, double &z_lab)
{
    // Position transform function from TCS to HCS

    static const char *const here = "TCS2HCS()";

    double cosang = cos(angle);
    double sinang = sin(angle);

    x_lab = y_tr * cosang + z_tr * sinang;
    y_lab = -x_tr;
    z_lab = z_tr * cosang - y_tr * sinang;

    if (fDebug > 3)
        Info(here, "%10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e\n", x_tr, y_tr, z_tr, x_lab, y_lab, z_lab);
}

void G2PAppBase::TCS2HCS(double t_tr, double p_tr, double angle, double &t_lab, double &p_lab)
{
    // Angle transform function from TCS to HCS

    static const char *const here = "TCS2HCS()";

    double x = tan(t_tr);
    double y = tan(p_tr);
    double z = 1.0;
    TCS2HCS(x, y, z, angle, x, y, z);
    t_lab = acos(z / sqrt(x * x + y * y + z * z));
    p_lab = atan2(y, x);

    if (fDebug > 3)
        Info(here, "%10.3e %10.3e -> %10.3e %10.3e\n", t_tr, p_tr, t_lab, p_lab);
}

void G2PAppBase::HCS2TCS(double x_lab, double y_lab, double z_lab, double angle, double &x_tr, double &y_tr, double &z_tr)
{
    // Position transform function from HCS to TCS

    static const char *const here = "HCS2TCS()";

    double cosang = cos(angle);
    double sinang = sin(angle);

    x_tr = -y_lab;
    y_tr = x_lab * cosang - z_lab * sinang;
    z_tr = x_lab * sinang + z_lab * cosang;

    if (fDebug > 3)
        Info(here, "%10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e\n", x_lab, y_lab, z_lab, x_tr, y_tr, z_tr);
}

void G2PAppBase::HCS2TCS(double t_lab, double p_lab, double angle, double &t_tr, double &p_tr)
{
    // Angle transform function from HCS to TCS

    static const char *const here = "HCS2TCS()";

    double x = sin(t_lab) * cos(p_lab);
    double y = sin(t_lab) * sin(p_lab);
    double z = cos(t_lab);
    HCS2TCS(x, y, z, angle, x, y, z);
    t_tr = atan2(x, z);
    p_tr = atan2(y, z);

    if (fDebug > 3)
        Info(here, "%10.3e %10.3e -> %10.3e %10.3e", t_lab, p_lab, t_tr, p_tr);
}

void G2PAppBase::Project(double x, double y, double z, double zout, double t, double p, double &xout, double &yout)
{
    // Project along z direction

    static const char *const here = "Project()";

    double xsave = x;
    double ysave = y;

    xout = xsave + (zout - z) * tan(t);
    yout = ysave + (zout - z) * tan(p);

    if (fDebug > 2)
        Info(here, "%10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e", xsave, ysave, z, xout, yout, zout);
}

void G2PAppBase::TRCS2FCS(const double *V5_tr, double angle, double *V5_fp)
{
    static const char *const here = "TRCS2FCS()";

    double x_tr = V5_tr[0];
    double t_tr = tan(V5_tr[1]);
    double y_tr = V5_tr[2];
    double p_tr = tan(V5_tr[3]);

    double x = x_tr;

    double x1 = x;
    double x2 = x1 * x;
    double x3 = x2 * x;

    const double *tMat;
    const double *yMat;
    const double *pMat;

    if (angle > 0) {
        tMat = kLT;
        yMat = kLY;
        pMat = kLP;
    } else {
        tMat = kRT;
        yMat = kRY;
        pMat = kRP;
    }

    double t_mat = tMat[0] + tMat[1] * x1 + tMat[2] * x2 + tMat[3] * x3;
    double y_mat = yMat[0] + yMat[1] * x1 + yMat[2] * x2 + yMat[3] * x3;
    double p_mat = pMat[0] + pMat[1] * x1 + pMat[2] * x2 + pMat[3] * x3;

    double tan_rho_0 = tMat[0];
    double cos_rho_0 = 1.0 / sqrt(1.0 + tan_rho_0 * tan_rho_0);

    double t_det = (t_tr - tan_rho_0) / (1 + t_tr * tan_rho_0);
    double p_det = p_tr * (1.0 - t_det * tan_rho_0) * cos_rho_0;

    double tan_rho = t_mat;
    double cos_rho = 1.0 / sqrt(1.0 + tan_rho * tan_rho);

    double y = y_tr - y_mat;
    double t = (t_det + tan_rho) / (1.0 - t_det * tan_rho);
    double p = (p_det - p_mat) / ((1.0 - t_det * tan_rho) * cos_rho);

    V5_fp[0] = x;
    V5_fp[1] = atan(t);
    V5_fp[2] = y;
    V5_fp[3] = atan(p);

    if (fDebug > 3)
        Info(here, "%10.3e %10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e %10.3e", V5_tr[0], V5_tr[1], V5_tr[2], V5_tr[3], V5_fp[0], V5_fp[1], V5_fp[2], V5_fp[3]);
}

void G2PAppBase::FCS2TRCS(const double *V5_fp, double angle, double *V5_tr)
{
    static const char *const here = "FCS2TRCS()";

    double x = V5_fp[0];
    double t = tan(V5_fp[1]);
    double y = V5_fp[2];
    double p = tan(V5_fp[3]);

    double x1 = x;
    double x2 = x1 * x;
    double x3 = x2 * x;

    const double *tMat;
    const double *yMat;
    const double *pMat;

    if (angle > 0) {
        tMat = kLT;
        yMat = kLY;
        pMat = kLP;
    } else {
        tMat = kRT;
        yMat = kRY;
        pMat = kRP;
    }

    double t_mat = tMat[0] + tMat[1] * x1 + tMat[2] * x2 + tMat[3] * x3;
    double y_mat = yMat[0] + yMat[1] * x1 + yMat[2] * x2 + yMat[3] * x3;
    double p_mat = pMat[0] + pMat[1] * x1 + pMat[2] * x2 + pMat[3] * x3;

    double tan_rho = t_mat;
    double cos_rho = 1.0 / sqrt(1.0 + tan_rho * tan_rho);

    double x_tr = x;
    double y_tr = y + y_mat;
    double t_det = (t - tan_rho) / (1.0 + t * tan_rho);
    double p_det = p * (1.0 - t_det * tan_rho) * cos_rho + p_mat;

    double tan_rho_0 = tMat[0];
    double cos_rho_0 = 1.0 / sqrt(1.0 + tan_rho_0 * tan_rho_0);

    double t_tr = (t_det + tan_rho_0) / (1.0 - t_det * tan_rho_0);
    double p_tr = p_det / (cos_rho_0 * (1.0 - t_det * tan_rho_0));

    V5_tr[0] = x_tr;
    V5_tr[1] = atan(t_tr);
    V5_tr[2] = y_tr;
    V5_tr[3] = atan(p_tr);

    if (fDebug > 3)
        Info(here, "%10.3e %10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e %10.3e", V5_fp[0], V5_fp[1], V5_fp[2], V5_fp[3], V5_tr[0], V5_tr[1], V5_tr[2], V5_tr[3]);
}

void G2PAppBase::TRCS2DCS(const double *V5_tr, double angle, double *V5_det)
{
    static const char *const here = "TRCS2DCS()";

    const double *tMat;

    if (angle > 0)
        tMat = kLT;
    else
        tMat = kRT;

    double tan_rho_0 = tMat[0];
    double cos_rho_0 = 1.0 / sqrt(1.0 + tan_rho_0 * tan_rho_0);

    double x_tr = V5_tr[0];
    double t_tr = tan(V5_tr[1]);
    double y_tr = V5_tr[2];
    double p_tr = tan(V5_tr[3]);

    double x_det = x_tr / (cos_rho_0 * (1 + t_tr * tan_rho_0));
    double y_det = y_tr - tan_rho_0 * cos_rho_0 * p_tr * x_det;
    double t_det = (t_tr - tan_rho_0) / (1 + t_tr * tan_rho_0);
    double p_det = p_tr * (1.0 - t_det * tan_rho_0) * cos_rho_0;

    V5_det[0] = x_det;
    V5_det[1] = atan(t_det);
    V5_det[2] = y_det;
    V5_det[3] = atan(p_det);

    if (fDebug > 3)
        Info(here, "%10.3e %10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e %10.3e", V5_tr[0], V5_tr[1], V5_tr[2], V5_tr[3], V5_det[0], V5_det[1], V5_det[2], V5_det[3]);
}

void G2PAppBase::DCS2TRCS(const double *V5_det, double angle, double *V5_tr)
{
    static const char *const here = "DCS2TRCS()";

    const double *tMat;

    if (angle > 0)
        tMat = kLT;
    else
        tMat = kRT;

    double tan_rho_0 = tMat[0];
    double cos_rho_0 = 1.0 / sqrt(1.0 + tan_rho_0 * tan_rho_0);

    double x_det = V5_det[0];
    double t_det = tan(V5_det[1]);
    double y_det = V5_det[2];
    double p_det = tan(V5_det[3]);

    double t_tr = (t_det + tan_rho_0) / (1.0 - t_det * tan_rho_0);
    double p_tr = p_det / (cos_rho_0 * (1.0 - t_det * tan_rho_0));
    double x_tr = x_det * cos_rho_0 * (1 + t_tr * tan_rho_0);
    double y_tr = y_det + tan_rho_0 * cos_rho_0 * p_tr * x_det;

    V5_tr[0] = x_tr;
    V5_tr[1] = atan(t_tr);
    V5_tr[2] = y_tr;
    V5_tr[3] = atan(p_tr);

    if (fDebug > 3)
        Info(here, "%10.3e %10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e %10.3e", V5_det[0], V5_det[1], V5_det[2], V5_det[3], V5_tr[0], V5_tr[1], V5_tr[2], V5_tr[3]);
}

void G2PAppBase::FCS2DCS(const double *V5_fp, double angle, double *V5_det)
{
    static const char *const here = "FCS2DCS()";

    double V5_tr[5];

    int save = fDebug;
    fDebug = 0;

    FCS2TRCS(V5_fp, angle, V5_tr);
    TRCS2DCS(V5_tr, angle, V5_det);

    fDebug = save;

    if (fDebug > 3)
        Info(here, "%10.3e %10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e %10.3e", V5_fp[0], V5_fp[1], V5_fp[2], V5_fp[3], V5_det[0], V5_det[1], V5_det[2], V5_det[3]);
}

void G2PAppBase::DCS2FCS(const double *V5_det, double angle, double *V5_fp)
{
    static const char *const here = "DCS2FCS()";

    double V5_tr[5];

    int save = fDebug;
    fDebug = 0;

    DCS2TRCS(V5_det, angle, V5_tr);
    TRCS2FCS(V5_tr, angle, V5_fp);

    fDebug = save;

    if (fDebug > 3)
        Info(here, "%10.3e %10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e %10.3e", V5_det[0], V5_det[1], V5_det[2], V5_det[3], V5_fp[0], V5_fp[1], V5_fp[2], V5_fp[3]);
}

int G2PAppBase::ConfigureFromList(const ConfDef *list, EMode mode)
{
    // Load configurations in the list from run manager

    static const char *const here = "LoadConfigFile()";

    if (!gG2PRun) {
        Error(here, "No run manager.");
        return ((mode == kREAD || mode == kTWOWAY) ? kINITERROR : kOK);
    }

    const ConfDef *item = list;

    if (mode == kREAD || mode == kTWOWAY) {
        while (item->name) {
            if (fConfigIsSet.find((unsigned long) item->var) != fConfigIsSet.end()) {
                if (mode == kTWOWAY)
                    gG2PRun->SetConfig(item, fPrefix);
            } else {
                if (gG2PRun->GetConfig(item, fPrefix))
                    fConfigIsSet.insert((unsigned long) item->var);
            }

            item++;
        }
    } else if (mode == kWRITE) {
        while (item->name) {
            gG2PRun->SetConfig(item, fPrefix);
            item++;
        }
    } else
        return kINITERROR;

    return kOK;
}

void G2PAppBase::MakePrefix(const char *basename)
{
    // Set up name prefix for global variables

    delete[] fPrefix;

    if (basename && *basename) {
        fPrefix = new char[strlen(basename) + 3];
        strcpy(fPrefix, basename);
        strcat(fPrefix, ".");
    } else {
        fPrefix = new char[2];
        *fPrefix = 0;
    }
}

ClassImp(G2PAppBase)
