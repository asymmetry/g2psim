// -*- C++ -*-

/* class G2PBPM
 * Calculate beam readout at BPM and target using kinematics from event generator.
 * Transport functions defined in G2PBPMTrans is used in this class.
 * Orbits are defined in G2PBPMTrans.
 *
 * Variables ending with "_bpm" are defined in a special coordinates.
 * TransBPM2Lab() will transform it to lab coordinates.
 * In output, these variables are labeled as "b_".
 */

// History:
//   Mar 2013, C. Gu, First public version.
//   Apr 2013, C. Gu, Add Pengjia's fitting result.
//   Jul 2013, C. Gu, Treat optics (no field) case specially.
//   Sep 2013, C. Gu, Rewrite it as a G2PProcBase class.
//

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"

#include "G2PAppBase.hh"
#include "G2PAppList.hh"
#include "G2PBPMTrans.hh"
#include "G2PDrift.hh"
#include "G2PGlobals.hh"
#include "G2PProcBase.hh"
#include "G2PRand.hh"
#include "G2PVar.hh"
#include "G2PVarDef.hh"
#include "G2PVarList.hh"

#include "G2PBPM.hh"

//#define USE_SINGLE_BPM 1

G2PBPM* G2PBPM::pG2PBPM = NULL;

G2PBPM::G2PBPM() :
fBeamEnergy(2.254), fFieldRatio(0.0), fBPMAX(0.0), fBPMAY(0.0), fBPMBX(0.0), fBPMBY(0.0), fBPMAZ(-957.6e-3), fBPMBZ(-692.0e-3), fBPMARes(0.3e-3), fBPMBRes(0.3e-3), pDrift(NULL), pfGetBPM(NULL)
{
    if (pG2PBPM) {
        Error("G2PBPM()", "Only one instance of G2PBPM allowed.");
        MakeZombie();
        return;
    }
    pG2PBPM = this;

    fPriority = 2;
    Clear();
}

G2PBPM::~G2PBPM()
{
    if (pG2PBPM == this) pG2PBPM = NULL;
}

int G2PBPM::Init()
{
    //static const char* const here = "Init()";

    if (G2PProcBase::Init() != 0) return fStatus;

    pDrift = static_cast<G2PDrift*> (gG2PApps->Find("G2PDrift"));
    if (!pDrift) {
        pDrift = new G2PDrift();
        gG2PApps->Add(pDrift);
    }

    return (fStatus = kOK);
}

int G2PBPM::Begin()
{
    //static const char* const here = "Begin()";

    if (G2PProcBase::Begin() != 0) return fStatus;

    SetBPMPos();

    return (fStatus = kOK);
}

int G2PBPM::Process()
{
    static const char* const here = "Process()";

    if (fDebug > 2) Info(here, " ");

    fV5beam_lab[0] = gG2PVars->FindSuffix("gun.beam.l_x")->GetValue();
    fV5beam_lab[1] = gG2PVars->FindSuffix("gun.beam.l_t")->GetValue();
    fV5beam_lab[2] = gG2PVars->FindSuffix("gun.beam.l_y")->GetValue();
    fV5beam_lab[3] = gG2PVars->FindSuffix("gun.beam.l_p")->GetValue();
    fV5beam_lab[4] = gG2PVars->FindSuffix("gun.beam.l_z")->GetValue();

    double V4[4];
    GetBPM(fV5beam_lab, fV5bpm_bpm, V4);
    fV2bpma_bpm[0] = V4[0];
    fV2bpma_bpm[1] = V4[1];
    fV2bpmb_bpm[0] = V4[2];
    fV2bpmb_bpm[1] = V4[3];

    TransBPM2Lab(fV5bpm_bpm, fV5bpm_lab);

    if (fDebug > 1) Info(here, "bpm_lab   : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5bpm_lab[0], fV5bpm_lab[1], fV5bpm_lab[2], fV5bpm_lab[3], fV5bpm_lab[4]);

    return 0;
}

void G2PBPM::Clear(Option_t* option)
{
    memset(fV5beam_lab, 0, sizeof (fV5beam_lab));
    memset(fV5bpm_bpm, 0, sizeof (fV5bpm_bpm));
    memset(fV5bpm_lab, 0, sizeof (fV5bpm_lab));
    memset(fV2bpma_bpm, 0, sizeof (fV2bpma_bpm));
    memset(fV2bpmb_bpm, 0, sizeof (fV2bpmb_bpm));

    G2PProcBase::Clear(option);
}

void G2PBPM::SetBPMRes(double a, double b)
{
    fBPMARes = a;
    fBPMBRes = b;

    fConfigIsSet.insert((unsigned long) &fBPMARes);
    fConfigIsSet.insert((unsigned long) &fBPMBRes);
}

void G2PBPM::GetBPM(const double* V5beam_lab, double* V5bpm_bpm, double* V4)
{
    static const char* const here = "GetBPM()";

    (this->*pfGetBPM)(V5beam_lab, V5bpm_bpm, V4);

    if (fDebug > 2) Info(here, "%10.3e %10.3e %10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e %10.3e %10.3e", V5beam_lab[0], V5beam_lab[1], V5beam_lab[2], V5beam_lab[3], V5beam_lab[4], V5bpm_bpm[0], V5bpm_bpm[1], V5bpm_bpm[2], V5bpm_bpm[3], V5bpm_bpm[4]);
}

void G2PBPM::TransBPM2Lab(const double* V5_bpm, double* V5_lab)
{
    static const char* const here = "TransBPM2Lab()";

    V5_lab[0] = V5_bpm[0];
    V5_lab[2] = V5_bpm[2];
    V5_lab[4] = V5_bpm[4];

    double p[3] = {tan(V5_bpm[3]), tan(V5_bpm[1]), 1.0};
    double pp = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
    V5_lab[1] = acos(1.0 / pp);
    V5_lab[3] = atan2(p[1], p[0]);

    if (fDebug > 3) Info(here, "%10.3e %10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e %10.3e", V5_bpm[0], V5_bpm[1], V5_bpm[2], V5_bpm[3], V5_lab[0], V5_lab[1], V5_lab[2], V5_lab[3]);
}

void G2PBPM::GetBPM0(const double* V5beam_lab, double* V5bpm_bpm, double* V4)
{
    // it is an estimation when there is target field
    // need to find better method

    double x[3] = {V5beam_lab[0], V5beam_lab[2], V5beam_lab[4]};
    double p[3] = {fBeamEnergy * sin(V5beam_lab[1]) * cos(V5beam_lab[3]), fBeamEnergy * sin(V5beam_lab[1]) * sin(V5beam_lab[3]), fBeamEnergy * cos(V5beam_lab[1])};

    pDrift->Drift(x, p, fBPMBZ, x, p);
    V4[2] = pRand->Gaus(x[0] - fBPMBX, fBPMBRes)*1e3;
    V4[3] = pRand->Gaus(x[1] - fBPMBY, fBPMBRes)*1e3;
    pDrift->Drift(x, p, fBPMAZ, x, p);
    V4[0] = pRand->Gaus(x[0] - fBPMAX, fBPMARes)*1e3;
    V4[1] = pRand->Gaus(x[1] - fBPMAY, fBPMARes)*1e3;

    double res = fBPMBRes / (fBPMBZ - fBPMAZ);
    // here theta phi are special coords defined by Pengjia
    // theta is dy/dz in hall coords
    // phi is dx/dz in hall coords
    double theta = pRand->Gaus(atan2(p[1], p[2]), res);
    double phi = pRand->Gaus(atan2(p[0], p[2]), res);
    p[1] = p[2] * tan(theta);
    p[0] = p[2] * tan(phi);
    double normF = fBeamEnergy / sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
    p[0] *= normF;
    p[1] *= normF;
    p[2] *= normF;
    x[0] = pRand->Gaus(x[0], fBPMARes);
    x[1] = pRand->Gaus(x[1], fBPMARes);
    pDrift->Drift(x, p, 0.0, x, p);
    V5bpm_bpm[0] = x[0];
    V5bpm_bpm[1] = atan2(p[1], p[2]);
    V5bpm_bpm[2] = x[1];
    V5bpm_bpm[3] = atan2(p[0], p[2]);
    V5bpm_bpm[4] = x[2];
}

void G2PBPM::GetBPM1(const double* V5beam_lab, double* V5bpm_bpm, double* V4)
{
    using namespace Orbit1;

    float x[4];
    GetBPMAB(V5beam_lab, x);
    for (int i = 0; i < 4; i++) V4[i] = x[i];

    V5bpm_bpm[0] = target_x(x, 4)*1e-3;
    V5bpm_bpm[1] = target_theta(x, 4);
    V5bpm_bpm[2] = target_y(x, 4)*1e-3;
    V5bpm_bpm[3] = target_phi(x, 4);
    V5bpm_bpm[4] = 0.0;
}

void G2PBPM::GetBPM4(const double* V5beam_lab, double* V5bpm_bpm, double* V4)
{
    using namespace Orbit4;

    float x[4];
    GetBPMAB(V5beam_lab, x);
    for (int i = 0; i < 4; i++) V4[i] = x[i];

    V5bpm_bpm[0] = target_x(x, 4)*1e-3;
    V5bpm_bpm[1] = target_theta(x, 4);
    V5bpm_bpm[2] = target_y(x, 4)*1e-3;
    V5bpm_bpm[3] = target_phi(x, 4);
    V5bpm_bpm[4] = 0.0;
}

void G2PBPM::GetBPM5(const double* V5beam_lab, double* V5bpm_bpm, double* V4)
{
    using namespace Orbit5;

    float x[4];
    GetBPMAB(V5beam_lab, x);
    for (int i = 0; i < 4; i++) V4[i] = x[i];

    V5bpm_bpm[0] = target_x(x, 4)*1e-3;
    V5bpm_bpm[1] = target_theta(x, 4);
    V5bpm_bpm[2] = target_y(x, 4)*1e-3;
    V5bpm_bpm[3] = target_phi(x, 4);
    V5bpm_bpm[4] = 0.0;
}

void G2PBPM::GetBPM7(const double* V5beam_lab, double* V5bpm_bpm, double* V4)
{
    using namespace Orbit7;

    float x[4];
    GetBPMAB(V5beam_lab, x);
    for (int i = 0; i < 4; i++) V4[i] = x[i];

    V5bpm_bpm[0] = target_x(x, 4)*1e-3;
    V5bpm_bpm[1] = target_theta(x, 4);
    V5bpm_bpm[2] = target_y(x, 4)*1e-3;
    V5bpm_bpm[3] = target_phi(x, 4);
    V5bpm_bpm[4] = 0.0;
}

void G2PBPM::GetBPM9(const double* V5beam_lab, double* V5bpm_bpm, double* V4)
{
    using namespace Orbit9;

    float x[4];
    GetBPMAB(V5beam_lab, x);
    for (int i = 0; i < 4; i++) V4[i] = x[i];

    V5bpm_bpm[0] = target_x(x, 4)*1e-3;
    V5bpm_bpm[1] = target_theta(x, 4);
    V5bpm_bpm[2] = target_y(x, 4)*1e-3;
    V5bpm_bpm[3] = target_phi(x, 4);
    V5bpm_bpm[4] = 0.0;
}

void G2PBPM::GetBPMO(const double* V5beam_lab, double* V5bpm_bpm, double* V4)
{
    // optics (no field) situation

    float x[4];
    GetBPMAB(V5beam_lab, x);
    for (int i = 0; i < 4; i++) V4[i] = x[i];

    double xbpm[4];
    xbpm[0] = x[0]*1e-3;
    xbpm[1] = x[1]*1e-3;
    xbpm[2] = x[2]*1e-3;
    xbpm[3] = x[3]*1e-3;

#ifdef USE_SINGLE_BPM
    V5bpm_bpm[0] = (xbpm[0] + xbpm[2]) / 2.0;
    V5bpm_bpm[1] = 0.0;
    V5bpm_bpm[2] = (xbpm[1] + xbpm[3]) / 2.0;
    V5bpm_bpm[3] = 0.0;
    V5bpm_bpm[4] = 0.0;
#else
    // here theta phi are special coords defined by Pengjia
    // theta is dy/dz in hall coords
    // phi is dx/dz in hall coords
    double theta = atan2(xbpm[3] - xbpm[1], fBPMBZ - fBPMAZ);
    double phi = atan2(xbpm[2] - xbpm[0], fBPMBZ - fBPMAZ);
    V5bpm_bpm[0] = xbpm[0] + (-fBPMAZ) * tan(phi);
    V5bpm_bpm[1] = theta;
    V5bpm_bpm[2] = xbpm[1] + (-fBPMAZ) * tan(theta);
    V5bpm_bpm[3] = phi;
    V5bpm_bpm[4] = 0.0;
#endif
}

void G2PBPM::GetBPMAB(const double* V5beam_lab, float* xout)
{
    static const char* const here = "GetBPMAB()";

    double x[3] = {V5beam_lab[0], V5beam_lab[2], V5beam_lab[4]};
    double p[3] = {fBeamEnergy * sin(V5beam_lab[1]) * cos(V5beam_lab[3]), fBeamEnergy * sin(V5beam_lab[1]) * sin(V5beam_lab[3]), fBeamEnergy * cos(V5beam_lab[1])};

    int save = pDrift->GetDebugLevel();
    if (fDebug <= 3) pDrift->SetDebugLevel(0);
    pDrift->Drift(x, p, fBPMBZ, x, p);
    xout[2] = pRand->Gaus(x[0] - fBPMBX, fBPMBRes)*1e3;
    xout[3] = pRand->Gaus(x[1] - fBPMBY, fBPMBRes)*1e3;
    pDrift->Drift(x, p, fBPMAZ, x, p);
    xout[0] = pRand->Gaus(x[0] - fBPMAX, fBPMARes)*1e3;
    xout[1] = pRand->Gaus(x[1] - fBPMAY, fBPMARes)*1e3;
    pDrift->SetDebugLevel(save);

    if (fDebug > 2) {
        Info(here, "%10.3e %10.3e", xout[0], xout[1]);
        Info(here, "%10.3e %10.3e", xout[2], xout[3]);
    }
}

void G2PBPM::SetBPMPos()
{
    static const char* const here = "SetBPMPos()";

    int orbit;
    if (fabs(fFieldRatio - 0.5) < 1e-8) {
        if (fabs(fBeamEnergy - 2.254) < 0.2) {
            fBPMAX = 0.743197468425e-3;
            fBPMAY = -97.9439291505e-3;
            fBPMAZ = -940.46966234e-3;
            fBPMBX = 1.04489437383e-3;
            fBPMBY = -69.5218412795e-3;
            fBPMBZ = -676.469580767e-3;
            pfGetBPM = &G2PBPM::GetBPM5;
            orbit = 5;
        } else if (fabs(fBeamEnergy - 1.706) < 0.2) {
            fBPMAX = 0.334417934854e-3;
            fBPMAY = -130.742151697e-3;
            fBPMAZ = -943.969646821e-3;
            fBPMBX = 0.436195607141e-3;
            fBPMBY = -93.1201715402e-3;
            fBPMBZ = -681.069555637e-3;
            pfGetBPM = &G2PBPM::GetBPM4;
            orbit = 4;
        } else if (fabs(fBeamEnergy - 1.159) < 0.2) {
            fBPMAX = 1.04109622687e-3;
            fBPMAY = -196.201718795e-3;
            fBPMAZ = -948.669889999e-3;
            fBPMBX = 1.24268316845e-3;
            fBPMBY = -143.579627734e-3;
            fBPMBZ = -688.370604604e-3;
            pfGetBPM = &G2PBPM::GetBPM1;
            orbit = 1;
        } else {
            pfGetBPM = &G2PBPM::GetBPM0;
            orbit = 0;
        }
    } else if (fabs(fFieldRatio - 1.0) < 1e-8) {
        if (fabs(fBeamEnergy - 2.254) < 0.2) {
            fBPMAX = 0.138592400799e-3;
            fBPMAY = -79.962134626e-3;
            fBPMAZ = -939.76959472e-3;
            fBPMBX = 0.442605162088e-3;
            fBPMBY = -55.6342835205e-3;
            fBPMBZ = -675.269596904e-3;
            pfGetBPM = &G2PBPM::GetBPM7;
            orbit = 7;
        } else if (fabs(fBeamEnergy - 3.355) < 0.2) {
            fBPMAX = 0.0609719900963e-3;
            fBPMAY = -53.5499078644e-3;
            fBPMAZ = -939.169701771e-3;
            fBPMBX = 0.662453421257e-3;
            fBPMBY = -37.0276045709e-3;
            fBPMBZ = -674.16965108e-3;
            pfGetBPM = &G2PBPM::GetBPM9;
            orbit = 9;
        } else {
            pfGetBPM = &G2PBPM::GetBPM0;
            orbit = 0;
        }
    } else if (fabs(fFieldRatio) < 1e-8) {
        pfGetBPM = &G2PBPM::GetBPMO;
        orbit = -1;
    } else {
        pfGetBPM = &G2PBPM::GetBPM0;
        orbit = 0;
    }

    if (fDebug > 0) {
        if (orbit >= 0) Info(here, "Using orbit %d.", orbit);
        else if (orbit == -1) Info(here, "Using optics orbit.");
    }
}

int G2PBPM::Configure(EMode mode)
{
    if (mode == kREAD || mode == kTWOWAY) {
        if (fIsInit) return 0;
        else fIsInit = true;
    }

    ConfDef confs[] = {
        {"run.e0", "Beam Energy", kDOUBLE, &fBeamEnergy},
        {"field.ratio", "Field Ratio", kDOUBLE, &fFieldRatio},
        {"a.res", "BPM A Resolution", kDOUBLE, &fBPMARes},
        {"b.res", "BPM B Resolution", kDOUBLE, &fBPMBRes},
        {0}
    };

    return ConfigureFromList(confs, mode);
}

int G2PBPM::DefineVariables(EMode mode)
{
    if (mode == kDEFINE && fIsSetup) return 0;
    fIsSetup = (mode == kDEFINE);

    VarDef vars[] = {
        {"b_t", "BPM T", kDOUBLE, &fV5bpm_bpm[1]},
        {"b_p", "BPM P", kDOUBLE, &fV5bpm_bpm[3]},
        {"l_x", "BPM X (lab)", kDOUBLE, &fV5bpm_lab[0]},
        {"l_t", "BPM T (lab)", kDOUBLE, &fV5bpm_lab[1]},
        {"l_y", "BPM Y (lab)", kDOUBLE, &fV5bpm_lab[2]},
        {"l_p", "BPM P (lab)", kDOUBLE, &fV5bpm_lab[3]},
        {"l_z", "BPM Z (lab)", kDOUBLE, &fV5bpm_lab[4]},
        {"a.x", "BPMA X", kDOUBLE, &fV2bpma_bpm[0]},
        {"a.y", "BPMA X", kDOUBLE, &fV2bpma_bpm[1]},
        {"b.x", "BPMB X", kDOUBLE, &fV2bpmb_bpm[0]},
        {"b.y", "BPMB X", kDOUBLE, &fV2bpmb_bpm[1]},
        {0}
    };

    return DefineVarsFromList(vars, mode);
}

void G2PBPM::MakePrefix()
{
    const char* base = "bpm";

    G2PAppBase::MakePrefix(base);
}

ClassImp(G2PBPM)
