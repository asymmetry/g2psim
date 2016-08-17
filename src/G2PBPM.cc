// -*- C++ -*-

/* class G2PBPM
 * Calculate beam readout at BPM and target using kinematics from event generator.
 * Transport functions defined in G2PBPMTrans is used in this class.
 * Orbits are defined in G2PBPMTrans.
 *
 * Variables ending with "_bpm" are defined in a special coordinates.
 * In output, these variables are labeled as "b_".
 */

// History:
//   Mar 2013, C. Gu, First public version.
//   Apr 2013, C. Gu, Add Pengjia's fitting result.
//   Jul 2013, C. Gu, Treat optics (no field) case specially.
//   Sep 2013, C. Gu, Rewrite it as a G2PProcBase class.
//   Jan 2015, C. Gu, Use new Drift() function in G2PProcBase.
//   May 2016, C. Gu, Only take an overall resolution.
//

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"

#include "G2PAppBase.hh"
#include "G2PAppList.hh"
#include "G2PGlobals.hh"
#include "G2PProcBase.hh"
#include "G2PRand.hh"
#include "G2PVar.hh"
#include "G2PVarDef.hh"
#include "G2PVarList.hh"

#include "G2PBPM.hh"
#include "G2PBPMTrans.hh"

G2PBPM *G2PBPM::pG2PBPM = NULL;

G2PBPM::G2PBPM() : fE0(0.0), fFieldType(0), fFieldRatio(0.0), fBPMResPos(0.0), fBPMResAngle(0.0), pfGetBPM(NULL)
{
    if (pG2PBPM) {
        Error("G2PBPM()", "Only one instance of G2PBPM allowed.");
        MakeZombie();
        return;
    }

    pG2PBPM = this;

    memset(fBPMA, 0, sizeof(fBPMA));
    memset(fBPMB, 0, sizeof(fBPMB));

    fPriority = 2;

    Clear();
}

G2PBPM::~G2PBPM()
{
    if (pG2PBPM == this)
        pG2PBPM = NULL;
}

int G2PBPM::Begin()
{
    //static const char* const here = "Begin()";

    if (G2PProcBase::Begin() != 0)
        return (fStatus = kBEGINERROR);

    SetBPMPos();
    Configure(kWRITE); // Write BPM positions to G2PRun

    return (fStatus = kOK);
}

int G2PBPM::Process()
{
    static const char *const here = "Process()";

    if (fDebug > 2)
        Info(here, " ");

    if (gG2PVars->FindSuffix("beam.l_x")) {
        fV5beam_lab[0] = gG2PVars->FindSuffix("beam.l_x")->GetValue();
        fV5beam_lab[1] = gG2PVars->FindSuffix("beam.l_t")->GetValue();
        fV5beam_lab[2] = gG2PVars->FindSuffix("beam.l_y")->GetValue();
        fV5beam_lab[3] = gG2PVars->FindSuffix("beam.l_p")->GetValue();
        fV5beam_lab[4] = gG2PVars->FindSuffix("beam.l_z")->GetValue();
    } else
        return -1;

    GetBPM(fV5beam_lab, fV5bpm_bpm, fV4bpmab_bpm);
    fV5bpm_bpm[0] = pRand->Gaus(fV5bpm_bpm[0], fBPMResPos);
    fV5bpm_bpm[1] = pRand->Gaus(fV5bpm_bpm[1], fBPMResAngle);
    fV5bpm_bpm[2] = pRand->Gaus(fV5bpm_bpm[2], fBPMResPos);
    fV5bpm_bpm[3] = pRand->Gaus(fV5bpm_bpm[3], fBPMResAngle);

    BPM2HCS(fV5bpm_bpm, fV5bpm_lab);
    HCS2TCS(fV5bpm_lab, fV5bpm_tr, fbpmz_tr);

    if (fDebug > 1)
        Info(here, "bpm_lab   : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5bpm_lab[0], fV5bpm_lab[1], fV5bpm_lab[2], fV5bpm_lab[3], fV5bpm_lab[4]);

    return 0;
}

void G2PBPM::Clear(Option_t *opt)
{
    memset(fV5beam_lab, 0, sizeof(fV5beam_lab));
    memset(fV5bpm_bpm, 0, sizeof(fV5bpm_bpm));
    memset(fV5bpm_lab, 0, sizeof(fV5bpm_lab));
    memset(fV5bpm_tr, 0, sizeof(fV5bpm_tr));
    memset(fV4bpmab_bpm, 0, sizeof(fV4bpmab_bpm));

    G2PProcBase::Clear(opt);
}

void G2PBPM::SetBPMRes(double pos, double angle)
{
    fBPMResPos = pos;
    fBPMResAngle = angle;

    fConfigIsSet.insert((unsigned long) &fBPMResPos);
    fConfigIsSet.insert((unsigned long) &fBPMResAngle);
}

void G2PBPM::GetBPM(const double *V5beam_lab, double *V5bpm_bpm, double *V4bpmab_bpm)
{
    static const char *const here = "GetBPM()";

    (this->*pfGetBPM)(V5beam_lab, V5bpm_bpm, V4bpmab_bpm);

    if (fDebug > 2)
        Info(here, "%10.3e %10.3e %10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e %10.3e %10.3e", V5beam_lab[0], V5beam_lab[1], V5beam_lab[2], V5beam_lab[3], V5beam_lab[4], V5bpm_bpm[0], V5bpm_bpm[1], V5bpm_bpm[2], V5bpm_bpm[3], V5bpm_bpm[4]);
}

void G2PBPM::SetBPMPos()
{
    static const char *const here = "SetBPMPos()";

    static const double bpmpos[10][6] = {
        {0, 0, -957.6e-3, 0, 0, -692.0e-3}, // orbit 0
        {1.0e-3, -196.1e-3, -964.1e-3, 1.2e-3, -143.4e-3, -703.8e-3}, // orbit 1
        { -0.4e-3, -33.6e-3, -955.5e-3, -0.6e-3, -17.4e-3, -690.5e-3}, // orbit 2
        {0.2e-3, -27.1e-3, -955.2e-3, 0.5e-3, -19.8e-3, -689.7e-3}, // orbit 3
        {0.3e-3, -130.8e-3, -959.4e-3, 0.4e-3, -93.1e-3, -696.5e-3}, // orbit 4
        {0.7e-3, -98.0e-3, -955.9e-3, 1.0e-3, -69.5e-3, -691.9e-3}, // orbit 5
        {0, 0, 0, 0, 0, 0},
        {0.1e-3, -80.0e-3, -955.2e-3, 0.4e-3, -55.6e-3, -690.7e-3}, // orbit 7
        {0.4e-3, -24.2e-3, -955.1e-3, 0.6e-3, -17.1e-3, -689.6e-3}, // orbit 8
        {0, -53.6e-3, -954.6e-3, 0.6e-3, -37.0e-3, -689.6e-3} // orbit 9
    };

    int orbit = 0;

    if (fFieldType == 10 || fFieldType == 11) {
        if (fabs(fFieldRatio - 0.5) < 1e-4) {
            if (fabs(fE0 - 2.254) < 0.2) {
                orbit = 5;
                pfGetBPM = &G2PBPM::GetBPM5;
            } else if (fabs(fE0 - 1.706) < 0.2) {
                orbit = 4;
                pfGetBPM = &G2PBPM::GetBPM4;
            } else if (fabs(fE0 - 1.159) < 0.2) {
                orbit = 1;
                pfGetBPM = &G2PBPM::GetBPM1;
            }
        } else if (fabs(fFieldRatio - 1.0) < 1e-4) {
            if (fabs(fE0 - 2.254) < 0.2) {
                orbit = 7;
                pfGetBPM = &G2PBPM::GetBPM7;
            } else if (fabs(fE0 - 3.355) < 0.2) {
                orbit = 9;
                pfGetBPM = &G2PBPM::GetBPM9;
            }
        }
    } else if (fFieldType == 30) {
        if (fabs(fE0 - 2.254) < 0.2) {
            orbit = 8;
            pfGetBPM = &G2PBPM::GetBPM8;
        } else if (fabs(fE0 - 1.706) < 0.2) {
            orbit = 3;
            pfGetBPM = &G2PBPM::GetBPM3;
        } else if (fabs(fE0 - 1.159) < 0.2) {
            orbit = 2;
            pfGetBPM = &G2PBPM::GetBPM2;
        }
    }

    if (orbit == 0)
        pfGetBPM = &G2PBPM::GetBPM0;

    for (int i = 0; i < 3; i++) {
        fBPMA[i] = bpmpos[orbit][i];
        fBPMB[i] = bpmpos[orbit][i + 3];
        fConfigIsSet.insert((unsigned long) &fBPMA[i]);
        fConfigIsSet.insert((unsigned long) &fBPMB[i]);
    }

    if (fDebug > 1)
        Info(here, "Using orbit %d.", orbit);
}

void G2PBPM::GetBPM0(const double *V5beam_lab, double *V5bpm_bpm, double *V4bpmab_bpm)
{
    // it is an estimation when there is target field
    // need to find better method

    GetBPMAB(V5beam_lab, V4bpmab_bpm);

    for (int i = 0; i < 4; i++)
        V4bpmab_bpm[i] *= 1e-3;

    double x[3] = {(V4bpmab_bpm[2] + V4bpmab_bpm[0]) / 2, (V4bpmab_bpm[3] + V4bpmab_bpm[1]) / 2, (fBPMB[2] + fBPMA[2]) / 2};
    double p[3] = {(V4bpmab_bpm[2] - V4bpmab_bpm[0]) / (fBPMB[2] - fBPMA[2]), (V4bpmab_bpm[3] - V4bpmab_bpm[1]) / (fBPMB[2] - fBPMA[2]), 1.0};
    double normE = fE0 / sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);

    for (int i = 0; i < 3; i++)
        p[i] *= normE;

    Drift("forward", x, p, 0.0, x, p);

    V5bpm_bpm[0] = x[0];
    V5bpm_bpm[1] = atan2(p[1], p[2]);
    V5bpm_bpm[2] = x[1];
    V5bpm_bpm[3] = atan2(p[0], p[2]);
    V5bpm_bpm[4] = x[2];
}

#define GETBPM(ORBITN) \
void G2PBPM::GetBPM##ORBITN(const double *V5beam_lab, double *V5bpm_bpm, double *V4bpmab_bpm) \
{ \
    using namespace Orbit##ORBITN; \
\
    GetBPMAB(V5beam_lab, V4bpmab_bpm); \
\
    V5bpm_bpm[0] = target_x(V4bpmab_bpm, 4) * 1e-3; \
    V5bpm_bpm[1] = target_theta(V4bpmab_bpm, 4); \
    V5bpm_bpm[2] = target_y(V4bpmab_bpm, 4) * 1e-3; \
    V5bpm_bpm[3] = target_phi(V4bpmab_bpm, 4); \
    V5bpm_bpm[4] = 0.0; \
}

GETBPM(1)
GETBPM(2)
GETBPM(3)
GETBPM(4)
GETBPM(5)
GETBPM(7)
GETBPM(8)
GETBPM(9)

#undef GETBPM

void G2PBPM::GetBPMAB(const double *V5beam_lab, double *V4bpmab_bpm)
{
    static const char *const here = "GetBPMAB()";

    double x[3] = {V5beam_lab[0], V5beam_lab[2], V5beam_lab[4]};
    double p[3] = {fE0 * sin(V5beam_lab[1]) *cos(V5beam_lab[3]), fE0 * sin(V5beam_lab[1]) *sin(V5beam_lab[3]), fE0 * cos(V5beam_lab[1])};

    int save = fDebug;
    fDebug = 0;
    Drift("backward", x, p, fBPMB[2], x, p);
    V4bpmab_bpm[2] = (x[0] - fBPMB[0]) * 1e3;
    V4bpmab_bpm[3] = (x[1] - fBPMB[1]) * 1e3;
    Drift("backward", x, p, fBPMA[2], x, p);
    V4bpmab_bpm[0] = (x[0] - fBPMA[0]) * 1e3;
    V4bpmab_bpm[1] = (x[1] - fBPMA[1]) * 1e3;
    fDebug = save;

    if (fDebug > 2) {
        Info(here, "%10.3e %10.3e", V4bpmab_bpm[0], V4bpmab_bpm[1]);
        Info(here, "%10.3e %10.3e", V4bpmab_bpm[2], V4bpmab_bpm[3]);
    }
}

int G2PBPM::Configure(EMode mode)
{
    if ((mode == kREAD || mode == kTWOWAY) && fConfigured)
        return 0;

    if (G2PProcBase::Configure(mode) != 0)
        return -1;

    ConfDef confs[] = {
        {"run.e0", "Beam Energy", kDOUBLE, &fE0},
        {"run.field", "Field Type", kINT, &fFieldType},
        {"field.ratio", "Field Ratio", kDOUBLE, &fFieldRatio},
        {"res.pos", "BPM Space Resolution", kDOUBLE, &fBPMResPos},
        {"res.angle", "BPM Angle Resolution", kDOUBLE, &fBPMResPos},
        {"a.x", "BPM A X Position", kDOUBLE, &fBPMA[0]},
        {"a.y", "BPM A Y Position", kDOUBLE, &fBPMA[1]},
        {"a.z", "BPM A Z Position", kDOUBLE, &fBPMA[2]},
        {"b.x", "BPM B X Position", kDOUBLE, &fBPMB[0]},
        {"b.y", "BPM B Y Position", kDOUBLE, &fBPMB[1]},
        {"b.z", "BPM B Z Position", kDOUBLE, &fBPMB[2]},
        {0}
    };

    return ConfigureFromList(confs, mode);
}

int G2PBPM::DefineVariables(EMode mode)
{
    if (mode == kDEFINE && fDefined)
        return 0;

    if (G2PProcBase::DefineVariables(mode) != 0)
        return -1;

    VarDef vars[] = {
        {"b_x", "BPM X (bpm)", kDOUBLE, &fV5bpm_bpm[0]},
        {"b_t", "BPM T (bpm)", kDOUBLE, &fV5bpm_bpm[1]},
        {"b_y", "BPM Y (bpm)", kDOUBLE, &fV5bpm_bpm[2]},
        {"b_p", "BPM P (bpm)", kDOUBLE, &fV5bpm_bpm[3]},
        {"b_z", "BPM Z (bpm)", kDOUBLE, &fV5bpm_bpm[4]},
        {"l_x", "BPM X (lab)", kDOUBLE, &fV5bpm_lab[0]},
        {"l_t", "BPM T (lab)", kDOUBLE, &fV5bpm_lab[1]},
        {"l_y", "BPM Y (lab)", kDOUBLE, &fV5bpm_lab[2]},
        {"l_p", "BPM P (lab)", kDOUBLE, &fV5bpm_lab[3]},
        {"l_z", "BPM Z (lab)", kDOUBLE, &fV5bpm_lab[4]},
        {"x", "BPM X", kDOUBLE, &fV5bpm_tr[0]},
        {"t", "BPM T", kDOUBLE, &fV5bpm_tr[1]},
        {"y", "BPM Y", kDOUBLE, &fV5bpm_tr[2]},
        {"p", "BPM P", kDOUBLE, &fV5bpm_tr[3]},
        {"z", "BPM Z", kDOUBLE, &fbpmz_tr},
        {"a.x", "BPMA X", kDOUBLE, &fV4bpmab_bpm[0]},
        {"a.y", "BPMA X", kDOUBLE, &fV4bpmab_bpm[1]},
        {"b.x", "BPMB X", kDOUBLE, &fV4bpmab_bpm[2]},
        {"b.y", "BPMB X", kDOUBLE, &fV4bpmab_bpm[3]},
        {0}
    };

    return DefineVarsFromList(vars, mode);
}

void G2PBPM::MakePrefix()
{
    const char *base = "bpm";

    G2PAppBase::MakePrefix(base);
}

ClassImp(G2PBPM)
