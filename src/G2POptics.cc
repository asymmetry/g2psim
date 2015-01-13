// -*- C++ -*-

/* class G2POptics
 * Special class to treat optics data.
 */

// History:
//   Jan 2014, C. Gu, First public version.
//

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <queue>
#include <vector>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"

#include "G2PAppBase.hh"
#include "G2PAppList.hh"
#include "G2PGlobals.hh"
#include "G2PProcBase.hh"
#include "G2PSieve.hh"
#include "G2PVar.hh"
#include "G2PVarDef.hh"
#include "G2PVarList.hh"

#include "G2POptics.hh"

using namespace std;

static double kPHI = (sqrt(5.0) - 1.0) / 2.0;

G2POptics *G2POptics::pG2POptics = NULL;

G2POptics::G2POptics()
{
    // Only for ROOT I/O
}

G2POptics::G2POptics(const char *filename) : fDataFile(filename), fTargetMass(0.0), fE0(0.0), fTiltAngle(0.0), fELoss(0.0), fNFoil(1), fHoleID(-1), pSieve(NULL)
{
    if (pG2POptics) {
        Error("G2POptics()", "Only one instance of G2POptics allowed.");
        MakeZombie();
        return;
    }

    pG2POptics = this;

    fPriority = 1;

    fHRSMomentumV.clear();
    fHRSMomentumV.push_back(0.0);
    fFoilZV.clear();
    fFoilZV.push_back(0.0);
    fTiltAngleV.clear();
    fTiltAngleV.push_back(0.0);
    fELossV.clear();
    fELossV.push_back(0.0);

    Clear();
}

G2POptics::~G2POptics()
{
    if (pG2POptics == this)
        pG2POptics = NULL;
}

int G2POptics::Begin()
{
    static const char *const here = "Begin()";

    if (G2PProcBase::Begin() != 0)
        return (fStatus = kBEGINERROR);

    pSieve = static_cast<G2PSieve *>(gG2PApps->Find("G2PSieve"));

    if (LoadData() != 0) {
        Error(here, "Cannot read data.");
        return (fStatus = kBEGINERROR);
    }

    return (fStatus = kOK);
}

int G2POptics::Process()
{
    static const char *const here = "Process()";

    sData tempdata;

    if (fData.empty())
        return -1;

    tempdata = fData.front();

    if (fDebug > 2)
        Info(here, "%04d %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e", tempdata.ind, tempdata.xf, tempdata.tf, tempdata.yf, tempdata.pf, tempdata.eb, tempdata.xb, tempdata.yb);

    fHoleID = tempdata.ind;
    int kineID = fHoleID / (pSieve->GetNRow() * pSieve->GetNCol() * fNFoil);
    int res = fHoleID % (pSieve->GetNRow() * pSieve->GetNCol() * fNFoil);
    int foilID = res / (pSieve->GetNRow() * pSieve->GetNCol());
    res = res % (pSieve->GetNRow() * pSieve->GetNCol());

    fE0 = tempdata.eb;

    fV3bpm_lab[0] = tempdata.xb;
    fV3bpm_lab[1] = tempdata.yb;

    fV5fp_det[0] = tempdata.xf;
    fV5fp_det[1] = tempdata.tf;
    fV5fp_det[2] = tempdata.yf;
    fV5fp_det[3] = tempdata.pf;

    fHRSMomentum = fHRSMomentumV[kineID];
    fTiltAngle = fTiltAngleV[kineID * fNFoil + foilID];
    fELoss = fELossV[foilID];
    fV3bpm_lab[2] = fFoilZV[foilID];

    if (fDebug > 1)
        Info(here, "bpm_lab   : %10.3e %10.3e %10.3e", fV3bpm_lab[0], fV3bpm_lab[1], fV3bpm_lab[2]);

    HCS2TCS(fV3bpm_lab[0], fV3bpm_lab[1], fV3bpm_lab[2], fV3bpm_tr[0], fV3bpm_tr[1], fV3bpm_tr[2]);

    double V3siv_tr[3];
    pSieve->GetPos(res, V3siv_tr);

    if (fDebug > 1)
        Info(here, "sieve_real: %10.3e %10.3e %10.3e", V3siv_tr[0], V3siv_tr[1], V3siv_tr[2]);

    double V3pd_tr[3];
    V3pd_tr[0] = V3siv_tr[0] - fV3bpm_tr[0];
    V3pd_tr[1] = V3siv_tr[1] - fV3bpm_tr[1];
    V3pd_tr[2] = V3siv_tr[2] - fV3bpm_tr[2];

    double ang[2];
    ang[0] = atan(V3pd_tr[0] / V3pd_tr[2]);
    ang[1] = atan(V3pd_tr[1] / V3pd_tr[2]);

    double pos[2];
    double dlast;
    CalPos(ang, pos);
    dlast = Distance(pos, V3siv_tr);

    while (dlast > 1e-5) {
        double step = dlast / pSieve->GetZ();

        double left[2], right[2];
        double testl[2], testr[2];

        left[0] = ang[0] - 2 * step;
        left[1] = ang[1];
        right[0] = ang[0] + 2 * step;
        right[1] = ang[1];
        testl[0] = right[0] - kPHI * (right[0] - left[0]);
        testl[1] = ang[1];
        testr[0] = left[0] + kPHI * (right[0] - left[0]);
        testr[1] = ang[1];

        double dl, dr;
        CalPos(testl, pos);
        dl = Distance(pos, V3siv_tr);
        CalPos(testr, pos);
        dr = Distance(pos, V3siv_tr);

        while (fabs(right[0] - left[0]) > 0.7e-5) {
            if (dl <= dr) {
                right[0] = testr[0];
                testr[0] = testl[0];
                dr = dl;
                testl[0] = right[0] - kPHI * (right[0] - left[0]);
                CalPos(testl, pos);
                dl = Distance(pos, V3siv_tr);
            } else {
                left[0] = testl[0];
                testl[0] = testr[0];
                dl = dr;
                testr[0] = left[0] + kPHI * (right[0] - left[0]);
                CalPos(testr, pos);
                dr = Distance(pos, V3siv_tr);
            }
        }

        ang[0] = (right[0] + left[0]) / 2;

        step = (dl + dr) / 2.0 / pSieve->GetZ();

        left[0] = ang[0];
        left[1] = ang[1] - 2 * step;
        right[0] = ang[0];
        right[1] = ang[1] + 2 * step;
        testl[0] = ang[0];
        testl[1] = right[1] - kPHI * (right[1] - left[1]);
        testr[0] = ang[0];
        testr[1] = left[1] + kPHI * (right[1] - left[1]);

        CalPos(testl, pos);
        dl = Distance(pos, V3siv_tr);
        CalPos(testr, pos);
        dr = Distance(pos, V3siv_tr);

        while (fabs(right[1] - left[1]) > 0.7e-5) {
            if (dl <= dr) {
                right[1] = testr[1];
                testr[1] = testl[1];
                dr = dl;
                testl[1] = right[1] - kPHI * (right[1] - left[1]);
                CalPos(testl, pos);
                dl = Distance(pos, V3siv_tr);
            } else {
                left[1] = testl[1];
                testl[1] = testr[1];
                dl = dr;
                testr[1] = left[1] + kPHI * (right[1] - left[1]);
                CalPos(testr, pos);
                dr = Distance(pos, V3siv_tr);
            }
        }

        ang[1] = (right[1] + left[1]) / 2;

        CalPos(ang, pos);
        dlast = Distance(pos, V3siv_tr);
    }

    if (fDebug > 1)
        Info(here, "sieve_tr  : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5sieve_tr[0], fV5sieve_tr[1], fV5sieve_tr[2], fV5sieve_tr[3], fV5sieve_tr[4]);

    Project(fV5sieve_tr, pSieve->GetZ(), 0.0, fV5tpproj_tr);

    if (fDebug > 1)
        Info(here, "tpproj_tr : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5tpproj_tr[0], fV5tpproj_tr[1], fV5tpproj_tr[2], fV5tpproj_tr[3], fV5tpproj_tr[4]);

    fData.pop();

    return 0;
}

void G2POptics::Clear(Option_t *option)
{
    fHoleID = -1;

    memset(fV3bpm_lab, 0, sizeof(fV3bpm_lab));
    memset(fV3bpm_tr, 0, sizeof(fV3bpm_tr));
    memset(fV5react_tr, 0, sizeof(fV5react_tr));
    memset(fV5sieve_tr, 0, sizeof(fV5sieve_tr));
    memset(fV5tpproj_tr, 0, sizeof(fV5tpproj_tr));
    memset(fV5fp_det, 0, sizeof(fV5fp_det));

    G2PProcBase::Clear(option);
}

void G2POptics::SetHRSMomentum(int n, double *value)
{
    fHRSMomentumV.clear();

    for (int i = 0; i < n; i++)
        fHRSMomentumV.push_back(value[i]);
}

void G2POptics::SetTiltAngle(int n, double *value)
{
    fTiltAngleV.clear();

    for (int i = 0; i < n; i++)
        fTiltAngleV.push_back(value[i]);
}

void G2POptics::SetFoilZ(int n, double *value)
{
    fNFoil = n;

    fFoilZV.clear();

    for (int i = 0; i < n; i++)
        fFoilZV.push_back(value[i]);
}

void G2POptics::SetEnergyLoss(int n, double *value)
{
    fNFoil = n;

    fELossV.clear();

    for (int i = 0; i < n; i++)
        fELossV.push_back(value[i]);
}

void G2POptics::CalPos(double *ang, double *pos)
{
    double V3pd_tr[3] = {tan(ang[0]), tan(ang[1]), 1.0};

    // Calculate delta based on angle
    double V3pd_lab[3];
    TCS2HCS(V3pd_tr[0], V3pd_tr[1], V3pd_tr[2], V3pd_lab[0], V3pd_lab[1], V3pd_lab[2]);

    double V3ed_lab[3] = {0, tan(fTiltAngle), 1.0};

    double cosscatangle = (V3ed_lab[0] * V3pd_lab[0] + V3ed_lab[1] * V3pd_lab[1] + V3ed_lab[2] * V3pd_lab[2]) / (sqrt(V3ed_lab[0] * V3ed_lab[0] + V3ed_lab[1] * V3ed_lab[1] + V3ed_lab[2] * V3ed_lab[2]) * sqrt(V3pd_lab[0] * V3pd_lab[0] + V3pd_lab[1] * V3pd_lab[1] + V3pd_lab[2] * V3pd_lab[2]));
    double scatmom = (fTargetMass * fE0) / (fTargetMass + fE0 - fE0 * cosscatangle);
    double delta = (scatmom - fHRSMomentum - fELoss) / fHRSMomentum;

    fV5react_tr[0] = fV3bpm_tr[0];
    fV5react_tr[1] = ang[0];
    fV5react_tr[2] = fV3bpm_tr[1];
    fV5react_tr[3] = ang[1];
    fV5react_tr[4] = delta;

    Drift("forward", fV5react_tr, fV3bpm_tr[2], pSieve->GetZ(), fV5sieve_tr);

    pos[0] = fV5sieve_tr[0];
    pos[1] = fV5sieve_tr[2];
}

double G2POptics::Distance(double *V2a, double *V2b)
{
    double V2diff[2] = {V2a[0] - V2b[0], V2a[1] - V2b[1]};
    return sqrt(V2diff[0] * V2diff[0] + V2diff[1] * V2diff[1]);
}

int G2POptics::LoadData()
{
    FILE *fp;

    if ((fp = fopen(fDataFile, "r")) == NULL)
        return -1;

    sData temp;
    fscanf(fp, "%d%lf%lf%lf%lf%lf%lf%lf", &temp.ind, &temp.xf, &temp.tf, &temp.yf, &temp.pf, &temp.eb, &temp.xb, &temp.yb);

    while (!feof(fp)) {
        fData.push(temp);
        fscanf(fp, "%d%lf%lf%lf%lf%lf%lf%lf", &temp.ind, &temp.xf, &temp.tf, &temp.yf, &temp.pf, &temp.eb, &temp.xb, &temp.yb);
    }

    fclose(fp);

    if (!fData.empty())
        return 0;
    else
        return -1;
}

int G2POptics::Configure(EMode mode)
{
    if ((mode == kREAD || mode == kTWOWAY) && fConfigured)
        return 0;

    if (G2PProcBase::Configure(mode) != 0)
        return -1;

    ConfDef confs[] = {
        {"run.target.mass", "Target Mass", kDOUBLE, &fTargetMass},
        {0}
    };

    return ConfigureFromList(confs, mode);
}

int G2POptics::DefineVariables(EMode mode)
{
    if (mode == kDEFINE && fDefined)
        return 0;

    if (G2PProcBase::DefineVariables(mode) != 0)
        return -1;

    VarDef vars[] = {
        {"id", "Hole ID", kINT, &fHoleID},
        {"energy", "Beam Energy", kDOUBLE, &fE0},
        {"bpm.l_x", "BPM X (lab)", kDOUBLE, &fV3bpm_lab[0]},
        {"bpm.l_y", "BPM Y (lab)", kDOUBLE, &fV3bpm_lab[1]},
        {"bpm.l_z", "BPM Z (lab)", kDOUBLE, &fV3bpm_lab[2]},
        {"bpm.x", "BPM X", kDOUBLE, &fV3bpm_tr[0]},
        {"bpm.y", "BPM Y", kDOUBLE, &fV3bpm_tr[1]},
        {"bpm.z", "BPM Z", kDOUBLE, &fV3bpm_tr[2]},
        {"react.x", "React Point X", kDOUBLE, &fV5react_tr[0]},
        {"react.t", "React Point T", kDOUBLE, &fV5react_tr[1]},
        {"react.y", "React Point Y", kDOUBLE, &fV5react_tr[2]},
        {"react.p", "React Point P", kDOUBLE, &fV5react_tr[3]},
        {"react.d", "React Point D", kDOUBLE, &fV5react_tr[4]},
        {"sieve.x", "Sieve X", kDOUBLE, &fV5sieve_tr[0]},
        {"sieve.t", "Sieve T", kDOUBLE, &fV5sieve_tr[1]},
        {"sieve.y", "Sieve Y", kDOUBLE, &fV5sieve_tr[2]},
        {"sieve.p", "Sieve P", kDOUBLE, &fV5sieve_tr[3]},
        {"sieve.d", "Sieve D", kDOUBLE, &fV5sieve_tr[4]},
        {"tp.proj.x", "Project to Target Plane X", kDOUBLE, &fV5tpproj_tr[0]},
        {"tp.proj.t", "Project to Target Plane T", kDOUBLE, &fV5tpproj_tr[1]},
        {"tp.proj.y", "Project to Target Plane Y", kDOUBLE, &fV5tpproj_tr[2]},
        {"tp.proj.p", "Project to Target Plane P", kDOUBLE, &fV5tpproj_tr[3]},
        {"tp.proj.d", "Project to Target Plane D", kDOUBLE, &fV5tpproj_tr[4]},
        {"fp.d_x", "Focus Plane X", kDOUBLE, &fV5fp_det[0]},
        {"fp.d_t", "Focus Plane T", kDOUBLE, &fV5fp_det[1]},
        {"fp.d_y", "Focus Plane Y", kDOUBLE, &fV5fp_det[2]},
        {"fp.d_p", "Focus Plane P", kDOUBLE, &fV5fp_det[3]},
        {0}
    };

    return DefineVarsFromList(vars, mode);
}

void G2POptics::MakePrefix()
{
    const char *basename = "optics";

    G2PAppBase::MakePrefix(basename);
}

ClassImp(G2POptics)
