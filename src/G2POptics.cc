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
#include <cstring>
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

G2POptics::G2POptics(const char *filename) : fDataFile(filename), fE0(0.0), fm(0.0), fM0(0.0), fELoss(0.0), fBPMZ(0), fNFoil(1), fHoleID(-1), pSieve(NULL)
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
    fM0V.clear();
    fM0V.push_back(0.0);
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

    if (fDebug > 2) {
        Info(here, "%04d %10.3e %10.3e %10.3e %10.3e %10.3e", tempdata.ind, tempdata.xf, tempdata.tf, tempdata.yf, tempdata.pf, tempdata.eb);
        Info(here, "%04d %10.3e %10.3e %10.3e %10.3e", tempdata.ind, tempdata.xb, tempdata.tb, tempdata.yb, tempdata.pb);
    }

    fHoleID = tempdata.ind;
    int kineID = fHoleID / (pSieve->GetNRow() * pSieve->GetNCol() * fNFoil);
    int res = fHoleID % (pSieve->GetNRow() * pSieve->GetNCol() * fNFoil);
    int foilID = res / (pSieve->GetNRow() * pSieve->GetNCol());
    res = res % (pSieve->GetNRow() * pSieve->GetNCol());

    fE0 = tempdata.eb / 1000.0;

    fV5bpm_bpm[0] = tempdata.xb / 1000.0;
    fV5bpm_bpm[1] = tempdata.tb;
    fV5bpm_bpm[2] = tempdata.yb / 1000.0;
    fV5bpm_bpm[3] = tempdata.pb;
    fV5bpm_bpm[4] = fBPMZ;

    fV5fp_det[0] = tempdata.xf;
    fV5fp_det[1] = atan(tempdata.tf);
    fV5fp_det[2] = tempdata.yf;
    fV5fp_det[3] = atan(tempdata.pf);
    fV5fp_det[4] = 0.0;

    BPM2HCS(fV5bpm_bpm, fV5bpm_lab);

    if (fDebug > 1)
        Info(here, "bpm_lab   : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5bpm_lab[0], fV5bpm_lab[1], fV5bpm_lab[2], fV5bpm_lab[3], fV5bpm_lab[4]);

    fFoilZ = fFoilZV[foilID];

    double x[3] = {fV5bpm_lab[0], fV5bpm_lab[2], fV5bpm_lab[4]};
    double p[3] = {fE0 * sin(fV5bpm_lab[1]) *cos(fV5bpm_lab[3]), fE0 * sin(fV5bpm_lab[1]) *sin(fV5bpm_lab[3]), fE0 * cos(fV5bpm_lab[1])};

    if (fV5bpm_lab[4] < fFoilZ)
        Drift("forward", x, p, fFoilZ, x, p);
    else if (fV5bpm_lab[4] > fFoilZ)
        Drift("backward", x, p, fFoilZ, x, p);

    fV5bpm_lab[0] = x[0];
    fV5bpm_lab[1] = acos(p[2] / fE0);
    fV5bpm_lab[2] = x[1];
    fV5bpm_lab[3] = atan2(p[1], p[0]);
    fV5bpm_lab[4] = x[2];

    DCS2TRCS(fV5fp_det, fV5fp_tr);
    TRCS2FCS(fV5fp_tr, fV5fp_rot);

    fHRSMomentum = fHRSMomentumV[kineID];
    fM0 = fM0V[foilID];
    fELoss = fELossV[foilID];

    if (fDebug > 1)
        Info(here, "bpm_lab   : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5bpm_lab[0], fV5bpm_lab[1], fV5bpm_lab[2], fV5bpm_lab[3], fV5bpm_lab[4]);

    HCS2TCS(fV5bpm_lab, fV5bpm_tr, fbpmz_tr);

    double V3siv_tr[3];
    pSieve->GetPos(res, V3siv_tr);

    if (fDebug > 1)
        Info(here, "sieve_real: %10.3e %10.3e %10.3e", V3siv_tr[0], V3siv_tr[1], V3siv_tr[2]);

    double V3pd_tr[3];
    V3pd_tr[0] = V3siv_tr[0] - fV5bpm_tr[0];
    V3pd_tr[1] = V3siv_tr[1] - fV5bpm_tr[2];
    V3pd_tr[2] = V3siv_tr[2] - fbpmz_tr;

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

void G2POptics::Clear(Option_t *opt)
{
    fHoleID = -1;
    fFoilZ = 0;
    fbpmz_tr = 0;

    memset(fV5bpm_bpm, 0, sizeof(fV5bpm_bpm));
    memset(fV5bpm_lab, 0, sizeof(fV5bpm_lab));
    memset(fV5bpm_tr, 0, sizeof(fV5bpm_tr));
    memset(fV5react_tr, 0, sizeof(fV5react_tr));
    memset(fV5sieve_tr, 0, sizeof(fV5sieve_tr));
    memset(fV5tpproj_tr, 0, sizeof(fV5tpproj_tr));
    memset(fV5fp_det, 0, sizeof(fV5fp_det));
    memset(fV5fp_tr, 0, sizeof(fV5fp_tr));
    memset(fV5fp_rot, 0, sizeof(fV5fp_rot));

    G2PProcBase::Clear(opt);
}

void G2POptics::SetHRSMomentum(int n, double *value)
{
    fHRSMomentumV.clear();

    for (int i = 0; i < n; i++)
        fHRSMomentumV.push_back(value[i]);
}

void G2POptics::SetBPMZ(double value)
{
    fBPMZ = value;
}

void G2POptics::SetFoilZ(int n, double *value)
{
    fNFoil = n;

    fFoilZV.clear();

    for (int i = 0; i < n; i++)
        fFoilZV.push_back(value[i]);
}

void G2POptics::SetTargetMass(int n, double *value)
{
    fNFoil = n;

    fM0V.clear();

    for (int i = 0; i < n; i++)
        fM0V.push_back(value[i]);
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

    double V3ed_lab[3] = {sin(fV5bpm_lab[1]) *cos(fV5bpm_lab[3]), sin(fV5bpm_lab[1]) *sin(fV5bpm_lab[3]), cos(fV5bpm_lab[1])};

    double cosang = (V3ed_lab[0] * V3pd_lab[0] + V3ed_lab[1] * V3pd_lab[1] + V3ed_lab[2] * V3pd_lab[2]) / sqrt(V3pd_lab[0] * V3pd_lab[0] + V3pd_lab[1] * V3pd_lab[1] + V3pd_lab[2] * V3pd_lab[2]);
    double P = sqrt(fE0 * fE0 - fm * fm);
    double scatmom = (P * fM0 / (fE0 + fM0 - P * cosang)) * (((fE0 + fM0) * sqrt(1 - (fm / fM0) * (fm / fM0) * (1 - cosang * cosang)) + (fE0 + (fm / fM0) * fm) * cosang) / (fE0 + fM0 + P * cosang));
    double delta = (scatmom - fHRSMomentum - fELoss) / fHRSMomentum;

    fV5react_tr[0] = fV5bpm_tr[0];
    fV5react_tr[1] = ang[0];
    fV5react_tr[2] = fV5bpm_tr[2];
    fV5react_tr[3] = ang[1];
    fV5react_tr[4] = delta;

    Drift("forward", fV5react_tr, fbpmz_tr, pSieve->GetZ(), fV5sieve_tr);

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
    fscanf(fp, "%d%lf%lf%lf%lf%lf%lf%lf%lf%lf", &temp.ind, &temp.xf, &temp.tf, &temp.yf, &temp.pf, &temp.eb, &temp.xb, &temp.tb, &temp.yb, &temp.pb);

    while (!feof(fp)) {
        fData.push(temp);
        fscanf(fp, "%d%lf%lf%lf%lf%lf%lf%lf%lf%lf", &temp.ind, &temp.xf, &temp.tf, &temp.yf, &temp.pf, &temp.eb, &temp.xb, &temp.tb, &temp.yb, &temp.pb);
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
        {"particle.mass", "Beam Particle Mass", kDOUBLE, &fm},
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
        {"bpm.b_x", "BPM X (bpm)", kDOUBLE, &fV5bpm_bpm[0]},
        {"bpm.b_t", "BPM T (bpm)", kDOUBLE, &fV5bpm_bpm[1]},
        {"bpm.b_y", "BPM Y (bpm)", kDOUBLE, &fV5bpm_bpm[2]},
        {"bpm.b_p", "BPM P (bpm)", kDOUBLE, &fV5bpm_bpm[3]},
        {"bpm.b_z", "BPM Z (bpm)", kDOUBLE, &fV5bpm_bpm[4]},
        {"bpm.l_x", "BPM X (lab)", kDOUBLE, &fV5bpm_lab[0]},
        {"bpm.l_t", "BPM T (lab)", kDOUBLE, &fV5bpm_lab[1]},
        {"bpm.l_y", "BPM Y (lab)", kDOUBLE, &fV5bpm_lab[2]},
        {"bpm.l_p", "BPM P (lab)", kDOUBLE, &fV5bpm_lab[3]},
        {"bpm.l_z", "BPM Z (lab)", kDOUBLE, &fV5bpm_lab[4]},
        {"bpm.x", "BPM X", kDOUBLE, &fV5bpm_tr[0]},
        {"bpm.t", "BPM T", kDOUBLE, &fV5bpm_tr[1]},
        {"bpm.y", "BPM Y", kDOUBLE, &fV5bpm_tr[2]},
        {"bpm.p", "BPM P", kDOUBLE, &fV5bpm_tr[3]},
        {"bpm.z", "BPM Z", kDOUBLE, &fbpmz_tr},
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
        {"fp.x", "FP X", kDOUBLE, &fV5fp_tr[0]},
        {"fp.t", "FP T", kDOUBLE, &fV5fp_tr[1]},
        {"fp.y", "FP Y", kDOUBLE, &fV5fp_tr[2]},
        {"fp.p", "FP P", kDOUBLE, &fV5fp_tr[3]},
        {"fp.r_x", "FP X (FCS)", kDOUBLE, &fV5fp_rot[0]},
        {"fp.r_t", "FP T (FCS)", kDOUBLE, &fV5fp_rot[1]},
        {"fp.r_y", "FP Y (FCS)", kDOUBLE, &fV5fp_rot[2]},
        {"fp.r_p", "FP P (FCS)", kDOUBLE, &fV5fp_rot[3]},
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
