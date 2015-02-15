// -*- C++ -*-

/* class G2PData
 * It loads real data to the simulation.
 */

// History:
//   Mar 2013, C. Gu, First public version.
//   Sep 2013, C. Gu, Rewrite G2PDataGun as a standalone G2PProcBase class.
//

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <queue>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"

#include "G2PAppBase.hh"
#include "G2PAppList.hh"
#include "G2PGlobals.hh"
#include "G2PProcBase.hh"
#include "G2PVarDef.hh"

#include "G2PData.hh"

using namespace std;

G2PData *G2PData::pG2PData = NULL;

G2PData::G2PData()
{
    // Only for ROOT I/O
}

G2PData::G2PData(const char *filename) : fDataFile(filename)
{
    if (pG2PData) {
        Error("G2PData()", "Only one instance of G2PData allowed.");
        MakeZombie();
        return;
    }

    pG2PData = this;

    fPriority = 1;

    Clear();
}

G2PData::~G2PData()
{
    if (pG2PData == this)
        pG2PData = NULL;
}

int G2PData::Begin()
{
    static const char *const here = "Begin()";

    if (G2PProcBase::Begin() != 0)
        return (fStatus = kBEGINERROR);

    if (LoadData() != 0) {
        Error(here, "Cannot read data.");
        return (fStatus = kBEGINERROR);
    }

    return (fStatus = kOK);
}

int G2PData::Process()
{
    static const char *const here = "Process()";

    if (fDebug > 2)
        Info(here, " ");

    sData tempdata;

    if (fData.empty())
        return -1;

    tempdata = fData.front();

    fV5bpm_bpm[0] = tempdata.xb / 1000.0;
    fV5bpm_bpm[1] = tempdata.tb;
    fV5bpm_bpm[2] = tempdata.yb / 1000.0;
    fV5bpm_bpm[3] = tempdata.pb;
    fV5bpm_bpm[4] = 0.0;

    fV5fp_det[0] = tempdata.xf;
    fV5fp_det[1] = atan(tempdata.tf);
    fV5fp_det[2] = tempdata.yf;
    fV5fp_det[3] = atan(tempdata.pf);
    fV5fp_det[4] = 0.0;

    DCS2TRCS(fV5fp_det, fV5fp_tr);
    TRCS2FCS(fV5fp_tr, fV5fp_rot);

    BPM2HCS(fV5bpm_bpm, fV5bpm_lab);
    HCS2TCS(fV5bpm_lab, fV5bpm_tr, fbpmz_tr);

    if (fDebug > 1) {
        Info(here, "bpm_lab   : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5bpm_lab[0], fV5bpm_lab[1], fV5bpm_lab[2], fV5bpm_lab[3], fV5bpm_lab[4]);
        Info(here, "fp_tr     : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5fp_tr[0], fV5fp_tr[1], fV5fp_tr[2], fV5fp_tr[3], fV5fp_tr[4]);
    }

    fData.pop();

    return 0;
}

void G2PData::Clear(Option_t *opt)
{
    fbpmz_tr = 0;

    memset(fV5bpm_bpm, 0, sizeof(fV5bpm_bpm));
    memset(fV5bpm_lab, 0, sizeof(fV5bpm_lab));
    memset(fV5bpm_tr, 0, sizeof(fV5bpm_tr));
    memset(fV5fp_det, 0, sizeof(fV5fp_det));
    memset(fV5fp_tr, 0, sizeof(fV5fp_tr));
    memset(fV5fp_rot, 0, sizeof(fV5fp_rot));

    G2PProcBase::Clear(opt);
}

int G2PData::LoadData()
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

int G2PData::DefineVariables(EMode mode)
{
    if (mode == kDEFINE && fDefined)
        return 0;

    if (G2PProcBase::DefineVariables(mode) != 0)
        return -1;

    VarDef vars[] = {
        {"fp.x", "FP X", kDOUBLE, &fV5fp_tr[0]},
        {"fp.t", "FP T", kDOUBLE, &fV5fp_tr[1]},
        {"fp.y", "FP Y", kDOUBLE, &fV5fp_tr[2]},
        {"fp.p", "FP P", kDOUBLE, &fV5fp_tr[3]},
        {"fp.r_x", "FP X (FCS)", kDOUBLE, &fV5fp_rot[0]},
        {"fp.r_t", "FP T (FCS)", kDOUBLE, &fV5fp_rot[1]},
        {"fp.r_y", "FP Y (FCS)", kDOUBLE, &fV5fp_rot[2]},
        {"fp.r_p", "FP P (FCS)", kDOUBLE, &fV5fp_rot[3]},
        {"bpm.b_x", "BPM X (bpm)", kDOUBLE, &fV5bpm_bpm[0]},
        {"bpm.b_t", "BPM T (bpm)", kDOUBLE, &fV5bpm_bpm[1]},
        {"bpm.b_y", "BPM Y (bpm)", kDOUBLE, &fV5bpm_bpm[2]},
        {"bpm.b_p", "BPM P (bpm)", kDOUBLE, &fV5bpm_bpm[3]},
        {"bpm.b_z", "BPM Z (bpm)", kDOUBLE, &fV5bpm_bpm[4]},
        {"bpm.l_x", "BPM X (bpm)", kDOUBLE, &fV5bpm_lab[0]},
        {"bpm.l_t", "BPM T (lab)", kDOUBLE, &fV5bpm_lab[1]},
        {"bpm.l_y", "BPM Y (lab)", kDOUBLE, &fV5bpm_lab[2]},
        {"bpm.l_p", "BPM P (lab)", kDOUBLE, &fV5bpm_lab[3]},
        {"bpm.l_z", "BPM Z (lab)", kDOUBLE, &fV5bpm_lab[4]},
        {"bpm.x", "BPM X", kDOUBLE, &fV5bpm_tr[0]},
        {"bpm.t", "BPM T", kDOUBLE, &fV5bpm_tr[1]},
        {"bpm.y", "BPM Y", kDOUBLE, &fV5bpm_tr[2]},
        {"bpm.p", "BPM P", kDOUBLE, &fV5bpm_tr[3]},
        {"bpm.z", "BPM Z", kDOUBLE, &fbpmz_tr},
        {0}
    };

    return DefineVarsFromList(vars, mode);
}

void G2PData::MakePrefix()
{
    const char *base = "data";

    G2PAppBase::MakePrefix(base);
}

ClassImp(G2PData)
