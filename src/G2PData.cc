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
#include <vector>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"

#include "G2PAppBase.hh"
#include "G2PAppList.hh"
#include "G2PGlobals.hh"
#include "G2PProcBase.hh"

#include "G2PData.hh"

using namespace std;

static const double kDEG = 3.14159265358979323846 / 180.0;

G2PData* G2PData::pG2PData = NULL;

G2PData::G2PData() : fDataFile(NULL)
{
    // Only for ROOT I/O
}

G2PData::G2PData(const char* filename) :
fDataFile(filename), fHRSAngle(5.767 * kDEG)
{
    if (pG2PData) {
        Error("G2PData()", "Only one instance of G2PData allowed.");
        MakeZombie();
        return;
    }
    pG2PData = this;

    fPriority = 1;
    fData.clear();
    Clear();
}

G2PData::~G2PData()
{
    if (pG2PData == this) pG2PData = NULL;

    fData.clear();
    Clear();
}

int G2PData::Begin()
{
    static const char* const here = "Begin()";

    if (G2PProcBase::Begin() != 0) return fStatus;

    if (LoadData() != 0) {
        Error(here, "Cannot read data.");
        return (fStatus = kINITERROR);
    }

    return (fStatus = kOK);
}

int G2PData::Process()
{
    static const char* const here = "Process()";

    if (fDebug > 2) Info(here, " ");

    sData tempdata;

    if (fData.empty()) return -1;
    tempdata = fData.back();
    fV5bpm_lab[0] = tempdata.xb;
    fV5bpm_lab[1] = tempdata.tb;
    fV5bpm_lab[2] = tempdata.yb;
    fV5bpm_lab[3] = tempdata.pb;
    fV5bpm_lab[4] = tempdata.zb;
    fV5fp_tr[0] = tempdata.xf;
    fV5fp_tr[1] = atan(tempdata.tf);
    fV5fp_tr[2] = tempdata.yf;
    fV5fp_tr[3] = atan(tempdata.pf);
    fV5fp_tr[4] = 0.0;

    TRCS2FCS(fV5fp_tr, fHRSAngle, fV5fp_rot);

    if (fDebug > 1) {
        Info(here, "bpm_lab   : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5bpm_lab[0], fV5bpm_lab[1], fV5bpm_lab[2], fV5bpm_lab[3], fV5bpm_lab[4]);
        Info(here, "fp_tr     : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5fp_tr[0], fV5fp_tr[1], fV5fp_tr[2], fV5fp_tr[3], fV5fp_tr[4]);
    }

    fData.pop_back();

    return 0;
}

void G2PData::Clear(Option_t* option)
{
    memset(fV5bpm_lab, 0, sizeof (fV5bpm_lab));
    memset(fV5fp_tr, 0, sizeof (fV5fp_tr));
    memset(fV5fp_rot, 0, sizeof (fV5fp_rot));

    G2PProcBase::Clear(option);
}

int G2PData::LoadData()
{
    FILE *fp;

    if ((fp = fopen(fDataFile, "r")) == NULL) return -1;

    sData temp;
    fscanf(fp, "%d%lf%lf%lf%lf%lf%lf%lf%lf%lf", &temp.ind, &temp.xb, &temp.tb, &temp.yb, &temp.pb, &temp.zb, &temp.xf, &temp.tf, &temp.yf, &temp.pf);
    while (!feof(fp)) {
        fData.push_back(temp);
        fscanf(fp, "%d%lf%lf%lf%lf%lf%lf%lf%lf%lf", &temp.ind, &temp.xb, &temp.tb, &temp.yb, &temp.pb, &temp.zb, &temp.xf, &temp.tf, &temp.yf, &temp.pf);
    }

    fclose(fp);

    if (!fData.empty()) return 0;
    else return -1;
}

int G2PData::Configure(EMode mode)
{
    if (mode == kREAD || mode == kTWOWAY) {
        if (fIsInit) return 0;
        else fIsInit = true;
    }

    ConfDef confs[] = {
        {"run.hrs.angle", "HRS Angle", kDOUBLE, &fHRSAngle},
        {0}
    };

    return ConfigureFromList(confs, mode);
}

int G2PData::DefineVariables(EMode mode)
{
    if (mode == kDEFINE && fIsSetup) return 0;
    fIsSetup = (mode == kDEFINE);

    VarDef vars[] = {
        {"fp.x", "FP X", kDOUBLE, &fV5fp_tr[0]},
        {"fp.t", "FP T", kDOUBLE, &fV5fp_tr[1]},
        {"fp.y", "FP Y", kDOUBLE, &fV5fp_tr[2]},
        {"fp.p", "FP P", kDOUBLE, &fV5fp_tr[3]},
        {"fp.r_x", "FP X (FCS)", kDOUBLE, &fV5fp_rot[0]},
        {"fp.r_t", "FP T (FCS)", kDOUBLE, &fV5fp_rot[1]},
        {"fp.r_y", "FP Y (FCS)", kDOUBLE, &fV5fp_rot[2]},
        {"fp.r_p", "FP P (FCS)", kDOUBLE, &fV5fp_rot[3]},
        {"bpm.l_x", "BPM X (lab)", kDOUBLE, &fV5bpm_lab[0]},
        {"bpm.l_t", "BPM T (lab)", kDOUBLE, &fV5bpm_lab[1]},
        {"bpm.l_y", "BPM Y (lab)", kDOUBLE, &fV5bpm_lab[2]},
        {"bpm.l_p", "BPM P (lab)", kDOUBLE, &fV5bpm_lab[3]},
        {"bpm.l_z", "BPM Z (lab)", kDOUBLE, &fV5bpm_lab[4]},
        {0}
    };

    return DefineVarsFromList(vars, mode);
}

void G2PData::MakePrefix()
{
    const char* base = "data";

    G2PAppBase::MakePrefix(base);
}

ClassImp(G2PData)
