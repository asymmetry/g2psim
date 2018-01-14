// -*- C++ -*-

/* class G2PBwdDB
 * Use the analyzer database to reconstruct target variables.
 * Several functions is developed from J. Huang's HRS optics class.
 */

// History:
//   Jun 2010, J. Huang, HRS optics matrix optimization class.
//   Jan 2013, C. Gu, First public version.
//   Sep 2013, C. Gu, Rewrite it as a G2PProcBase class.
//   Nov 2013, C. Gu, Change back to a G2PAppBase class.
//   Apr 2014, C. Gu, Apply effective bpm fitting.
//

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TString.h"

#include "G2PAppBase.hh"
#include "G2PAppList.hh"
#include "G2PGlobals.hh"
#include "G2PProcBase.hh"
#include "G2PSieve.hh"
#include "G2PVar.hh"
#include "G2PVarDef.hh"
#include "G2PVarList.hh"

#include "G2PBwdDB.hh"

using namespace std;

G2PBwdDB *G2PBwdDB::pG2PBwdDB = NULL;

G2PBwdDB::G2PBwdDB()
{
    // Only for ROOT I/O
}

G2PBwdDB::G2PBwdDB(const char *name) : fDBFile(name), fE(0.0), fFieldRatio(0.0), frecz_lab(0.0)
{
    if (pG2PBwdDB) {
        Error("G2PBwdDB()", "Only one instance of G2PBwdDB allowed.");
        MakeZombie();
        return;
    }

    pG2PBwdDB = this;

    memset(fFitPars, 0, sizeof(fFitPars));
    memset(fCorT, 0, sizeof(fCorT));
    memset(fCorP, 0, sizeof(fCorP));
    memset(fCorD, 0, sizeof(fCorD));

    fPriority = 5;

    Clear();
}

G2PBwdDB::~G2PBwdDB()
{
    if (pG2PBwdDB == this)
        pG2PBwdDB = NULL;
}

int G2PBwdDB::Begin()
{
    static const char *const here = "Begin()";

    if (G2PProcBase::Begin() != 0)
        return (fStatus = kBEGINERROR);

    pSieve = static_cast<G2PSieve *>(gG2PApps->Find("G2PSieve"));

    const char *dbprefix = "L";

    if (fHRSAngle > 0) {
        if (fDBFile == NULL)
            fDBFile = "db_L.vdc.dat";
    } else {
        dbprefix = "R";

        if (fDBFile == NULL)
            fDBFile = "db_R.vdc.dat";
    }

    // Read VDC database
    ifstream ifs(fDBFile, ios_base::in);

    if (!ifs.good()) {
        Error(here, "Cannot initialize, database file \"%s\" does not exist.", fDBFile);
        return (fStatus = kBEGINERROR);
    }

    TString tag(dbprefix);
    Ssiz_t pos = tag.Index(".");

    if (pos != kNPOS)
        tag = tag(0, pos + 1);
    else
        tag.Append(".");

    tag.Prepend("[");
    tag.Append("global]");
    tag.ToLower();

    const int LEN = 300;
    char buff[LEN];
    TString tmpline;

    // Locate the matrix label in the database file
    bool found = false;

    while ((!found) && (ifs.getline(buff, LEN))) {
        tmpline = ::Compress(buff); //strip blanks

        if (tmpline.EndsWith("\n"))
            tmpline.Chop();

        tmpline.ToLower();

        if (tag == tmpline)
            found = true;
    }

    if (!found) {
        Error(here, "Cannot initialize, no matrix in database file \"%s\".", fDBFile);
        ifs.close();
        return (fStatus = kBEGINERROR);
    }

    ifs.getline(buff, LEN);
    ifs.getline(buff, LEN); // skip 2 comment line

    ftMatrixElems.clear();
    fyMatrixElems.clear();
    fpMatrixElems.clear();
    fTMatrixElems.clear();
    fDMatrixElems.clear();
    fPMatrixElems.clear();
    fPTAMatrixElems.clear();
    fYMatrixElems.clear();
    fYTAMatrixElems.clear();

    // Each matrix element is assigned a integer vector to indicate its power
    // This map is the length of this vector
    typedef vector<string>::size_type vsiz_t;
    map<string, vsiz_t> power;
    power["t"] = 3; // transport to focal-plane tensors
    power["y"] = 3;
    power["p"] = 3;
    power["D"] = 3; // focal-plane to target tensors
    power["T"] = 3;
    power["Y"] = 3;
    power["YTA"] = 4;
    power["P"] = 3;
    power["PTA"] = 4;
    power["L"] = 4; // pathlength from z=0 (target) to focal plane (meters)
    power["XF"] = 5; // forward: target to focal-plane (I think)
    power["TF"] = 5;
    power["PF"] = 5;
    power["YF"] = 5;

    map<string, vector<THaMatrixElement>*> matrix_map;
    matrix_map["t"] = &ftMatrixElems;
    matrix_map["y"] = &fyMatrixElems;
    matrix_map["p"] = &fpMatrixElems;
    matrix_map["D"] = &fDMatrixElems;
    matrix_map["T"] = &fTMatrixElems;
    matrix_map["Y"] = &fYMatrixElems;
    matrix_map["YTA"] = &fYTAMatrixElems;
    matrix_map["P"] = &fPMatrixElems;
    matrix_map["PTA"] = &fPTAMatrixElems;

    if (fDebug > 0)
        Info(here, "Loading matrix from %s ...", fDBFile);

    // Read matrix elements line by line
    while (ifs.getline(buff, LEN)) {
        TString tmpline(buff);

        if (tmpline.EndsWith("\n"))
            tmpline.Chop();

        istringstream ist(tmpline.Data());
        string tmpstr;
        vector<string> line_spl;

        while (ist >> tmpstr)
            line_spl.push_back(tmpstr);

        if (line_spl.empty()) {
            continue;    // ignore empty lines
        }

        const char *w = line_spl[0].c_str();
        vsiz_t npow = power[w];

        if (npow == 0) {
            break;    // stop if the line does not start with a string referring to a known type of matrix element
        }

        THaMatrixElement ME;
        ME.fPower.resize(npow);
        ME.fIsZero = true;
        ME.fOrder = 0;

        vsiz_t pos;

        for (pos = 1; (pos <= npow) && (pos < line_spl.size()); pos++)
            ME.fPower[pos - 1] = atoi(line_spl[pos].c_str());

        vsiz_t p_cnt;

        for (p_cnt = 0; (pos < line_spl.size()) && (p_cnt < kPORDER) && (pos <= npow + kPORDER); pos++, p_cnt++) {
            ME.fPoly[p_cnt] = atof(line_spl[pos].c_str());

            if (ME.fPoly[p_cnt] != 0.0) {
                ME.fIsZero = false;
                ME.fOrder = p_cnt + 1;
            }
        }

        if (p_cnt < 1) {
            Error(here, "Cannot Initialize, matrix element %s%d%d%d has error.", w, ME.fPower[0], ME.fPower[1], ME.fPower[2]);
            ifs.close();
            return (fStatus = kBEGINERROR);
        }

        if (ME.fIsZero)
            continue;

        vector<THaMatrixElement> *mat = matrix_map[w];

        if (mat) {
            bool match = false;

            for (vector<THaMatrixElement>::iterator it = mat->begin(); (it != mat->end()) && (!(match = it->IsMatch(ME))); it++);

            if (match)
                Warning(here, "Duplicate definition of matrix element %s%d%d%d.", w, ME.fPower[0], ME.fPower[1], ME.fPower[2]);
            else
                mat->push_back(ME);
        }
    }

    ifs.close();

    return (fStatus = kOK);
}

int G2PBwdDB::Process()
{
    static const char *const here = "Process()";

    if (fDebug > 2)
        Info(here, " ");

    if (gG2PVars->FindSuffix("bpm.x") && gG2PVars->FindSuffix("fp.x")) {
        fE = gG2PVars->FindSuffix("phys.e")->GetValue();

        fV5bpm_tr[0] = gG2PVars->FindSuffix("bpm.x")->GetValue();
        fV5bpm_tr[1] = gG2PVars->FindSuffix("bpm.t")->GetValue();
        fV5bpm_tr[2] = gG2PVars->FindSuffix("bpm.y")->GetValue();
        fV5bpm_tr[3] = gG2PVars->FindSuffix("bpm.p")->GetValue();
        fbpmz_tr = gG2PVars->FindSuffix("bpm.z")->GetValue();

        fV5fp_tr[0] = gG2PVars->FindSuffix("fp.x")->GetValue();
        fV5fp_tr[1] = gG2PVars->FindSuffix("fp.t")->GetValue();
        fV5fp_tr[2] = gG2PVars->FindSuffix("fp.y")->GetValue();
        fV5fp_tr[3] = gG2PVars->FindSuffix("fp.p")->GetValue();
    } else
        return -1;

    fV5bpm_tr[4] = fE / fHRSMomentum - 1;

    int save = fDebug;
    fDebug = 0;

    if (fbpmz_tr < 0.0)
        Drift("forward", fV5bpm_tr, fbpmz_tr, 0.0, fV5bpm_tr); // Drift to target plane (z_tr = 0)
    else
        Drift("backward", fV5bpm_tr, fbpmz_tr, 0.0, fV5bpm_tr);

    fDebug = save;

    double bpm_temp[5];
    fV5fp_tr[4] = fV5bpm_tr[0];

    if (!Backward(fV5fp_tr, fV5tpmat_tr))
        return -1;

    if (fDebug > 1)
        Info(here, "tpmat_tr  : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5tpmat_tr[0], fV5tpmat_tr[1], fV5tpmat_tr[2], fV5tpmat_tr[3], fV5tpmat_tr[4]);

    // first iteration
    GetEffBPM(fV5tpmat_tr, fV5bpm_tr, bpm_temp);
    Correct(bpm_temp, fV5tpmat_tr, fV5tpcorr_tr);

    // second iteration
    GetEffBPM(fV5tpcorr_tr, fV5bpm_tr, bpm_temp);
    Correct(bpm_temp, fV5tpmat_tr, fV5tpcorr_tr);

    // third iteration
    GetEffBPM(fV5tpcorr_tr, fV5bpm_tr, bpm_temp);
    Correct(bpm_temp, fV5tpmat_tr, fV5tpcorr_tr);

    GetEffBPM(fV5tpcorr_tr, fV5bpm_tr, bpm_temp);
    fV5tpcorr_tr[0] = bpm_temp[0];
    fV5tpcorr_tr[2] = bpm_temp[2];

    if (fDebug > 1)
        Info(here, "tpcorr_tr : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5tpcorr_tr[0], fV5tpcorr_tr[1], fV5tpcorr_tr[2], fV5tpcorr_tr[3], fV5tpcorr_tr[4]);

    Project(fV5tpcorr_tr, 0.0, pSieve->GetZ(), fV5sieveproj_tr);

    if (fDebug > 1)
        Info(here, "sivproj_tr: %10.3e %10.3e %10.3e %10.3e %10.3e", fV5sieveproj_tr[0], fV5sieveproj_tr[1], fV5sieveproj_tr[2], fV5sieveproj_tr[3], fV5sieveproj_tr[4]);

    Drift("backward", fV5sieveproj_tr, pSieve->GetZ(), 0.0, fV5tprec_tr);
    TCS2HCS(fV5tprec_tr, 0.0, fV5tprec_lab);

    if (fDebug > 1) {
        Info(here, "tprec_tr  : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5tprec_tr[0], fV5tprec_tr[1], fV5tprec_tr[2], fV5tprec_tr[3], fV5tprec_tr[4]);
        Info(here, "tprec_lab : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5tprec_lab[0], fV5tprec_lab[1], fV5tprec_lab[2], fV5tprec_lab[3], fV5tprec_lab[4]);
    }

    for (int i = 0; i < 5; i++) {
        fV5rec_tr[i] = fV5tprec_tr[i];
        fV5rec_lab[i] = fV5tprec_lab[i];
    }

    if (fabs(fV5tprec_lab[4]) > 1.0e-5) {
        double z_tr;

        if (fV5tprec_lab[4] < frecz_lab)
            Drift("forward", fV5tprec_tr, 0.0, frecz_lab, fV5rec_tr, z_tr);
        else
            Drift("backward", fV5tprec_tr, 0.0, frecz_lab, fV5rec_tr, z_tr);

        TCS2HCS(fV5rec_tr, z_tr, fV5rec_lab);
    }

    if (fDebug > 1)
        Info(here, "rec_tr    : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5rec_tr[0], fV5rec_tr[1], fV5rec_tr[2], fV5rec_tr[3], fV5rec_tr[4]);

    return 0;
}

void G2PBwdDB::Clear(Option_t *opt)
{
    fbpmz_tr = 0;

    memset(fV5bpm_tr, 0, sizeof(fV5bpm_tr));
    memset(fV5fp_tr, 0, sizeof(fV5fp_tr));
    memset(fV5tpmat_tr, 0, sizeof(fV5tpmat_tr));
    memset(fV5tpcorr_tr, 0, sizeof(fV5tpcorr_tr));
    memset(fV5sieveproj_tr, 0, sizeof(fV5sieveproj_tr));
    memset(fV5tprec_tr, 0, sizeof(fV5tprec_tr));
    memset(fV5tprec_lab, 0, sizeof(fV5tprec_lab));
    memset(fV5rec_tr, 0, sizeof(fV5rec_tr));
    memset(fV5rec_lab, 0, sizeof(fV5rec_lab));

    G2PProcBase::Clear(opt);
}

void G2PBwdDB::SetParsX(const double *pars)
{
    fFitPars[0][0] = pars[0];
    fFitPars[0][1] = pars[1];
    fFitPars[0][2] = pars[2];

    fConfigIsSet.insert((unsigned long) &fFitPars[0][0]);
    fConfigIsSet.insert((unsigned long) &fFitPars[0][1]);
    fConfigIsSet.insert((unsigned long) &fFitPars[0][2]);
}

void G2PBwdDB::SetParsY(const double *pars)
{
    fFitPars[1][0] = pars[0];
    fFitPars[1][1] = pars[1];
    fFitPars[1][2] = pars[2];

    fConfigIsSet.insert((unsigned long) &fFitPars[1][0]);
    fConfigIsSet.insert((unsigned long) &fFitPars[1][1]);
    fConfigIsSet.insert((unsigned long) &fFitPars[1][2]);
}

void G2PBwdDB::SetRecZ(double z)
{
    frecz_lab = z;

    fConfigIsSet.insert((unsigned long) &frecz_lab);
}

void G2PBwdDB::GetEffBPM(const double *V5tp_tr, const double *V5bpm_tr, double *V5bpmeff_tr)
{
    static const char *const here = "GetEffBPM()";

    V5bpmeff_tr[0] = V5bpm_tr[0];
    V5bpmeff_tr[1] = V5bpm_tr[1];
    V5bpmeff_tr[2] = V5bpm_tr[2];
    V5bpmeff_tr[3] = V5bpm_tr[3];
    V5bpmeff_tr[4] = V5bpm_tr[4];

    if (fFieldRatio > 1e-5) {
        // Fit:
        // (Xbpm_tr-Xeffbpm_tr) vs P
        // ([0]+[1]/x)
        double p = (1 + V5tp_tr[4]) * fHRSMomentum;
        V5bpmeff_tr[0] = V5bpm_tr[0] - (fFitPars[0][0] + (fFitPars[0][1] + fFitPars[0][2] * V5bpm_tr[2]) / p) / 1000;
        V5bpmeff_tr[2] = V5bpm_tr[2] - (fFitPars[1][0] + (fFitPars[1][1] + fFitPars[1][2] * V5bpm_tr[0]) / p) / 1000;
    }

    if (fDebug > 2)
        Info(here, "effbpm_tr : %10.3e %10.3e", V5bpmeff_tr[0], V5bpmeff_tr[2]);
}

void G2PBwdDB::Correct(const double *V5bpm_tr, const double *V5tp_tr, double *V5corr_tr)
{
    V5corr_tr[0] = V5tp_tr[0];
    V5corr_tr[1] = V5tp_tr[1] + fCorT[0] + fCorT[1] * V5bpm_tr[0] + fCorT[2] * V5bpm_tr[2];
    V5corr_tr[2] = V5tp_tr[2];
    V5corr_tr[3] = V5tp_tr[3] + fCorP[0] + fCorP[1] * V5bpm_tr[0] + fCorP[2] * V5bpm_tr[2];
    V5corr_tr[4] = V5tp_tr[4] + fCorD[0] + fCorD[1] * V5bpm_tr[0] + fCorD[2] * V5bpm_tr[2];
}

double G2PBwdDB::CalcVar(const double powers[][5], vector<THaMatrixElement> &matrix)
{
    // Calculate the value of a variable at the target
    // Must already have x values for the matrix elements
    double value = 0.0;
    double v = 0.0;

    for (vector<THaMatrixElement>::iterator it = matrix.begin(); it != matrix.end(); it++) {
        if (it->fValue != 0.0) {
            v = it->fValue;
            vector<int>::size_type np = it->fPower.size();

            for (vector<int>::size_type i = 0; i < np; i++)
                v *= powers[it->fPower[i]][i + 1];

            value += v;
        }
    }

    return value;
}

void G2PBwdDB::CalcMatrix(const double x, vector<THaMatrixElement> &matrix)
{
    // Calculate the value of a matrix element for a given x
    double value = 0.0;

    for (vector<THaMatrixElement>::iterator it = matrix.begin(); it != matrix.end(); it++) {
        value = 0.0;

        if (it->fOrder > 0) {
            for (int i = it->fOrder - 1; i >= 1; i--)
                value = x * (value + it->fPoly[i]);

            value += it->fPoly[0];
        }

        it->fValue = value;
    }
}

bool G2PBwdDB::Backward(const double *V5fp_tr, double *V5tp_tr)
{
    static const char *const here = "Backward()";

    double V5fp_rot[5] = {0};

    TRCS2FCS(V5fp_tr, V5fp_rot);

    // Calculate target coordinates from focal plane coordinates
    double x_fp, y_fp, t_fp, p_fp;
    double powers[kNUM_PRECOMP_POW][5];

    x_fp = V5fp_rot[0];
    t_fp = tan(V5fp_rot[1]);
    y_fp = V5fp_rot[2];
    p_fp = tan(V5fp_rot[3]);

    for (int i = 0; i < kNUM_PRECOMP_POW; i++) {
        powers[i][0] = pow(x_fp, i);
        powers[i][1] = pow(t_fp, i);
        powers[i][2] = pow(y_fp, i);
        powers[i][3] = pow(p_fp, i);
        powers[i][4] = pow(fabs(t_fp), i);
    }

    CalcMatrix(x_fp, fDMatrixElems);
    CalcMatrix(x_fp, fTMatrixElems);
    CalcMatrix(x_fp, fPMatrixElems);
    CalcMatrix(x_fp, fYMatrixElems);
    CalcMatrix(x_fp, fPTAMatrixElems);
    CalcMatrix(x_fp, fYTAMatrixElems);

    V5tp_tr[0] = 0.0;
    V5tp_tr[1] = atan(CalcVar(powers, fTMatrixElems));
    V5tp_tr[2] = CalcVar(powers, fYMatrixElems);
    V5tp_tr[3] = atan(CalcVar(powers, fPMatrixElems));
    V5tp_tr[4] = CalcVar(powers, fDMatrixElems);

    if (fDebug > 2)
        Info(here, "%10.3e %10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e %10.3e %10.3e", V5fp_tr[0], V5fp_tr[1], V5fp_tr[2], V5fp_tr[3], V5tp_tr[0], V5tp_tr[1], V5tp_tr[2], V5tp_tr[3], V5tp_tr[4]);

    return true;
}

int G2PBwdDB::Configure(EMode mode)
{
    if ((mode == kREAD || mode == kTWOWAY) && fConfigured)
        return 0;

    if (G2PProcBase::Configure(mode) != 0)
        return -1;

    ConfDef confs[] = {
        {"field.ratio", "Field Ratio", kDOUBLE, &fFieldRatio},
        {"x.p0", "Effective X p0", kDOUBLE, &fFitPars[0][0]},
        {"x.p1", "Effective X p1", kDOUBLE, &fFitPars[0][1]},
        {"x.p2", "Effective X p2", kDOUBLE, &fFitPars[0][2]},
        {"y.p0", "Effective Y p0", kDOUBLE, &fFitPars[1][0]},
        {"y.p1", "Effective Y p1", kDOUBLE, &fFitPars[1][1]},
        {"y.p2", "Effective Y p2", kDOUBLE, &fFitPars[1][2]},
        {"t.p0", "T Const Correction ", kDOUBLE, &fCorT[0]},
        {"t.px", "T Correction vs Target X", kDOUBLE, &fCorT[1]},
        {"t.py", "T Correction vs Target Y", kDOUBLE, &fCorT[2]},
        {"p.p0", "P Const Correction", kDOUBLE, &fCorP[0]},
        {"p.px", "P Correction vs Target X", kDOUBLE, &fCorP[1]},
        {"p.py", "P Correction vs Target Y", kDOUBLE, &fCorP[2]},
        {"d.p0", "D Const Correction", kDOUBLE, &fCorD[0]},
        {"d.px", "D Correction vs Target X", kDOUBLE, &fCorD[1]},
        {"d.py", "D Correction vs Target Y", kDOUBLE, &fCorD[2]},
        {"z", "Rec Z (lab)", kDOUBLE, &frecz_lab},
        {0}
    };

    return ConfigureFromList("rec.", confs, mode);
}

int G2PBwdDB::DefineVariables(EMode mode)
{
    if (mode == kDEFINE && fDefined)
        return 0;

    if (G2PProcBase::DefineVariables(mode) != 0)
        return -1;

    VarDef gvars[] = {
        {"x", "Rec X", kDOUBLE, &fV5rec_tr[0]},
        {"t", "Rec T", kDOUBLE, &fV5rec_tr[1]},
        {"y", "Rec Y", kDOUBLE, &fV5rec_tr[2]},
        {"p", "Rec P", kDOUBLE, &fV5rec_tr[3]},
        {"d", "Rec D", kDOUBLE, &fV5rec_tr[4]},
        {"l_x", "Rec X (lab)", kDOUBLE, &fV5rec_lab[0]},
        {"l_t", "Rec T (lab)", kDOUBLE, &fV5rec_lab[1]},
        {"l_y", "Rec Y (lab)", kDOUBLE, &fV5rec_lab[2]},
        {"l_p", "Rec P (lab)", kDOUBLE, &fV5rec_lab[3]},
        {"l_z", "Rec Z (lab)", kDOUBLE, &fV5rec_lab[4]},
        {0}
    };

    if (DefineVarsFromList("rec.", gvars, mode) != 0)
        return -1;

    VarDef vars[] = {
        {"tp.mat.x", "Matrix Rec to Target Plane X", kDOUBLE, &fV5tpmat_tr[0]},
        {"tp.mat.t", "Matrix Rec to Target Plane T", kDOUBLE, &fV5tpmat_tr[1]},
        {"tp.mat.y", "Matrix Rec to Target Plane Y", kDOUBLE, &fV5tpmat_tr[2]},
        {"tp.mat.p", "Matrix Rec to Target Plane P", kDOUBLE, &fV5tpmat_tr[3]},
        {"tp.mat.d", "Matrix Rec to Target Plane D", kDOUBLE, &fV5tpmat_tr[4]},
        {"tp.corr.x", "After Correction X", kDOUBLE, &fV5tpcorr_tr[0]},
        {"tp.corr.t", "After Correction T", kDOUBLE, &fV5tpcorr_tr[1]},
        {"tp.corr.y", "After Correction Y", kDOUBLE, &fV5tpcorr_tr[2]},
        {"tp.corr.p", "After Correction P", kDOUBLE, &fV5tpcorr_tr[3]},
        {"tp.corr.d", "After Correction D", kDOUBLE, &fV5tpcorr_tr[4]},
        {"sieve.proj.x", "Project to Sieve X", kDOUBLE, &fV5sieveproj_tr[0]},
        {"sieve.proj.t", "Project to Sieve T", kDOUBLE, &fV5sieveproj_tr[1]},
        {"sieve.proj.y", "Project to Sieve Y", kDOUBLE, &fV5sieveproj_tr[2]},
        {"sieve.proj.p", "Project to Sieve P", kDOUBLE, &fV5sieveproj_tr[3]},
        {"sieve.proj.d", "Project to Sieve D", kDOUBLE, &fV5sieveproj_tr[3]},
        {"tp.rec.x", "Rec to Target Plane X", kDOUBLE, &fV5tprec_tr[0]},
        {"tp.rec.t", "Rec to Target Plane T", kDOUBLE, &fV5tprec_tr[1]},
        {"tp.rec.y", "Rec to Target Plane Y", kDOUBLE, &fV5tprec_tr[2]},
        {"tp.rec.p", "Rec to Target Plane P", kDOUBLE, &fV5tprec_tr[3]},
        {"tp.rec.d", "Rec to Target Plane D", kDOUBLE, &fV5tprec_tr[4]},
        {0}
    };

    return DefineVarsFromList(vars, mode);
}

void G2PBwdDB::MakePrefix()
{
    const char *base = "bwd";

    G2PAppBase::MakePrefix(base);
}

G2PBwdDB::THaMatrixElement::THaMatrixElement() :
    fIsZero(true), fPower(3), fOrder(0), fPoly(G2PBwdDB::kPORDER), fValue(0)
{
    // Nothing to do
}

G2PBwdDB::THaMatrixElement::~THaMatrixElement()
{
    // Nothing to do
}

bool G2PBwdDB::THaMatrixElement::IsMatch(const THaMatrixElement &rhs) const
{
    // Compare coefficients of this matrix element to another
    if (fPower.size() != rhs.fPower.size())
        return false;

    for (vector<int>::size_type i = 0; i < fPower.size(); i++) {
        if (fPower[i] != rhs.fPower[i])
            return false;
    }

    return true;
}

void G2PBwdDB::THaMatrixElement::SkimPoly()
{
    if (fIsZero)
        return;

    while ((!fPoly[fOrder - 1]) && (fOrder > 0)) {
        fPoly.pop_back();
        fOrder = fOrder - 1;
    }

    if (fOrder == 0)
        fIsZero = true;
}

ClassImp(G2PBwdDB)
