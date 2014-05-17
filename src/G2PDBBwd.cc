// -*- C++ -*-

/* class G2PDBBwd
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
#include "G2PDrift.hh"
#include "G2PGlobals.hh"
#include "G2PProcBase.hh"
#include "G2PSieve.hh"
#include "G2PVar.hh"
#include "G2PVarDef.hh"
#include "G2PVarList.hh"

#include "G2PDBBwd.hh"

using namespace std;

G2PDBBwd* G2PDBBwd::pG2PDBBwd = NULL;

G2PDBBwd::G2PDBBwd()
{
    // Only for ROOT I/O
}

G2PDBBwd::G2PDBBwd(const char* name) :
fHRSAngle(0.0), fHRSMomentum(0.0), fFieldRatio(0.0), frecz_lab(0.0), pDrift(NULL), pSieve(NULL), fDBPrefix(NULL), fDBFile(name)
{
    if (pG2PDBBwd) {
        Error("G2PDBBwd()", "Only one instance of G2PDBBwd allowed.");
        MakeZombie();
        return;
    }
    pG2PDBBwd = this;

    memset(fFitPars, 0, sizeof (fFitPars));

    fPriority = 5;

    Clear();
}

G2PDBBwd::~G2PDBBwd()
{
    if (pG2PDBBwd == this) pG2PDBBwd = NULL;
}

int G2PDBBwd::Init()
{
    //static const char* const here = "Init()";

    if (G2PProcBase::Init() != 0) return fStatus;

    pDrift = static_cast<G2PDrift*> (gG2PApps->Find("G2PDrift"));
    if (!pDrift) {
        pDrift = new G2PDrift();
        gG2PApps->Add(pDrift);
    }

    pSieve = static_cast<G2PSieve*> (gG2PApps->Find("G2PSieve"));
    if (!pSieve) {
        pSieve = new G2PSieve();
        gG2PApps->Add(pSieve);
    }

    return (fStatus = kOK);
}

int G2PDBBwd::Begin()
{
    static const char* const here = "Begin()";

    if (G2PAppBase::Begin() != 0) return fStatus;

    if (fHRSAngle > 0) {
        fDBPrefix = "L";
        if (fDBFile == NULL) fDBFile = "db_L.vdc.dat";
    } else {
        fDBPrefix = "R";
        if (fDBFile == NULL) fDBFile = "db_R.vdc.dat";
    }

    // Read VDC database
    ifstream ifs(fDBFile, ios_base::in);
    if (!ifs.good()) {
        Error(here, "Cannot initialize, database file \"%s\" does not exist.", fDBFile);
        return (fStatus = kINITERROR);
    }

    TString tag(fDBPrefix);
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
    while ((!found)&&(ifs.getline(buff, LEN) != 0)) {
        tmpline = ::Compress(buff); //strip blanks

        if (tmpline.EndsWith("\n")) tmpline.Chop();
        tmpline.ToLower();

        if (tag == tmpline) found = true;
    }

    if (!found) {
        Error(here, "Cannot initialize, no matrix in database file \"%s\".", fDBFile);
        ifs.close();
        return (fStatus = kINITERROR);
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

    if (fDebug > 0) Info(here, "Loading matrix from %s ...", fDBFile);

    // Read matrix elements line by line
    while (ifs.getline(buff, LEN) != 0) {
        TString tmpline(buff);

        if (tmpline.EndsWith("\n")) tmpline.Chop();

        istringstream ist(tmpline.Data());
        string tmpstr;
        vector<string> line_spl;
        while (ist >> tmpstr) {
            line_spl.push_back(tmpstr);
        }

        if (line_spl.empty()) continue; // ignore empty lines
        const char* w = line_spl[0].c_str();
        vsiz_t npow = power[w];
        if (npow == 0) break; // stop if the line does not start with a string referring to a known type of matrix element

        THaMatrixElement ME;
        ME.fPower.resize(npow);
        ME.fIsZero = true;
        ME.fOrder = 0;

        vsiz_t pos;
        for (pos = 1; (pos <= npow)&&(pos < line_spl.size()); pos++)
            ME.fPower[pos - 1] = atoi(line_spl[pos].c_str());
        vsiz_t p_cnt;
        for (p_cnt = 0; (pos < line_spl.size())&&(p_cnt < kPORDER)&&(pos <= npow + kPORDER); pos++, p_cnt++) {
            ME.fPoly[p_cnt] = atof(line_spl[pos].c_str());
            if (ME.fPoly[p_cnt] != 0.0) {
                ME.fIsZero = false;
                ME.fOrder = p_cnt + 1;
            }
        }
        if (p_cnt < 1) {
            Error(here, "Cannot Initialize, matrix element %s%d%d%d has error.", w, ME.fPower[0], ME.fPower[1], ME.fPower[2]);
            ifs.close();
            return (fStatus = kINITERROR);
        }

        if (ME.fIsZero) continue;

        vector<THaMatrixElement>* mat = matrix_map[w];
        if (mat) {
            bool match = false;
            for (vector<THaMatrixElement>::iterator it = mat->begin(); (it != mat->end())&&(!(match = it->IsMatch(ME))); it++);
            if (match)
                Warning(here, "Duplicate definition of matrix element %s%d%d%d.", w, ME.fPower[0], ME.fPower[1], ME.fPower[2]);
            else
                mat->push_back(ME);
        }
    }

    fStatus = kOK;
    ifs.close();

    return (fStatus = kOK);
}

int G2PDBBwd::Process()
{
    static const char* const here = "Process()";

    if (fDebug > 2) Info(here, " ");

    fV5bpm_lab[0] = gG2PVars->FindSuffix("bpm.l_x")->GetValue();
    fV5bpm_lab[1] = gG2PVars->FindSuffix("bpm.l_t")->GetValue();
    fV5bpm_lab[2] = gG2PVars->FindSuffix("bpm.l_y")->GetValue();
    fV5bpm_lab[3] = gG2PVars->FindSuffix("bpm.l_p")->GetValue();
    fV5bpm_lab[4] = gG2PVars->FindSuffix("bpm.l_z")->GetValue();

    fV5fp_tr[0] = gG2PVars->FindSuffix("fp.x")->GetValue();
    fV5fp_tr[1] = gG2PVars->FindSuffix("fp.t")->GetValue();
    fV5fp_tr[2] = gG2PVars->FindSuffix("fp.y")->GetValue();
    fV5fp_tr[3] = gG2PVars->FindSuffix("fp.p")->GetValue();

    HCS2TCS(fV5bpm_lab[0], fV5bpm_lab[2], fV5bpm_lab[4], fHRSAngle, fV5bpm_tr[0], fV5bpm_tr[2], fV5bpm_tr[4]);
    HCS2TCS(fV5bpm_lab[1], fV5bpm_lab[3], fHRSAngle, fV5bpm_tr[1], fV5bpm_tr[3]);

    fV5fp_tr[4] = fV5bpm_tr[0];

    if (!Backward(fV5fp_tr, fV5tpmat_tr)) return -1;

    fV5tpmat_tr[0] = GetEffBPM(0); // Get effective bpm x
    fV5tpmat_tr[2] = GetEffBPM(1); // Get effective bpm y

    if (fDebug > 1) {
        Info(here, "tpsnake_tr: %10.3e %10.3e %10.3e %10.3e %10.3e", fV5tpmat_tr[0], fV5tpmat_tr[1], fV5tpmat_tr[2], fV5tpmat_tr[3], fV5tpmat_tr[4]);
    }

    Project(fV5tpmat_tr[0], fV5tpmat_tr[2], 0.0, pSieve->GetZ(), fV5tpmat_tr[1], fV5tpmat_tr[3], fV5sieveproj_tr[0], fV5sieveproj_tr[2]);
    fV5sieveproj_tr[1] = fV5tpmat_tr[1];
    fV5sieveproj_tr[3] = fV5tpmat_tr[3];
    fV5sieveproj_tr[4] = fV5tpmat_tr[4];

    if (fDebug > 1) {
        Info(here, "sivproj_tr: %10.3e %10.3e %10.3e %10.3e %10.3e", fV5sieveproj_tr[0], fV5sieveproj_tr[1], fV5sieveproj_tr[2], fV5sieveproj_tr[3], fV5sieveproj_tr[4]);
    }

    pDrift->Drift(fV5sieveproj_tr, pSieve->GetZ(), fHRSMomentum, fHRSAngle, 0.0, fV5tprec_tr);
    TCS2HCS(fV5tprec_tr[0], fV5tprec_tr[2], 0.0, fHRSAngle, fV5tprec_lab[0], fV5tprec_lab[2], fV5tprec_lab[4]);
    TCS2HCS(fV5tprec_tr[1], fV5tprec_tr[3], fHRSAngle, fV5tprec_lab[1], fV5tprec_lab[3]);

    if (fabs(frecz_lab) > 1.0e-5) {
        double x[3] = {fV5tprec_lab[0], fV5tprec_lab[2], fV5tprec_lab[4]};
        double p[3] = {fHRSMomentum * (1 + fV5tprec_tr[4]) * sin(fV5tprec_lab[1]) * cos(fV5tprec_lab[3]),
            fHRSMomentum * (1 + fV5tprec_tr[4]) * sin(fV5tprec_lab[1]) * sin(fV5tprec_lab[3]),
            fHRSMomentum * (1 + fV5tprec_tr[4]) * cos(fV5tprec_lab[1])};
        pDrift->Drift(x, p, frecz_lab, x, p);
        double z_tr;
        fV5tprec_lab[0] = x[0];
        fV5tprec_lab[1] = acos(p[2] / (fHRSMomentum * (1 + fV5tprec_tr[4])));
        fV5tprec_lab[2] = x[1];
        fV5tprec_lab[3] = atan2(p[1], p[0]);
        fV5tprec_lab[4] = x[2];
        HCS2TCS(x[0], x[1], x[2], fHRSAngle, fV5tprec_tr[0], fV5tprec_tr[2], z_tr);
        HCS2TCS(fV5tprec_lab[1], fV5tprec_lab[2], fHRSAngle, fV5tprec_tr[1], fV5tprec_tr[3]);
    }

    if (fDebug > 1) {
        Info(here, "tprec_tr  : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5tprec_tr[0], fV5tprec_tr[1], fV5tprec_tr[2], fV5tprec_tr[3], fV5tprec_tr[4]);
    }

    return 0;
}

void G2PDBBwd::Clear(Option_t* option)
{
    memset(fV5bpm_lab, 0, sizeof (fV5bpm_lab));
    memset(fV5bpm_tr, 0, sizeof (fV5bpm_tr));
    memset(fV5fp_tr, 0, sizeof (fV5fp_tr));
    memset(fV5tpmat_tr, 0, sizeof (fV5tpmat_tr));
    memset(fV5sieveproj_tr, 0, sizeof (fV5sieveproj_tr));
    memset(fV5tprec_tr, 0, sizeof (fV5tprec_tr));
    memset(fV5tprec_lab, 0, sizeof (fV5tprec_lab));

    G2PProcBase::Clear(option);
}

void G2PDBBwd::SetParsX(const double* pars)
{
    fFitPars[0][0] = pars[0];
    fFitPars[0][1] = pars[1];
    fFitPars[0][2] = pars[2];

    fConfigIsSet.insert((unsigned long) &fFitPars[0][0]);
    fConfigIsSet.insert((unsigned long) &fFitPars[0][1]);
    fConfigIsSet.insert((unsigned long) &fFitPars[0][2]);
}

void G2PDBBwd::SetParsY(const double* pars)
{
    fFitPars[1][0] = pars[0];
    fFitPars[1][1] = pars[1];
    fFitPars[1][2] = pars[2];

    fConfigIsSet.insert((unsigned long) &fFitPars[1][0]);
    fConfigIsSet.insert((unsigned long) &fFitPars[1][1]);
    fConfigIsSet.insert((unsigned long) &fFitPars[1][2]);
}

void G2PDBBwd::SetRecZ(double z)
{
    frecz_lab = z;

    fConfigIsSet.insert((unsigned long) &frecz_lab);
}

double G2PDBBwd::GetEffBPM(int axis)
{
    static const char* const here = "GetEffBPM()";

    double xbpm_tr = fV5bpm_tr[0];
    double ybpm_tr = fV5bpm_tr[2];

    double effbpm_tr;
    if (axis == 0)
        effbpm_tr = xbpm_tr;
    else if (axis == 1)
        effbpm_tr = ybpm_tr;
    else return 1e38;

    if (fFieldRatio > 1e-5) {
        // Fit:
        // (Xbpm_tr-Xeffbpm_tr) vs P
        // ([0]+[1]/x)
        double p = (1 + fV5tpmat_tr[4]) * fHRSMomentum;
        if (axis == 0)
            effbpm_tr = xbpm_tr - (fFitPars[0][0] + (fFitPars[0][1] + fFitPars[0][2] * ybpm_tr) / p) / 1000;
        else if (axis == 1)
            effbpm_tr = ybpm_tr - (fFitPars[0][0] + (fFitPars[0][1] + fFitPars[0][2] * xbpm_tr) / p) / 1000;
    }

    if (fDebug > 2) {
        Info(here, "effbpm_tr :%10.3e", effbpm_tr);
    }

    return effbpm_tr;
}

bool G2PDBBwd::Backward(const double* V5fp_tr, double* V5tp_tr)
{
    static const char* const here = "Backward()";

    double V5fp_rot[5] = {0};

    TRCS2FCS(V5fp_tr, fHRSAngle, V5fp_rot);

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

    if (fDebug > 2) Info(here, "%10.3e %10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e %10.3e %10.3e", V5fp_tr[0], V5fp_tr[1], V5fp_tr[2], V5fp_tr[3], V5tp_tr[0], V5tp_tr[1], V5tp_tr[2], V5tp_tr[3], V5tp_tr[4]);

    return true;
}

int G2PDBBwd::Configure(EMode mode)
{
    if (mode == kREAD || mode == kTWOWAY) {
        if (fIsInit) return 0;
        else fIsInit = true;
    }

    ConfDef confs[] = {
        {"run.hrs.angle", "HRS Angle", kDOUBLE, &fHRSAngle},
        {"run.hrs.p0", "HRS Momentum", kDOUBLE, &fHRSMomentum},
        {"field.ratio", "Field Ratio", kDOUBLE, &fFieldRatio},
        {"fit.x.p0", "Effective X p0", kDOUBLE, &fFitPars[0][0]},
        {"fit.x.p1", "Effective X p1", kDOUBLE, &fFitPars[0][1]},
        {"fit.x.p2", "Effective X p2", kDOUBLE, &fFitPars[0][2]},
        {"fit.y.p0", "Effective Y p0", kDOUBLE, &fFitPars[1][0]},
        {"fit.y.p1", "Effective Y p1", kDOUBLE, &fFitPars[1][1]},
        {"fit.y.p2", "Effective Y p2", kDOUBLE, &fFitPars[1][2]},
        {"l_z", "Rec Z (lab)", kDOUBLE, &frecz_lab},
        {0}
    };

    return ConfigureFromList(confs, mode);
}

int G2PDBBwd::DefineVariables(EMode mode)
{
    if (mode == kDEFINE && fIsSetup) return 0;
    fIsSetup = (mode == kDEFINE);

    VarDef vars[] = {
        {"tp.mat.x", "MATRIX rec to Target Plane X", kDOUBLE, &fV5tpmat_tr[0]},
        {"tp.mat.t", "MATRIX rec to Target Plane T", kDOUBLE, &fV5tpmat_tr[1]},
        {"tp.mat.y", "MATRIX rec to Target Plane Y", kDOUBLE, &fV5tpmat_tr[2]},
        {"tp.mat.p", "MATRIX rec to Target Plane P", kDOUBLE, &fV5tpmat_tr[3]},
        {"tp.mat.d", "MATRIX rec to Target Plane D", kDOUBLE, &fV5tpmat_tr[4]},
        {"sieve.proj.x", "Project to Sieve X", kDOUBLE, &fV5sieveproj_tr[0]},
        {"sieve.proj.t", "Project to Sieve T", kDOUBLE, &fV5sieveproj_tr[1]},
        {"sieve.proj.y", "Project to Sieve Y", kDOUBLE, &fV5sieveproj_tr[2]},
        {"sieve.proj.p", "Project to Sieve P", kDOUBLE, &fV5sieveproj_tr[3]},
        {"sieve.proj.d", "Project to Sieve D", kDOUBLE, &fV5sieveproj_tr[3]},
        {"tp.rec.x", "Rec X", kDOUBLE, &fV5tprec_tr[0]},
        {"tp.rec.t", "Rec T", kDOUBLE, &fV5tprec_tr[1]},
        {"tp.rec.y", "Rec Y", kDOUBLE, &fV5tprec_tr[2]},
        {"tp.rec.p", "Rec P", kDOUBLE, &fV5tprec_tr[3]},
        {"tp.rec.d", "Rec D", kDOUBLE, &fV5tprec_tr[4]},
        {"tp.rec.l_x", "Rec X (lab)", kDOUBLE, &fV5tprec_lab[0]},
        {"tp.rec.l_t", "Rec T (lab)", kDOUBLE, &fV5tprec_lab[1]},
        {"tp.rec.l_y", "Rec Y (lab)", kDOUBLE, &fV5tprec_lab[2]},
        {"tp.rec.l_p", "Rec P (lab)", kDOUBLE, &fV5tprec_lab[3]},
        {"tp.rec.l_z", "Rec Z (lab)", kDOUBLE, &fV5tprec_lab[4]},
        {0}
    };

    return DefineVarsFromList(vars, mode);
}

void G2PDBBwd::MakePrefix()
{
    const char* base = "bwd";

    G2PAppBase::MakePrefix(base);
}

double G2PDBBwd::CalcVar(const double powers[][5], vector<THaMatrixElement> &matrix)
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

void G2PDBBwd::CalcMatrix(const double x, vector<THaMatrixElement> &matrix)
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

G2PDBBwd::THaMatrixElement::THaMatrixElement() :
fIsZero(true), fPower(3), fOrder(0), fPoly(G2PDBBwd::kPORDER), fValue(0)
{
    // Nothing to do
}

G2PDBBwd::THaMatrixElement::~THaMatrixElement()
{
    // Nothing to do
}

bool G2PDBBwd::THaMatrixElement::IsMatch(const THaMatrixElement& rhs) const
{
    // Compare coefficients of this matrix element to another
    if (fPower.size() != rhs.fPower.size()) return false;
    for (vector<int>::size_type i = 0; i < fPower.size(); i++) {
        if (fPower[i] != rhs.fPower[i]) return false;
    }
    return true;
}

void G2PDBBwd::THaMatrixElement::SkimPoly()
{
    if (fIsZero) return;

    while ((!fPoly[fOrder - 1])&&(fOrder > 0)) {
        fPoly.pop_back();
        fOrder = fOrder - 1;
    }
    if (fOrder == 0) fIsZero = true;
}

ClassImp(G2PDBBwd)
