// -*- C++ -*-

/* class G2PDBRec
 * Use the analyzer database to reconstruct target variables.
 * Several functions is developed from J. Huang's HRS optics class.
 */

// History:
//   Jun 2010, J. Huang, HRS optics matrix optimization class.
//   Jan 2013, C. Gu, First public version.
//   Sep 2013, C. Gu, Rewrite it as a G2PProcBase class.
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
#include "G2PVar.hh"
#include "G2PVarDef.hh"
#include "G2PVarList.hh"

#include "G2PDBRec.hh"

using namespace std;

static const double kDEG = 3.14159265358979323846 / 180.0;

G2PDBRec* G2PDBRec::pG2PDBRec = NULL;

G2PDBRec::G2PDBRec() :
fDBPrefix(NULL), fDBFile(NULL), fHRSAngle(5.767 * kDEG) {
    if (pG2PDBRec) {
        Error("G2PDBRec()", "Only one instance of G2PDBRec allowed.");
        MakeZombie();
        return;
    }
    pG2PDBRec = this;

    fPriority = 5;
    Clear();
}

G2PDBRec::~G2PDBRec() {
    if (pG2PDBRec == this) pG2PDBRec = NULL;
}

int G2PDBRec::Begin() {
    static const char* const here = "Begin()";

    if (G2PProcBase::Begin() != 0) return fStatus;

    if (fHRSAngle > 0) {
        fDBPrefix = "L";
        fDBFile = "db_L.vdc.dat";
    }
    else {
        fDBPrefix = "R";
        fDBFile = "db_R.vdc.dat";
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

int G2PDBRec::Process() {
    static const char* const here = "Process()";

    if (fDebug > 2) Info(here, " ");

    fV5fp_tr[0] = gG2PVars->FindSuffix("fp.x")->GetValue();
    fV5fp_tr[1] = gG2PVars->FindSuffix("fp.t")->GetValue();
    fV5fp_tr[2] = gG2PVars->FindSuffix("fp.y")->GetValue();
    fV5fp_tr[3] = gG2PVars->FindSuffix("fp.p")->GetValue();

    TRCS2DCS(fV5fp_tr, fHRSAngle, fV5fp_rot);

    // Calculate target coordinates from focal plane coordinates
    double x_fp, y_fp, t_fp, p_fp;
    double powers[kNUM_PRECOMP_POW][5];

    x_fp = fV5fp_rot[0];
    t_fp = tan(fV5fp_rot[1]);
    y_fp = fV5fp_rot[2];
    p_fp = tan(fV5fp_rot[3]);

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

    fV5rec_tr[0] = 0.0;
    fV5rec_tr[1] = atan(CalcVar(powers, fTMatrixElems));
    fV5rec_tr[2] = CalcVar(powers, fYMatrixElems);
    fV5rec_tr[3] = atan(CalcVar(powers, fPMatrixElems));
    fV5rec_tr[4] = CalcVar(powers, fDMatrixElems);

    TCS2HCS(fV5rec_tr[0], fV5rec_tr[2], 0.0, fHRSAngle, fV5rec_lab[0], fV5rec_lab[2], fV5rec_lab[4]);
    TCS2HCS(fV5rec_tr[1], fV5rec_tr[3], fHRSAngle, fV5rec_lab[1], fV5rec_lab[3]);

    if (fDebug > 1) Info(here, "dbrec_tr  : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5rec_tr[0], fV5rec_tr[1], fV5rec_tr[2], fV5rec_tr[3], fV5rec_tr[4]);

    return 0;
}

void G2PDBRec::Clear(Option_t* option) {
    memset(fV5fp_tr, 0, sizeof (fV5fp_tr));
    memset(fV5fp_rot, 0, sizeof (fV5fp_rot));
    memset(fV5rec_tr, 0, sizeof (fV5rec_tr));
    memset(fV5rec_lab, 0, sizeof (fV5rec_lab));
    
    G2PProcBase::Clear(option);
}

void G2PDBRec::PrintDataBase() {
    static const char* const here = "TransDet2Rot()";

    map<int, vector<THaMatrixElement>*> matrix_map;
    matrix_map[1] = &ftMatrixElems;
    matrix_map[2] = &fyMatrixElems;
    matrix_map[3] = &fpMatrixElems;
    matrix_map[4] = &fDMatrixElems;
    matrix_map[5] = &fTMatrixElems;
    matrix_map[6] = &fPMatrixElems;
    matrix_map[7] = &fYMatrixElems;
    matrix_map[8] = &fPTAMatrixElems;
    matrix_map[9] = &fYTAMatrixElems;

    map<int, const char*> name_map;
    name_map[1] = "t";
    name_map[2] = "y";
    name_map[3] = "p";
    name_map[4] = "D";
    name_map[5] = "T";
    name_map[6] = "P";
    name_map[7] = "Y";
    name_map[8] = "PTA";
    name_map[9] = "YTA";

    Info(here, "Reconstruction database:");

    for (int i = 1; i <= 9; i++) {
        vector<THaMatrixElement>* mat = matrix_map[i];
        if (mat->size() > 0) {
            for (vector<THaMatrixElement>::iterator it = mat->begin(); it != mat->end(); it++) {
                cout << name_map[i] << " ";
                it->Print();
            }
        }
    }
}

double G2PDBRec::CalcVar(const double powers[][5], vector<THaMatrixElement> &matrix) {
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

void G2PDBRec::CalcMatrix(const double x, vector<THaMatrixElement> &matrix) {
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

int G2PDBRec::Configure(EMode mode) {
    if (mode == kREAD || mode == kTWOWAY) {
        if (fIsInit) return 0;
        else fIsInit = true;
    }

    ConfDef confs[] = {
        {"run.hrs.angle", "Beam Energy", kDOUBLE, &fHRSAngle},
        {0}
    };

    return ConfigureFromList(confs, mode);
}

int G2PDBRec::DefineVariables(EMode mode) {
    if (mode == kDEFINE && fIsSetup) return 0;
    fIsSetup = (mode == kDEFINE);

    VarDef vars[] = {
        {"rec.x", "Rec TP X", kDOUBLE, &fV5rec_tr[0]},
        {"rec.t", "Rec TP T", kDOUBLE, &fV5rec_tr[1]},
        {"rec.y", "Rec TP Y", kDOUBLE, &fV5rec_tr[2]},
        {"rec.p", "Rec TP P", kDOUBLE, &fV5rec_tr[3]},
        {"rec.d", "Rec TP D", kDOUBLE, &fV5rec_tr[4]},
        {"rec.l_x", "Rec TP X (lab)", kDOUBLE, &fV5rec_lab[0]},
        {"rec.l_t", "Rec TP T (lab)", kDOUBLE, &fV5rec_lab[1]},
        {"rec.l_y", "Rec TP Y (lab)", kDOUBLE, &fV5rec_lab[2]},
        {"rec.l_p", "Rec TP P (lab)", kDOUBLE, &fV5rec_lab[3]},
        {"rec.l_z", "Rec TP Z (lab)", kDOUBLE, &fV5rec_lab[4]},
        {0}
    };

    return DefineVarsFromList(vars, mode);
}

void G2PDBRec::MakePrefix() {
    const char* base = "db";

    G2PAppBase::MakePrefix(base);
}

G2PDBRec::THaMatrixElement::THaMatrixElement()
: fIsZero(true), fPower(3), fOrder(0), fPoly(G2PDBRec::kPORDER),
fValue(0) {
    // Nothing to do
}

G2PDBRec::THaMatrixElement::~THaMatrixElement() {
    // Nothing to do
}

bool G2PDBRec::THaMatrixElement::IsMatch(const THaMatrixElement& rhs) const {
    // Compare coefficients of this matrix element to another
    if (fPower.size() != rhs.fPower.size()) return false;
    for (vector<int>::size_type i = 0; i < fPower.size(); i++) {
        if (fPower[i] != rhs.fPower[i]) return false;
    }
    return true;
}

void G2PDBRec::THaMatrixElement::SkimPoly() {
    if (fIsZero) return;

    while ((!fPoly[fOrder - 1])&&(fOrder > 0)) {
        fPoly.pop_back();
        fOrder = fOrder - 1;
    }
    if (fOrder == 0) fIsZero = true;
}

void G2PDBRec::THaMatrixElement::Print() {
    if (fIsZero) {
        cout << "This element is zero" << endl;
        return;
    }

    for (vector<int>::size_type i = 0; i < fPower.size(); i++) cout << fPower[i] << " ";
    cout << "\t";
    cout.setf(ios::scientific, ios::floatfield);
    cout.precision(6);
    for (vector<double>::size_type i = 0; i < fPoly.size(); i++) cout << fPoly[i] << "\t";
    cout << endl;
    cout.unsetf(ios::floatfield);
}

ClassImp(G2PDBRec)
