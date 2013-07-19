// -*- C++ -*-

/* class G2PDBRec
 * This file defines a class G2PDBRec.
 * It use the analyzer database to reconstruct target variables.
 * It also provides transform functions among 3 different coordinates on focus plane.
 * It calculates the beam position at BPM and target using kinematics from event generator.
 * Several functions is developed from J. Huang's HRS optics class.
 * G2PProcBase classes will call CalcTargetCoords() to get target variables.
 */

// History:
//   Jun 2010, J. Huang, HRS optics matrix optimization class.
//   Jan 2013, C. Gu, First public version.
//

#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>

#include "TROOT.h"
#include "TObject.h"
#include "TError.h"
#include "TString.h"

#include "G2PAppBase.hh"
#include "G2PGlobals.hh"
#include "G2PRunBase.hh"

#include "G2PDBRec.hh"

using namespace std;

G2PDBRec* G2PDBRec::pG2PDBRec = NULL;

G2PDBRec::G2PDBRec() :
pPrefix(NULL), pDBName(NULL) {
    if (pG2PDBRec) {
        Error("G2PDBRec()", "Only one instance of G2PDBRec allowed.");
        MakeZombie();
        return;
    }
    pG2PDBRec = this;
}

G2PDBRec::~G2PDBRec() {
    if (pG2PDBRec == this) pG2PDBRec = NULL;
}

int G2PDBRec::Init() {
    static const char* const here = "Init()";

    if (G2PAppBase::Init() != 0) return fStatus;

    if (gG2PRun->GetHRSAngle() > 0) {
        pPrefix = "L";
        pDBName = "db_L.vdc.dat";
    }
    else {
        pPrefix = "R";
        pDBName = "db_R.vdc.dat";
    }

    // Read VDC database
    ifstream ifs(pDBName, ios_base::in);
    if (!ifs.good()) {
        Error(here, "Cannot initialize, database file \"%s\" does not exist.", pDBName);
        return (fStatus = kINITERROR);
    }

    TString tag(pPrefix);
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
    while ((!found)&&(ifs.getline(buff, LEN) != NULL)) {
        tmpline = ::Compress(buff); //strip blanks

        if (tmpline.EndsWith("\n")) tmpline.Chop();
        tmpline.ToLower();

        if (tag == tmpline) found = true;
    }

    if (!found) {
        Error(here, "Cannot initialize, no matrix in database file \"%s\".", pDBName);
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

    if (fDebug > 0) Info(here, "Loading matrix from %s ...", pDBName);

    // Read matrix elements line by line
    while (ifs.getline(buff, LEN) != NULL) {
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
        ME.iPower.resize(npow);
        ME.bIsZero = true;
        ME.iOrder = 0;

        vsiz_t pos;
        for (pos = 1; (pos <= npow)&&(pos < line_spl.size()); pos++)
            ME.iPower[pos - 1] = atoi(line_spl[pos].c_str());
        vsiz_t p_cnt;
        for (p_cnt = 0; (pos < line_spl.size())&&(p_cnt < kPORDER)&&(pos <= npow + kPORDER); pos++, p_cnt++) {
            ME.fPoly[p_cnt] = atof(line_spl[pos].c_str());
            if (ME.fPoly[p_cnt] != 0.0) {
                ME.bIsZero = false;
                ME.iOrder = p_cnt + 1;
            }
        }
        if (p_cnt < 1) {
            Error(here, "Cannot Initialize, matrix element %s%d%d%d has error.", w, ME.iPower[0], ME.iPower[1], ME.iPower[2]);
            ifs.close();
            return (fStatus = kINITERROR);
        }

        if (ME.bIsZero) continue;

        vector<THaMatrixElement>* mat = matrix_map[w];
        if (mat) {
            bool match = false;
            for (vector<THaMatrixElement>::iterator it = mat->begin(); (it != mat->end())&&(!(match = it->IsMatch(ME))); it++);
            if (match)
                Warning(here, "Duplicate definition of matrix element %s%d%d%d.", w, ME.iPower[0], ME.iPower[1], ME.iPower[2]);
            else
                mat->push_back(ME);
        }
    }

    fStatus = kOK;
    ifs.close();

    return (fStatus = kOK);
}

void G2PDBRec::CalcTargetCoords(const double* V5fp_rot, double* V5tg_tr) {
    static const char* const here = "CalcTargetCoords()";

    // Calculate target coordinates from focal plane coordinates
    double x_fp, y_fp, th_fp, ph_fp;
    double powers[kNUM_PRECOMP_POW][5];
    double x, y, theta, phi, dp;

    x_fp = V5fp_rot[0];
    th_fp = tan(V5fp_rot[1]);
    y_fp = V5fp_rot[2];
    ph_fp = tan(V5fp_rot[3]);

    for (int i = 0; i < kNUM_PRECOMP_POW; i++) {
        powers[i][0] = pow(x_fp, i);
        powers[i][1] = pow(th_fp, i);
        powers[i][2] = pow(y_fp, i);
        powers[i][3] = pow(ph_fp, i);
        powers[i][4] = pow(fabs(th_fp), i);
    }

    CalcMatrix(x_fp, fDMatrixElems);
    CalcMatrix(x_fp, fTMatrixElems);
    CalcMatrix(x_fp, fPMatrixElems);
    CalcMatrix(x_fp, fYMatrixElems);
    CalcMatrix(x_fp, fPTAMatrixElems);
    CalcMatrix(x_fp, fYTAMatrixElems);

    theta = CalcVar(powers, fTMatrixElems);
    phi = CalcVar(powers, fPMatrixElems);
    y = CalcVar(powers, fYMatrixElems);
    dp = CalcVar(powers, fDMatrixElems);

    x = 0.0;

    V5tg_tr[0] = x;
    V5tg_tr[1] = atan(theta);
    V5tg_tr[2] = y;
    V5tg_tr[3] = atan(phi);
    V5tg_tr[4] = dp;

    if (fDebug > 2) Info(here, "%10.3e %10.3e %10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e %10.3e %10.3e", V5fp_rot[0], V5fp_rot[1], V5fp_rot[2], V5fp_rot[3], V5fp_rot[4], V5tg_tr[0], V5tg_tr[1], V5tg_tr[2], V5tg_tr[3], V5tg_tr[4]);
}

void G2PDBRec::TransTr2Rot(const double* V5fp_tr, double* V5fp_rot) {
    static const char* const here = "TransTr2Rot()";

    double x, y, theta, phi;

    double V5fp_det[5];

    TransTr2Det(V5fp_tr, V5fp_det);

    double th_det = tan(V5fp_det[1]);
    double ph_det = tan(V5fp_det[3]);

    double x_tr = V5fp_tr[0];
    double y_tr = V5fp_tr[2];

    x = x_tr;

    CalcMatrix(x, ftMatrixElems);
    CalcMatrix(x, fyMatrixElems);
    CalcMatrix(x, fpMatrixElems);

    double tan_rho = ftMatrixElems[0].fValue;
    double cos_rho = 1.0 / sqrt(1.0 + tan_rho * tan_rho);

    y = y_tr - fyMatrixElems[0].fValue;
    theta = (th_det + tan_rho) / (1.0 - th_det * tan_rho);
    phi = (ph_det - fpMatrixElems[0].fValue) / ((1.0 - th_det * tan_rho) * cos_rho);

    V5fp_rot[0] = x;
    V5fp_rot[1] = atan(theta);
    V5fp_rot[2] = y;
    V5fp_rot[3] = atan(phi);

    if (fDebug > 3) Info(here, "%10.3e %10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e %10.3e", V5fp_tr[0], V5fp_tr[1], V5fp_tr[2], V5fp_tr[3], V5fp_rot[0], V5fp_rot[1], V5fp_rot[2], V5fp_rot[3]);
}

void G2PDBRec::TransRot2Tr(const double* V5fp_rot, double* V5fp_tr) {
    static const char* const here = "TransRot2Tr()";

    double x = V5fp_rot[0];
    double t = tan(V5fp_rot[1]);
    double y = V5fp_rot[2];
    double p = tan(V5fp_rot[3]);

    CalcMatrix(x, ftMatrixElems);
    CalcMatrix(x, fyMatrixElems);
    CalcMatrix(x, fpMatrixElems);

    double tan_rho = ftMatrixElems[0].fValue;
    double cos_rho = 1.0 / sqrt(1.0 + tan_rho * tan_rho);

    double x_tr = x;
    double y_tr = y + fyMatrixElems[0].fValue;
    double t_det = (t - tan_rho) / (1.0 + t * tan_rho);
    double p_det = p * (1.0 - t_det * tan_rho) * cos_rho + fpMatrixElems[0].fValue;

    double tan_rho_0 = ftMatrixElems[0].fPoly[0];
    double cos_rho_0 = 1.0 / sqrt(1.0 + tan_rho_0 * tan_rho_0);

    double t_tr = (t_det + tan_rho_0) / (1.0 - t_det * tan_rho_0);
    double p_tr = p_det / (cos_rho_0 * (1.0 - t_det * tan_rho_0));

    V5fp_tr[0] = x_tr;
    V5fp_tr[1] = atan(t_tr);
    V5fp_tr[2] = y_tr;
    V5fp_tr[3] = atan(p_tr);

    if (fDebug > 3) Info(here, "%10.3e %10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e %10.3e", V5fp_rot[0], V5fp_rot[1], V5fp_rot[2], V5fp_rot[3], V5fp_tr[0], V5fp_tr[1], V5fp_tr[2], V5fp_tr[3]);
}

void G2PDBRec::TransTr2Det(const double* V5fp_tr, double* V5fp_det) {
    static const char* const here = "TransTr2Det()";

    double tan_rho_0 = ftMatrixElems[0].fPoly[0];
    double cos_rho_0 = 1.0 / sqrt(1.0 + tan_rho_0 * tan_rho_0);

    double x_tr = V5fp_tr[0];
    double th_tr = tan(V5fp_tr[1]);
    double y_tr = V5fp_tr[2];
    double ph_tr = tan(V5fp_tr[3]);

    double x_det = x_tr / (cos_rho_0 * (1 + th_tr * tan_rho_0));
    double y_det = y_tr - tan_rho_0 * cos_rho_0 * ph_tr*x_det;
    double th_det = (th_tr - tan_rho_0) / (1 + th_tr * tan_rho_0);
    double ph_det = ph_tr * (1.0 - th_det * tan_rho_0) * cos_rho_0;

    V5fp_det[0] = x_det;
    V5fp_det[1] = atan(th_det);
    V5fp_det[2] = y_det;
    V5fp_det[3] = atan(ph_det);

    if (fDebug > 3) Info(here, "%10.3e %10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e %10.3e", V5fp_tr[0], V5fp_tr[1], V5fp_tr[2], V5fp_tr[3], V5fp_det[0], V5fp_det[1], V5fp_det[2], V5fp_det[3]);
}

void G2PDBRec::TransDet2Tr(const double* V5fp_det, double* V5fp_tr) {
    static const char* const here = "TransDet2Tr()";

    double tan_rho_0 = ftMatrixElems[0].fPoly[0];
    double cos_rho_0 = 1.0 / sqrt(1.0 + tan_rho_0 * tan_rho_0);

    double x_det = V5fp_det[0];
    double th_det = tan(V5fp_det[1]);
    double y_det = V5fp_det[2];
    double ph_det = tan(V5fp_det[3]);

    double th_tr = (th_det + tan_rho_0) / (1.0 - th_det * tan_rho_0);
    double ph_tr = ph_det / (cos_rho_0 * (1.0 - th_det * tan_rho_0));
    double x_tr = x_det * cos_rho_0 * (1 + th_tr * tan_rho_0);
    double y_tr = y_det + tan_rho_0 * cos_rho_0 * ph_tr*x_det;

    V5fp_tr[0] = x_tr;
    V5fp_tr[1] = atan(th_tr);
    V5fp_tr[2] = y_tr;
    V5fp_tr[3] = atan(ph_tr);

    if (fDebug > 3) Info(here, "%10.3e %10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e %10.3e", V5fp_det[0], V5fp_det[1], V5fp_det[2], V5fp_det[3], V5fp_tr[0], V5fp_tr[1], V5fp_tr[2], V5fp_tr[3]);
}

void G2PDBRec::TransRot2Det(const double* V5fp_rot, double* V5fp_det) {
    static const char* const here = "TransRot2Det()";

    double V5fp_tr[5];

    TransRot2Tr(V5fp_rot, V5fp_tr);
    TransTr2Det(V5fp_tr, V5fp_det);

    if (fDebug > 3) Info(here, "%10.3e %10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e %10.3e", V5fp_rot[0], V5fp_rot[1], V5fp_rot[2], V5fp_rot[3], V5fp_det[0], V5fp_det[1], V5fp_det[2], V5fp_det[3]);
}

void G2PDBRec::TransDet2Rot(const double* V5fp_det, double* V5fp_rot) {
    static const char* const here = "TransDet2Rot()";

    double V5fp_tr[5];

    TransDet2Tr(V5fp_det, V5fp_tr);
    TransTr2Rot(V5fp_tr, V5fp_rot);

    if (fDebug > 3) Info(here, "%10.3e %10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e %10.3e", V5fp_det[0], V5fp_det[1], V5fp_det[2], V5fp_det[3], V5fp_rot[0], V5fp_rot[1], V5fp_rot[2], V5fp_rot[3]);
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
            vector<int>::size_type np = it->iPower.size();
            for (vector<int>::size_type i = 0; i < np; i++)
                v *= powers[it->iPower[i]][i + 1];
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
        if (it->iOrder > 0) {
            for (int i = it->iOrder - 1; i >= 1; i--)
                value = x * (value + it->fPoly[i]);
            value += it->fPoly[0];
        }
        it->fValue = value;
    }
}

G2PDBRec::THaMatrixElement::THaMatrixElement()
: bIsZero(true), iPower(3), iOrder(0), fPoly(G2PDBRec::kPORDER),
fValue(0) {
    // Nothing to do
}

G2PDBRec::THaMatrixElement::~THaMatrixElement() {
    // Nothing to do
}

bool G2PDBRec::THaMatrixElement::IsMatch(const THaMatrixElement& rhs) const {
    // Compare coefficients of this matrix element to another
    if (iPower.size() != rhs.iPower.size()) return false;
    for (vector<int>::size_type i = 0; i < iPower.size(); i++) {
        if (iPower[i] != rhs.iPower[i]) return false;
    }
    return true;
}

void G2PDBRec::THaMatrixElement::SkimPoly() {
    if (bIsZero) return;

    while ((!fPoly[iOrder - 1])&&(iOrder > 0)) {
        fPoly.pop_back();
        iOrder = iOrder - 1;
    }
    if (iOrder == 0) bIsZero = true;
}

void G2PDBRec::THaMatrixElement::Print() {
    if (bIsZero) {
        cout << "This element is zero" << endl;
        return;
    }

    for (vector<int>::size_type i = 0; i < iPower.size(); i++) cout << iPower[i] << " ";
    cout << "\t";
    cout.setf(ios::scientific, ios::floatfield);
    cout.precision(6);
    for (vector<double>::size_type i = 0; i < fPoly.size(); i++) cout << fPoly[i] << "\t";
    cout << endl;
    cout.unsetf(ios::floatfield);
}

ClassImp(G2PDBRec)
