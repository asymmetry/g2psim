#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <cmath>

#include "TROOT.h"
#include "TObject.h"
#include "TString.h"
#include "TMath.h"

#include "HRSRecUseDB.hh"

//#define RECUSEDB_DEBUG 1

using namespace std;

ClassImp(HRSRecUseDB);

HRSRecUseDB::HRSRecUseDB()
    :fIsInit(false)
{
    // Nothing to do
}

HRSRecUseDB::HRSRecUseDB(const char* prefix, const char* dbname)
    :fIsInit(false)
{
    SetPrefix(prefix);
    SetDBName(dbname);
    LoadDataBase();
}

HRSRecUseDB::~HRSRecUseDB()
{
    // Nothing to do
}

int HRSRecUseDB::LoadDataBase()
{
    // Read VDC database
    ifstream ifs(fDBName, ios_base::in);
    if (!ifs.good()) {
        cout << "***Error: Can not open database file \"" << fDBName << "\" !" << endl;
        return 1;
    }
    else
        cout << "Reading database file \"" << fDBName << "\"" << endl;

    TString tag(fPrefix);
    Ssiz_t pos = tag.Index(".");
    if (pos != kNPOS)
        tag = tag(0, pos+1);
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
    while (!found && (ifs.getline(buff, LEN) != NULL)) {
		tmpline = ::Compress(buff);  //strip blanks
        
		if (tmpline.EndsWith("\n")) tmpline.Chop();
        tmpline.ToLower();

        if (tag == tmpline) found = true;
    }

    if (!found){
        cout << "***Error: Can not find matrix in database file \"" << fDBName << "\" !" << endl;
        return 1;
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
	map<string,vsiz_t> power;
    power["t"] = 3;  // transport to focal-plane tensors
	power["y"] = 3;
	power["p"] = 3;
	power["D"] = 3;  // focal-plane to target tensors
	power["T"] = 3;
	power["Y"] = 3;
	power["YTA"] = 4;
	power["P"] = 3;
	power["PTA"] = 4;
	power["L"] = 4;  // pathlength from z=0 (target) to focal plane (meters)
	power["XF"] = 5; // forward: target to focal-plane (I think)
	power["TF"] = 5;
	power["PF"] = 5;
	power["YF"] = 5;
    
	map<string,vector<THaMatrixElement>*> matrix_map;
    matrix_map["t"] = &ftMatrixElems;
	matrix_map["y"] = &fyMatrixElems;
	matrix_map["p"] = &fpMatrixElems;
	matrix_map["D"] = &fDMatrixElems;
    matrix_map["T"] = &fTMatrixElems;
	matrix_map["Y"] = &fYMatrixElems;
	matrix_map["YTA"] = &fYTAMatrixElems;
	matrix_map["P"] = &fPMatrixElems;
	matrix_map["PTA"] = &fPTAMatrixElems;

    cout << "Loading matrix ..." << endl;

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
        for (pos = 1; pos <= npow && pos < line_spl.size(); pos++)
            ME.iPower[pos-1] = atoi(line_spl[pos].c_str());
        vsiz_t p_cnt;
        for (p_cnt = 0; pos < line_spl.size() && p_cnt < kPORDER && pos <= npow + kPORDER; pos++, p_cnt++) {
            ME.fPoly[p_cnt] = atof(line_spl[pos].c_str());
            if (ME.fPoly[p_cnt] != 0.0) {
                ME.bIsZero = false;
                ME.iOrder = p_cnt + 1;
            }
        }
        if (p_cnt < 1) {
            cout << "***Error: Could not read in Matrix Element " << w << ME.iPower[0] << ME.iPower[1] << ME.iPower[2] << " !" << endl;
            cout << "          Line looks like : " << tmpline.Data() << endl;
            ifs.close();
            return 1;
        }

        if (ME.bIsZero) continue;

        vector<THaMatrixElement> *mat = matrix_map[w];
        if (mat) {
            bool match = false;
            for (vector<THaMatrixElement>::iterator it = mat->begin(); it != mat->end() && !(match = it->IsMatch(ME)); it++);
            if (match)
                cout << "***Warning: Duplicate definition of matrix element " << w << ME.iPower[0] << ME.iPower[1] << ME.iPower[2] << ". Using first definition." << endl;
            else
                mat->push_back(ME);
        }
    }

    fIsInit = true;
    ifs.close();

#ifdef RECUSEDB_DEBUG
    PrintDataBase();
#endif
    
    return 0;
}

void HRSRecUseDB::PrintDataBase()
{
    map<int,vector<THaMatrixElement>*> matrix_map;
    matrix_map[1] = &ftMatrixElems;
	matrix_map[2] = &fyMatrixElems;
	matrix_map[3] = &fpMatrixElems;
	matrix_map[4] = &fDMatrixElems;
	matrix_map[5] = &fTMatrixElems;
	matrix_map[6] = &fPMatrixElems;
	matrix_map[7] = &fYMatrixElems;
	matrix_map[8] = &fPTAMatrixElems;
	matrix_map[9] = &fYTAMatrixElems;

    map<int,const char*> name_map;
    name_map[1] = "t";
    name_map[2] = "y";
    name_map[3] = "p";
    name_map[4] = "D";
    name_map[5] = "T";
    name_map[6] = "P";
    name_map[7] = "Y";
    name_map[8] = "PTA";
    name_map[9] = "YTA";

    cout << "Reconstruction database: " << endl;
    
    for (int i = 1; i <= 9; i++) {
        vector<THaMatrixElement> *mat = matrix_map[i];
        if (mat->size() > 0) {
            for (vector<THaMatrixElement>::iterator it = mat->begin(); it != mat->end(); it++){
                cout << name_map[i] << " ";
                it->Print();
            }
        }
    }
}

void HRSRecUseDB::CalcTargetCoords(const double *V5fp_rot, double *V5tg_tr)
{
    // Calculate target coordinates from focal plane coordinates
    double x_fp, y_fp, th_fp, ph_fp;
    double powers[kNUM_PRECOMP_POW][5];
    double x, y, theta, phi, dp;

    x_fp = V5fp_rot[0];
    th_fp = V5fp_rot[1];
    y_fp = V5fp_rot[2];
    ph_fp = V5fp_rot[3];

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
    phi   = CalcVar(powers, fPMatrixElems);
    y     = CalcVar(powers, fYMatrixElems);
    dp    = CalcVar(powers, fDMatrixElems);

    x = 0.0;
    
    V5tg_tr[0] = x;
    V5tg_tr[1] = theta;
    V5tg_tr[2] = y;
    V5tg_tr[3] = phi;
    V5tg_tr[4] = dp;
    
#ifdef RECUSEDB_DEBUG
    double V5fp_tr[5];
    TransRot2Tr(V5fp_rot, V5fp_tr);
    printf("%e\t%e\t%e\t%e\n", V5fp_tr[0], V5fp_tr[1], V5fp_tr[2], V5fp_tr[3]);
    double V5fp_det[5];
    TransTr2Det(V5fp_tr, V5fp_det);
    printf("%e\t%e\t%e\t%e\n", V5fp_det[0], V5fp_det[1], V5fp_det[2], V5fp_det[3]);
    printf("%e\t%e\t%e\t%e\n\n", V5fp_rot[0], V5fp_rot[1], V5fp_rot[2], V5fp_rot[3]);
#endif
}

void HRSRecUseDB::TransTr2Rot(const double *V5fp_tr, double *V5fp_rot)
{
    double x, y, theta, phi;
    
    double V5fp_det[5];

    TransTr2Det(V5fp_tr, V5fp_det);
    
    double th_det = V5fp_det[1];
    double ph_det = V5fp_det[3];

    double x_tr = V5fp_tr[0];
    double y_tr = V5fp_tr[2];

    x = x_tr;

    CalcMatrix(x, ftMatrixElems);
    CalcMatrix(x, fyMatrixElems);
    CalcMatrix(x, fpMatrixElems);

    double tan_rho = ftMatrixElems[0].fValue;
    double cos_rho = 1.0/sqrt(1.0+tan_rho*tan_rho);

    y = y_tr-fyMatrixElems[0].fValue;
    theta = (th_det+tan_rho)/(1.0-th_det*tan_rho);
    phi = (ph_det-fpMatrixElems[0].fValue)/((1.0-th_det*tan_rho)*cos_rho);

    V5fp_rot[0] = x;
    V5fp_rot[1] = theta;
    V5fp_rot[2] = y;
    V5fp_rot[3] = phi;
}

void HRSRecUseDB::TransRot2Tr(const double *V5fp_rot, double *V5fp_tr)
{
    double x = V5fp_rot[0];
    double t = V5fp_rot[1];
    double y = V5fp_rot[2];
    double p = V5fp_rot[3];

    CalcMatrix(x, ftMatrixElems);
    CalcMatrix(x, fyMatrixElems);
    CalcMatrix(x, fpMatrixElems);

    double tan_rho = ftMatrixElems[0].fValue;
    double cos_rho = 1.0/sqrt(1.0+tan_rho*tan_rho);
    
    double x_tr = x;
    double y_tr = y+fyMatrixElems[0].fValue;
    double t_det = (t-tan_rho)/(1.0+t*tan_rho);
    double p_det = p*(1.0-t_det*tan_rho)*cos_rho+fpMatrixElems[0].fValue;

    double tan_rho_0 = ftMatrixElems[0].fPoly[0];
    double cos_rho_0 = 1.0/sqrt(1.0+tan_rho_0*tan_rho_0);

    double t_tr = (t_det+tan_rho_0)/(1.0-t_det*tan_rho_0);
    double p_tr = p_det/(cos_rho_0*(1.0-t_det*tan_rho_0));
    
    // double t_tr =
    //     -.004791569093294506*x*x*x +.02738386993285832*t*x*x
    //     -.005985192096309318*x*x   +.1654807057996422*t*t*x
    //     +.1654807236128732*x       +.9999999874969879*t;
    // double p_tr =
    //     +6.0069583028972e-4*x*x*x  -.005906841761909994*t*x*x
    //     +.01369193472004655*p*x*x  +.006232486773147315*x*x
    //     +3.35069502985833e-4*t*t*x +.1654807212124184*p*t*x
    //     +.001746885379624117*t*x   -0.00140983544384478*x
    //     +.002024824822383892*t     +p
    //     -.002022529283042846;

    V5fp_tr[0] = x_tr;
    V5fp_tr[1] = t_tr;
    V5fp_tr[2] = y_tr;
    V5fp_tr[3] = p_tr;
}

void HRSRecUseDB::TransTr2Det(const double *V5fp_tr, double *V5fp_det)
{
    double tan_rho_0 = ftMatrixElems[0].fPoly[0];
    double cos_rho_0 = 1.0/sqrt(1.0+tan_rho_0*tan_rho_0);

    double x_tr = V5fp_tr[0];
    double th_tr = V5fp_tr[1];
    double y_tr = V5fp_tr[2];
    double ph_tr = V5fp_tr[3];

    double x_det = x_tr/(cos_rho_0*(1+th_tr*tan_rho_0));
    double y_det = y_tr-tan_rho_0*cos_rho_0*ph_tr*x_det;
    double th_det = (th_tr-tan_rho_0)/(1+th_tr*tan_rho_0);
    double ph_det = ph_tr*(1.0-th_det*tan_rho_0)*cos_rho_0;

    V5fp_det[0] = x_det;
    V5fp_det[1] = th_det;
    V5fp_det[2] = y_det;
    V5fp_det[3] = ph_det;
}

void HRSRecUseDB::TransDet2Tr(const double *V5fp_det, double *V5fp_tr)
{
    double tan_rho_0 = ftMatrixElems[0].fPoly[0];
    double cos_rho_0 = 1.0/sqrt(1.0+tan_rho_0*tan_rho_0);

    double x_det = V5fp_det[0];
    double th_det = V5fp_det[1];
    double y_det = V5fp_det[2];
    double ph_det = V5fp_det[3];

    double th_tr = (th_det+tan_rho_0)/(1.0-th_det*tan_rho_0);
    double ph_tr = ph_det/(cos_rho_0*(1.0-th_det*tan_rho_0));
    double x_tr = x_det*cos_rho_0*(1+th_tr*tan_rho_0);
    double y_tr = y_det+tan_rho_0*cos_rho_0*ph_tr*x_det;

    V5fp_tr[0] = x_tr;
    V5fp_tr[1] = th_tr;
    V5fp_tr[2] = y_tr;
    V5fp_tr[3] = ph_tr;
}

double HRSRecUseDB::CalcVar(const double powers[][5], const vector<THaMatrixElement> &matrix)
{
    // Calculate the value of a variable at the target
    // Must already have x values for the matrix elements
    double value = 0.0;
    double v = 0.0;
    for (vector<THaMatrixElement>::const_iterator it = matrix.begin(); it != matrix.end(); it++) {
        if (it->fValue != 0.0) {
            v = it->fValue;
            vector<int>::size_type np = it->iPower.size();
            for (vector<int>::size_type i = 0; i < np; i++)
                v *= powers[it->iPower[i]][i+1];
            value += v;
        }
    }
    return value;
}


void HRSRecUseDB::CalcMatrix(const double x, vector<THaMatrixElement> &matrix)
{
    // Calculate the value of a matrix element for a given x
    double value = 0.0;
    
    for (vector<THaMatrixElement>::iterator it = matrix.begin(); it != matrix.end(); it++) {
        value = 0.0;
        if (it->iOrder > 0) {
            for (int i = it->iOrder-1; i>=1; i--)
                value = x * (value + it->fPoly[i]);
            value+= it->fPoly[0];
        }
        it->fValue = value;
    }   
}

HRSRecUseDB::THaMatrixElement::THaMatrixElement()
    :bIsZero(true), iPower(3), iOrder(0), fPoly(HRSRecUseDB::kPORDER),
     fValue(0)
{
    // Nothing to do
}

HRSRecUseDB::THaMatrixElement::~THaMatrixElement()
{
    // Nothing to do
}

bool HRSRecUseDB::THaMatrixElement::IsMatch(const THaMatrixElement& rhs) const
{
    // Compare coefficients of this matrix element to another
    if (iPower.size() != rhs.iPower.size()) return false;
    for (vector<int>::size_type i = 0; i < iPower.size(); i++) {
        if (iPower[i] != rhs.iPower[i] ) return false;
    }
    return true;
}

void HRSRecUseDB::THaMatrixElement::SkimPoly()
{
    if (bIsZero) return;

    while (!fPoly[iOrder-1] && iOrder > 0) {
        fPoly.pop_back();
        iOrder = iOrder - 1;
    }
    if (iOrder == 0) bIsZero = true;
}

void HRSRecUseDB::THaMatrixElement::Print()
{
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
