#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <cmath>

#include "TROOT.h"
#include "TString.h"

#include "HRSRecUseDB.hh"

using namespace std;

HRSRecUseDB::HRSRecUseDB()
    : fIsInit(false)
{
}

HRSRecUseDB::HRSRecUseDB(const char* prefix, const char* dbname)
    : fIsInit(false)
{
    SetPrefix(prefix);
    SetDBName(dbname);
    LoadDataBase();
}

HRSRecUseDB::~HRSRecUseDB()
{
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
        ME.pw.resize(npow);
        ME.iszero = true;
        ME.order = 0;

        vsiz_t pos;
        for (pos = 1; pos <= npow && pos < line_spl.size(); pos++)
            ME.pw[pos-1] = atoi(line_spl[pos].c_str());
        vsiz_t p_cnt;
        for (p_cnt = 0; pos < line_spl.size() && p_cnt < kPORDER && pos <= npow + kPORDER; pos++, p_cnt++) {
            ME.poly[p_cnt] = atof(line_spl[pos].c_str());
            if (ME.poly[p_cnt] != 0.0) {
                ME.iszero = false;
                ME.order = p_cnt + 1;
            }
        }
        if (p_cnt < 1) {
            cout << "***Error: Could not read in Matrix Element " << w << ME.pw[0] << ME.pw[1] << ME.pw[2] << " !" << endl;
            cout << "          Line looks like : " << tmpline.Data() << endl;
            ifs.close();
            return 1;
        }

        if (ME.iszero) continue;

        vector<THaMatrixElement> *mat = matrix_map[w];
        if (mat) {
            bool match = false;
            for (vector<THaMatrixElement>::iterator it = mat->begin(); it != mat->end() && !(match = it->match(ME)); it++);
            if (match)
                cout << "***Warning: Duplicate definition of matrix element " << w << ME.pw[0] << ME.pw[1] << ME.pw[2] << ". Using first definition." << endl;
            else
                mat->push_back(ME);
        }
    }

    fIsInit = true;
    ifs.close();
    return 0;
}

void HRSRecUseDB::PrintDataBase()
{
    map<int,vector<THaMatrixElement>*> matrix_map;
	matrix_map[1] = &fDMatrixElems;
	matrix_map[2] = &fTMatrixElems;
	matrix_map[3] = &fPMatrixElems;
	matrix_map[4] = &fYMatrixElems;
	matrix_map[5] = &fPTAMatrixElems;
	matrix_map[6] = &fYTAMatrixElems;

    map<int,const char*> name_map;
    name_map[1] = "D";
    name_map[2] = "T";
    name_map[3] = "P";
    name_map[4] = "Y";
    name_map[5] = "PTA";
    name_map[6] = "YTA";

    cout << "Reconstruction database: " << endl;
    
    for (int i = 1; i <= 6; i++) {
        vector<THaMatrixElement> *mat = matrix_map[i];
        if (mat->size() > 0) {
            for (vector<THaMatrixElement>::iterator it = mat->begin(); it != mat->end(); it++){
                cout << name_map[i] << " ";
                it->Print();
            }
        }
    }
}

void HRSRecUseDB::CalcTargetCoords(const double *fpcoords, double *tgcoords)
{
    // Calculate target coordinates from focal plane coordinates
    double x_fp, y_fp, th_fp, ph_fp;
    double powers[kNUM_PRECOMP_POW][5];
    double x, y, theta, phi, dp;

    x_fp = fpcoords[0];
    th_fp = fpcoords[1];
    y_fp = fpcoords[2];
    ph_fp = fpcoords[3];

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

    theta = CalcTargetVar(powers, fTMatrixElems);
    phi   = CalcTargetVar(powers, fPMatrixElems);
    y     = CalcTargetVar(powers, fYMatrixElems);
    dp    = CalcTargetVar(powers, fDMatrixElems);

    x = 0.0;
    
    tgcoords[0] = x;
    tgcoords[1] = theta;
    tgcoords[2] = y;
    tgcoords[3] = phi;
    tgcoords[4] = dp;
}

double HRSRecUseDB::CalcTargetVar(const double powers[][5], const vector<THaMatrixElement> &matrix)
{
    // Calculate the value of a variable at the target
    // Must already have x values for the matrix elements
    double value = 0.0;
    double v = 0.0;
    for (vector<THaMatrixElement>::const_iterator it = matrix.begin(); it != matrix.end(); it++) {
        if (it->v != 0.0) {
            v = it->v;
            vector<int>::size_type np = it->pw.size();
            for (vector<int>::size_type i = 0; i < np; i++)
                v *= powers[it->pw[i]][i+1];
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
        if (it->order > 0) {
            for (int i = it->order-1; i>=1; i--)
                value = x * (value + it->poly[i]);
            value+= it->poly[0];
        }
        it->v = value;
    }   
}

bool THaMatrixElement::match(const THaMatrixElement& rhs) const
{
    // Compare coefficients of this matrix element to another
    if (pw.size() != rhs.pw.size()) return false;
    for (vector<int>::size_type i = 0; i < pw.size(); i++) {
        if (pw[i] != rhs.pw[i] ) return false;
    }
    return true;
}

void THaMatrixElement::SkimPoly()
{
    if (iszero) return;

    while (!poly[order-1] && order > 0) {
        poly.pop_back();
        order = order - 1;
    }
    if (order == 0) iszero = true;
}

void THaMatrixElement::Print()
{
    if (iszero) {
        cout << "This element is zero" << endl;
        return;
    }

    for (vector<int>::size_type i = 0; i < pw.size(); i++) cout << pw[i] << " ";
    cout << "\t";
    //cout.setf(ios::left, ios::adjustfield);
    cout.setf(ios::scientific, ios::floatfield);
    //cout.width(10);
    cout.precision(6);
    for (vector<double>::size_type i = 0; i < poly.size(); i++) cout << poly[i] << "\t";
    cout << endl;
    //cout.unsetf(ios::adjustfield);
    cout.unsetf(ios::floatfield);
}
