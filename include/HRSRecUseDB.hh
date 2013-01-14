#ifndef HRS_RECUSEDB_H
#define HRS_RECUSEDB_H

#include <vector>

using namespace std;

class THaMatrixElement;

class HRSRecUseDB
{
public:
    HRSRecUseDB();
    HRSRecUseDB(const char* prefix, const char* dbname);
    ~HRSRecUseDB();

    int LoadDataBase();
    void PrintDataBase();

    void CalcTargetCoords(const double *pV5fp_rot, double *pV5tg_tr);
    void TransTr2Rot(const double *pV5fp_tr, double *pV5fp_rot);
    void TransRot2Tr(const double *pV5fp_rot, double *pV5fp_tr);
    void TransTr2Det(const double *pV5fp_tr, double *pV5fp_det);
    void TransDet2Tr(const double *pV5fp_det, double *pV5fp_tr);
    
    bool IsInit() { return fIsInit; }
    void SetPrefix(const char* prefix) { fPrefix = prefix; }
    void SetDBName(const char* dbname) { fDBName = dbname; }

    enum { kPORDER = 7, kNUM_PRECOMP_POW = 10 }; // constants

private:
    const char* fPrefix;
    const char* fDBName;

    bool fIsInit;

    void CalcMatrix(const double x, vector<THaMatrixElement> &matrix);
    double CalcVar(const double powers[][5], const vector<THaMatrixElement> &matrix);

    friend class THaMatrixElement;
    vector<THaMatrixElement> *fCurrentMatrixElems;
    vector<THaMatrixElement> ftMatrixElems;
    vector<THaMatrixElement> fyMatrixElems;
    vector<THaMatrixElement> fpMatrixElems;
    vector<THaMatrixElement> fTMatrixElems;
    vector<THaMatrixElement> fDMatrixElems;
    vector<THaMatrixElement> fPMatrixElems;
    vector<THaMatrixElement> fPTAMatrixElems; // involves abs(theta_fp)
    vector<THaMatrixElement> fYMatrixElems;
    vector<THaMatrixElement> fYTAMatrixElems; // involves abs(theta_fp)
};

// Class for storing matrix element data
class THaMatrixElement
{
public:
    THaMatrixElement() : iszero(true), pw(3), order(0), poly(HRSRecUseDB::kPORDER) {}

    void Print();
    
    bool iszero;         // whether the element is zero
    vector<int> pw;      // exponents of matrix element, e.g. D100 = { 1, 0, 0 }
    bool match( const THaMatrixElement& rhs ) const;

    int order;
    vector<double> poly; // the associated polynomial
    void SkimPoly();     // reduce order to highest non-zero poly

    double v;            // the final value once x is given
};

#endif
