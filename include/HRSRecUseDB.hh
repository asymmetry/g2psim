#ifndef HRSRecUseDB_h
#define HRSRecUseDB_h 1

#include "TString.h"

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
    void CalcRotateCoords(const double *pV5fp_tr, double *pV5fp_rot);
    void CalcTransCoords(const double *pV5fp_det, double *pV5fp_tr);
    void CalcDetectorCoords(const double *pV5fp_tr, double *pV5fp_det);
    
    bool IsInit() {return fIsInit;}
    int SetPrefix(const char* prefix) {fPrefix = prefix; return 0;}
    int SetDBName(const char* dbname) {fDBName = dbname; return 0;}

    enum {kPORDER = 7, kNUM_PRECOMP_POW = 10}; // constants

private:
    const char* fPrefix;
    const char* fDBName;

    bool fIsInit;

    void CalcMatrix(const double x, std::vector<THaMatrixElement> &matrix);
    double CalcVar(const double powers[][5], const std::vector<THaMatrixElement> &matrix);

    friend class THaMatrixElement;
    std::vector<THaMatrixElement> *fCurrentMatrixElems;
    std::vector<THaMatrixElement> ftMatrixElems;
    std::vector<THaMatrixElement> fyMatrixElems;
    std::vector<THaMatrixElement> fpMatrixElems;
    std::vector<THaMatrixElement> fTMatrixElems;
    std::vector<THaMatrixElement> fDMatrixElems;
    std::vector<THaMatrixElement> fPMatrixElems;
    std::vector<THaMatrixElement> fPTAMatrixElems; // involves abs(theta_fp)
    std::vector<THaMatrixElement> fYMatrixElems;
    std::vector<THaMatrixElement> fYTAMatrixElems; // involves abs(theta_fp)
};

// class for storing matrix element data
class THaMatrixElement
{
public:
    THaMatrixElement() : iszero(true), pw(3), order(0), poly(HRSRecUseDB::kPORDER) {}

    void Print();
    
    bool iszero;              // whether the element is zero
    std::vector<int> pw;      // exponents of matrix element, e.g. D100 = { 1, 0, 0 }
    bool match( const THaMatrixElement& rhs ) const;

    int  order;
    std::vector<double> poly; // the associated polynomial
    void SkimPoly();          // reduce order to highest non-zero poly

    double v;                 // the final value of this matrix element once x is given
};

#endif
