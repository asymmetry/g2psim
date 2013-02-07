// This file defines a class HRSRecUseDB.
// This class is used in G2PSim class together with HRSTransport.
// It reads replay database and reconstruct target variables using this
//+database.
// It also provides transform functions among 3 different coordinates on focus
//+plane.
// Several functions is developed from J. Huang's HRS optics class.
//
// History:
//   Jun 2010, J. Huang, HRS optics matrix optimization class.
//   Jan 2013, C. Gu, First public version.
//

#ifndef HRS_RECUSEDB_H
#define HRS_RECUSEDB_H

#include <vector>

#include "TROOT.h"
#include "TObject.h"

using namespace std;

class HRSRecUseDB : public TObject
{
public:
    HRSRecUseDB();
    HRSRecUseDB(const char* prefix, const char* dbname);
    ~HRSRecUseDB();
    
    bool IsInit() { return bIsInit; }
    void SetPrefix(const char* prefix) { pPrefix = prefix; }
    void SetDBName(const char* dbname) { pDBName = dbname; }

    int LoadDataBase();
    void PrintDataBase();

    void CalcTargetCoords(const double* V5fp_rot, double* V5tg_tr);
    
    void TransTr2Rot(const double* V5fp_tr, double* V5fp_rot);
    void TransRot2Tr(const double* V5fp_rot, double* V5fp_tr);
    void TransTr2Det(const double* V5fp_tr, double* V5fp_det);
    void TransDet2Tr(const double* V5fp_det, double* V5fp_tr);
    void TransRot2Det(const double* V5fp_rot, double* V5fp_det);
    void TransDet2Rot(const double* V5fp_det, double* V5fp_rot);

    enum { kPORDER = 7, kNUM_PRECOMP_POW = 10 }; // constants

private:
    const char* pPrefix;
    const char* pDBName;

    bool bIsInit;

    // Class for storing matrix element data
    class THaMatrixElement
    {
    public:
        THaMatrixElement();
        ~THaMatrixElement();

        void SkimPoly();         // reduce order to highest non-zero poly

        bool IsMatch(const THaMatrixElement& rhs) const;
        void Print();

        bool bIsZero;            // whether the element is zero
        vector<int> iPower;      // exponents of matrix element, e.g. D100 = {
                                 // 1, 0, 0 }
        
        int iOrder;
        vector<double> fPoly;    // the associated polynomial

        double fValue;           // the final value once x is given
    };

    void CalcMatrix(const double x, vector<THaMatrixElement> &matrix);
    double CalcVar(const double powers[][5], vector<THaMatrixElement> &matrix);

    vector<THaMatrixElement>* fCurrentMatrixElems;
    vector<THaMatrixElement> ftMatrixElems;
    vector<THaMatrixElement> fyMatrixElems;
    vector<THaMatrixElement> fpMatrixElems;
    vector<THaMatrixElement> fTMatrixElems;
    vector<THaMatrixElement> fDMatrixElems;
    vector<THaMatrixElement> fPMatrixElems;
    vector<THaMatrixElement> fPTAMatrixElems; // involves abs(theta_fp)
    vector<THaMatrixElement> fYMatrixElems;
    vector<THaMatrixElement> fYTAMatrixElems; // involves abs(theta_fp)

    ClassDef(HRSRecUseDB,1);
};

#endif
