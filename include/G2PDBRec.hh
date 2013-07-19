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

#ifndef G2P_RECUSEDB_H
#define G2P_RECUSEDB_H

#include <vector>

#include "G2PAppBase.hh"

using namespace std;

class G2PDBRec : public G2PAppBase {
public:
    G2PDBRec();
    ~G2PDBRec();

    int Init();

    void CalcTargetCoords(const double* V5fp_rot, double* V5tg_tr);

    void TransTr2Rot(const double* V5fp_tr, double* V5fp_rot);
    void TransRot2Tr(const double* V5fp_rot, double* V5fp_tr);
    void TransTr2Det(const double* V5fp_tr, double* V5fp_det);
    void TransDet2Tr(const double* V5fp_det, double* V5fp_tr);
    void TransRot2Det(const double* V5fp_rot, double* V5fp_det);
    void TransDet2Rot(const double* V5fp_det, double* V5fp_rot);

    static G2PDBRec* GetInstance() {
        return pG2PDBRec;
    }

protected:

    enum {
        kPORDER = 7, kNUM_PRECOMP_POW = 10
    }; // constants

    void PrintDataBase();

    // Class for storing matrix element data

    class THaMatrixElement {
    public:
        THaMatrixElement();
        ~THaMatrixElement();

        void SkimPoly(); // reduce order to highest non-zero poly

        bool IsMatch(const THaMatrixElement& rhs) const;
        void Print();

        bool bIsZero; // whether the element is zero
        vector<int> iPower; // exponents of matrix element, e.g. D100 = {1, 0, 0}

        int iOrder;
        vector<double> fPoly; // the associated polynomial

        double fValue; // the final value once x is given
    };

    void CalcMatrix(const double x, vector<THaMatrixElement> &matrix);
    double CalcVar(const double powers[][5], vector<THaMatrixElement> &matrix);

    const char* pPrefix;
    const char* pDBName;

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

private:
    static G2PDBRec* pG2PDBRec;

    ClassDef(G2PDBRec, 1)
};

#endif
