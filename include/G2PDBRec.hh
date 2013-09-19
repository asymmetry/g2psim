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

#ifndef G2P_DBREC_H
#define G2P_DBREC_H

#include <vector>

#include "G2PProcBase.hh"

using namespace std;

class G2PDBRec : public G2PProcBase {
public:
    G2PDBRec();
    virtual ~G2PDBRec();

    virtual int Begin();
    virtual int Process();
    virtual void Clear();

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

        bool fIsZero; // whether the element is zero
        vector<int> fPower; // exponents of matrix element, e.g. D100 = {1, 0, 0}

        int fOrder;
        vector<double> fPoly; // the associated polynomial

        double fValue; // the final value once x is given
    };

    void CalcMatrix(const double x, vector<THaMatrixElement> &matrix);
    double CalcVar(const double powers[][5], vector<THaMatrixElement> &matrix);

    virtual int Configure(EMode mode = kTWOWAY);
    virtual int DefineVariables(EMode mode = kDEFINE);
    virtual void MakePrefix();

    const char* fDBPrefix;
    const char* fDBFile;

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

    double fHRSAngle;

    double fV5fp_tr[5];
    double fV5fp_rot[5];
    double fV5rec_tr[5];
    double fV5rec_lab[5];

private:
    static G2PDBRec* pG2PDBRec;

    ClassDef(G2PDBRec, 1)
};

#endif
