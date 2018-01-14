// -*- C++ -*-

/* class G2PBwdDB
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

#ifndef G2P_BWDDB_H
#define G2P_BWDDB_H

#include <vector>

#include "G2PProcBase.hh"

using namespace std;

class G2PSieve;

class G2PBwdDB : public G2PProcBase
{
public:
    G2PBwdDB(const char *name);
    virtual ~G2PBwdDB();

    virtual int Begin();
    virtual int Process();
    virtual void Clear(Option_t *opt = "");

    // Gets

    // Sets
    void SetParsX(const double *pars);
    void SetParsY(const double *pars);
    void SetRecZ(double z);

    // Class for storing matrix element data
    class THaMatrixElement
    {
    public:
        THaMatrixElement();
        ~THaMatrixElement();

        void SkimPoly(); // reduce order to highest non-zero poly

        bool IsMatch(const THaMatrixElement &rhs) const;

        bool fIsZero; // whether the element is zero
        vector<int> fPower; // exponents of matrix element, e.g. D100 = {1, 0, 0}

        int fOrder;
        vector<double> fPoly; // the associated polynomial

        double fValue; // the final value once x is given
    };

protected:
    G2PBwdDB(); // Only for ROOT I/O

    enum {
        kPORDER = 7, kNUM_PRECOMP_POW = 10
    }; // constants

    void GetEffBPM(const double *V5tp_tr, const double *V5bpm_tr, double *V5bpmeff_tr);
    void Correct(const double *V5bpm_tr, const double *V5tp_tr, double *V5corr_tr);

    void CalcMatrix(const double x, vector<THaMatrixElement> &matrix);
    double CalcVar(const double powers[][5], vector<THaMatrixElement> &matrix);

    bool Backward(const double *V5fp_tr, double *V5tp_tr);

    virtual int Configure(EMode mode = kTWOWAY);
    virtual int DefineVariables(EMode mode = kDEFINE);
    virtual void MakePrefix();

    const char *fDBFile;

    double fE;
    double fFieldRatio;

    double fFitPars[2][3];
    double fCorT[3];
    double fCorP[3];
    double fCorD[3];
    double frecz_lab;

    double fV5bpm_tr[5];
    double fbpmz_tr;

    double fV5fp_tr[5];

    double fV5tpmat_tr[5];
    double fV5tpcorr_tr[5];
    double fV5sieveproj_tr[5];

    double fV5tprec_tr[5];
    double fV5tprec_lab[5];

    double fV5rec_tr[5];
    double fV5rec_lab[5];

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

    G2PSieve *pSieve;

private:
    static G2PBwdDB *pG2PBwdDB;

    ClassDef(G2PBwdDB, 1)
};

#endif
