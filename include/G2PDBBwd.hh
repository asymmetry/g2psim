// -*- C++ -*-

/* class G2PDBBwd
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

#ifndef G2P_DBBWD_H
#define G2P_DBBWD_H

#include <vector>

#include "G2PProcBase.hh"

using namespace std;

class G2PDrift;
class G2PGeoSieve;

class G2PDBBwd : public G2PProcBase
{
public:
    G2PDBBwd(const char *name);
    virtual ~G2PDBBwd();

    virtual int Init();
    virtual int Begin();
    virtual int Process();
    virtual void Clear(Option_t * /*option*/ = "");

    // Gets

    // Sets
    void SetParsX(const double *pars);
    void SetParsY(const double *pars);
    void SetRecZ(double z);

protected:
    G2PDBBwd(); // Only for ROOT I/O

    double GetEffBPM(int axis);

    bool Backward(const double *V5fp_tr, double *V5tp_tr);

    virtual int Configure(EMode mode = kTWOWAY);
    virtual int DefineVariables(EMode mode = kDEFINE);
    virtual void MakePrefix();

    double fHRSMomentum;
    double fFieldRatio;

    double fFitPars[2][3];
    double frecz_lab;

    double fV5bpm_lab[5];
    double fV5bpm_tr[5];

    double fV5fp_tr[5];

    double fV5tpmat_tr[5];
    double fV5sieveproj_tr[5];

    double fV5tprec_tr[5];
    double fV5tprec_lab[5];

    G2PDrift *pDrift;
    G2PGeoSieve *pSieve;

private:

    enum {
        kPORDER = 7, kNUM_PRECOMP_POW = 10
    }; // constants

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

    void CalcMatrix(const double x, vector<THaMatrixElement> &matrix);
    double CalcVar(const double powers[][5], vector<THaMatrixElement> &matrix);

    const char *fDBPrefix;
    const char *fDBFile;

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

    static G2PDBBwd *pG2PDBBwd;

    ClassDef(G2PDBBwd, 1)
};

#endif
