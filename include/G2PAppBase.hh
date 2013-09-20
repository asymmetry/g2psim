// -*- C++ -*-

/* class G2PAppBase
 * Abstract base class for g2p simulation tools.
 * It provides fundamental functions like coordinates transport.
 * No instance allowed for this class.
 * Many functions are modified from THaAnalysisObject. Thanks to O. Hansen.
 */

// History:
//   Mar 2013, C. Gu, First public version.
//   Sep 2013, C. Gu, Add configure functions.
//

#ifndef G2P_APPBASE_H
#define G2P_APPBASE_H

#include <set>

#include "TObject.h"

#include "G2PVarDef.hh"

using namespace std;

class G2PRand;

class G2PAppBase : public TObject {
public:
    virtual ~G2PAppBase();

    enum EStatus {
        kOK = 0, kNOTINIT, kINITERROR
    };

    enum EMode {
        kREAD = 0, kWRITE, kTWOWAY, kDEFINE = 0, kDELETE
    };

    //static const double kLARGE;

    // General processes
    virtual int Init();
    virtual int Begin();
    virtual int End();
    virtual void Clear();

    // Gets
    int GetDebugLevel() const;
    const char* GetPrefix() const;
    bool IsInit() const;
    bool IsOK() const;
    EStatus Status() const;
    int GetPriority() const;

    // Sets
    void SetDebugLevel(int level);

    static void SetSeed(unsigned n);

protected:
    G2PAppBase(); // No instance allowed for this class

    // Geometry utility functions
    virtual void TCS2HCS(double x_tr, double y_tr, double z_tr, double angle, double &x_lab, double &y_lab, double &z_lab);
    virtual void TCS2HCS(double t_tr, double p_tr, double angle, double &t_lab, double &p_lab);
    virtual void HCS2TCS(double x_lab, double y_lab, double z_lab, double angle, double &x_tr, double &y_tr, double &z_tr);
    virtual void HCS2TCS(double t_lab, double p_lab, double angle, double &t_tr, double &p_tr);
    virtual void Project(double x, double y, double z, double zout, double t, double p, double &xout, double &yout);

    virtual void TRCS2FCS(const double* V5_tr, double angle, double* V5_fp);
    virtual void FCS2TRCS(const double* V5_fp, double angle, double* V5_tr);
    virtual void TRCS2DCS(const double* V5_tr, double angle, double* V5_det);
    virtual void DCS2TRCS(const double* V5_det, double angle, double* V5_tr);
    virtual void FCS2DCS(const double* V5_fp, double angle, double* V5_det);
    virtual void DCS2FCS(const double* V5_det, double angle, double* V5_fp);

    // Configure functions
    virtual int Configure(EMode mode = kTWOWAY) = 0;
    int ConfigureFromList(const ConfDef* list, EMode mode = kTWOWAY);
    virtual int WriteConfs();

    // Make Prefix
    virtual void MakePrefix() = 0;
    void MakePrefix(const char* basename);

    // General status variables
    char* fPrefix;
    EStatus fStatus;
    int fDebug;
    bool fIsInit;
    bool fIsSetup;

    // FIXME :
    // I understand that to use a priority variable is a bad idea.
    // If anyone has a better idea to organize the processing order of processes, I will appreciate that.
    int fPriority;

    set<unsigned long> fConfigIsSet;

    // Random number generator
    static G2PRand* pRand;

private:
    // Prevent default copy and assignment function
    G2PAppBase(const G2PAppBase&);
    G2PAppBase& operator=(const G2PAppBase&);

    ClassDef(G2PAppBase, 1)
};

#endif
