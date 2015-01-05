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
//   Nov 2014, C. Gu, Set random seed in G2PRun class.
//   Nov 2014, C. Gu, Add fHRSAngle and several interface function for coordinate transform.
//

#ifndef G2P_APPBASE_H
#define G2P_APPBASE_H

#include <set>

#include "TObject.h"

#include "G2PVarDef.hh"

using namespace std;

class G2PRand;

class G2PAppBase : public TObject
{
public:
    virtual ~G2PAppBase();

    enum EStatus {
        kOK = 0, kINITERROR, kBEGINERROR
    };

    enum EMode {
        kREAD = 0, kWRITE, kTWOWAY, kDEFINE = 0, kDELETE
    };

    // General processes
    virtual int Init();
    virtual int Begin();
    virtual int End();
    virtual void Clear(Option_t * /*option*/ = "");

    // Gets
    EStatus Status() const;
    bool IsInit() const;
    int GetDebugLevel() const;
    int GetPriority() const;

    // Sets
    void SetDebugLevel(int level);

protected:
    G2PAppBase(); // No instance allowed for this class

    // Geometry utility functions
    virtual void TCS2HCS(const double *V5_tr, double z_tr, double *V5_lab);
    virtual void TCS2HCS(double x_tr, double y_tr, double z_tr, double &x_lab, double &y_lab, double &z_lab);
    virtual void TCS2HCS(double t_tr, double p_tr, double &t_lab, double &p_lab);

    virtual void HCS2TCS(const double *V5_lab, double *V5_tr, double &z_tr);
    virtual void HCS2TCS(double x_lab, double y_lab, double z_lab, double &x_tr, double &y_tr, double &z_tr);
    virtual void HCS2TCS(double t_lab, double p_lab, double &t_tr, double &p_tr);

    virtual void TRCS2FCS(const double *V5_tr, double *V5_fp);
    virtual void FCS2TRCS(const double *V5_fp, double *V5_tr);
    virtual void TRCS2DCS(const double *V5_tr, double *V5_det);
    virtual void DCS2TRCS(const double *V5_det, double *V5_tr);
    virtual void FCS2DCS(const double *V5_fp, double *V5_det);
    virtual void DCS2FCS(const double *V5_det, double *V5_fp);

    // Configure functions
    virtual int Configure(EMode mode = kTWOWAY);
    int ConfigureFromList(const ConfDef *list, EMode mode = kTWOWAY);

    // Make Prefix
    virtual void MakePrefix() = 0;
    void MakePrefix(const char *basename);

    char *fPrefix;

    // General status variables
    EStatus fStatus;
    bool fConfigured; // Configure flag

    int fDebug; // Debug level

    // FIXME :
    // I understand that to use a priority variable is a bad idea.
    // If anyone has a better idea to organize the processing order of processes, I will appreciate that.
    int fPriority;

    set<unsigned long> fConfigIsSet;

    double fHRSAngle;

    // Random number generator
    static G2PRand *pRand;

private:
    ClassDef(G2PAppBase, 1)
};

#endif
