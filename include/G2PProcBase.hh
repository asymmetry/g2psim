// -*- C++ -*-

/* class G2PProcBase
 * Abstract base class for g2p simulation processes.
 * It provides fundamental functions like variable registration.
 * No instance allowed for this class.
 * Derived class must set its own internal variables and register them to the global variable list.
 */

// History:
//   Apr 2013, C. Gu, First public version.
//   Sep 2013, C. Gu, Move DefineVariables() function from G2PAppBase to here.
//   Dec 2014, C. Gu, Add drifting functions and projection functions.
//

#ifndef G2P_PROCBASE_H
#define G2P_PROCBASE_H

#include "G2PAppBase.hh"

class G2PDrift;
class G2PGeoBase;

class G2PProcBase : public G2PAppBase
{
public:
    virtual ~G2PProcBase();

    enum EStage {
        kREADY = 0, kDONE, kSTOP
    };

    virtual int Init();
    virtual int Begin();
    virtual int Process() = 0;

    // Gets
    EStage GetStage();

    // Sets
    void SetStage(EStage stage);

protected:
    G2PProcBase(); // No instance allowed for this class

    int ArrayCopy(double *out, const double *in, int length);

    virtual double Drift(const char *dir, const double *x, const double *p, double zf, double *xout, double *pout); // HCS
    virtual double Drift(const char *dir, const double *V5_tr, double z_tr, double zf_tr, double *V5out_tr); // TCS
    virtual double Drift(const char *dir, const double *V5_tr, double z_tr, double zf_lab, double *V5out_tr, double &zout_tr); // TCS
    virtual double Drift(const char *dir, const double *x, const double *p, G2PGeoBase *geo, double *xout, double *pout); // HCS with geometry
    virtual double Drift(const char *dir, const double *V5_tr, double z_tr, G2PGeoBase *geo, double *V5out_tr, double &zout_tr); // TCS with Geometry

    virtual double Project(const double *x, const double *p, double zf, double *xout, double *pout); // HCS
    virtual double Project(const double *V5_tr, double z_tr, double zf_tr, double *V5out_tr); // TCS
    virtual double Project(double x, double y, double z, double zf, double t, double p, double &xout, double &yout);

    virtual int Configure(EMode mode = kTWOWAY);

    // Global variable functions
    virtual int DefineVariables(EMode mode = kDEFINE);
    int DefineVarsFromList(const VarDef *list, EMode mode = kDEFINE) const;
    virtual int RemoveVariables();

    EStage fStage;

    bool fDefined;

    double fHRSMomentum;

    G2PDrift *pDrift;

private:
    ClassDef(G2PProcBase, 1)
};

#endif
