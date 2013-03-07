#ifndef G2P_APPSBASE_H
#define G2P_APPSBASE_H

#include "TObject.h"

class G2PRand;
class TList;
class TTree;

class G2PAppsBase : public TObject
{
public:
    virtual ~G2PAppsBase();

    enum EStatus { kOK, kNOTINIT, kINITERROR };

    void SetDebug(int level) { fDebug = level; }

    virtual EStatus Init();
    virtual int Begin() { return 0; }
    virtual int End() { return 0; }
    virtual void Clear() { }

    bool IsInit() const { return IsOK(); }
    bool IsOK() const { return (fStatus==kOK); }
    EStatus Status() const { return fStatus; }
    int GetDebug() const { return fDebug; }

    virtual int RegisterModel() { return 0; }
    virtual int DefineVariables(TTree *t) { return 0; }

    static void SetSeed(int n);

protected:
    G2PAppsBase(); // No instance allowed for this class

    // Geometry utility functions
    virtual void TCS2HCS(double x_tr, double y_tr, double z_tr, double angle, double &x_lab, double &y_lab, double &z_lab);
    virtual void TCS2HCS(double t_tr, double p_tr, double angle, double &t_lab, double &p_lab);
    virtual void HCS2TCS(double x_lab, double y_lab, double z_lab, double angle, double &x_tr, double &y_tr, double &z_tr);
    virtual void HCS2TCS(double t_lab, double p_lab, double angle, double &t_tr, double &p_tr);

    G2PAppsBase* FindModule(const char* classname);
    void FindModule(const char* classname, TList* list);

    EStatus fStatus;
    int fDebug;

    TList* fApps;

    static G2PRand* pRand;

private:
    static TList* pAppsBase;

    ClassDef(G2PAppsBase, 1)
};

#endif
