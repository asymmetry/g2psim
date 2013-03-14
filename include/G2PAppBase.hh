#ifndef G2P_APPBASE_H
#define G2P_APPBASE_H

#include "TObject.h"

#include "G2PRand.hh"
#include "G2PVarDef.hh"

class TList;
class TTree;

class G2PAppBase : public TObject
{
public:
    virtual ~G2PAppBase();

    enum EStatus { kOK = 0, kNOTINIT, kINITERROR, kERROR };
    enum EMode { kDefine = 0, kDelete };

    const double kLARGE;

    void SetDebug(int level) { fDebug = level; }

    virtual int Init();
    virtual int Begin();
    virtual int End();
    virtual void Clear() { }

    bool IsInit() const { return IsOK(); }
    bool IsOK() const { return (fStatus==kOK); }
    EStatus Status() const { return fStatus; }

    static void SetSeed(int n) { pRand->SetSeed(n); }

protected:
    G2PAppBase(); // No instance allowed for this class

    // Geometry utility functions
    virtual void TCS2HCS(double x_tr, double y_tr, double z_tr, double angle, double &x_lab, double &y_lab, double &z_lab);
    virtual void TCS2HCS(double t_tr, double p_tr, double angle, double &t_lab, double &p_lab);
    virtual void HCS2TCS(double x_lab, double y_lab, double z_lab, double angle, double &x_tr, double &y_tr, double &z_tr);
    virtual void HCS2TCS(double t_lab, double p_lab, double angle, double &t_tr, double &p_tr);
    virtual void Project(double x, double y, double z, double z_out, double t, double p, double &xout, double &yout);

    virtual int DefineVariables(EMode mode = kDefine) { return 0; }
    int DefineVarsFromList(const VarDef* list, EMode mode = kDefine) const;

    virtual void MakePrefix() { }
    void MakePrefix(const char* basename);

    G2PAppBase* FindModule(const char* classname) const;
    void FindModule(const char* classname, TList* list) const;

    bool bIsSetup;

    EStatus fStatus;
    char* fPrefix;

    int fDebug;

    static G2PRand* pRand;

private:
    static TList* fgAppBase;

    ClassDef(G2PAppBase, 1)
};

#endif
