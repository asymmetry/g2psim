#include <cstdlib>
#include <cstring>
#include <cmath>

#include "TROOT.h"
#include "TObject.h"
#include "TError.h"
#include "TClass.h"
#include "TList.h"

#include "G2PGlobals.hh"
#include "G2PRand.hh"
#include "G2PVarDef.hh"
#include "G2PVarList.hh"

#include "G2PAppBase.hh"

G2PRand* G2PAppBase::pRand = G2PRand::GetInstance();
TList* G2PAppBase::fgAppBase = NULL;

G2PAppBase::G2PAppBase() :
    kLARGE(1.0e38), bIsSetup(false), fStatus(kNOTINIT),
    fPrefix(NULL), fDebug(0)
{
    if (!fgAppBase) fgAppBase = new TList;
    fgAppBase->Add(this);
}

G2PAppBase::~G2PAppBase()
{
    if (fgAppBase) {
        fgAppBase->Remove(this);
        if (fgAppBase->GetSize()==0) {
            delete fgAppBase;
            fgAppBase = NULL;
        }
    }

    delete [] fPrefix; fPrefix = 0;
}

int G2PAppBase::Init()
{
    static const char* const here = "Init()";

    if (fDebug>0) Info(here, "Initializing ...");

    if (IsZombie()) return (fStatus = kNOTINIT);

    EStatus status = kOK;

    MakePrefix();

    if (DefineVariables(kDefine)) status = kINITERROR;

    return (fStatus = status);
}

int G2PAppBase::Begin()
{
    static const char* const here = "Begin()";

    if (fDebug>0) Info(here, "Beginning ...");

    if (!IsInit()) return fStatus;

    return (fStatus = kOK);
}

int G2PAppBase::End()
{
    static const char* const here = "End()";

    if (fDebug>0) Info(here, "Ending ...");

    fStatus = kNOTINIT;

    if (DefineVariables(kDelete)) fStatus = kERROR;

    return (fStatus==kNOTINIT)?0:fStatus;
}

void G2PAppBase::TCS2HCS(double x_tr, double y_tr, double z_tr, double angle, double &x_lab, double &y_lab, double &z_lab)
{
    static const char* const here = "TCS2HCS()";

    double cosang = cos(angle);
    double sinang = sin(angle);

    x_lab =  y_tr*cosang+z_tr*sinang;
    y_lab = -x_tr;
    z_lab =  z_tr*cosang-y_tr*sinang;

    if (fDebug>3) Info(here, "%10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e\n", x_tr, y_tr, z_tr, x_lab, y_lab, z_lab);
}

void G2PAppBase::TCS2HCS(double t_tr, double p_tr, double angle, double &t_lab, double &p_lab)
{
    static const char* const here = "TCS2HCS()";

    double x = tan(t_tr);
    double y = tan(p_tr);
    double z = 1.0;
    TCS2HCS(x, y, z, angle, x, y, z);
    t_lab = acos(z/sqrt(x*x+y*y+z*z));
    p_lab = atan2(y, x);

    if (fDebug>3) Info(here, "%10.3e %10.3e -> %10.3e %10.3e\n", t_tr, p_tr, t_lab, p_lab);
}

void G2PAppBase::HCS2TCS(double x_lab, double y_lab, double z_lab, double angle, double &x_tr, double &y_tr, double &z_tr)
{
    static const char* const here = "HCS2TCS()";

    double cosang = cos(angle);
    double sinang = sin(angle);

    x_tr = -y_lab;
    y_tr =  x_lab*cosang-z_lab*sinang;
	z_tr =  x_lab*sinang+z_lab*cosang;

    if (fDebug>3) Info(here, "%10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e\n", x_lab, y_lab, z_lab, x_tr, y_tr, z_tr);
}

void G2PAppBase::HCS2TCS(double t_lab, double p_lab, double angle, double &t_tr, double &p_tr)
{
    static const char* const here = "HCS2TCS()";

    double x = sin(t_lab)*cos(p_lab);
    double y = sin(t_lab)*sin(p_lab);
    double z = cos(t_lab);
    HCS2TCS(x, y, z, angle, x, y, z);
    t_tr = atan2(x, z);
    p_tr = atan2(y, z);

    if (fDebug>3) Info(here, "%10.3e %10.3e -> %10.3e %10.3e\n", t_lab, p_lab, t_tr, p_tr);
}

void G2PAppBase::Project(double x, double y, double z, double z_out, double t, double p, double &xout, double &yout)
{
    static const char* const here = "Project()";

    double xsave = x;
    double ysave = y;

    xout = xsave + (z_out-z)*tan(t);
    yout = ysave + (z_out-z)*tan(p);

    if (fDebug>3) Info(here, "%10.3e %10.3e -> %10.3e %10.3e\n", xsave, ysave, xout, yout);
}

int G2PAppBase::DefineVarsFromList(const VarDef* list, EMode mode) const
{
    static const char* const here = "DefineVarsFromList()";

    if (!gG2PVars) {
        Warning(here, "No global variable list found.");
        return (mode==kDefine?-1:0);
    }

    if (mode==kDefine) {
        gG2PVars->DefineVariables(list, fPrefix);
    }
    else if (mode==kDelete) {
        const VarDef* item;
        while ((item = list++)&&item->name) {
            gG2PVars->RemoveName(Form("%s%s", fPrefix, item->name));
        }
    }

    return 0;
}

void G2PAppBase::MakePrefix(const char* basename)
{
    delete [] fPrefix;
    if (basename&&*basename){
        fPrefix = new char[strlen(basename)+3];
        strcpy(fPrefix, basename);
        strcat(fPrefix, ".");
    }
    else {
        fPrefix = new char[2];
        *fPrefix = 0;
    }
}

G2PAppBase* G2PAppBase::FindModule(const char* classname) const
{
    static const char* const here = "FindModule()";
    static const char* const g2papp = "G2PAppBase";

    if (fDebug>0) Info(here, "Searching %s ...", classname);

    TObjLink* lnk = fgAppBase->FirstLink();
    while (lnk) {
        TObject* obj = lnk->GetObject();
        if ((obj->IsA()->InheritsFrom(classname))&&(obj->IsA()->InheritsFrom(g2papp))) {
            G2PAppBase* aobj = static_cast<G2PAppBase*>(obj);
            return aobj;
        }
        lnk = lnk->Next();
    }
    
    if (fDebug>0) {
        Info(here, "Module %s does not exist.", classname);
    }

    return NULL;
}

void G2PAppBase::FindModule(const char* classname, TList* list) const
{
    static const char* const here = "FindModule()";

    if (fDebug>0) Info(here, "Searching %s ...", classname);

    list->Clear();
    TObjLink* lnk = fgAppBase->FirstLink();
    int n = 0;
    while (lnk) {
        TObject* obj = lnk->GetObject();
        if (obj->IsA()->InheritsFrom(classname)) {
            G2PAppBase* aobj = static_cast<G2PAppBase*>(obj);
            list->Add(aobj);
            n++;
        }
        lnk = lnk->Next();
    }

    if (fDebug>0) Info(here, "Found %d module %s.", n, classname);
}

ClassImp(G2PAppBase)
