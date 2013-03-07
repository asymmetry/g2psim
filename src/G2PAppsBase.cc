#include <cstdlib>
#include <cmath>

#include "TROOT.h"
#include "TObject.h"
#include "TError.h"
#include "TClass.h"
#include "TList.h"
#include "TTree.h"

#include "G2PGlobals.hh"
#include "G2PRand.hh"

#include "G2PAppsBase.hh"

G2PRand* G2PAppsBase::pRand = G2PRand::GetInstance();
TList* G2PAppsBase::pAppsBase = NULL;

G2PAppsBase::G2PAppsBase() :
    fStatus(kNOTINIT), fDebug(0)
{
    if (!pAppsBase) pAppsBase = new TList;
    pAppsBase->Add(this);

    fApps = new TList;
}

G2PAppsBase::~G2PAppsBase()
{
    TIter next(fApps);
    while (TObject* obj = next()) fApps->Remove(obj);
    delete fApps;

    if (pAppsBase) {
        pAppsBase->Remove(this);
        if (pAppsBase->GetSize()==0) {
            delete pAppsBase;
            pAppsBase = NULL;
        }
    }
}

G2PAppsBase::EStatus G2PAppsBase::Init()
{
    static const char* const here = "Init()";

    if (fDebug>0) Info(here, "Initializing ...");

    TIter next(fApps);
    while (G2PAppsBase* aobj = static_cast<G2PAppsBase*>(next())) {
        if (!aobj->IsInit()) {
            if (aobj->Init()) return (fStatus = kINITERROR);
        }
    }

    return (fStatus = kOK);
}

void G2PAppsBase::SetSeed(int n)
{
    pRand->SetSeed(n);
}

void G2PAppsBase::TCS2HCS(double x_tr, double y_tr, double z_tr, double angle, double &x_lab, double &y_lab, double &z_lab)
{
    static const char* const here = "TCS2HCS()";

    double cosang = cos(angle);
    double sinang = sin(angle);

    x_lab =  y_tr*cosang+z_tr*sinang;
    y_lab = -x_tr;
    z_lab =  z_tr*cosang-y_tr*sinang;

    if (fDebug>3) Info(here, "%10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e\n", x_tr, y_tr, z_tr, x_lab, y_lab, z_lab);
}

void G2PAppsBase::TCS2HCS(double t_tr, double p_tr, double angle, double &t_lab, double &p_lab)
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

void G2PAppsBase::HCS2TCS(double x_lab, double y_lab, double z_lab, double angle, double &x_tr, double &y_tr, double &z_tr)
{
    static const char* const here = "HCS2TCS()";

    double cosang = cos(angle);
    double sinang = sin(angle);

    x_tr = -y_lab;
    y_tr =  x_lab*cosang-z_lab*sinang;
	z_tr =  x_lab*sinang+z_lab*cosang;

    if (fDebug>3) Info(here, "%10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e\n", x_lab, y_lab, z_lab, x_tr, y_tr, z_tr);
}

void G2PAppsBase::HCS2TCS(double t_lab, double p_lab, double angle, double &t_tr, double &p_tr)
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

G2PAppsBase* G2PAppsBase::FindModule(const char* classname)
{
    static const char* const here = "FindModule()";
    static const char* const g2papp = "G2PAppsBase";

    if (fDebug>0) Info(here, "Searching %s ...", classname);

    TObjLink* lnk = pAppsBase->FirstLink();
    while (lnk) {
        TObject* obj = lnk->GetObject();
        if ((obj->IsA()->InheritsFrom(classname))&&(obj->IsA()->InheritsFrom(g2papp))) {
            G2PAppsBase* aobj = static_cast<G2PAppsBase*>(obj);
            return aobj;
        }
        lnk = lnk->Next();
    }
    
    if (fDebug>0) {
        Info(here, "Module %s does not exist.", classname);
    }

    return NULL;
}

void G2PAppsBase::FindModule(const char* classname, TList* list)
{
    static const char* const here = "FindModule()";

    if (fDebug>0) Info(here, "Searching %s ...", classname);

    list->Clear();
    TObjLink* lnk = pAppsBase->FirstLink();
    int n = 0;
    while (lnk) {
        TObject* obj = lnk->GetObject();
        if (obj->IsA()->InheritsFrom(classname)) {
            G2PAppsBase* aobj = static_cast<G2PAppsBase*>(obj);
            list->Add(aobj);
            n++;
        }
        lnk = lnk->Next();
    }

    if (fDebug>0) Info(here, "Found %d module %s.", n, classname);
}

ClassImp(G2PAppsBase)
