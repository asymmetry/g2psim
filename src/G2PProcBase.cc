// -*- C++ -*-

/* class G2PProcBase
 * This file defines a class G2PProcBase.
 * It is the base class of g2p process classes.
 * It provides fundamental functions like pass variables between processes.
 */

// History:
//   Apr 2013, C. Gu, First public version.
//

#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <map>

#include "TROOT.h"
#include "TObject.h"
#include "TError.h"
#include "TList.h"

#include "G2PAppBase.hh"
#include "G2PBPM.hh"
#include "G2PDBRec.hh"
#include "G2PDrift.hh"
#include "G2PGlobals.hh"
#include "G2PGun.hh"
#include "G2PPointGun.hh"
#include "G2PHRSTrans.hh"
#include "G2PPhys.hh"

#include "G2PProcBase.hh"

using namespace std;

G2PProcBase::G2PProcBase() {
    mName.clear();
    mLength.clear();

    fApps = new TList;
    fAppsList.clear();

    Clear();
}

G2PProcBase::~G2PProcBase() {
    TIter next(fApps);
    while (TObject * obj = next()) fApps->Remove(obj);
    delete fApps;
}

int G2PProcBase::Init() {
    //static const char* const here = "Init()";

    if (G2PAppBase::Init() != 0) return fStatus;

    for (vector<const char*>::iterator it = fAppsList.begin(); it != fAppsList.end(); it++) {
        if (FindModule(*it) == NULL) {
            const char* name = *it;
            if (Add(name) != 0) return (fStatus = kINITERROR);
        }
    }

    return (fStatus = kOK);
}

int G2PProcBase::Begin() {
    static const char* const here = "Begin()";

    if (G2PAppBase::Begin() != 0) return fStatus;

    if (fDebug > 0) Info(here, "Beginning ...");

    EStatus status = kOK;
    TIter next(fApps);
    while (G2PProcBase * aobj = static_cast<G2PProcBase*> (next())) {
        if (!aobj->IsInit()) status = kERROR;
    }

    return (fStatus = status);
}

int G2PProcBase::SetValue(const char* name, double* value) {
    if (mName.find(name) == mName.end()) return -1;

    return ArrayCopy(mName[name], value, mLength[name]);
}

int G2PProcBase::GetValue(const char* name, double* value) {
    if (mName.find(name) == mName.end()) return -1;

    return ArrayCopy(value, mName[name], mLength[name]);
}

int G2PProcBase::Add(const char* name) {
    static const char* const here = "Add()";

    map<string, int> temp;
    temp["G2PBPM"] = 0;
    temp["G2PDBRec"] = 1;
    temp["G2PDrift"] = 2;
    temp["G2PGun"] = 3;
    temp["G2PHRSTrans"] = 4;
    temp["G2PPhys"] = 5;

    switch (temp[name]) {
    case 0:
        gG2PApps->Add(new G2PBPM());
        break;
    case 1:
        gG2PApps->Add(new G2PDBRec());
        break;
    case 2:
        gG2PApps->Add(new G2PDrift());
        break;
    case 3:
        gG2PApps->Add(new G2PPointGun());
        break;
    case 4:
        gG2PApps->Add(new G2PHRSTrans("484816"));
        break;
    case 5:
        gG2PApps->Add(new G2PPhys("elastic"));
        break;
    default:
        Error(here, "Bad app name.");
        return -1;
    }

    return 0;
}

int G2PProcBase::ArrayCopy(double* out, const double* in, int length) {
    if (!(out && in)) return -1;

    for (int i = 0; i < length; i++) out[i] = in[i];

    return 0;
}

ClassImp(G2PProcBase)
