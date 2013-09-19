// -*- C++ -*-

/* class G2PVarList
 * A collection of G2PVar variables.
 */

// History:
//   Apr 2013, C. Gu, First public version.
//

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TList.h"
#include "TRegexp.h"
#include "TString.h"

#include "G2PVar.hh"
#include "G2PVarDef.hh"

#include "G2PVarList.hh"

using namespace std;

G2PVarList::G2PVarList() {
    // Nothing to do
}

G2PVarList::~G2PVarList() {
    Clear();
}

G2PVar* G2PVarList::DefineByType(const char* name, const char* descript, const void* var, VarType type) {
    // Define a variable in the list with given type

    static const char* const here = "DefineByType()";

    // Check duplicate names
    G2PVar* ptr = Find(name);
    if (ptr) {
        Warning(here, "Variable %s already exists. Not redefined.", ptr->GetName());
        return NULL;
    }

    switch (type) {
    case kBOOL:
        ptr = new G2PVar(name, descript, static_cast<const bool*> (var));
        break;
    case kCHAR:
        ptr = new G2PVar(name, descript, static_cast<const char*> (var));
        break;
    case kINT:
        ptr = new G2PVar(name, descript, static_cast<const int*> (var));
        break;
    case kSHORT:
        ptr = new G2PVar(name, descript, static_cast<const short*> (var));
        break;
    case kLONG:
        ptr = new G2PVar(name, descript, static_cast<const long*> (var));
        break;
    case kFLOAT:
        ptr = new G2PVar(name, descript, static_cast<const float*> (var));
        break;
    case kDOUBLE:
        ptr = new G2PVar(name, descript, static_cast<const double*> (var));
        break;
    }

    if (ptr)
        AddLast(ptr);
    else
        Error(here, "Error allocating new variable %s. No variable defined.", name);

    return ptr;
}

int G2PVarList::DefineVariables(const VarDef* list, const char* prefix) {
    // Add variables in "list" to the list

    static const char* const here = "DefineVariables()";

    if (!list) {
        Warning(here, "No input variable list.");
        return -1;
    }

    const VarDef* item = list;
    int ndef = 0;
    while (item->name) {
        const char* description = item->desc;
        if (!description || !*description) description = item->name;
        G2PVar* var = DefineByType(Form("%s%s", prefix, item->name), description, item->loc, item->type);
        if (var) ndef++;
        item++;
    }

    return ndef;
}

G2PVar* G2PVarList::Find(const char* name) const {
    return static_cast<G2PVar*> (FindObject(name));
}

G2PVar* G2PVarList::FindSuffix(const char* suf) const {
    static const char* const here = "FindSuffix()";

    G2PVar* p = NULL;
    int nfind = 0;
    TString name;
    TIter next(this);
    while (G2PVar * ptr = static_cast<G2PVar*> (next())) {
        name = ptr->GetName();
        if (name.EndsWith(suf)) {
            if (nfind == 0) p = ptr;
            nfind++;
        }
    }

    if (nfind > 1) {
        Warning(here, "Find %d variables name contain \"%s\", check definition.", nfind, suf);
    }

    return p;
}

G2PVar* G2PVarList::FindRegexp(const char* expr) const {
    static const char* const here = "FindRegexp()";

    G2PVar* p = NULL;
    TRegexp re(expr, kTRUE);
    if (re.Status()) return NULL;

    int nfind = 0;
    TString name;
    TIter next(this);
    while (G2PVar * ptr = static_cast<G2PVar*> (next())) {
        name = ptr->GetName();
        if (name.Index(re) != kNPOS) {
            if (nfind == 0) p = ptr;
            nfind++;
        }
    }

    if (nfind > 1) {
        Warning(here, "Find %d variables name satisfy \"%s\", check definition.", nfind, expr);
    }

    return p;
}

int G2PVarList::RemoveName(const char* name) {
    G2PVar* ptr = Find(name);
    if (ptr) {
        G2PVar* p = static_cast<G2PVar*> (TList::Remove(ptr));
        delete p;
        return 1;
    }
    else return 0;
}

int G2PVarList::RemoveRegexp(const char* expr) {
    TRegexp re(expr, kTRUE);
    if (re.Status()) return -1;

    int ndel = 0;
    TString name;
    TIter next(this);
    while (G2PVar * ptr = static_cast<G2PVar*> (next())) {
        name = ptr->GetName();
        if (name.Index(re) != kNPOS) {
            G2PVar* p = static_cast<G2PVar*> (TList::Remove(ptr));
            delete p;
            ndel++;
        }
    }

    return ndel;
}

ClassImp(G2PVarList)
