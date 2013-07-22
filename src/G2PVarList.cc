// -*- C++ -*-

/* class G2PVarList
 * This file defines a class G2PVarList.
 * It defines a collection of G2PVar variables.
 */

// History:
//   Apr 2013, C. Gu, First public version.
//

#include <cstdlib>
#include <cmath>

#include "TROOT.h"
#include "TList.h"
#include "TError.h"
#include "TRegexp.h"

#include "G2PVar.hh"
#include "G2PVarDef.hh"

#include "G2PVarList.hh"

G2PVarList::G2PVarList() :
TList() {
    // Nothing to do
}

G2PVarList::~G2PVarList() {
    Clear();
}

void G2PVarList::Clear() {
    while (fFirst) {
        G2PVar* obj = (G2PVar*) TList::Remove(fFirst);
        delete obj;
    }
}

G2PVar* G2PVarList::DefineByType(const char* name, const char* descript, const void* var, VarType type) {
    static const char* const here = "DefineByType()";

    G2PVar* ptr = Find(name);
    if (ptr) {
        Warning(here, "Variable %s already exists. Not redefined.", ptr->GetName());
        return NULL;
    }

    switch (type) {
    case kChar:
        ptr = new G2PVar(name, descript, static_cast<const char*> (var));
        break;
    case kInt:
        ptr = new G2PVar(name, descript, static_cast<const int*> (var));
        break;
    case kShort:
        ptr = new G2PVar(name, descript, static_cast<const short*> (var));
        break;
    case kLong:
        ptr = new G2PVar(name, descript, static_cast<const long*> (var));
        break;
    case kFloat:
        ptr = new G2PVar(name, descript, static_cast<const float*> (var));
        break;
    case kDouble:
        ptr = new G2PVar(name, descript, static_cast<const double*> (var));
        break;
    case kBool:
        ptr = new G2PVar(name, descript, static_cast<const bool*> (var));
        break;
    }

    if (ptr)
        AddLast(ptr);
    else
        Error(here, "Error allocating new variable %s. No variable defined.", name);

    return ptr;
}

int G2PVarList::DefineVariables(const VarDef* list, const char* prefix) {
    if (!list) return -1;

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
