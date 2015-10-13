// -*- C++ -*-

/* class G2PVarList
 * A collection of G2PVar variables.
 */

// History:
//   Apr 2013, C. Gu, First public version.
//

#ifndef G2P_VARLIST_H
#define G2P_VARLIST_H

#include "TList.h"

#include "G2PVar.hh"
#include "G2PVarDef.hh"

class G2PVarList : public TList
{
public:
    G2PVarList();
    virtual ~G2PVarList();

    // Define() with reference to variable
    G2PVar *Define(const char *name, const char *descript, const bool &var);
    G2PVar *Define(const char *name, const char *descript, const char &var);
    G2PVar *Define(const char *name, const char *descript, const int &var);
    G2PVar *Define(const char *name, const char *descript, const long &var);
    G2PVar *Define(const char *name, const char *descript, const short &var);
    G2PVar *Define(const char *name, const char *descript, const float &var);
    G2PVar *Define(const char *name, const char *descript, const double &var);

    // Define() with pointer to variable
    G2PVar *Define(const char *name, const char *descript, const bool *var);
    G2PVar *Define(const char *name, const char *descript, const char *var);
    G2PVar *Define(const char *name, const char *descript, const int *var);
    G2PVar *Define(const char *name, const char *descript, const long *var);
    G2PVar *Define(const char *name, const char *descript, const short *var);
    G2PVar *Define(const char *name, const char *descript, const float *var);
    G2PVar *Define(const char *name, const char *descript, const double *var);

    G2PVar *DefineByType(const char *name, const char *desc, const void *var, VarType type);
    int DefineVariables(const VarDef *list, const char *prefix = "");

    G2PVar *Find(const char *name) const;
    G2PVar *FindSuffix(const char *suf) const;
    G2PVar *FindRegexp(const char *expr) const;
    int RemoveName(const char *name);
    int RemoveRegexp(const char *expr);

private:
    ClassDef(G2PVarList, 1)
};

#endif
