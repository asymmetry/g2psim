#ifndef G2P_VARLIST_H
#define G2P_VARLIST_H

#include "TList.h"

#include "G2PVar.hh"
#include "G2PVarDef.hh"

class G2PVarList : public TList 
{
public:
    G2PVarList();
    ~G2PVarList();

    G2PVar* Define(const char* name, const char* descript, const char& var)
        { return DefineByType(name, descript, &var, kChar); }
    G2PVar* Define(const char* name, const char* descript, const int& var)
        { return DefineByType(name, descript, &var, kInt); }
    G2PVar* Define(const char* name, const char* descript, const long& var)
        { return DefineByType(name, descript, &var, kLong); }
    G2PVar* Define(const char* name, const char* descript, const short& var)
        { return DefineByType(name, descript, &var, kShort); }
    G2PVar* Define(const char* name, const char* descript, const float& var)
        { return DefineByType(name, descript, &var, kFloat); }
    G2PVar* Define(const char* name, const char* descript, const double& var)
        { return DefineByType(name, descript, &var, kDouble); }
    G2PVar* Define(const char* name, const char* descript, const bool& var)
        { return DefineByType(name, descript, &var, kBool); }

    G2PVar* Define(const char* name, const char* descript, const char* var)
        { return DefineByType(name, descript, var, kChar); }
    G2PVar* Define(const char* name, const char* descript, const int* var)
        { return DefineByType(name, descript, var, kInt); }
    G2PVar* Define(const char* name, const char* descript, const long* var)
        { return DefineByType(name, descript, var, kLong); }
    G2PVar* Define(const char* name, const char* descript, const short* var)
        { return DefineByType(name, descript, var, kShort); }
    G2PVar* Define(const char* name, const char* descript, const float* var)
        { return DefineByType(name, descript, var, kFloat); }
    G2PVar* Define(const char* name, const char* descript, const double* var)
        { return DefineByType(name, descript, var, kDouble); }
    G2PVar* Define(const char* name, const char* descript, const bool* var)
        { return DefineByType(name, descript, var, kBool); }

    void Clear();
    G2PVar* DefineByType(const char* name, const char* desc, const void* var, VarType type);
    int DefineVariables(const VarDef* list, const char* prefix = "");
    G2PVar* Find(const char* name) const { return static_cast<G2PVar*>(FindObject(name)); }
    int RemoveName(const char* name);
    int RemoveRegexp(const char* expr);

private:
    ClassDef(G2PVarList, 2)
};

#endif
