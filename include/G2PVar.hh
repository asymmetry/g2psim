// -*- C++ -*-

/* class G2PVar
 * Global variables in the simulation.
 * It can be used to retrieve data from an object.
 */

// History:
//   Apr 2013, C. Gu, First public version.
//

#ifndef G2P_VAR_H
#define G2P_VAR_H

#include "TNamed.h"

#include "G2PVarDef.hh"

class G2PVar : public TNamed
{
public:
    G2PVar();
    G2PVar(const G2PVar &rhs);
    G2PVar &operator=(const G2PVar &);
    virtual ~G2PVar();

    // Constructors by type
    G2PVar(const char *name, const char *descript, const bool *var);
    G2PVar(const char *name, const char *descript, const char *var);
    G2PVar(const char *name, const char *descript, const int *var);
    G2PVar(const char *name, const char *descript, const short *var);
    G2PVar(const char *name, const char *descript, const long *var);
    G2PVar(const char *name, const char *descript, const float *var);
    G2PVar(const char *name, const char *descript, const double *var);

    // Gets
    VarType GetType() const;
    const char *GetTypeName() const;
    int GetTypeSize() const;
    double GetValue() const;

    // Sets
    void SetVar(const bool &var);
    void SetVar(const char &var);
    void SetVar(const int &var);
    void SetVar(const short &var);
    void SetVar(const long &var);
    void SetVar(const float &var);
    void SetVar(const double &var);

protected:

    union {
        const bool *fValueB;
        const char *fValueC;
        const int *fValueI;
        const short *fValueS;
        const long *fValueL;
        const float *fValueF;
        const double *fValueD;
    };

    VarType fType;

private:
    ClassDef(G2PVar, 1)
};

#endif
