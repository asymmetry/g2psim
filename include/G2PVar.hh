// -*- C++ -*-

/* class G2PVar
 * This file defines a class G2PVar.
 * It defines the data structure of output variables in the simulation.
 *
 * Use this class, different types of data can be put in one TList.
 */

// History:
//   Apr 2013, C. Gu, First public version.
//

#ifndef G2P_VAR_H
#define G2P_VAR_H

#include "TNamed.h"

#include "G2PVarDef.hh"

class G2PVar : public TNamed {
public:
    G2PVar();
    G2PVar(const G2PVar& rhs);
    G2PVar& operator=(const G2PVar&);
    ~G2PVar();

    G2PVar(const char* name, const char* descript, const char* var) :
    TNamed(name, descript), fValueC(var), fType(kChar) { }

    G2PVar(const char* name, const char* descript, const int* var) :
    TNamed(name, descript), fValueI(var), fType(kInt) { }

    G2PVar(const char* name, const char* descript, const short* var) :
    TNamed(name, descript), fValueS(var), fType(kShort) { }

    G2PVar(const char* name, const char* descript, const long* var) :
    TNamed(name, descript), fValueL(var), fType(kLong) { }

    G2PVar(const char* name, const char* descript, const float* var) :
    TNamed(name, descript), fValueF(var), fType(kFloat) { }

    G2PVar(const char* name, const char* descript, const double* var) :
    TNamed(name, descript), fValueD(var), fType(kDouble) { }

    G2PVar(const char* name, const char* descript, const bool* var) :
    TNamed(name, descript), fValueB(var), fType(kBool) { }

    void SetVar(const char& var) {
        fValueC = &var;
        fType = kChar;
    }

    void SetVar(const int& var) {
        fValueI = &var;
        fType = kInt;
    }

    void SetVar(const short& var) {
        fValueS = &var;
        fType = kShort;
    }

    void SetVar(const long& var) {
        fValueL = &var;
        fType = kLong;
    }

    void SetVar(const float& var) {
        fValueF = &var;
        fType = kFloat;
    }

    void SetVar(const double& var) {
        fValueD = &var;
        fType = kDouble;
    }

    void SetVar(const bool& var) {
        fValueB = &var;
        fType = kBool;
    }

    VarType GetType() const {
        return fType;
    }
    double GetValue() const;

protected:

    union {
        const char* fValueC;
        const int* fValueI;
        const short* fValueS;
        const long* fValueL;
        const float* fValueF;
        const double* fValueD;
        const bool* fValueB;
    };

    VarType fType;

private:
    ClassDef(G2PVar, 0)
};

#endif
