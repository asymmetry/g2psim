// -*- C++ -*-

/* class G2PVar
 * Global variables in the simulation.
 * It can be used to retrieve data from an object.
 */

// History:
//   Apr 2013, C. Gu, First public version.
//

#include <cstdlib>
#include <cstdio>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TNamed.h"

#include "G2PVarDef.hh"

#include "G2PVar.hh"

using namespace std;

G2PVar::G2PVar() : fValueD(NULL), fType(kDOUBLE)
{
    // Nothing to do
}

G2PVar::G2PVar(const G2PVar &rhs) : TNamed(rhs), fValueD(rhs.fValueD), fType(rhs.fType)
{
    // Nothing to do
}

G2PVar &G2PVar::operator=(const G2PVar &rhs)
{
    if (this != &rhs) {
        TNamed::operator=(rhs);
        fValueD = rhs.fValueD;
        fType = rhs.fType;
    }

    return *this;
}

G2PVar::~G2PVar()
{
    // Nothing to do
}

G2PVar::G2PVar(const char *name, const char *descript, const bool *var) : TNamed(name, descript), fValueB(var), fType(kBOOL)
{
    // Nothing to do
}

G2PVar::G2PVar(const char *name, const char *descript, const char *var) : TNamed(name, descript), fValueC(var), fType(kCHAR)
{
    // Nothing to do
}

G2PVar::G2PVar(const char *name, const char *descript, const int *var) : TNamed(name, descript), fValueI(var), fType(kINT)
{
    // Nothing to do
}

G2PVar::G2PVar(const char *name, const char *descript, const short *var) : TNamed(name, descript), fValueS(var), fType(kSHORT)
{
    // Nothing to do
}

G2PVar::G2PVar(const char *name, const char *descript, const long *var) : TNamed(name, descript), fValueL(var), fType(kLONG)
{
    // Nothing to do
}

G2PVar::G2PVar(const char *name, const char *descript, const float *var) : TNamed(name, descript), fValueF(var), fType(kFLOAT)
{
    // Nothing to do
}

G2PVar::G2PVar(const char *name, const char *descript, const double *var) : TNamed(name, descript), fValueD(var), fType(kDOUBLE)
{
    // Nothing to do
}

const char *G2PVar::GetTypeName() const
{
    static const char *const type[] = {"Bool", "Char", "Int", "Short", "Long", "Float", "Double"};

    return type[fType];
}

int G2PVar::GetTypeSize() const
{
    static const int size[] = {sizeof(bool), sizeof(char), sizeof(int), sizeof(short), sizeof(long), sizeof(float), sizeof(double)};

    return size[fType];
}

double G2PVar::GetValue() const
{
    switch (fType) {
    case kBOOL:
        return static_cast<double>(*fValueB);

    case kCHAR:
        return static_cast<double>(*fValueC);

    case kINT:
        return static_cast<double>(*fValueI);

    case kSHORT:
        return static_cast<double>(*fValueS);

    case kLONG:
        return static_cast<double>(*fValueL);

    case kFLOAT:
        return static_cast<double>(*fValueF);

    case kDOUBLE:
        return static_cast<double>(*fValueD);
    }

    return 1.0e38;
}

VarType G2PVar::GetType() const
{
    return fType;
}

void G2PVar::SetVar(const bool &var)
{
    fValueB = &var;
    fType = kBOOL;
}

void G2PVar::SetVar(const char &var)
{
    fValueC = &var;
    fType = kCHAR;
}

void G2PVar::SetVar(const int &var)
{
    fValueI = &var;
    fType = kINT;
}

void G2PVar::SetVar(const short &var)
{
    fValueS = &var;
    fType = kSHORT;
}

void G2PVar::SetVar(const long &var)
{
    fValueL = &var;
    fType = kLONG;
}

void G2PVar::SetVar(const float &var)
{
    fValueF = &var;
    fType = kFLOAT;
}

void G2PVar::SetVar(const double &var)
{
    fValueD = &var;
    fType = kDOUBLE;
}

ClassImp(G2PVar)
