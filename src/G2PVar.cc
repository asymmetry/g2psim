#include "TROOT.h"
#include "TNamed.h"
#include "TError.h"

#include "G2PVarDef.hh"

#include "G2PVar.hh"

G2PVar::G2PVar() :
    fValueD(0), fType(kDouble)
{
    // Nothing to do
}

G2PVar::G2PVar(const G2PVar& rhs) :
    TNamed(rhs), fValueD(rhs.fValueD), fType(rhs.fType)
{
    // Nothing to do
}

G2PVar& G2PVar::operator=(const G2PVar& rhs) 
{
    if (this!=&rhs) {
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

double G2PVar::GetValue() const
{
    switch (fType) {
    case kChar:
        return static_cast<double>(*fValueC);
    case kInt:
        return static_cast<double>(*fValueI);
    case kShort:
        return static_cast<double>(*fValueS);
    case kLong:
        return static_cast<double>(*fValueL);
    case kFloat:
        return static_cast<double>(*fValueF);
    case kDouble:
        return static_cast<double>(*fValueD);
    case kBool:
        return static_cast<double>(*fValueB);
    }

    return 1.0e38;
}

ClassImp(G2PVar)
