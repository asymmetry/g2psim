#include <cstdlib>
#include <cstring>
#include <cmath>
#include <map>

#include "TROOT.h"
#include "TObject.h"
#include "TError.h"
#include "TList.h"

#include "G2PAppBase.hh"
#include "G2PGlobals.hh"

#include "G2PProcBase.hh"

using namespace std;

G2PProcBase::G2PProcBase()
{
    mName.clear();
    mLength.clear();

    fApps = new TList;

    Clear();
}

G2PProcBase::~G2PProcBase()
{
    TIter next(fApps);
    while (TObject* obj = next()) fApps->Remove(obj);
    delete fApps;
}

int G2PProcBase::Init()
{
    //static const char* const here = "Init()";

    if (G2PAppBase::Init()!=0) return fStatus;

    return (fStatus = kOK);
}

int G2PProcBase::Begin()
{
    //static const char* const here = "Begin()";

    if (G2PAppBase::Begin()!=0) return fStatus;

    EStatus status = kOK;
    TIter next(fApps);
    while (G2PProcBase* aobj = static_cast<G2PProcBase*>(next())) {
        if (!aobj->IsInit()) status = kERROR;
    }

    return (fStatus = status);
}

void G2PProcBase::Clear()
{
    G2PAppBase::Clear();

    bIsGood = false;
}

int G2PProcBase::SetValue(const char* name, double* value)
{
    if (mName.find(name)==mName.end()) return -1;
    
    return ArrayCopy(mName[name], value, mLength[name]);
}

int G2PProcBase::GetValue(const char* name, double* value)
{
    if (mName.find(name)==mName.end()) return -1;

    return ArrayCopy(value, mName[name], mLength[name]);
}

int G2PProcBase::ArrayCopy(double* out, const double* in, int length)
{
    if (!(out&&in)) return -1;
    
    for (int i = 0; i<length; i++) out[i] = in[i];

    return 0;
}

ClassImp(G2PProcBase)
