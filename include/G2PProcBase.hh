#ifndef G2P_PROCBASE_H
#define G2P_PROCBASE_H

#include <cstring>
#include <map>

#include "G2PAppBase.hh"

using namespace std;

class TList;

class G2PProcBase : public G2PAppBase
{
public:
    virtual ~G2PProcBase();

    virtual int Init();
    virtual int Begin();
    virtual int Process() = 0;
    virtual void Clear();

    virtual int SetValue(const char* name, double* value);
    virtual int GetValue(const char* name, double* value);

protected:
    G2PProcBase(); // No instance allowed for this class

    int ArrayCopy(double* out, const double* in, int length);

    virtual int DefineVariables(EMode mode = kDefine) { return 0; }    

    virtual void MakePrefix() { }

    bool bIsGood;

    map<string, double*> mName;
    map<string, int> mLength;

    TList* fApps;

private:
    ClassDef(G2PProcBase, 1)
};

#endif
