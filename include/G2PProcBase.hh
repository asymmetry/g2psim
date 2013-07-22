// -*- C++ -*-

/* class G2PProcBase
 * This file defines a class G2PProcBase.
 * It is the base class of g2p process classes.
 * It provides fundamental functions like passing variables between processes.
 */

// History:
//   Apr 2013, C. Gu, First public version.
//

#ifndef G2P_PROCBASE_H
#define G2P_PROCBASE_H

#include <cstring>
#include <vector>
#include <map>

#include "G2PAppBase.hh"

using namespace std;

class TList;

class G2PProcBase : public G2PAppBase {
public:
    virtual ~G2PProcBase();

    virtual int Init();
    virtual int Begin();
    virtual int Process() = 0;

    virtual void Clear() { }

    virtual int SetValue(const char* name, double* value);
    virtual int GetValue(const char* name, double* value);

protected:
    G2PProcBase(); // No instance allowed for this class

    int Add(const char* name);

    int ArrayCopy(double* out, const double* in, int length);

    virtual int DefineVariables(EMode mode = kDefine) {
        return 0;
    }

    virtual void MakePrefix() { }

    map<string, double*> mName;
    map<string, int> mLength;

    TList* fApps;
    vector<const char*> fAppsList;

private:
    ClassDef(G2PProcBase, 1)
};

#endif
