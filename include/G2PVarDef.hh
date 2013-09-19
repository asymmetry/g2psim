// -*- C++ -*-

/* struct VarDef
 * Type definition of G2PVar and G2POutput.
 *
 * struct ConfDef
 * Type definition of configurations of G2PAppBase.
 */

// History:
//   Apr 2013, C. Gu, First public version.
//

#ifndef G2P_VARDEF_H
#define G2P_VARDEF_H

enum VarType {
    kBOOL = 0, kCHAR, kINT, kSHORT, kLONG, kFLOAT, kDOUBLE
};

struct VarDef {
    const char* name;
    const char* desc;
    VarType type;
    const void* loc;
};

struct ConfDef {
    const char* name;
    const char* desc;
    VarType type;
    void* var;
};

#endif
