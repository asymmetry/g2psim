#ifndef G2P_VARDEF_H
#define G2P_VARDEF_H

enum VarType { kChar = 0, kInt, kShort, kLong, kFloat, kDouble, kBool };

struct VarDef 
{
    const char* name;
    const char* desc;
    VarType type;
    const void* loc;
};

#endif
