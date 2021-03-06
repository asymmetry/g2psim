#ifndef ROOT_VarDef
#define ROOT_VarDef

//////////////////////////////////////////////////////////////////////////
//
// VarDef
//
// Structures for calls to THaVarList::DefineVariables() and derived
// functions.
// 
// Also defines structure for THaAnalysisModule::LoadDB().
//
//////////////////////////////////////////////////////////////////////////

#include "VarType.h"

struct VarDef {
  const char*      name;     // Variable name
  const char*      desc;     // Variable description
  Int_t            type;     // Variable data type (see VarType.h)
  Int_t            size;     // Size of array (0/1 = scalar)
  const void*      loc;      // Location of data
  const Int_t*     count;    // Optional: Actual size of variable size array
};

struct RVarDef {
  const char*      name;     // Variable name
  const char*      desc;     // Variable description
  const char*      def;      // Definition of data (data member or method name)
};

struct TagDef {
  const char*      name;     // Key name
  void*            var;      // Destination of result (default to Double*)
  Int_t            fatal;    // optional: Error to return if tag not found (0=ignore)
  UInt_t           expected; // optional: number of elements to write/read (0/1 scaler)
  VarType          type;     // optional(default kDouble): data type (kInt, etc.)
};

struct DBRequest {
  const char*      name;     // Key name
  void*            var;      // Pointer to data (default to Double*)
  VarType          type;     // (opt) data type (see VarType.h, default Double_t)
  UInt_t           nelem;    // (opt) number of array elements (0/1 = 1 or auto)
  Bool_t           optional; // (opt) If true, missing key is ok
  Int_t            search;   // (opt) Search for key along name tree
  const char*      descript; // (opt) Key description (if 0, same as name)
};

#endif
