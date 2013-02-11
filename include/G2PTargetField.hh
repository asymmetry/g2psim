#ifndef G2P_TARGETFIELD_H
#define G2P_TARGETFIELD_H

#include "TROOT.h"
#include "TObject.h"

#include "G2PHallBField.hh"

class G2PTargetField : public TObject
{
public:
    G2PTargetField();
    G2PTargetField(const char* mapfile);
    ~G2PTargetField();
    
private:
    void ReadMap();
    
    G2PHallBField GetHallBField; // Use this as a function
    
    ClassDef(G2PTargetField, 1);
};

#endif
