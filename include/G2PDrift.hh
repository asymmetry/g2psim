#ifndef G2P_DRIFT_H
#define G2P_DRIFT_H

#include "TROOT.h"
#include "TObject.h"

class G2PDrift : public TObject
{
public:
    G2PDrift();
    ~G2PDrift();

    void SetArm(bool isleftarm) { bIsLeftArm = isleftarm; }
    void SetHRSAngle(double value) { fHRSAngle = value; }
    
private:
    bool bIsLeftArm;
    double fHRSAngle;

    bool Forward();
    bool Backward();
    
    ClassDef(G2PDrift, 1);
};

#endif
