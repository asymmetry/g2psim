// -*- C++ -*-

/* class G2PGun
 * Abstract base class of g2p event generator classes.
 * Use function Shoot() to get reaction point kinematics.
 * Shoot() is a pure virtual method so each derived class should define its own implement.
 */

// History:
//   Mar 2013, C. Gu, First public version.
//   Sep 2013, C. Gu, Rewrite it as a G2PProcBase class.
//

#ifndef G2P_GUN_H
#define G2P_GUN_H

#include "G2PProcBase.hh"

using namespace std;

class G2PDrift;

class G2PGun : public G2PProcBase
{
public:
    G2PGun();
    virtual ~G2PGun();

    virtual int Init();
    virtual int Begin();
    virtual int Process();
    virtual void Clear(Option_t * /*option*/ = "");

    // Gets

    // Sets
    void SetBeamPos(double x, double y, double z);
    void SetTiltAngle(double theta, double phi);
    void SetReactZ(double low, double high);
    void SetRasterSize(double val);
    void SetTargetTh(double low, double high);
    void SetTargetPh(double low, double high);
    void SetDelta(double low, double high);
    void SetDelta(const char *elastic);

protected:
    virtual int Shoot(double *V51, double *V52) = 0;

    void SetTiltAngle();
    void GetReactPoint(double x, double y, double reactz, double *V5);

    virtual int Configure(EMode mode = kTWOWAY);
    virtual int DefineVariables(EMode mode = kDEFINE);
    virtual void MakePrefix();

    double fHRSMomentum;
    double fBeamEnergy;
    double fParticleMass;
    double fTargetMass;
    double fFieldRatio;

    bool fForceElastic;

    double fBeamX_lab, fBeamY_lab, fBeamZ_lab;
    double fBeamR_lab;

    double fReactZLow_lab;
    double fReactZHigh_lab;

    double fTargetThLow_tr;
    double fTargetThHigh_tr;
    double fTargetPhLow_tr;
    double fTargetPhHigh_tr;

    double fDeltaLow; // in the unit of delta
    double fDeltaHigh;

    double fTiltTheta_bpm;
    double fTiltPhi_bpm;

    double fV5beam_lab[5];
    double fV5react_tr[5];
    double freactz_tr;
    double fV5react_lab[5];

    double fV5tp_tr[5];

    G2PDrift *pDrift;

private:
    static G2PGun *pG2PGun;

    ClassDef(G2PGun, 1)
};

#endif
