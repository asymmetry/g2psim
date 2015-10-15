// -*- C++ -*-

/* class G2PEvGen
 * Abstract base class of g2p event generator classes.
 * Use function Shoot() to get reaction point kinematics.
 * Shoot() is a pure virtual method so each derived class should define its own implement.
 */

// History:
//   Mar 2013, C. Gu, First public version.
//   Sep 2013, C. Gu, Rewrite it as a G2PProcBase class.
//   Nov 2014, J. Liu, Add calculation of the energy loss before the scattering.
//   Dec 2014, C. Gu, Rewrite with new G2P geometry classes.
//   Dec 2014, C. Gu, Merge G2PFlatGun into this class and rename this class to G2PEvGen.
//

#ifndef G2P_EVGEN_H
#define G2P_EVGEN_H

#include "G2PProcBase.hh"

class G2PEvGen : public G2PProcBase
{
public:
    G2PEvGen();
    virtual ~G2PEvGen();

    virtual int Begin();
    virtual int Process();
    virtual void Clear(Option_t *opt = "");

    // Gets

    // Sets
    void SetBeamPos(double x, double y, double z);
    void SetTiltAngle(double theta, double phi);
    void SetReactZ(double low, double high);
    void SetFastRasterSize(double sizex, double sizey);
    void SetSlowRasterSize(double sizex, double sizey);
    void SetTargetTh(double low, double high);
    void SetTargetPh(double low, double high);
    void SetDelta(double low, double high);
    void SetDelta(const char *elastic);
    void SetCoords(const char *coords);

protected:
    void SetTiltAngle();

    void GetReactPoint(double x, double y, double reactz, double *V5);

    virtual int Configure(EMode mode = kTWOWAY);
    virtual int DefineVariables(EMode mode = kDEFINE);
    virtual void MakePrefix();

    bool fUseTrans;

    double fE0;
    double fm, fM0;
    double fFieldRatio;

    bool fForceElastic;

    double fBeamX_bpm, fBeamT_bpm, fBeamY_bpm, fBeamP_bpm, fBeamZ_bpm;

    double fBeamFastRx_lab, fBeamFastRy_lab;
    double fBeamSlowRx_lab, fBeamSlowRy_lab;

    double fE; // Beam energy after energy loss
    double fELoss;
    double fTb;

    double fReactZLow_lab;
    double fReactZHigh_lab;

    double fTargetThLow_tr;
    double fTargetThHigh_tr;
    double fTargetPhLow_tr;
    double fTargetPhHigh_tr;

    double fDeltaLow; // in the unit of delta
    double fDeltaHigh;

    double fV5beam_lab[5];
    double fV5react_tr[5];
    double freactz_tr;
    double fV5react_lab[5];

    double fV5tp_tr[5];

private:
    static G2PEvGen *pG2PEvGen;

    ClassDef(G2PEvGen, 1)
};

#endif
