// -*- C++ -*-

/* class G2PEvGen
 * Abstract base class of g2p event generator classes.
 * Use function Process() to get reaction point kinematics.
 * It will calculate the energy loss before the scattering point.
 */

// History:
//   Mar 2013, C. Gu, First public version.
//   Sep 2013, C. Gu, Rewrite it as a G2PProcBase class.
//   Nov 2014, J. Liu, Add calculation of the energy loss before the scattering.
//   Dec 2014, C. Gu, Rewrite with new G2P geometry classes.
//   Dec 2014, C. Gu, Merge G2PFlatGun into this class and rename this class to G2PEvGen.
//

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"

#include "G2PAppBase.hh"
#include "G2PAppList.hh"
#include "G2PGeoBase.hh"
#include "G2PGlobals.hh"
#include "G2PMaterial.hh"
#include "G2PProcBase.hh"
#include "G2PRand.hh"
#include "G2PVarDef.hh"

#include "G2PEvGen.hh"

using namespace std;

static const double kDEG = 3.14159265358979323846 / 180.0;

G2PEvGen *G2PEvGen::pG2PEvGen = NULL;

G2PEvGen::G2PEvGen() : fUseTrans(true), fE0(0.0), fm(0.0), fM0(0.0), fFieldRatio(0.0), fForceElastic(false), fBeamX_lab(0.0), fBeamY_lab(0.0), fBeamZ_lab(0.0), fBeamR_lab(0.0), fE(0.0), fELoss(0.0), fTb(0.0), fReactZLow_lab(0.0), fReactZHigh_lab(0.0), fTargetThLow_tr(0.0), fTargetThHigh_tr(0.0), fTargetPhLow_tr(0.0), fTargetPhHigh_tr(0.0), fDeltaLow(0.0), fDeltaHigh(0.0), fTiltTheta_bpm(0.0), fTiltPhi_bpm(0.0)
{
    if (pG2PEvGen) {
        Error("G2PEvGen()", "Only one instance of G2PEvGen allowed.");
        MakeZombie();
        return;
    }

    pG2PEvGen = this;

    fPriority = 1;
    Clear();
}

G2PEvGen::~G2PEvGen()
{
    if (pG2PEvGen == this)
        pG2PEvGen = NULL;
}

int G2PEvGen::Begin()
{
    //static const char* const here = "Begin()";

    if (G2PProcBase::Begin() != 0)
        return (fStatus = kBEGINERROR);

    SetTiltAngle();

    return (fStatus = kOK);
}

int G2PEvGen::Process()
{
    static const char *const here = "Process()";

    if (fDebug > 2)
        Info(here, " ");

    double X_lab, Y_lab;

    if (fBeamR_lab > 1e-5) {
        do {
            X_lab = pRand->Uniform(-fBeamR_lab, fBeamR_lab);
            Y_lab = pRand->Uniform(-fBeamR_lab, fBeamR_lab);
        } while (X_lab * X_lab + Y_lab * Y_lab > fBeamR_lab * fBeamR_lab);
    } else {
        X_lab = 0.0;
        Y_lab = 0.0;
    }

    X_lab += fBeamX_lab;
    Y_lab += fBeamY_lab;
    double ReactZ_lab = pRand->Uniform(fReactZLow_lab, fReactZHigh_lab);

    GetReactPoint(X_lab, Y_lab, ReactZ_lab, fV5beam_lab);

    if (fDebug > 1)
        Info(here, "beam_lab  : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5beam_lab[0], fV5beam_lab[1], fV5beam_lab[2], fV5beam_lab[3], fV5beam_lab[4]);

    HCS2TCS(fV5beam_lab[0], fV5beam_lab[2], fV5beam_lab[4], fV5react_tr[0], fV5react_tr[2], freactz_tr);
    fV5react_lab[0] = fV5beam_lab[0];
    fV5react_lab[2] = fV5beam_lab[2];
    fV5react_lab[4] = fV5beam_lab[4];

    if (fUseTrans) {
        fV5react_tr[1] = pRand->Uniform(fTargetThLow_tr, fTargetThHigh_tr);
        fV5react_tr[3] = pRand->Uniform(fTargetPhLow_tr, fTargetPhHigh_tr);
        TCS2HCS(fV5react_tr[1], fV5react_tr[3], fV5react_lab[1], fV5react_lab[3]);
    } else {
        double cos_theta_high = cos(fTargetThLow_tr);
        double cos_theta_low = cos(fTargetThHigh_tr); // In lab coordinates, the scattering angle is always larger than 0
        double cos_theta = pRand->Uniform(cos_theta_low, cos_theta_high);
        fV5react_lab[1] = acos(cos_theta);
        fV5react_lab[3] = pRand->Uniform(fTargetPhLow_tr, fTargetPhHigh_tr);
        HCS2TCS(fV5react_lab[1], fV5react_lab[3], fV5react_tr[1], fV5react_tr[3]);
    }

    // Calculate the beam energy after energy loss
    double x[3] = {fV5beam_lab[0], fV5beam_lab[2], fV5beam_lab[4]};
    double p[3] = {fE0 * sin(fV5beam_lab[1]) *cos(fV5beam_lab[3]), fE0 * sin(fV5beam_lab[1]) *sin(fV5beam_lab[3]), fE0 * cos(fV5beam_lab[1])};

    double l, eloss;
    TIter next(gG2PGeos);
    fE = fE0;

    while (G2PGeoBase *geo = static_cast<G2PGeoBase *>(next())) {
        if (geo->IsInside(x)) {
            l = Drift("backward", x, p, geo, x, p);
            G2PMaterial *mat = geo->GetMaterial();

            if (mat) {
                eloss = mat->EnergyLoss(fE, l);
                fELoss += eloss;

                for (int i = 0; i < 3; i++)
                    p[i] *= (fE - eloss) / fE;

                fE -= eloss;
                fTb += (l * 100) / (mat->GetRadLen() / mat->GetDensity());
            }

            next.Reset();
        }
    }

    // Additional energy loss (Be window in this case)
    static G2PMaterial Be("Be", 4, 9.0122, 65.19, 1.85, 63.7, 2.7847);
    eloss = Be.EnergyLoss(fE, 0.0381e-2); // Be window 0.0381cm
    fELoss += eloss;
    fE -= eloss;
    fTb += 0.0381 / (Be.GetRadLen() / Be.GetDensity());

    // Calculate Internal Bremsstrahlung
    double Pi[3] = {sin(fV5beam_lab[1]) *cos(fV5beam_lab[3]), sin(fV5beam_lab[1]) *sin(fV5beam_lab[3]), cos(fV5beam_lab[1])};
    double Pf[3] = {sin(fV5react_lab[1]) *cos(fV5react_lab[3]), sin(fV5react_lab[1]) *sin(fV5react_lab[3]), cos(fV5react_lab[1])};
    double cosang = Pi[0] * Pf[0] + Pi[1] * Pf[1] + Pi[2] * Pf[2];

    double scatangle = acos(cosang);
    eloss = InterBremsstrahlung(fE, scatangle);
    fELoss += eloss;
    fE -= eloss;

    if (fDebug > 2)
        Info("EnergyLoss()", "%10.3e %10.3e", fELoss, fTb);

    // calculate elastic scattering momentum
    double P = sqrt(fE * fE - fm * fm);
    double scatmom = (P * fM0 / (fE + fM0 - P * cosang)) * (((fE + fM0) * sqrt(1 - (fm / fM0) * (fm / fM0) * (1 - cosang * cosang)) + (fE + (fm / fM0) * fm) * cosang) / (fE + fM0 + P * cosang));
    double elasticd = scatmom / fHRSMomentum - 1;

    if (fForceElastic)
        fV5react_tr[4] = elasticd;
    else { // no larger than elastic momentum
        double tempd = pRand->Uniform(fDeltaLow, fDeltaHigh);

        if (tempd > elasticd)
            return -1;
        else fV5react_tr[4] = tempd;
    }

    if (freactz_tr < 0.0)
        Drift("forward", fV5react_tr, freactz_tr, 0.0, fV5tp_tr); // Drift to target plane (z_tr = 0)
    else
        Drift("backward", fV5react_tr, freactz_tr, 0.0, fV5tp_tr);

    if (fDebug > 1)
        Info(here, "tp_tr     : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5tp_tr[0], fV5tp_tr[1], fV5tp_tr[2], fV5tp_tr[3], fV5tp_tr[4]);

    return 0;
}

void G2PEvGen::Clear(Option_t *opt)
{
    fE = 0;
    fELoss = 0;
    fTb = 0;

    freactz_tr = 0;

    memset(fV5beam_lab, 0, sizeof(fV5beam_lab));
    memset(fV5react_tr, 0, sizeof(fV5react_tr));
    memset(fV5react_lab, 0, sizeof(fV5react_lab));
    memset(fV5tp_tr, 0, sizeof(fV5tp_tr));

    G2PProcBase::Clear(opt);
}

void G2PEvGen::SetBeamPos(double x, double y, double z)
{
    fBeamX_lab = x;
    fBeamY_lab = y;
    fBeamZ_lab = z;

    fConfigIsSet.insert((unsigned long) &fBeamX_lab);
    fConfigIsSet.insert((unsigned long) &fBeamY_lab);
    fConfigIsSet.insert((unsigned long) &fBeamZ_lab);
}

void G2PEvGen::SetTiltAngle(double theta, double phi)
{
    fTiltTheta_bpm = theta;
    fTiltPhi_bpm = phi;

    fConfigIsSet.insert((unsigned long) &fTiltTheta_bpm);
    fConfigIsSet.insert((unsigned long) &fTiltPhi_bpm);
}

void G2PEvGen::SetReactZ(double low, double high)
{
    fReactZLow_lab = low;
    fReactZHigh_lab = high;

    fConfigIsSet.insert((unsigned long) &fReactZLow_lab);
    fConfigIsSet.insert((unsigned long) &fReactZHigh_lab);
}

void G2PEvGen::SetRasterSize(double val)
{
    fBeamR_lab = val;

    fConfigIsSet.insert((unsigned long) &fBeamR_lab);
}

void G2PEvGen::SetTargetTh(double low, double high)
{
    fTargetThLow_tr = low;
    fTargetThHigh_tr = high;

    fConfigIsSet.insert((unsigned long) &fTargetThLow_tr);
    fConfigIsSet.insert((unsigned long) &fTargetThHigh_tr);
}

void G2PEvGen::SetTargetPh(double low, double high)
{
    fTargetPhLow_tr = low;
    fTargetPhHigh_tr = high;

    fConfigIsSet.insert((unsigned long) &fTargetPhLow_tr);
    fConfigIsSet.insert((unsigned long) &fTargetPhHigh_tr);
}

void G2PEvGen::SetDelta(double low, double high)
{
    fDeltaLow = low;
    fDeltaHigh = high;

    fConfigIsSet.insert((unsigned long) &fDeltaLow);
    fConfigIsSet.insert((unsigned long) &fDeltaHigh);
}

void G2PEvGen::SetDelta(const char *elastic)
{
    fForceElastic = true;

    fConfigIsSet.insert((unsigned long) &fForceElastic);
}

void G2PEvGen::SetCoords(const char *coords)
{
    string s = coords;

    if (s == "lab")
        fUseTrans = false;

    fConfigIsSet.insert((unsigned long) &fUseTrans);
}

void G2PEvGen::SetTiltAngle()
{
    static const char *const here = "SetTiltAngle()";

    if ((fConfigIsSet.find((unsigned long) &fTiltTheta_bpm) == fConfigIsSet.end())
            && (fConfigIsSet.find((unsigned long) &fTiltPhi_bpm) == fConfigIsSet.end())) {
        // Default values
        if (fabs(fFieldRatio - 0.5) < 1e-8) {
            if (fabs(fE0 - 2.254) < 0.2)
                fTiltTheta_bpm = 3.31 * kDEG;
            else if (fabs(fE0 - 1.706) < 0.2)
                fTiltTheta_bpm = 4.03 * kDEG;
            else if (fabs(fE0 - 1.158) < 0.2)
                fTiltTheta_bpm = 5.97 * kDEG;
            else
                fTiltTheta_bpm = 0.0;
        } else if (fabs(fFieldRatio - 1.0) < 1e-8) {
            if (fabs(fE0 - 2.254) < 0.2)
                fTiltTheta_bpm = 0.0;
            else if (fabs(fE0 - 3.355) < 0.2)
                fTiltTheta_bpm = 0.0;
            else
                fTiltTheta_bpm = 0.0;
        } else
            fTiltTheta_bpm = 0.0;

        fTiltPhi_bpm = 0.0;

        if (fDebug > 0)
            Info(here, "Default set to %10.3e %10.3e", fTiltTheta_bpm / kDEG, fTiltPhi_bpm / kDEG);
    } else {
        if (fDebug > 0)
            Info(here, "Manually set to %10.3e %10.3e", fTiltTheta_bpm / kDEG, fTiltPhi_bpm / kDEG);
    }
}

void G2PEvGen::GetReactPoint(double x, double y, double reactz, double *V5)
{
    static const char *const here = "GetReactPoint()";

    double xb[3] = {x, y, fBeamZ_lab};
    double p[3] = {tan(fTiltPhi_bpm), tan(fTiltTheta_bpm), 1.0};
    double pp = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
    double pb[3] = {fE0 *p[0] / pp, fE0 *p[1] / pp, fE0 *p[2] / pp};    // convert tilt angle from bpm coordinate system to lab system

    if (fBeamZ_lab < reactz)
        Drift("forward", xb, pb, reactz, xb, pb);
    else
        Drift("backward", xb, pb, reactz, xb, pb);

    V5[0] = xb[0];

    if (fabs(pb[2] - fE0) < 1e-8) {
        V5[1] = acos(1.0);    // pb[2] may be a bit larger than fBeam Energy because of round-off error
    } else
        V5[1] = acos(pb[2] / fE0);

    V5[2] = xb[1];
    V5[3] = atan2(pb[1], pb[0]);
    V5[4] = xb[2];

    if (fDebug > 2)
        Info(here, "%10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e %10.3e %10.3e", x, y, fBeamZ_lab, V5[0], V5[1], V5[2], V5[3], V5[4]);
}

int G2PEvGen::Configure(EMode mode)
{
    if ((mode == kREAD || mode == kTWOWAY) && fConfigured)
        return 0;

    if (G2PProcBase::Configure(mode) != 0)
        return -1;

    ConfDef confs[] = {
        {"run.e0", "Beam Energy", kDOUBLE, &fE0},
        {"run.particle.mass", "Beam Particle Mass", kDOUBLE, &fm},
        {"run.target.mass", "Target Mass", kDOUBLE, &fM0},
        {"field.ratio", "Field Ratio", kDOUBLE, &fFieldRatio},
        {"beam.l_x", "Beam X", kDOUBLE, &fBeamX_lab},
        {"beam.l_y", "Beam Y", kDOUBLE, &fBeamY_lab},
        {"beam.l_z", "Beam Z (set by BPM)", kDOUBLE, &fBeamZ_lab},
        {"beam.tilt.t", "Beam Tilt Angle Theta", kDOUBLE, &fTiltTheta_bpm},
        {"beam.tilt.p", "Beam Tilt Angle Phi", kDOUBLE, &fTiltPhi_bpm},
        {"raster.size", "Raster Size", kDOUBLE, &fBeamR_lab},
        {"react.t.min", "Theta Min", kDOUBLE, &fTargetThLow_tr},
        {"react.t.max", "Theta Max", kDOUBLE, &fTargetThHigh_tr},
        {"react.p.min", "Phi Min", kDOUBLE, &fTargetPhLow_tr},
        {"react.p.max", "Phi Max", kDOUBLE, &fTargetPhHigh_tr},
        {"react.d.min", "Delta Min", kDOUBLE, &fDeltaLow},
        {"react.d.max", "Delta Max", kDOUBLE, &fDeltaHigh},
        {"react.d.elastic", "Delta Elastic", kBOOL, &fForceElastic},
        {"react.l_z.min", "React Z Min", kDOUBLE, &fReactZLow_lab},
        {"react.l_z.max", "React Z Max", kDOUBLE, &fReactZHigh_lab},
        {0}
    };
    return ConfigureFromList(confs, mode);
}

int G2PEvGen::DefineVariables(EMode mode)
{
    if (mode == kDEFINE && fDefined)
        return 0;

    if (G2PProcBase::DefineVariables(mode) != 0)
        return -1;

    VarDef gvars[] = {
        {"e", "Beam Energy after Energy Loss", kDOUBLE, &fE},
        {"eloss.b", "Energy Loss before Scattering", kDOUBLE, &fELoss},
        {"tb", "Relative Thickness before Scattering", kDOUBLE, &fTb},
        {0}
    };

    if (DefineVarsFromList("phys.", gvars, mode) != 0)
        return -1;

    VarDef vars[] = {
        {"beam.l_x", "Beam X (lab)", kDOUBLE, &fV5beam_lab[0]},
        {"beam.l_t", "Beam T (lab)", kDOUBLE, &fV5beam_lab[1]},
        {"beam.l_y", "Beam Y (lab)", kDOUBLE, &fV5beam_lab[2]},
        {"beam.l_p", "Beam P (lab)", kDOUBLE, &fV5beam_lab[3]},
        {"beam.l_z", "Beam Z (lab)", kDOUBLE, &fV5beam_lab[4]},
        {"react.x", "React Point X", kDOUBLE, &fV5react_tr[0]},
        {"react.t", "React Point T", kDOUBLE, &fV5react_tr[1]},
        {"react.y", "React Point Y", kDOUBLE, &fV5react_tr[2]},
        {"react.p", "React Point P", kDOUBLE, &fV5react_tr[3]},
        {"react.z", "React Point Z", kDOUBLE, &freactz_tr},
        {"react.d", "React Point D", kDOUBLE, &fV5react_tr[4]},
        {"react.l_x", "React Point X (lab)", kDOUBLE, &fV5react_lab[0]},
        {"react.l_t", "React Point T (lab)", kDOUBLE, &fV5react_lab[1]},
        {"react.l_y", "React Point Y (lab)", kDOUBLE, &fV5react_lab[2]},
        {"react.l_p", "React Point P (lab)", kDOUBLE, &fV5react_lab[3]},
        {"react.l_z", "React Point Z (lab)", kDOUBLE, &fV5react_lab[4]},
        {"tp.x", "Target Plane X", kDOUBLE, &fV5tp_tr[0]},
        {"tp.t", "Target Plane T", kDOUBLE, &fV5tp_tr[1]},
        {"tp.y", "Target Plane Y", kDOUBLE, &fV5tp_tr[2]},
        {"tp.p", "Target Plane P", kDOUBLE, &fV5tp_tr[3]},
        {"tp.d", "Target Plane D", kDOUBLE, &fV5tp_tr[4]},
        {0}
    };

    return DefineVarsFromList(vars, mode);
}

void G2PEvGen::MakePrefix()
{
    const char *base = "gen";

    G2PAppBase::MakePrefix(base);
}

ClassImp(G2PEvGen)
