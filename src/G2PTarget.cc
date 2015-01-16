// -*- C++ -*-

/* class G2PTarget
 * Fake class to initialize and configure the geometries used in the simulation.
 * The geometries are stored in a G2PAppList.
 */

// History:
//   Dec 2014, C. Gu, Add this class for g2p geometries.
//

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"

#include "G2PAppBase.hh"
#include "G2PAppList.hh"
#include "G2PGeoBase.hh"
#include "G2PGeoTube.hh"
#include "G2PGeoSub.hh"
#include "G2PGlobals.hh"
#include "G2PMaterial.hh"
#include "G2PVarDef.hh"

#include "G2PTarget.hh"

static double kDEG = 3.14159265358979323846 / 180.0;

G2PTarget::G2PTarget()
{
    fMats = new G2PAppList();
    fGeos = new G2PAppList();
}

G2PTarget::~G2PTarget()
{
    TIter mat_iter(fMats);

    while (G2PAppBase *aobj = static_cast<G2PAppBase *>(mat_iter())) {
        fMats->Remove(aobj);
        aobj->Delete();
    }

    TIter geo_iter(fGeos);

    while (G2PAppBase *aobj = static_cast<G2PAppBase *>(geo_iter())) {
        fGeos->Remove(aobj);
        aobj->Delete();
    }

    delete fMats;
    delete fGeos;
}

int G2PTarget::Begin()
{
    if (G2PAppBase::Begin() != 0)
        return (fStatus = kBEGINERROR);

    // Define Materials
    G2PMaterial *C = new G2PMaterial("C", 6, 12.0107, 42.7, 2.0, 78.0, 2.9925);
    G2PMaterial *Al = new G2PMaterial("Al", 13, 26.9815, 24.01, 2.70, 166.0, 4.2395);
    G2PMaterial *LHe = new G2PMaterial("LHe", 2, 4.0026, 94.32, 0.145, 41.8, 4.36875);
    G2PMaterial *PCTFE = new G2PMaterial("PCTFE", 9.33, 19.411, 28.186, 2.10, 120.7, 3.8551);

    double target_Z = (10.0 * fPF * 0.817 / 17.0305 + 2.0 * (1 - fPF) * 0.145 / 4.0026) / (fPF * 0.817 / 17.0305 * 4.0 + (1 - fPF) * 0.145 / 4.0026);
    double target_A = (fPF * 0.817 + (1 - fPF) * 0.145) / (fPF * 0.817 / 17.0305 * 4.0 + (1 - fPF) * 0.145 / 4.0026);
    double target_x0 = (fPF * 0.817 + (1 - fPF) * 0.145) / (fPF * 0.817 / 40.8739 + (1 - fPF) * 0.145 / 94.32);
    double target_density = fPF * 0.817 + (1 - fPF) * 0.145;
    double target_ion = 10.0 * fPF * 0.817 * log(53.047) / 17.0305 + 2 * (1 - fPF) * 0.145 * log(41.8) / 4.0026;
    target_ion = target_ion / (fPF * 0.817 / 17.0305 * 4.0 + (1 - fPF) * 0.145 / 4.0026);
    target_ion = exp(target_ion);
    double hnup = 28.816 * sqrt(target_density * target_Z / target_A);
    double target_cor = 2.0 * log(target_ion / hnup) + 1;
    G2PMaterial *target = new G2PMaterial("target", target_Z, target_A, target_x0, target_density, target_ion, target_cor);

    fMats->Add(C);
    fMats->Add(Al);
    fMats->Add(LHe);
    fMats->Add(PCTFE);
    fMats->Add(target);

    // Define Geometries
    double min = -1.0e-8;
    double target_r = 13.6144e-3; // 0.536"
    double target_l = 28.2702e-3; // 1.113"
    double cap_thick = 0.01778e-3; // 0.7mil
    double carbon_thick = 0.0;

    if (fTargetType == 20 || fTargetType == 21)
        carbon_thick = 1.016e-3;
    else if (fTargetType == 22 || fTargetType == 23)
        carbon_thick = 3.175e-3;

    double cell_thick = 0.889e-3; // 0.035"
    double nose_r = 21.0058e-3; // 0.827"
    double nose_window_thick = 0.1016e-3; // 4mil
    double shield_4k_r = 38.1e-3; // 1.5"
    double shield_4k_thick = 0.0381e-3; // 1.5mil
    double shield_LN2_r = 419.1e-3; // 16.5"
    double shield_LN2_thick = 0.0381e-3; // 1.5mil
    double chamber_r = 479.425e-3; // 18.875"
    double chamber_h = 1.3589; // 53.5"
    double window_beam_thick;

    if ((fFieldType == 11) || (fFieldType == 20)) // longitudinal or gep
        window_beam_thick = 0.508e-3; // 20mil
    else
        window_beam_thick = 0.1778e-3; // 7mil

    double window_out_thick = 0.508e-3; // 20mil

    // Define Geometries
    G2PGeoBase *production = new G2PGeoTube(min, target_r, target_l - 2 * cap_thick);
    production->SetMaterial(target);

    G2PGeoBase *optics = new G2PGeoTube(min, target_r, carbon_thick);
    optics->SetOrigin(0, 0, -(target_l - carbon_thick) / 2);
    optics->SetMaterial(C);
    G2PGeoBase *optics_LHe = new G2PGeoTube(min, target_r, target_l - carbon_thick);
    optics_LHe->SetOrigin(0, 0, carbon_thick / 2);
    optics_LHe->SetMaterial(LHe);

    G2PGeoBase *dummy = new G2PGeoTube(min, target_r, target_l - 2 * cap_thick);
    dummy->SetMaterial(LHe);

    G2PGeoBase *dummy_nocap = new G2PGeoTube(min, target_r, target_l);
    dummy_nocap->SetMaterial(LHe);

    G2PGeoBase *empty = new G2PGeoTube(min, target_r + cell_thick, target_l);
    empty->SetMaterial(LHe);

    G2PGeoBase *cell = new G2PGeoTube(target_r, target_r + cell_thick, target_l);
    cell->SetMaterial(PCTFE);

    G2PGeoBase *cap_u = new G2PGeoTube(min, target_r, cap_thick);
    cap_u->SetOrigin(0, 0, -(target_l - cap_thick) / 2);
    cap_u->SetMaterial(Al);
    G2PGeoBase *cap_d = new G2PGeoTube(min, target_r, cap_thick);
    cap_d->SetOrigin(0, 0, (target_l - cap_thick) / 2);
    cap_d->SetMaterial(Al);

    G2PGeoSub *nose = new G2PGeoSub(new G2PGeoTube(min, nose_r, chamber_h));
    nose->SetEulerAngle(0.0, -90.0 * kDEG, 0.0); // vertical
    nose->SetMaterial(LHe);

    G2PGeoBase *window_nose = new G2PGeoTube(nose_r, nose_r + nose_window_thick, chamber_h);
    window_nose->SetEulerAngle(0.0, -90.0 * kDEG, 0.0); // vertical
    window_nose->SetMaterial(Al);

    G2PGeoBase *vacuum_4k = new G2PGeoTube(nose_r + nose_window_thick, shield_4k_r, chamber_h);
    vacuum_4k->SetEulerAngle(0.0, -90.0 * kDEG, 0.0); // vertical

    G2PGeoBase *shield_4k = new G2PGeoTube(shield_4k_r, shield_4k_r + shield_4k_thick, chamber_h);
    shield_4k->SetEulerAngle(0.0, -90.0 * kDEG, 0.0); // vertical
    shield_4k->SetMaterial(Al);

    G2PGeoBase *vacuum_LN2 = new G2PGeoTube(shield_4k_r + shield_4k_thick, shield_LN2_r, chamber_h);
    vacuum_LN2->SetEulerAngle(0.0, -90.0 * kDEG, 0.0); // vertical

    G2PGeoBase *shield_LN2 = new G2PGeoTube(shield_LN2_r, shield_LN2_r + shield_LN2_thick, chamber_h);
    shield_LN2->SetEulerAngle(0.0, -90.0 * kDEG, 0.0); // vertical
    shield_LN2->SetMaterial(Al);

    G2PGeoBase *chamber = new G2PGeoTube(shield_LN2_r + shield_LN2_thick, chamber_r, chamber_h);
    chamber->SetEulerAngle(0.0, -90.0 * kDEG, 0.0); // vertical

    G2PGeoBase *window_chamber_u = new G2PGeoTube(chamber_r, chamber_r + window_beam_thick, chamber_h, 0 * kDEG, 180 * kDEG); // beam // 0.1778mm = 7mil
    window_chamber_u->SetEulerAngle(0.0, -90.0 * kDEG, 0.0);
    window_chamber_u->SetMaterial(Al);
    G2PGeoBase *window_chamber_d = new G2PGeoTube(chamber_r, chamber_r + window_out_thick, chamber_h, 180 * kDEG, 180 * kDEG); //scattered e- // 0.508mm = 20mil
    window_chamber_d->SetEulerAngle(0.0, -90.0 * kDEG, 0.0);
    window_chamber_d->SetMaterial(Al);

    // Add target geometries to fGeos
    switch (fTargetType) {
    case 10: // production target
        fGeos->Add(production);
        fGeos->Add(cell);
        fGeos->Add(cap_u);
        fGeos->Add(cap_d);
        break;

    case 20: // optics target, 40mil carbon, no LHe
    case 22: // optics target, 125mil carbon, no LHe
        nose->SetMaterial(NULL);
        optics_LHe->SetMaterial(NULL);

    case 21: // optics target, 40mil carbon, with LHe
    case 23: // optics target, 40mil carbon, with LHe
        fGeos->Add(optics);
        fGeos->Add(optics_LHe);
        fGeos->Add(cell);
        break;

    case 30: // dummy target
        fGeos->Add(dummy);
        fGeos->Add(cell);
        fGeos->Add(cap_u);
        fGeos->Add(cap_d);
        break;

    case 31: // dummy target without end cap
        fGeos->Add(dummy_nocap);
        fGeos->Add(cell);
        break;

    case 40:
        fGeos->Add(empty);
    }

    // Subtract target geometries from nose
    TIter sub_iter(fGeos);

    while (G2PGeoBase *geo = static_cast<G2PGeoBase *>(sub_iter()))
        nose->Subtract(geo);

    fGeos->Add(nose);
    fGeos->Add(window_nose);
    fGeos->Add(vacuum_4k);
    fGeos->Add(shield_4k);
    fGeos->Add(vacuum_LN2);
    fGeos->Add(shield_LN2);
    fGeos->Add(chamber);
    fGeos->Add(window_chamber_u);
    fGeos->Add(window_chamber_d);

    TIter mat_iter(fMats);

    while (G2PMaterial *mat = static_cast<G2PMaterial *>(mat_iter())) {
        gG2PApps->Add(mat);

        if (mat->Begin() != 0)
            return (fStatus = kBEGINERROR);
    }

    TIter geo_iter(fGeos);

    while (G2PGeoBase *geo = static_cast<G2PGeoBase *>(geo_iter())) {
        gG2PApps->Add(geo);

        if (!geo->IsInit()) {
            geo->SetDebugLevel(fDebug);

            if (geo->Begin() != 0)
                return (fStatus = kBEGINERROR);
        }
    }

    gG2PGeos = fGeos;

    return (fStatus = kOK);
}

int G2PTarget::Configure(EMode mode)
{
    if ((mode == kREAD || mode == kTWOWAY) && fConfigured)
        return 0;

    if (G2PAppBase::Configure(mode) != 0)
        return -1;

    ConfDef confs[] = {
        {"run.target.type", "Run Type", kINT, &fTargetType},
        {"run.target.production.pf", "Packing Fraction", kDOUBLE, &fPF},
        {"field.type", "GEP", kINT, &fFieldType},
        {0}
    };

    return ConfigureFromList(confs, mode);
}

void G2PTarget::MakePrefix()
{
    const char *base = "target";

    G2PAppBase::MakePrefix(base);
}

ClassImp(G2PTarget)
