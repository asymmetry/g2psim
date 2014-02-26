// -*- C++ -*-

/* class G2PFwdProc
 * It simulates the movement of the scatted particles in the spectrometers.
 * G2PDrift, G2PHRS, G2PMaterial and G2PSieve are used in this class.
 * Input variables: fV5tp_tr, fV5react_lab (register in gG2PVars).
 */

// History:
//   Apr 2013, C. Gu, First public version.
//   Oct 2013, J. Liu, Add Energy loss and Multiple scattering.
//

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"

#include "G2PAppBase.hh"
#include "G2PAppList.hh"
#include "G2PDrift.hh"
#include "G2PGlobals.hh"
#include "G2PHRS.hh"
#include "G2PMaterial.hh"
#include "G2PProcBase.hh"
#include "G2PRand.hh"
#include "G2PSieve.hh"
#include "G2PVar.hh"
#include "G2PVarDef.hh"
#include "G2PVarList.hh"

#include "G2PFwdProc.hh"

using namespace std;

static double kPI = 3.14159265358979323846;

G2PFwdProc* G2PFwdProc::pG2PFwdProc = NULL;

G2PFwdProc::G2PFwdProc() :
fHRSAngle(0.0), fHRSMomentum(0.0), fSieveOn(false), fHoleID(-1), pDrift(NULL), pHRS(NULL), pSieve(NULL)
{
    if (pG2PFwdProc) {
        Error("G2PFwdProc()", "Only one instance of G2PFwdProc allowed.");
        MakeZombie();
        return;
    }
    pG2PFwdProc = this;

    fPriority = 3;

    Clear();
}

G2PFwdProc::~G2PFwdProc()
{
    if (pG2PFwdProc == this) pG2PFwdProc = NULL;
}

int G2PFwdProc::Init()
{
    //static const char* const here = "Init()";

    if (G2PProcBase::Init() != 0) return fStatus;

    pDrift = static_cast<G2PDrift*> (gG2PApps->Find("G2PDrift"));
    if (!pDrift) {
        pDrift = new G2PDrift();
        gG2PApps->Add(pDrift);
    }

    pSieve = static_cast<G2PSieve*> (gG2PApps->Find("G2PSieve"));
    if (!pSieve) {
        pSieve = new G2PSieve();
        gG2PApps->Add(pSieve);
    }

    pHRS = static_cast<G2PHRS*> (gG2PApps->Find("G2PHRS"));
    if (!pHRS) return (fStatus == kINITERROR);

    return (fStatus = kOK);
}

int G2PFwdProc::Begin()
{
    //static const char* const here = "Begin()";

    if (G2PProcBase::Begin() != 0) return fStatus;

    return (fStatus = kOK);
}

int G2PFwdProc::Process()
{
    static const char* const here = "Process()";

    if (fDebug > 2) Info(here, " ");

    double V5react_lab[5], V5react_tr[5];

    V5react_lab[0] = gG2PVars->FindSuffix("react.l_x")->GetValue();
    V5react_lab[1] = gG2PVars->FindSuffix("react.l_t")->GetValue();
    V5react_lab[2] = gG2PVars->FindSuffix("react.l_y")->GetValue();
    V5react_lab[3] = gG2PVars->FindSuffix("react.l_p")->GetValue();
    V5react_lab[4] = gG2PVars->FindSuffix("react.l_z")->GetValue();

    V5react_tr[0] = gG2PVars->FindSuffix("react.x")->GetValue();
    V5react_tr[1] = gG2PVars->FindSuffix("react.t")->GetValue();
    V5react_tr[2] = gG2PVars->FindSuffix("react.y")->GetValue();
    V5react_tr[3] = gG2PVars->FindSuffix("react.p")->GetValue();
    V5react_tr[4] = gG2PVars->FindSuffix("react.d")->GetValue();

    double V5troj[5]; // scatted e- trajectory

    // Define materials
    // static double packing_fraction = 0.6;
    // static double target_fZ = (10 * packing_fraction * 0.817 / 17.0305 + 2 * (1 - packing_fraction)*0.145 / 4.0026) / (packing_fraction * 0.817 / 17.0305 + 2 * (1 - packing_fraction)*0.145 / 4.0026);
    // static double target_fA = (packing_fraction * 0.817 + (1 - packing_fraction)*0.145) / (packing_fraction * 0.817 / 17.0305 + 2 * (1 - packing_fraction)*0.145 / 4.0026);
    // static double target_x0 = (packing_fraction * 0.817 + (1 - packing_fraction)*0.145) / (packing_fraction * 0.817 / 24.603003 + (1 - packing_fraction)*0.145 / 94.32);
    // static double target_density = packing_fraction * 0.817 + (1 - packing_fraction)*0.145;

    static G2PMaterial *pHe = new G2PMaterial("He", 2, 4.0026, 94.32, 0.00016);
    static G2PMaterial *pAl = new G2PMaterial("Al", 13, 26.9815, 24.01, 2.70);
    static G2PMaterial *pTi = new G2PMaterial("Ti", 22, 47.867, 16.16, 4.54);
    static G2PMaterial *pKapton = new G2PMaterial("Kapton", 5.02, 9.80, 40.56, 1.42);
    // static G2PMaterial *pTarget = new G2PMaterial("target", target_fZ, target_fA, target_x0, target_density);

    double driftlength, E, eloss, angle, rot;
    double dlentot = 0.0, elosstot = 0.0;
    double z_tr;

    HCS2TCS(V5react_lab[0], V5react_lab[2], V5react_lab[4], fHRSAngle, V5react_tr[0], V5react_tr[2], z_tr);

    // Drift to target wall or end-cap
    switch (fRunType) {
    case 10: // production target
        break;
    case 20: // 40 mil carbon target, without LHe
        RunType20(1.016e-3, V5react_tr, z_tr, V5troj, dlentot, elosstot);
        break;
    case 21: // 125 mil carbon target, without LHe
        RunType20(3.175e-3, V5react_tr, z_tr, V5troj, dlentot, elosstot);
        break;
    case 22: // 40 mil carbon target, with LHe
        RunType21(1.016e-3, V5react_tr, z_tr, V5troj, dlentot, elosstot);
        break;
    case 23: // 125 mil carbon target, with LHe
        RunType21(3.175e-3, V5react_tr, z_tr, V5troj, dlentot, elosstot);
        break;
    case 31: // pure LHe
        RunType31(V5react_tr, z_tr, V5troj, dlentot, elosstot);
    default:
        break;
    }

    // Drift to target nose window
    driftlength = pDrift->Drift(V5troj, z_tr, fHRSMomentum, fHRSAngle, 21.133e-3, V5troj, z_tr); // vertical cylinder
    E = fHRSMomentum * (1 + V5troj[4]);
    eloss = pAl->EnergyLoss(E, driftlength);
    V5troj[4] = V5troj[4] - eloss / fHRSMomentum;
    angle = pAl->MultiScattering(E, driftlength);
    rot = pRand->Uniform(0, 2 * kPI);
    V5troj[1] += angle * cos(rot);
    V5troj[3] += angle * sin(rot);
    dlentot += driftlength;
    elosstot += eloss;

    // Drift in vacuum
    driftlength = pDrift->Drift(V5troj, z_tr, fHRSMomentum, fHRSAngle, 38.1000e-3, V5troj, z_tr); // vertical cylinder
    dlentot += driftlength;

    // Drift in 4k shield
    driftlength = pDrift->Drift(V5troj, z_tr, fHRSMomentum, fHRSAngle, 38.1127e-3, V5troj, z_tr); // vertical cylinder
    E = fHRSMomentum * (1 + V5troj[4]);
    eloss = pAl->EnergyLoss(E, driftlength);
    V5troj[4] = V5troj[4] - eloss / fHRSMomentum;
    angle = pAl->MultiScattering(E, driftlength);
    rot = pRand->Uniform(0, 2 * kPI);
    V5troj[1] += angle * cos(rot);
    V5troj[3] += angle * sin(rot);
    dlentot += driftlength;
    elosstot += eloss;

    // Drift in vacuum
    driftlength = pDrift->Drift(V5troj, z_tr, fHRSMomentum, fHRSAngle, 419.100e-3, V5troj, z_tr); // vertical cylinder
    dlentot += driftlength;

    // Drift to LN2 window
    driftlength = pDrift->Drift(V5troj, z_tr, fHRSMomentum, fHRSAngle, 419.138e-3, V5troj, z_tr); // vertical cylinder
    E = fHRSMomentum * (1 + V5troj[4]);
    eloss = pAl->EnergyLoss(E, driftlength);
    V5troj[4] = V5troj[4] - eloss / fHRSMomentum;
    angle = pAl->MultiScattering(E, driftlength);
    rot = pRand->Uniform(0, 2 * kPI);
    V5troj[1] += angle * cos(rot);
    V5troj[3] += angle * sin(rot);
    dlentot += driftlength;
    elosstot += eloss;

    // Drift in vacuum
    driftlength = pDrift->Drift(V5troj, z_tr, fHRSMomentum, fHRSAngle, 479.425e-3, V5troj, z_tr); // vertical cylinder
    dlentot += driftlength;

    // Drift to chamber exit window
    driftlength = pDrift->Drift(V5troj, z_tr, fHRSMomentum, fHRSAngle, 479.933e-3, V5troj, z_tr); // vertical cylinder
    E = fHRSMomentum * (1 + V5troj[4]);
    eloss = pAl->EnergyLoss(E, driftlength);
    V5troj[4] = V5troj[4] - eloss / fHRSMomentum;
    angle = pAl->MultiScattering(E, driftlength);
    rot = pRand->Uniform(0, 2 * kPI);
    V5troj[1] += angle * cos(rot);
    V5troj[3] += angle * sin(rot);
    dlentot += driftlength;
    elosstot += eloss;

    // Local dump front face
    double x[3] = {0};
    TCS2HCS(V5troj[0], V5troj[2], z_tr, fHRSAngle, x[0], x[1], x[2]);
    double V5_lab[5];
    TCS2HCS(V5troj[1], V5troj[3], fHRSAngle, V5_lab[1], V5_lab[3]);
    double pp = fHRSMomentum * (1 + V5troj[4]);
    double p[3] = {pp * sin(V5_lab[1]) * cos(V5_lab[3]), pp * sin(V5_lab[1]) * sin(V5_lab[3]), pp * cos(V5_lab[1])};

    driftlength = 0;
    driftlength += pDrift->Drift(x, p, 640.0e-3, x, p); // along z direction in lab coordinate
    if ((fabs(x[0]) < 46.0e-3) || (fabs(x[0]) > 87.0e-3)) return -1;
    if ((x[1] < -43.0e-3) || (x[1] > 50.0e-3)) return -1;

    // Local dump back face
    driftlength += pDrift->Drift(x, p, 790.0e-3, x, p);
    if ((fabs(x[0]) < 58.0e-3) || (fabs(x[0]) > 106.0e-3)) return -1;
    if ((x[1] < -53.0e-3) || (x[1] > 58.0e-3)) return -1;

    HCS2TCS(x[0], x[1], x[2], fHRSAngle, V5troj[0], V5troj[2], z_tr);
    HCS2TCS(acos(p[2] / pp), atan2(p[1], p[0]), fHRSAngle, V5troj[1], V5troj[3]);

    // He4 energy loss
    driftlength += pDrift->Drift(V5troj, z_tr, fHRSMomentum, fHRSAngle, pSieve->GetZ(), fV5sieve_tr);
    E = fHRSMomentum * (1 + fV5sieve_tr[4]);
    eloss = pHe->EnergyLoss(E, driftlength);
    fV5sieve_tr[4] = fV5sieve_tr[4] - eloss / fHRSMomentum;
    angle = pHe->MultiScattering(E, driftlength);
    rot = pRand->Uniform(0, 2 * kPI);
    fV5sieve_tr[1] += angle * cos(rot);
    fV5sieve_tr[3] += angle * sin(rot);
    dlentot += driftlength;
    elosstot += eloss;

    // He bag
    E = fHRSMomentum * (1 + fV5sieve_tr[4]);
    driftlength = (170.9e-2 - pSieve->GetZ()) / cos(fV5sieve_tr[1]);
    eloss = pHe->EnergyLoss(E, driftlength);
    fV5sieve_tr[4] = fV5sieve_tr[4] - eloss / fHRSMomentum;
    elosstot += eloss;

    // HRS entrance
    E = fHRSMomentum * (1 + fV5sieve_tr[4]);
    driftlength = 0.0254e-2 / cos(fV5sieve_tr[1]);
    eloss = pKapton->EnergyLoss(E, driftlength);
    fV5sieve_tr[4] = fV5sieve_tr[4] - eloss / fHRSMomentum;
    elosstot += eloss;

    // HRS exit
    E = fHRSMomentum * (1 + fV5sieve_tr[4]);
    driftlength = 0.01016e-2;
    eloss = pTi->EnergyLoss(E, driftlength);
    fV5sieve_tr[4] = fV5sieve_tr[4] - eloss / fHRSMomentum;
    elosstot += eloss;

    if (fDebug > 2) {
        Info("EnergyLoss()", "%10.3e %10.3e", dlentot, elosstot);
    }

    if (fDebug > 1) {
        Info(here, "sieve_tr  : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5sieve_tr[0], fV5sieve_tr[1], fV5sieve_tr[2], fV5sieve_tr[3], fV5sieve_tr[4]);
    }

    if (fSieveOn) {
        fHoleID = pSieve->CanPass(fV5sieve_tr);
        if (fHoleID < 0) return -1;
    }

    Project(fV5sieve_tr[0], fV5sieve_tr[2], pSieve->GetZ(), 0.0, fV5sieve_tr[1], fV5sieve_tr[3], fV5tpproj_tr[0], fV5tpproj_tr[2]);
    fV5tpproj_tr[1] = fV5sieve_tr[1];
    fV5tpproj_tr[3] = fV5sieve_tr[3];
    fV5tpproj_tr[4] = fV5sieve_tr[4];

    if (fDebug > 1) {
        Info(here, "tpproj_tr : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5tpproj_tr[0], fV5tpproj_tr[1], fV5tpproj_tr[2], fV5tpproj_tr[3], fV5tpproj_tr[4]);
    }

    if (!pHRS->Forward(fV5tpproj_tr, fV5fp_tr)) return -1;
    ApplyVDCRes(fV5fp_tr);
    TRCS2FCS(fV5fp_tr, fHRSAngle, fV5fp_rot);

    if (fDebug > 1) {
        Info(here, "fp_tr     : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5fp_tr[0], fV5fp_tr[1], fV5fp_tr[2], fV5fp_tr[3], fV5fp_tr[4]);
    }

    return 0;
}

void G2PFwdProc::Clear(Option_t* option)
{
    fHoleID = -1;

    memset(fV5sieve_tr, 0, sizeof (fV5sieve_tr));
    memset(fV5tpproj_tr, 0, sizeof (fV5tpproj_tr));
    memset(fV5fp_tr, 0, sizeof (fV5fp_tr));
    memset(fV5fp_rot, 0, sizeof (fV5fp_rot));

    G2PProcBase::Clear(option);
}

void G2PFwdProc::RunType20(double thickness, double* V5react_tr, double& z_tr, double* V5troj, double& dlentot, double& elosstot)
{
    // carbon target, without LHe

    double driftlength, E, eloss, angle, rot;
    int surf;

    // Define materials
    static G2PMaterial *pC = new G2PMaterial("C", 6, 12.0107, 42.7, 2.0);
    static G2PMaterial *pPCTFE = new G2PMaterial("PCTFE", 9.33, 19.411, 28.186, 2.12);

    driftlength = pDrift->Drift(V5react_tr, z_tr, fHRSMomentum, fHRSAngle, 13.6144e-3, -14.1351e-3 + thickness, V5troj, z_tr, surf); // longitudinal cylinder
    E = fHRSMomentum * (1 + V5troj[4]);
    eloss = pC->EnergyLoss(E, driftlength);
    V5troj[4] = V5troj[4] - eloss / fHRSMomentum;
    angle = pC->MultiScattering(E, driftlength);
    rot = pRand->Uniform(0, 2 * kPI);
    V5troj[1] += angle * cos(rot);
    V5troj[3] += angle * sin(rot);
    dlentot += driftlength;
    elosstot += eloss;

    if (surf == 0) {
        driftlength = pDrift->Drift(V5troj, z_tr, fHRSMomentum, fHRSAngle, 13.6144e-3, 14.1351e-3, V5troj, z_tr, surf); // longitudinal cylinder
        dlentot += driftlength;
    }

    if (surf == 1) {
        driftlength = pDrift->Drift(V5troj, z_tr, fHRSMomentum, fHRSAngle, 14.5034e-3, 14.1351e-3, V5troj, z_tr, surf); // longitudinal cylinder
        E = fHRSMomentum * (1 + V5troj[4]);
        eloss = pPCTFE->EnergyLoss(E, driftlength);
        V5troj[4] = V5troj[4] - eloss / fHRSMomentum;
        angle = pPCTFE->MultiScattering(E, driftlength);
        rot = pRand->Uniform(0, 2 * kPI);
        V5troj[1] += angle * cos(rot);
        V5troj[3] += angle * sin(rot);
        dlentot += driftlength;
        elosstot += eloss;
    }

    // Drift in vacuum
    driftlength = pDrift->Drift(V5troj, z_tr, fHRSMomentum, fHRSAngle, 21.0058e-3, V5troj, z_tr); // vertical cylinder
    dlentot += driftlength;
}

void G2PFwdProc::RunType21(double thickness, double* V5react_tr, double& z_tr, double* V5troj, double& dlentot, double& elosstot)
{
    // carbon target, with LHe

    double driftlength, E, eloss, angle, rot;
    int surf;

    // Define materials
    static G2PMaterial *pC = new G2PMaterial("C", 6, 12.0107, 42.7, 2.0);
    static G2PMaterial *pLHe = new G2PMaterial("LHe", 2, 4.0026, 94.32, 0.145);
    static G2PMaterial *pPCTFE = new G2PMaterial("PCTFE", 9.33, 19.411, 28.186, 2.12);

    driftlength = pDrift->Drift(V5react_tr, z_tr, fHRSMomentum, fHRSAngle, 13.6144e-3, -14.1351e-3 + thickness, V5troj, z_tr, surf); // longitudinal cylinder
    E = fHRSMomentum * (1 + V5troj[4]);
    eloss = pC->EnergyLoss(E, driftlength);
    V5troj[4] = V5troj[4] - eloss / fHRSMomentum;
    angle = pC->MultiScattering(E, driftlength);
    rot = pRand->Uniform(0, 2 * kPI);
    V5troj[1] += angle * cos(rot);
    V5troj[3] += angle * sin(rot);
    dlentot += driftlength;
    elosstot += eloss;

    if (surf == 0) {
        driftlength = pDrift->Drift(V5troj, z_tr, fHRSMomentum, fHRSAngle, 13.6144e-3, 14.1351e-3, V5troj, z_tr, surf); // longitudinal cylinder
        E = fHRSMomentum * (1 + V5troj[4]);
        eloss = pLHe->EnergyLoss(E, driftlength);
        V5troj[4] = V5troj[4] - eloss / fHRSMomentum;
        angle = pLHe->MultiScattering(E, driftlength);
        rot = pRand->Uniform(0, 2 * kPI);
        V5troj[1] += angle * cos(rot);
        V5troj[3] += angle * sin(rot);
        dlentot += driftlength;
        elosstot += eloss;
    }

    if (surf == 1) {
        driftlength = pDrift->Drift(V5troj, z_tr, fHRSMomentum, fHRSAngle, 14.5034e-3, 14.1351e-3, V5troj, z_tr, surf); // longitudinal cylinder
        E = fHRSMomentum * (1 + V5troj[4]);
        eloss = pPCTFE->EnergyLoss(E, driftlength);
        V5troj[4] = V5troj[4] - eloss / fHRSMomentum;
        angle = pPCTFE->MultiScattering(E, driftlength);
        rot = pRand->Uniform(0, 2 * kPI);
        V5troj[1] += angle * cos(rot);
        V5troj[3] += angle * sin(rot);
        dlentot += driftlength;
        elosstot += eloss;
    }

    // Drift in LHe
    driftlength = pDrift->Drift(V5troj, z_tr, fHRSMomentum, fHRSAngle, 21.0058e-3, V5troj, z_tr); // vertical cylinder
    E = fHRSMomentum * (1 + V5troj[4]);
    eloss = pLHe->EnergyLoss(E, driftlength);
    V5troj[4] = V5troj[4] - eloss / fHRSMomentum;
    angle = pLHe->MultiScattering(E, driftlength);
    rot = pRand->Uniform(0, 2 * kPI);
    V5troj[1] += angle * cos(rot);
    V5troj[3] += angle * sin(rot);
    dlentot += driftlength;
    elosstot += eloss;
}

void G2PFwdProc::RunType31(double* V5react_tr, double& z_tr, double* V5troj, double& dlentot, double& elosstot)
{
    static G2PMaterial *pLHe = new G2PMaterial("LHe", 2, 4.0026, 94.32, 0.145);

    double driftlength, E, eloss, angle, rot;

    // Drift in LHe
    driftlength = pDrift->Drift(V5react_tr, z_tr, fHRSMomentum, fHRSAngle, 21.0058e-3, V5troj, z_tr); // vertical cylinder
    E = fHRSMomentum * (1 + V5troj[4]);
    eloss = pLHe->EnergyLoss(E, driftlength);
    V5troj[4] = V5troj[4] - eloss / fHRSMomentum;
    angle = pLHe->MultiScattering(E, driftlength);
    rot = pRand->Uniform(0, 2 * kPI);
    V5troj[1] += angle * cos(rot);
    V5troj[3] += angle * sin(rot);
    dlentot += driftlength;
    elosstot += eloss;
}

void G2PFwdProc::ApplyVDCRes(double* V5fp)
{
    // VDC Res set to 0.1mm (pos) and 0.5mrad (angle), HallA NIM ch3.3

    double WireChamberResX = 0.0001; //m;
    double WireChamberResY = 0.0001; //m;
    double WireChamberResT = 0.0005; //rad;
    double WireChamberResP = 0.0005; //rad;

    V5fp[0] = pRand->Gaus(V5fp[0], WireChamberResX);
    V5fp[1] = pRand->Gaus(V5fp[1], WireChamberResT);
    V5fp[2] = pRand->Gaus(V5fp[2], WireChamberResY);
    V5fp[3] = pRand->Gaus(V5fp[3], WireChamberResP);
}

int G2PFwdProc::Configure(EMode mode)
{
    if (mode == kREAD || mode == kTWOWAY) {
        if (fIsInit) return 0;
        else fIsInit = true;
    }

    ConfDef confs[] = {
        {"run.type", "Run Type", kINT, &fRunType},
        {"run.hrs.angle", "HRS Angle", kDOUBLE, &fHRSAngle},
        {"run.hrs.p0", "HRS Momentum", kDOUBLE, &fHRSMomentum},
        {"run.sieveon", "Sieve On", kBOOL, &fSieveOn},
        {0}
    };

    return ConfigureFromList(confs, mode);
}

int G2PFwdProc::DefineVariables(EMode mode)
{
    if (mode == kDEFINE && fIsSetup) return 0;
    fIsSetup = (mode == kDEFINE);

    VarDef vars[] = {
        {"id", "Hole ID", kINT, &fHoleID},
        {"sieve.x", "Sieve X", kDOUBLE, &fV5sieve_tr[0]},
        {"sieve.t", "Sieve T", kDOUBLE, &fV5sieve_tr[1]},
        {"sieve.y", "Sieve Y", kDOUBLE, &fV5sieve_tr[2]},
        {"sieve.p", "Sieve P", kDOUBLE, &fV5sieve_tr[3]},
        {"sieve.d", "Sieve D", kDOUBLE, &fV5sieve_tr[4]},
        {"tp.proj.x", "Project to Target Plane X", kDOUBLE, &fV5tpproj_tr[0]},
        {"tp.proj.t", "Project to Target Plane T", kDOUBLE, &fV5tpproj_tr[1]},
        {"tp.proj.y", "Project to Target Plane Y", kDOUBLE, &fV5tpproj_tr[2]},
        {"tp.proj.p", "Project to Target Plane P", kDOUBLE, &fV5tpproj_tr[3]},
        {"tp.proj.d", "Project to Target Plane D", kDOUBLE, &fV5tpproj_tr[4]},
        {"fp.x", "Focus Plane X", kDOUBLE, &fV5fp_tr[0]},
        {"fp.t", "Focus Plane T", kDOUBLE, &fV5fp_tr[1]},
        {"fp.y", "Focus Plane Y", kDOUBLE, &fV5fp_tr[2]},
        {"fp.p", "Focus Plane P", kDOUBLE, &fV5fp_tr[3]},
        {"fp.r_x", "Focus Plane X (rot)", kDOUBLE, &fV5fp_rot[0]},
        {"fp.r_t", "Focus Plane T (rot)", kDOUBLE, &fV5fp_rot[1]},
        {"fp.r_y", "Focus Plane Y (rot)", kDOUBLE, &fV5fp_rot[2]},
        {"fp.r_p", "Focus Plane P (rot)", kDOUBLE, &fV5fp_rot[3]},
        {0}
    };

    return DefineVarsFromList(vars, mode);
}

void G2PFwdProc::MakePrefix()
{
    const char* basename = "fwd";

    G2PAppBase::MakePrefix(basename);
}

ClassImp(G2PFwdProc)
