// -*- C++ -*-

/* class G2PFwdHRS
 * It simulates the movement of the scatted particles in the spectrometers.
 * G2PDrift, G2PMaterial and G2PGeoSieve are used in this class.
 * Input variables: fV5tp_tr, fV5react_lab (register in gG2PVars).
 */

// History:
//   Apr 2013, C. Gu, First public version.
//   Oct 2013, J. Liu, Add Energy loss and Multiple scattering.
//   Jan 2015, C. Gu, Rewrite with geometry classes.
//

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <map>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"

#include "HRSTransBase.hh"
#include "G2PTrans400016/G2PTrans400016.hh"
#include "G2PTrans400016OLD/G2PTrans400016OLD.hh"
#include "G2PTrans484816/G2PTrans484816.hh"
#include "G2PTrans484816OLD/G2PTrans484816OLD.hh"
#include "GDHTransLargeX0/GDHTransLargeX0.hh"
#include "GDHTransSTD/GDHTransSTD.hh"
#include "HRSTransSTD/HRSTransSTD.hh"
#include "G2PAppBase.hh"
#include "G2PAppList.hh"
#include "G2PGeoBase.hh"
#include "G2PGlobals.hh"
#include "G2PMaterial.hh"
#include "G2PProcBase.hh"
#include "G2PRand.hh"
#include "G2PSieve.hh"
#include "G2PVar.hh"
#include "G2PVarDef.hh"
#include "G2PVarList.hh"

#include "G2PFwdHRS.hh"

using namespace std;

static const double kPI = 3.14159265358979323846;

G2PFwdHRS *G2PFwdHRS::pG2PFwdHRS = NULL;

G2PFwdHRS::G2PFwdHRS()
{
    //Only for ROOT I/O
}

G2PFwdHRS::G2PFwdHRS(const char *name) : fSetting(1), fSieveOn(false), fHoleID(-1), fEndPlane(0), pSieve(NULL), pModel(NULL)
{
    if (pG2PFwdHRS) {
        Error("G2PFwdHRS()", "Only one instance of G2PFwdHRS allowed.");
        MakeZombie();
        return;
    }

    pG2PFwdHRS = this;

    map<string, int> model_map;
    model_map["484816"] = 1;
    model_map["403216"] = 2;
    model_map["400016"] = 3;
    model_map["gdhLargeX0"] = 4;
    model_map["gdhSTD"] = 5;
    model_map["hrsSTD"] = 6;
    model_map["484816OLD"] = 11;
    model_map["400016OLD"] = 21;

    fSetting = model_map[name];
    fConfigIsSet.insert((unsigned long) &fSetting);

    memset(fVDCRes, 0, sizeof(fVDCRes));

    fPriority = 3;

    Clear();
}

G2PFwdHRS::~G2PFwdHRS()
{
    if (pModel) {
        delete pModel;
        pModel = NULL;
    }

    if (pG2PFwdHRS == this)
        pG2PFwdHRS = NULL;
}

int G2PFwdHRS::Begin()
{
    static const char *const here = "Begin()";

    if (G2PProcBase::Begin() != 0)
        return (fStatus = kBEGINERROR);

    pSieve = static_cast<G2PSieve *>(gG2PApps->Find("G2PSieve"));

    switch (fSetting) {
    case 1:
        pModel = new G2PTrans484816();
        break;

    case 2:
        pModel = new G2PTrans484816(); // FIXME: should be 403216 here
        break;

    case 3:
        pModel = new G2PTrans400016();
        break;

    case 4:
        pModel = new GDHTransLargeX0();
        break;

    case 5:
        pModel = new GDHTransSTD(); 
        break;

    case 6:
        pModel = new HRSTransSTD();
        break;

    case 11:
        pModel = new G2PTrans484816OLD();
        break;

    case 21:
        pModel = new G2PTrans400016OLD();
        break;

    default:
        Error(here, "Invalid setting.");
        return (fStatus = kBEGINERROR);
        break;
    }

    return (fStatus = kOK);
}

int G2PFwdHRS::Process()
{
    static const char *const here = "Process()";

    if (fDebug > 2)
        Info(here, " ");

    if (gG2PVars->FindSuffix("beam.l_x") && gG2PVars->FindSuffix("react.l_x") && gG2PVars->FindSuffix("react.x") && gG2PVars->FindSuffix("react.z")) {
        fV5beam_lab[1] = gG2PVars->FindSuffix("beam.l_t")->GetValue();
        fV5beam_lab[3] = gG2PVars->FindSuffix("beam.l_p")->GetValue();

        fV5react_lab[1] = gG2PVars->FindSuffix("react.l_t")->GetValue();
        fV5react_lab[3] = gG2PVars->FindSuffix("react.l_p")->GetValue();

        fV5react_tr[0] = gG2PVars->FindSuffix("react.x")->GetValue();
        fV5react_tr[1] = gG2PVars->FindSuffix("react.t")->GetValue();
        fV5react_tr[2] = gG2PVars->FindSuffix("react.y")->GetValue();
        fV5react_tr[3] = gG2PVars->FindSuffix("react.p")->GetValue();
        fV5react_tr[4] = gG2PVars->FindSuffix("react.d")->GetValue();

        freactz_tr = gG2PVars->FindSuffix("react.z")->GetValue();
    }

    double V5troj[5] = {fV5react_tr[0], fV5react_tr[1], fV5react_tr[2], fV5react_tr[3], fV5react_tr[4]};
    double ztroj = freactz_tr;

    double l, eloss, angle, rot;
    double E = (1 + fV5react_tr[4]) * fHRSMomentum;

    // Calculate Internal Bremsstrahlung
    double Pi[3] = {sin(fV5beam_lab[1]) *cos(fV5beam_lab[3]), sin(fV5beam_lab[1]) *sin(fV5beam_lab[3]), cos(fV5beam_lab[1])};
    double Pf[3] = {sin(fV5react_lab[1]) *cos(fV5react_lab[3]), sin(fV5react_lab[1]) *sin(fV5react_lab[3]), cos(fV5react_lab[1])};
    double cosang = Pi[0] * Pf[0] + Pi[1] * Pf[1] + Pi[2] * Pf[2];

    double scatangle = acos(cosang);
    eloss = InterBremsstrahlung(E, scatangle);
    fELoss += eloss;
    E -= eloss;
    V5troj[4] = V5troj[4] - eloss / fHRSMomentum;
    TIter next(gG2PGeos);

    while (G2PGeoBase *geo = static_cast<G2PGeoBase *>(next())) {
        if (geo->IsInside(V5troj, ztroj)) {
            l = Drift("forward", V5troj, ztroj, geo, V5troj, ztroj);
            G2PMaterial *mat = geo->GetMaterial();

            if (mat) {
                eloss = mat->EnergyLoss(E, l);
                V5troj[4] = V5troj[4] - eloss / fHRSMomentum;
                fELoss += eloss;
                E -= eloss;
                fTa += (l * 100) / (mat->GetRadLen() / mat->GetDensity());

                angle = mat->MultiScattering(E, l);
                rot = pRand->Uniform(0, 2 * kPI);
                V5troj[1] += angle * cos(rot);
                V5troj[3] += angle * sin(rot);
            }

            next.Reset();
        }
    }

    // Local dump front face
    double x[3];

    Drift("forward", V5troj, ztroj, 640.0e-3, V5troj, ztroj);
    TCS2HCS(V5troj[0], V5troj[2], ztroj, x[0], x[1], x[2]);
    fDumpFront[0] = x[0];
    fDumpFront[1] = x[1];
    if ((fabs(x[0]) < 46.0e-3) || (fabs(x[0]) > 87.0e-3) || (x[1] < -43.0e-3) || (x[1] > 50.0e-3))
      {
        fEndPlane = -2;
        return -1;
      }
    // Local dump back face
    Drift("forward", V5troj, ztroj, 790.0e-3, V5troj, ztroj);
    TCS2HCS(V5troj[0], V5troj[2], ztroj, x[0], x[1], x[2]);
    fDumpBack[0] = x[0];
    fDumpBack[1] = x[1];
    if ((fabs(x[0]) < 58.0e-3) || (fabs(x[0]) > 106.0e-3) || (x[1] < -53.0e-3) || (x[1] > 58.0e-3))
      {
        fEndPlane = -1;
        return -1;
      }
    // Sieve plane
    static G2PMaterial He("He", 2, 4.0026, 94.32, 0.00016, 41.8, 11.139393);
    l = Drift("forward", V5troj, ztroj, pSieve->GetZ(), fV5sieve_tr);
    eloss = He.EnergyLoss(E, l);
    fV5sieve_tr[4] = fV5sieve_tr[4] - eloss / fHRSMomentum;
    fELoss += eloss;
    E -= eloss;
    fTa += (l * 100) / (He.GetRadLen() / He.GetDensity());

    angle = He.MultiScattering(E, l);
    rot = pRand->Uniform(0, 2 * kPI);
    fV5sieve_tr[1] += angle * cos(rot);
    fV5sieve_tr[3] += angle * sin(rot);

    // He bag
    l = (170.9e-2 - pSieve->GetZ()) / cos(fV5sieve_tr[1]);
    eloss = He.EnergyLoss(E, l);
    fV5sieve_tr[4] = fV5sieve_tr[4] - eloss / fHRSMomentum;
    fELoss += eloss;
    E -= eloss;
    fTa += (l * 100) / (He.GetRadLen() / He.GetDensity());

    // HRS entrance
    static G2PMaterial Kapton("Kapton", 5.02, 9.80, 40.56, 1.42, 79.6, 3.3497);
    l = 0.0254e-2 / cos(fV5sieve_tr[1]);
    eloss = Kapton.EnergyLoss(E, l);
    fV5sieve_tr[4] = fV5sieve_tr[4] - eloss / fHRSMomentum;
    fELoss += eloss;
    E -= eloss;
    fTa += (l * 100) / (Kapton.GetRadLen() / Kapton.GetDensity());

    /*
    // HRS exit
    static G2PMaterial Ti("Ti", 22, 47.867, 16.16, 4.54, 230.0, 4.4450);
    l = 0.01016e-2;
    eloss = Ti.EnergyLoss(E, l);
    fV5sieve_tr[4] = fV5sieve_tr[4] - eloss / fHRSMomentum;
    fELoss += eloss;
    E -= eloss;
    fTa += (l * 100) / (Ti.GetRadLen() / Ti.GetDensity());
    */

    if (fDebug > 2)
        Info("EnergyLoss()", "%10.3e %10.3e", fELoss, fTa);

    if (fDebug > 1)
        Info(here, "sieve_tr  : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5sieve_tr[0], fV5sieve_tr[1], fV5sieve_tr[2], fV5sieve_tr[3], fV5sieve_tr[4]);

    if (fSieveOn) {
        if (!pSieve->CanPass(fV5sieve_tr, fHoleID))
            return -1;
    }

    Project(fV5sieve_tr, pSieve->GetZ(), 0.0, fV5tpproj_tr);

    if (fDebug > 1)
        Info(here, "tpproj_tr : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5tpproj_tr[0], fV5tpproj_tr[1], fV5tpproj_tr[2], fV5tpproj_tr[3], fV5tpproj_tr[4]);

    if (!Forward(fV5tpproj_tr, fV5fp_tr))
        return -1;

    ApplyVDCRes(fV5fp_tr);
    TRCS2FCS(fV5fp_tr, fV5fp_rot);

    if (fDebug > 1)
        Info(here, "fp_tr     : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5fp_tr[0], fV5fp_tr[1], fV5fp_tr[2], fV5fp_tr[3], fV5fp_tr[4]);

    return 0;
}

void G2PFwdHRS::Clear(Option_t *opt)
{
    fHoleID = -1;
    fELoss = 0;
    fTa = 0;
    freactz_tr = 0;

    memset(fV5beam_lab, 0, sizeof(fV5beam_lab));
    memset(fV5react_lab, 0, sizeof(fV5react_lab));
    memset(fV5react_tr, 0, sizeof(fV5react_tr));
    memset(fV5sieve_tr, 0, sizeof(fV5sieve_tr));
    memset(fV5tpproj_tr, 0, sizeof(fV5tpproj_tr));
    memset(fV5fp_tr, 0, sizeof(fV5fp_tr));
    memset(fV5fp_rot, 0, sizeof(fV5fp_rot));
    memset(fPlanePosX, 0, sizeof(fPlanePosX));
    memset(fPlanePosY, 0, sizeof(fPlanePosY));
    memset(fDumpFront, 0, sizeof(fDumpFront));
    memset(fDumpBack, 0, sizeof(fDumpBack));

    G2PProcBase::Clear(opt);
}

void G2PFwdHRS::SetSieve(const char *opt)
{
    string str = opt;

    if (str == "in") {
        fSieveOn = true;
        fConfigIsSet.insert((unsigned long) &fSieveOn);
    } else {
        fSieveOn = false;
        fConfigIsSet.insert((unsigned long) &fSieveOn);
    }
}

void G2PFwdHRS::SetVDCRes(double x, double t, double y, double p)
{
    fVDCRes[0] = x;
    fVDCRes[1] = t;
    fVDCRes[2] = y;
    fVDCRes[3] = p;
}

bool G2PFwdHRS::Forward(const double *V5tp_tr, double *V5fp_tr)
{
    static const char *const here = "Forward()";

    // Definition of variables
    // V5tp_tr = {x_tp, theta_tp, y_tp, phi_tp, delta@tp};
    // V5fp_tr = {x_fp, theta_fp, y_fp, phi_fp, delta@tp};
    // delta does not change

    double V5[5];

    V5[0] = V5tp_tr[0];
    V5[1] = tan(V5tp_tr[1]);
    V5[2] = V5tp_tr[2];
    V5[3] = tan(V5tp_tr[3]);
    V5[4] = V5tp_tr[4];

    bool isgood = false;

    if (fHRSAngle > 0) {
        //pModel->CoordsCorrection(fHRSAngle-fModelAngle, V5);
      fEndPlane = pModel->TransLeftHRS(V5, fPlanePosX, fPlanePosY);

        if (!fEndPlane)
            isgood = true;

        //pModel->FPCorrLeft(V5tp_tr, V5);
    } else {
        //pModel->CoordsCorrection(fHRSAngle+fModelAngle, V5);
        fEndPlane = pModel->TransRightHRS(V5, fPlanePosX, fPlanePosY);

        if (!fEndPlane)
            isgood = true;

        //pModel->FPCorrRight(V5tp_tr, V5);
    }

    V5fp_tr[0] = V5[0];
    V5fp_tr[1] = atan(V5[1]);
    V5fp_tr[2] = V5[2];
    V5fp_tr[3] = atan(V5[3]);
    V5fp_tr[4] = V5[4];

    if (fDebug > 2)
        Info(here, "%10.3e %10.3e %10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e %10.3e %10.3e", V5tp_tr[0], V5tp_tr[1], V5tp_tr[2], V5tp_tr[3], V5tp_tr[4], V5fp_tr[0], V5fp_tr[1], V5fp_tr[2], V5fp_tr[3], V5fp_tr[4]);


    return isgood;
}

void G2PFwdHRS::ApplyVDCRes(double *V5fp_tr)
{
    for (int i = 0; i < 4; i++)
        V5fp_tr[i] = pRand->Gaus(V5fp_tr[i], fVDCRes[i]);
}

int G2PFwdHRS::Configure(EMode mode)
{
    if ((mode == kREAD || mode == kTWOWAY) && fConfigured)
        return 0;

    if (G2PProcBase::Configure(mode) != 0)
        return -1;

    ConfDef confs[] = {
        {"vdc.res.x", "VDC Resolution X", kDOUBLE, &fVDCRes[0]},
        {"vdc.res.t", "VDC Resolution T", kDOUBLE, &fVDCRes[1]},
        {"vdc.res.y", "VDC Resolution Y", kDOUBLE, &fVDCRes[2]},
        {"vdc.res.p", "VDC Resolution P", kDOUBLE, &fVDCRes[3]},
        {0}
    };

    return ConfigureFromList(confs, mode);
}

int G2PFwdHRS::DefineVariables(EMode mode)
{
    if (mode == kDEFINE && fDefined)
        return 0;

    if (G2PProcBase::DefineVariables(mode) != 0)
        return -1;

    VarDef gvars[] = {
        {"eloss.a", "Energy Loss after Scattering", kDOUBLE, &fELoss},
        {"ta", "Relative Thickness after Scattering", kDOUBLE, &fTa},
        {0}
    };

    if (DefineVarsFromList("phys.", gvars, mode) != 0)
        return -1;

    VarDef vars[] = {
        {"id.hole", "Hole ID", kINT, &fHoleID},
        {"id.plane", "End Plane Pass ID", kINT, &fEndPlane},
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
        {"dump.en.x", "Local Dump Entrance X", kDOUBLE, &fDumpFront[0]},
        {"dump.en.y", "Local Dump Entrance Y", kDOUBLE, &fDumpFront[1]},
        {"dump.ex.x", "Local Dump Exit X", kDOUBLE, &fDumpBack[0]},
        {"dump.ex.y", "Local Dump Exit Y", kDOUBLE, &fDumpBack[1]},
        {"septum.en.x", "Sepump Entrance X", kDOUBLE, &fPlanePosX[5]},
        {"septum.en.y", "Sepump Entrance Y", kDOUBLE, &fPlanePosY[5]},
        {"septum.ex.x", "Sepump Exit X", kDOUBLE, &fPlanePosX[9]},
        {"septum.ex.y", "Sepump Exit Y", kDOUBLE, &fPlanePosY[9]},
        {"q1.en.x", "q1 Entrance X", kDOUBLE, &fPlanePosX[10]},
        {"q1.en.y", "q1 Entrance Y", kDOUBLE, &fPlanePosY[10]},
        {"q1.ex.x", "q1 Exit X", kDOUBLE, &fPlanePosX[13]},
        {"q1.ex.y", "q1 Exit Y", kDOUBLE, &fPlanePosY[13]},
        {"q2.en.x", "q2 Entrance X", kDOUBLE, &fPlanePosX[17]},
        {"q2.en.y", "q2 Entrance Y", kDOUBLE, &fPlanePosY[17]},
        {"q2.ex.x", "q2 Exit X", kDOUBLE, &fPlanePosX[20]},
        {"q2.ex.y", "q2 Exit Y", kDOUBLE, &fPlanePosY[20]},
        {"dipole.en.x", "dipole Entrance X", kDOUBLE, &fPlanePosX[23]},
        {"dipole.en.y", "dipole Entrance Y", kDOUBLE, &fPlanePosY[23]},
        {"dipole.ex.x", "dipole Exit X", kDOUBLE, &fPlanePosX[24]},
        {"dipole.ex.y", "dipole Exit Y", kDOUBLE, &fPlanePosY[24]},
        {"q3.en.x", "q3 Entrance X", kDOUBLE, &fPlanePosX[26]},
        {"q3.en.y", "q3 Entrance Y", kDOUBLE, &fPlanePosY[26]},
        {"q3.ex.x", "q3 Exit X", kDOUBLE, &fPlanePosX[29]},
        {"q3.ex.y", "q3 Exit Y", kDOUBLE, &fPlanePosY[29]},

        {0}
    };

    return DefineVarsFromList(vars, mode);
}

void G2PFwdHRS::MakePrefix()
{
    const char *base = "fwd";

    G2PAppBase::MakePrefix(base);
}

ClassImp(G2PFwdHRS)
