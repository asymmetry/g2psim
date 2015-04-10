// -*- C++ -*-

/* class G2PPhys
 * It will calculate cross sections at reaction point.
 *
 * Elastic cross section models (G2PPhysEl):
 * * All: Form factors from K. C. Stansfield et al., Phys. Rev. C, 3(1971)1448
 * * H1 : Form factors from S. Venkat et al., Phys. Rev. C, 83(2011)015203 (global fit, with TPE correction)
 *                          J. Arrington et al., Phys. Rev. C 76(2007)035201 (low Q2, with/without TPE correction)
 * * 4He: Charge and magnetization densities from De Jager, At. Data Nucl. Data Tables, 14(1974)
 * * 12C: Charge distribution from L. S. Cardman et al., Phys. Lett. B, 91(1970)203
 * * 14N: Charge and magnetization densities from De Jager, At. Data Nucl. Data Tables, 14(1974)
 *
 * Inelastic cross section models:
 * * G2PPhysEPC: J. W. Lightbody et al, Computers in Physics, 2(1988)57
 * * G2PPhysPB: P. E. Bosted et al, Phys. Rev. C, 78(2008)015202 and arXiv:1203.2262
 * * G2PPhysQFS: J. W. Lightbody et al, Computers in Physics, 2(1988)57
 * * G2PPhysWISER: D. E. Wiser, Ph.D. Thesis
 *
 * Radiative correction added for P. Bosted model and QFS.
 *
 * Meaning of parameters:
 * fPID: incident particle ID, following the PDG definition:
 *       2212 for p        ;   2112 for n     ;   211 for pi+   ;
 *       -211 for pi-      ;   111  for pi0   ;   11  for e-    ;
 *       22   for photon   ;
 * fZ, fA: proton and mass number of the nucleus.
 *
 * Please also read headers of QFS, PBosted, EPC and WISER models. They contains very important usage information!
 *
 * Unit of elastic cross section is ub/sr.
 * Unit of inelastic cross section is ub/MeV-sr.
 */

// History:
//   Mar 2013, C. Gu, First public version.
//   Apr 2013, C. Gu, Add P. Bosted's model.
//   Apr 2013, C. Gu, Add WISER model.
//   May 2013, C. Gu, Add L. Cardman's C12 elastic model.
//   Oct 2013, C. Gu, Add H, He, N form factors.
//   Feb 2014, C. Gu, Add EPC model.
//   Apr 2014, C. Gu, Update H form factors.
//

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <map>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"

#include "G2PPhysBase.hh"
#include "G2PPhysEl/G2PPhysEl.hh"
#include "G2PPhysEPC/G2PPhysEPC.hh"
#include "G2PPhysPB/G2PPhysPB.hh"
#include "G2PPhysQFS/G2PPhysQFS.hh"
#include "G2PPhysWISER/G2PPhysWISER.hh"

#include "G2PAppBase.hh"
#include "G2PAppList.hh"
#include "G2PGlobals.hh"
#include "G2PProcBase.hh"
#include "G2PVar.hh"
#include "G2PVarDef.hh"
#include "G2PVarList.hh"

#include "G2PPhys.hh"

using namespace std;

static const double kDEG = 3.14159265358979323846 / 180.0;

G2PPhys *G2PPhys::pG2PPhys = NULL;

G2PPhys::G2PPhys()
{
    // Only for ROOT I/O
}

G2PPhys::G2PPhys(const char *model) : fSetting(1), fPID(11), fZ(1), fA(1), fTargetMass(0.0), fParticleMass(0.0), fPars(NULL), fNPars(0), fHRSMomentum(2.251),  pModel(NULL)
{
    if (pG2PPhys) {
        Error("G2PPhys()", "Only one instance of G2PPhys allowed.");
        MakeZombie();
        return;
    }

    pG2PPhys = this;

    fPriority = 7;
    map<string, int> model_map;
    model_map["elastic"] = 1;
    model_map["epc"] = 21;
    model_map["pbosted"] = 11;
    model_map["qfs"] = 12;
    model_map["wiser"] = 22;

    fSetting = model_map[model];
}

G2PPhys::~G2PPhys()
{
    if (pModel) {
        delete pModel;
        pModel = NULL;
    }

    if (pG2PPhys == this)
        pG2PPhys = NULL;
}

int G2PPhys::Begin()
{
    static const char *const here = "Begin()";

    if (G2PProcBase::Begin() != 0)
        return (fStatus = kBEGINERROR);

    switch (fSetting) {
    case 1:
        pModel = new G2PPhysEl();
        break;

    case 11:
        pModel = new G2PPhysPB();
        break;

    case 12:
        pModel = new G2PPhysQFS();
        break;

    case 21:
        pModel = new G2PPhysEPC();
        break;

    case 22:
        pModel = new G2PPhysWISER();
        break;

    default:
        Error(here, "Cannot initialize, invalid setting.");
        return (fStatus = kBEGINERROR);
        break;
    }

    pModel->SetTarget(fZ, fA);
    pModel->SetParticle(fPID);
    pModel->SetPars(fPars, fNPars);

    return (fStatus = kOK);
}

int G2PPhys::Process()
{
    static const char *const here = "Process()";

    if (fDebug > 2)
        Info(here, " ");

    double V51[5], V52[5];

    if (gG2PVars->FindSuffix("phys.e"))
        fE = gG2PVars->FindSuffix("phys.e")->GetValue();

    if (gG2PVars->FindSuffix("phys.tb") && gG2PVars->FindSuffix("phys.ta")) {
        fTb = gG2PVars->FindSuffix("phys.tb")->GetValue();
        fTa = gG2PVars->FindSuffix("phys.ta")->GetValue();
    }

    if (gG2PVars->FindSuffix("beam.l_x") && gG2PVars->FindSuffix("react.x")) {
        V51[0] = gG2PVars->FindSuffix("beam.l_x")->GetValue();
        V51[1] = gG2PVars->FindSuffix("beam.l_t")->GetValue();
        V51[2] = gG2PVars->FindSuffix("beam.l_y")->GetValue();
        V51[3] = gG2PVars->FindSuffix("beam.l_p")->GetValue();
        V51[4] = gG2PVars->FindSuffix("beam.l_z")->GetValue();

        V52[0] = gG2PVars->FindSuffix("react.x")->GetValue();
        V52[1] = gG2PVars->FindSuffix("react.t")->GetValue();
        V52[2] = gG2PVars->FindSuffix("react.y")->GetValue();
        V52[3] = gG2PVars->FindSuffix("react.p")->GetValue();
        V52[4] = gG2PVars->FindSuffix("react.d")->GetValue();

        fXSreact = CalXS(V51, V52, fTHreact);
    } else
        fXSreact = 0;

    if (gG2PVars->FindSuffix("bpm.l_x") && gG2PVars->FindSuffix("tp.rec.x")) {
        V51[0] = gG2PVars->FindSuffix("bpm.l_x")->GetValue();
        V51[1] = gG2PVars->FindSuffix("bpm.l_t")->GetValue();
        V51[2] = gG2PVars->FindSuffix("bpm.l_y")->GetValue();
        V51[3] = gG2PVars->FindSuffix("bpm.l_p")->GetValue();
        V51[4] = gG2PVars->FindSuffix("bpm.l_z")->GetValue();

        V52[0] = gG2PVars->FindSuffix("tp.rec.x")->GetValue();
        V52[1] = gG2PVars->FindSuffix("tp.rec.t")->GetValue();
        V52[2] = gG2PVars->FindSuffix("tp.rec.y")->GetValue();
        V52[3] = gG2PVars->FindSuffix("tp.rec.p")->GetValue();
        V52[4] = gG2PVars->FindSuffix("tp.rec.d")->GetValue();

        fXSrec = CalXS(V51, V52, fTHrec);
    } else
        fXSrec = 0;

    if (fDebug > 1) {
        Info(here, "phys_react: %10.3e %10.3e", fTHreact / kDEG, fXSreact);
        Info(here, "phys_rec  : %10.3e %10.3e", fTHrec / kDEG, fXSrec);
    }

    return 0;
}

void G2PPhys::Clear(Option_t *opt)
{
    fE = 0;
    fTb = 0;
    fTa = 0;

    fTHreact = 0.0;
    fXSreact = 0.0;
    fTHrec = 0.0;
    fXSrec = 0.0;

    G2PProcBase::Clear(opt);
}

void G2PPhys::SetPars(double *array, int n)
{
    fPars = array;
    fNPars = n;
}

double G2PPhys::CalXS(const double *V5lab, const double *V5tr, double &scatangle)
{
    static const char *const here = "CalXS()";

    double Eb[3] = {sin(V5lab[1]) *cos(V5lab[3]), sin(V5lab[1]) *sin(V5lab[3]), cos(V5lab[1])};

    double theta, phi;
    TCS2HCS(V5tr[1], V5tr[3], theta, phi);

    double Ef[3] = {sin(theta) *cos(phi), sin(theta) *sin(phi), cos(theta)};

    scatangle = acos(Eb[0] * Ef[0] + Eb[1] * Ef[1] + Eb[2] * Ef[2]);

    double Ebval = fE;
    double Efval = (1 + V5tr[4]) * fHRSMomentum;

    if (fDebug > 2)
        Info(here, "%10.3e %10.3e %10.3e", Ebval, Efval, scatangle / kDEG);

    double m = fParticleMass;
    double QQ = 2.0 * Ebval * Efval * (1- cos(scatangle));
    double radcor = (9/28.0 - 13/6.0*log(QQ/(m * m)));
    radcor -=  TDiLog((Efval - Ebval)/Efval);
    // radcor =0;
    double bornxs = pModel->GetXS(Ebval, Efval, scatangle);
    return exp(- 0.0023234* radcor) * bornxs;
}

double G2PPhys::TDiLog( double x)
{
   // The DiLogarithm function
   // Code translated by R.Brun from CERNLIB DILOG function C332

   const Double_t hf  = 0.5;
   const Double_t pi  = 3.141592627;
   const Double_t pi2 = pi*pi;
   const Double_t pi3 = pi2/3;
   const Double_t pi6 = pi2/6;
   const Double_t pi12 = pi2/12;
   const Double_t c[20] = {0.42996693560813697, 0.40975987533077105,
     -0.01858843665014592, 0.00145751084062268,-0.00014304184442340,
      0.00001588415541880,-0.00000190784959387, 0.00000024195180854,
     -0.00000003193341274, 0.00000000434545063,-0.00000000060578480,
      0.00000000008612098,-0.00000000001244332, 0.00000000000182256,
     -0.00000000000027007, 0.00000000000004042,-0.00000000000000610,
      0.00000000000000093,-0.00000000000000014, 0.00000000000000002};

   Double_t t,h,y,s,a,alfa,b1,b2,b0;

   if (x == 1) {
      h = pi6;
   } else if (x == -1) {
      h = -pi12;
   } else {
      t = -x;
      if (t <= -2) {
         y = -1/(1+t);
         s = 1;
         b1= log(-t);
         b2= log(1+1/t);
         a = -pi3+hf*(b1*b1-b2*b2);
      } else if (t < -1) {
         y = -1-t;
         s = -1;
         a = log(-t);
         a = -pi6+a*(a+log(1+1/t));
      } else if (t <= -0.5) {
         y = -(1+t)/t;
         s = 1;
         a = log(-t);
         a = -pi6+a*(-hf*a+log(1+t));
      } else if (t < 0) {
         y = -t/(1+t);
         s = -1;
         b1= log(1+t);
         a = hf*b1*b1;
      } else if (t <= 1) {
         y = t;
         s = 1;
         a = 0;
      } else {
         y = 1/t;
         s = -1;
         b1= log(t);
         a = pi6+hf*b1*b1;
      }
      h    = y+y-1;
      alfa = h+h;
      b1   = 0;
      b2   = 0;
      for (Int_t i=19;i>=0;i--){
         b0 = c[i] + alfa*b1-b2;
         b2 = b1;
         b1 = b0;
      }
      h = -(s*(b0-h*b2)+a);
   }
   return h;
}


int G2PPhys::Configure(EMode mode)
{
    if ((mode == kREAD || mode == kTWOWAY) && fConfigured)
        return 0;

    if (G2PProcBase::Configure(mode) != 0)
        return -1;

    ConfDef confs[] = {
        {"run.particle.id", "Particle ID", kINT, &fPID},
        {"run.target.z", "Target Z", kINT, &fZ},
        {"run.target.a", "Target A", kINT, &fA},
        {"run.target.mass", "Target Mass", kDOUBLE, &fTargetMass},
        {"run.hrs.p0", "HRS Momentum", kDOUBLE, &fHRSMomentum},
        {"run.particle.mass", "Beam Particle Mass", kDOUBLE, &fParticleMass},
        {0}
    };

    return ConfigureFromList(confs, mode);
}

int G2PPhys::DefineVariables(EMode mode)
{
    if (mode == kDEFINE && fDefined)
        return 0;

    if (G2PProcBase::DefineVariables(mode) != 0)
        return -1;

    VarDef vars[] = {
        {"react.angle", "Real scattering angle", kDOUBLE, &fTHreact},
        {"react.xs", "Cross section with real kins", kDOUBLE, &fXSreact},
        {"rec.angle", "Rec scattering angle", kDOUBLE, &fTHrec},
        {"rec.xs", "Cross section with rec kins", kDOUBLE, &fXSrec},
        {0}
    };

    return DefineVarsFromList(vars, mode);
}

void G2PPhys::MakePrefix()
{
    const char *base = "phys";

    G2PAppBase::MakePrefix(base);
}

ClassImp(G2PPhys)
