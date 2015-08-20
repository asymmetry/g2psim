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

#ifndef G2P_PHYS_H
#define G2P_PHYS_H

#include "G2PProcBase.hh"

class G2PPhysBase;

class G2PPhys : public G2PProcBase
{
public:
    G2PPhys(const char *name);
    virtual ~G2PPhys();

    virtual int Begin();
    virtual int Process();
    virtual void Clear(Option_t *opt = "");

    // Gets

    // Sets
    void SetPars(double *array, int n);

protected:
    G2PPhys(); // Only for ROOT I/O

    double CalXS(const double *V5lab, const double *V5tr, double &scatangle);
    double TDiLog(double x);
    virtual int Configure(EMode mode = kTWOWAY);
    virtual int DefineVariables(EMode mode = kDEFINE);
    virtual void MakePrefix();

    int fSetting;

    int fPID;
    int fZ, fA; // Define Target
    double fTargetMass;
    double fParticleMass;
    double *fPars;
    int fNPars;

    double fHRSMomentum;

    double fE;
    double fTb, fTa;

    double fXSreact;
    double fTHreact;
    double fXSrec;
    double fTHrec;

    G2PPhysBase *pModel;

private:
    static G2PPhys *pG2PPhys;

    ClassDef(G2PPhys, 1)
};

#endif
