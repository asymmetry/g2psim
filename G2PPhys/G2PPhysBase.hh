// -*- C++ -*-

/* class G2PPhysBase
 * Abstract base class of G2PPhys classes.
 * It provides interface functions.
 *
 * Elastic cross section models (G2PPhysEl):
 * * All: Form factors from K. C. Stansfield et al., Phys. Rev. C, 3(1971)1448
 * * H1 : Form factors from S. Venkat et al., Phys. Rev. C, 83(2011)015203 (global fit, with TPE correction)
 *                          J. Arrington et al., Phys. Rev. C 76(2007)035201 (low Q2, with/without TPE correction)
 * * He4: Charge and magnetization densities from De Jager, At. Data Nucl. Data Tables, 14(1974)
 * * C12: Charge and magnetization densities from De Jager, At. Data Nucl. Data Tables, 14(1974)
 *        Charge distribution from L. S. Cardman et al., Phys. Lett. B, 91(1970)203
 * * N14: Charge and magnetization densities from De Jager, At. Data Nucl. Data Tables, 14(1974)
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
 * Unit of cross section is ub/MeV-sr.
 */

// History:
//   Mar 2013, C. Gu, First public version.
//   Nov 2016, C. Gu, Rewrite parameter parser.
//

#ifndef G2P_PHYSBASE_H
#define G2P_PHYSBASE_H

#include <vector>

using namespace std;

class G2PPhysBase
{
public:
    G2PPhysBase();
    virtual ~G2PPhysBase();

    void SetTarget(int Z, int A);
    void SetTargetMass(double value);
    void SetParticle(int pid);

    virtual void SetPar(int id, double value) = 0;
    virtual double GetXS(double Ei, double Ef, double theta) = 0;

protected:
    void SetTargetMass();

    int fZ, fA; // Define Target
    double fTargetMass;

    int fPID; // Define particle

    vector<double> fPars;
};

#endif
