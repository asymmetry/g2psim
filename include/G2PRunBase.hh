#ifndef G2P_RUNBASE_H
#define G2P_RUNBASE_H

#include <cstring>

#include "G2PAppsBase.hh"

class G2PBPM;
class G2PDrift;
class G2PGunBase;
class G2PHRSTrans;
class G2PPhys;
class G2PRecUseDB;
class TTree;

class G2PRunBase : public G2PAppsBase
{
public:
    G2PRunBase();
    virtual ~G2PRunBase();

    void SetHRSAngle(double angle) { fHRSAngle = angle; }
    void SetHRSMomentum(double momentum) { fHRSMomentum = momentum; }
    void SetBeamEnergy(double value) { fBeamEnergy = value; }
    void SetTarget(int Z, int A) { iTargetZ = Z; iTargetA = A; }
    void SetTargetMass(double M) { fTargetMass = M; }
    void SetEnergyLoss(double E) { fEnergyLoss = E; }
    void SetParticlePID(int pid) { iParticlePID = pid; }
    void SetParticleMass(double M0) { fParticleM0 = M0; }
    void SetParticleCharge(double Q) { fParticleQ = Q; }

    virtual EStatus Init();
    virtual int Begin() { return 0; }
    virtual int End() { return 0; }
    virtual void Clear() { }

    virtual int Run() = 0;

    double GetHRSAngle() { return fHRSAngle; }
    double GetHRSMomentum() { return fHRSMomentum; }
    double GetBeamEnergy() { return fBeamEnergy; }
    int GetTargetZ() { return iTargetZ; }
    int GetTargetA() { return iTargetA; }
    double GetTargetMass() { return fTargetMass; }
    double GetEnergyLoss() { return fEnergyLoss; }
    int GetParticlePID() { return iParticlePID; }
    double GetParticleMass() { return fParticleM0; }
    double GetParticleCharge() { return fParticleQ; }

    virtual int RegisterModel();
    virtual int DefineVariables(TTree* t) { return 0; }

    static G2PRunBase* GetInstance() { return pG2PRunBase; }

protected:
    virtual void Project(double x, double y, double z, double z_out, double t, double p, double &xout, double &yout);

    double fHRSAngle;
    double fHRSMomentum;
    double fBeamEnergy;
    int iTargetZ, iTargetA;
    double fTargetMass;
    double fEnergyLoss;
    int iParticlePID;
    double fParticleM0, fParticleQ;

    G2PBPM* pBPM;
    G2PDrift* pDrift;
    G2PGunBase* pGun;
    G2PHRSTrans* pHRS;
    G2PPhys* pPhys;
    G2PRecUseDB* pRecUseDB;
    
private:
    static G2PRunBase* pG2PRunBase;

    ClassDef(G2PRunBase, 1)
};

#endif
