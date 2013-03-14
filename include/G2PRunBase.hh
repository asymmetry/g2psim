#ifndef G2P_RUNBASE_H
#define G2P_RUNBASE_H

#include <vector>

#include "TObject.h"

using namespace std;

class TList;

class G2PRunBase : public TObject
{
public:
    virtual ~G2PRunBase();
    
    enum EStatus { kOK = 0, kNOTINIT, kINITERROR, kERROR };

    void SetHRSAngle(double angle) { fHRSAngle = angle; }
    void SetHRSMomentum(double momentum) { fHRSMomentum = momentum; }
    void SetBeamEnergy(double value) { fBeamEnergy = value; }
    void SetTarget(int Z, int A) { iTargetZ = Z; iTargetA = A; }
    void SetTargetMass(double M) { fTargetMass = M; }
    void SetEnergyLoss(double E) { fEnergyLoss = E; }
    void SetParticlePID(int pid) { iParticlePID = pid; }
    void SetParticleMass(double M0) { fParticleM0 = M0; }
    void SetParticleCharge(double Q) { fParticleQ = Q; }

    void SetDebug(int n) { fDebug = n; }

    virtual int Init();
    virtual int Begin();
    virtual int Process();
    virtual int End() { return 0; }
    virtual void Clear();

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

    static G2PRunBase* GetInstance() { return pG2PRunBase; }

protected:
    G2PRunBase();

    EStatus fStatus;
    int fDebug;

    double fBeamEnergy;
    double fHRSAngle;
    double fHRSMomentum;

    int iTargetZ, iTargetA;
    double fTargetMass;
    double fEnergyLoss;

    int iParticlePID;
    double fParticleM0, fParticleQ;

    TList* fProcs;
    vector<vector<const char*> > fProcReqs;

private:
    static G2PRunBase* pG2PRunBase;

    ClassDef(G2PRunBase, 1)
};

#endif
