// -*- C++ -*-

/* class G2PRun
 * Run manager for g2p simulation.
 * Parse the configuration file and store all run parameters.
 * It will allocate G2PFwdProc and G2PBwdProc automatically.
 */

// History:
//   Jan 2013, C. Gu, First public version.
//   May 2013, C. Gu, Add G2PProcBase classes, G2PRunBase is more general.
//   Sep 2013, C. Gu, Combine G2PRunBase and G2PRun to be the new run manager.
//

#ifndef G2P_RUN_H
#define G2P_RUN_H

#include <cstring>
#include <map>

#include "TObject.h"

#include "G2PVarDef.hh"

using namespace std;

class G2PRun : public TObject {
public:
    G2PRun();
    virtual ~G2PRun();

    virtual int Init();
    virtual int Begin();
    virtual int End();
    virtual void Clear();

    int GetConfig(const ConfDef* item, const char* prefix);
    int SetConfig(const ConfDef* item, const char* prefix);

    //Gets
    int GetDebugLevel();
    double GetHRSAngle();
    double GetHRSMomentum();
    int GetParticleID();
    double GetParticleMass();
    double GetParticleCharge();
    double GetBeamEnergy();
    int GetTargetZ();
    int GetTargetA();
    double GetTargetMass();
    double GetFieldRatio();

    //Sets
    void SetConfigFile(const char* file);
    void SetDebugLevel(int n);
    void SetSeed(unsigned n);
    void SetHRSAngle(double angle);
    void SetHRSMomentum(double P0);
    void SetParticleID(int pid);
    void SetParticleMass(double M0);
    void SetParticleCharge(double Q);
    void SetBeamEnergy(double value);
    void SetTarget(int Z, int A);
    void SetTargetMass(double M);
    void SetFieldRatio(double ratio);
    void SetSieve();

protected:
    const char* fConfigFile;

    map<string, double> fConfig;

private:
    static G2PRun* pG2PRun;

    ClassDef(G2PRun, 1)
};

#endif
