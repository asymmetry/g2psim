// -*- C++ -*-

/* class G2PRun
 * Run manager for g2p simulation.
 * Parse the configuration file and store all run parameters.
 * It use libconfig, a 3rd party package, to do the parsing.
 */

// History:
//   Jan 2013, C. Gu, First public version.
//   May 2013, C. Gu, Add G2PProcBase classes, G2PRunBase is more general.
//   Sep 2013, C. Gu, Combine G2PRunBase and G2PRun to be the new run manager.
//   Nov 2014, J. Liu, Add target type and field type.
//   Oct 2015, C. Gu, Use G2PRun to handle all parameters.
//

#ifndef G2P_RUN_H
#define G2P_RUN_H

#include <cstring>
#include <map>
#include <set>
#include <sstream>

#include "TObject.h"

#include "libconfig.h"

#include "G2PVarDef.hh"

using namespace std;

class G2PRand;

class G2PRun : public TObject
{
public:
    G2PRun();
    virtual ~G2PRun();

    virtual int Begin();
    virtual int End();

    virtual void Print(Option_t *opt = "") const;

    int GetConfig(const ConfDef *item, const char *prefix);
    int SetConfig(const ConfDef *item, const char *prefix);

    //Gets
    int GetConfigList(ConfDef *&conf)const;

    //Sets
    void SetConfigFile(const char *file);

    void SetNEvent(int n);
    void SetDebugLevel(int n);
    void SetSeed(unsigned n);

    void SetRunType(const char *type);
    void SetTargetType(const char *type);
    void SetTargetType(int id);
    void SetFieldType(const char *type);
    void SetFieldType(int id);

    void SetParticle(int id);
    void SetBeamEnergy(double value);

    void SetHRSAngle(double angle);
    void SetHRSMomentum(double P0);

    void SetTarget(int Z, int A);
    void SetTargetMass(double M);

    void SetTargetOffset(double x, double y, double z);
    void SetTargetPF(double pf);

    void SetFieldRatio(double ratio);
    void SetFieldOrigin(double x, double y, double z);
    void SetFieldAngle(double alpha, double beta, double gamma);

protected:
    int ParseConfigFile();
    int ParseSetting(const char *prefix, const config_setting_t *setting);
    double GetValue(const config_setting_t *setting);

    const char *fConfigFile;

    map<string, double> fConfig;
    set<string> fConfigIsSet;

private:
    // Random number generator
    static G2PRand *pRand;
    static void StaticSetSeed(unsigned n);

    static G2PRun *pG2PRun;

    ClassDef(G2PRun, 1)
};

#endif
