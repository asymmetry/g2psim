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
//

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <map>
#include <set>
#include <sstream>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TTree.h"

#include "libconfig.h"

#include "G2PAppBase.hh"
#include "G2PAppList.hh"
#include "G2PDrift.hh"
#include "G2PField.hh"
#include "G2PGeoBase.hh"
#include "G2PRand.hh"
#include "G2PSieve.hh"
#include "G2PTarget.hh"
#include "G2PVar.hh"
#include "G2PVarDef.hh"
#include "G2PVarList.hh"

#include "G2PRun.hh"

using namespace std;

static const double e = 1.60217656535e-19;
static const double kDEG = 3.14159265358979323846 / 180.0;
static const double kU = 0.931494028;

G2PAppList *gG2PApps = NULL;
G2PAppList *gG2PGeos = NULL;
G2PRun *gG2PRun = NULL;
G2PVarList *gG2PVars = NULL;

G2PRand *G2PRun::pRand = G2PRand::GetInstance();
G2PRun *G2PRun::pG2PRun = NULL;

G2PRun::G2PRun() : fConfigFile("sim.cfg")
{
    if (pG2PRun) {
        Error("G2PRun()", "Only one instance of G2PRun allowed.");
        MakeZombie();
        return;
    }

    pG2PRun = this;

    fConfig.clear();
    fConfig["run.n"] = 50000;
    fConfig["run.debuglevel"] = 1;
    fConfig["run.seed"] = 0;
    fConfig["run.e0"] = 2.2535;
    fConfig["run.particle.id"] = 11;
    fConfig["run.particle.mass"] = 0.51099892811e-3;
    fConfig["run.particle.charge"] = -1 * e;
    fConfig["run.target.type"] = 10;
    fConfig["run.target.z"] = 1;
    fConfig["run.target.a"] = 1;
    fConfig["run.target.mass"] = 1.00794 * kU;
    fConfig["run.target.pf"] = 0.55;
    fConfig["run.hrs.angle"] = 5.767 * kDEG;
    fConfig["run.hrs.p0"] = 2.24949;
    fConfig["field.type"] = 0;
    fConfig["field.ratio"] = 0.0;
    fConfigIsSet.clear();

    gG2PApps = new G2PAppList();
    gG2PGeos = new G2PAppList();
    gG2PVars = new G2PVarList();

    gG2PRun = this;
}

G2PRun::~G2PRun()
{
    if (pG2PRun == this)
        pG2PRun = NULL;

    TIter next(gG2PApps);

    while (G2PAppBase *aobj = static_cast<G2PAppBase *>(next())) {
        if (aobj != NULL) {
            gG2PApps->Remove(aobj);
            delete aobj;
        }
    }

    if (gG2PApps) {
        delete gG2PApps;
        gG2PApps = NULL;
    }

    TIter geo_iter(gG2PGeos);

    while (G2PGeoBase *aobj = static_cast<G2PGeoBase *>(geo_iter())) {
        if (aobj != NULL) {
            gG2PGeos->Remove(aobj);

            if (aobj->GetNUsed() == 0)
                delete aobj;
            else
                aobj->SetNUsed(aobj->GetNUsed() - 1);
        }
    }

    if (gG2PGeos) {
        delete gG2PGeos;
        gG2PGeos = NULL;
    }

    TIter var_iter(gG2PVars);

    while (G2PVar *aobj = static_cast<G2PVar *>(var_iter())) {
        if (aobj != NULL) {
            gG2PVars->Remove(aobj);
            delete aobj;
        }
    }

    if (gG2PVars) {
        delete gG2PVars;
        gG2PVars = NULL;
    }

    gG2PRun = NULL;
}

int G2PRun::Begin()
{
    //static const char* const here = "Begin()";

    if (fConfigFile && *fConfigFile) {
        if (ParseConfigFile() != 0)
            return -1;
    }

    // Set random seed
    G2PRun::StaticSetSeed((int)(fConfig["run.seed"]));

    // Set target
    if (gG2PApps->Find("G2PTarget") == NULL) {
        G2PTarget *target = new G2PTarget();
        gG2PApps->Add(target);
    }

    // Set sieve
    if (gG2PApps->Find("G2PSieve") == NULL) {
        G2PSieve *sieve = new G2PSieve();
        gG2PApps->Add(sieve);
    }

    // Set field
    if (gG2PApps->Find("G2PField") == NULL) {
        G2PField *field = new G2PField();
        gG2PApps->Add(field);
    }

    // Set Drift
    if (gG2PApps->Find("G2PDrift") == NULL) {
        G2PDrift *drift = new G2PDrift();
        gG2PApps->Add(drift);
    }

    return 0;
}

int G2PRun::End()
{
    if (fConfig["run.debuglevel"] > 1)
        Print();

    return 0;
}

void G2PRun::Print(Option_t *opt) const
{
    static const char *const here = "Print()";

    Info(here, "The configuration is :");

    map<string, double>::const_iterator it = fConfig.begin();

    while (it != fConfig.end()) {
        ostringstream ostr; // String stream has a better output format
        ostr << it->first << " = " << it->second;
        Info(here, "%s", ostr.str().c_str());
        it++;
    }
}

int G2PRun::GetConfig(const ConfDef *item, const char *prefix)
{
    // Get value of item

    static const char *const here = "GetConfig()";

    if (!item) {
        Error(here, "Bad variable.");
        return -1;
    }

    if (item->var) {
        double dval;
        string keystr(prefix);
        keystr.append(item->name);
        string key(item->name);

        if (fConfig.count(key) > 0)
            dval = fConfig[key];
        else if (fConfig.count(keystr) > 0)
            dval = fConfig[keystr];
        else
            return -1;

        switch (item->type) {
        case kBOOL:
            *((bool *) item->var) = (bool) dval;
            break;

        case kCHAR:
            *((char *) item->var) = (char) dval;
            break;

        case kINT:
            *((int *) item->var) = (int) dval;
            break;

        case kSHORT:
            *((short *) item->var) = (short) dval;
            break;

        case kLONG:
            *((long *) item->var) = (long) dval;
            break;

        case kFLOAT:
            *((float *) item->var) = (float) dval;
            break;

        case kDOUBLE:
            *((double *) item->var) = dval;
            break;
        }
    }

    return 0;
}

int G2PRun::SetConfig(const ConfDef *item, const char *prefix)
{
    // Set value of item

    static const char *const here = "SetConfig()";

    if (!item) {
        Error(here, "Bad variable.");
        return -1;
    }

    if (item->var) {
        double dval = 0;

        switch (item->type) {
        case kBOOL:
            dval = (double)(*((bool *) item->var));
            break;

        case kCHAR:
            dval = (double)(*((char *) item->var));
            break;

        case kINT:
            dval = (double)(*((int *) item->var));
            break;

        case kSHORT:
            dval = (double)(*((short *) item->var));
            break;

        case kLONG:
            dval = (double)(*((long *) item->var));
            break;

        case kFLOAT:
            dval = (double)(*((float *) item->var));
            break;

        case kDOUBLE:
            dval = *((double *) item->var);
            break;
        }

        string keystr(prefix);
        keystr.append(item->name);
        string key(item->name);

        if (fConfig.count(keystr) > 0)
            fConfig[keystr] = dval;
        else if (fConfig.count(key) > 0)
            fConfig[key] = dval;
        else
            fConfig[keystr] = dval;
    }

    return 0;
}

int G2PRun::GetConfigList(ConfDef *&conf) const
{
    map<string, double>::size_type n = fConfig.size();

    ConfDef *result = new ConfDef[n];
    int i = 0;

    for (map<string, double>::const_iterator it = fConfig.begin(); it != fConfig.end(); ++it) {
        result[i].desc = it->first.c_str();
        result[i].name = it->first.c_str();
        result[i].type = kDOUBLE;
        result[i].var = (void *) new double(it->second);
        i++;
    }

    conf = result;

    return ((int) n);
}

void G2PRun::SetConfigFile(const char *file)
{
    fConfigFile = file;
}

void G2PRun::SetNEvent(int n)
{
    fConfig["run.n"] = (double) n;

    fConfigIsSet.insert("run.n");
}

void G2PRun::SetDebugLevel(int n)
{
    fConfig["run.debuglevel"] = (double) n;

    fConfigIsSet.insert("run.debuglevel");
}

void G2PRun::SetSeed(unsigned n)
{
    fConfig["run.seed"] = (double) n;

    fConfigIsSet.insert("run.seed");
}

void G2PRun::SetRunType(const char *type)
{
    SetTargetType(type);
}

void G2PRun::SetTargetType(const char *type)
{
    map<string, int> tempmap;
    tempmap["prod"] = 10;
    tempmap["production"] = 10;

    tempmap["sprod"] = 11;
    tempmap["short"] = 11;
    tempmap["shortproduction"] = 11;

    tempmap["optics"] = 20;
    tempmap["optics20"] = 20;
    tempmap["optics_C40"] = 20;

    tempmap["optics21"] = 21;
    tempmap["optics_C40He"] = 21;
    tempmap["dilu"] = 21;
    tempmap["dilution"] = 21;
    tempmap["dilu_C40"] = 21;
    tempmap["dilution_C40"] = 21;

    tempmap["optics22"] = 22;
    tempmap["optics_C125"] = 22;

    tempmap["optics23"] = 23;
    tempmap["optics_C125He"] = 23;
    tempmap["dilu_C125"] = 23;
    tempmap["dilution_C125"] = 23;

    tempmap["dummy"] = 30;
    tempmap["dummy30"] = 30;

    tempmap["dummy31"] = 31;
    tempmap["dummy_nocap"] = 31;

    tempmap["empty"] = 40;

    if (tempmap.count(type) > 0)
        SetTargetType(tempmap[type]);
    else
        SetTargetType(10);
}

void G2PRun::SetTargetType(int id)
{
    fConfig["run.target.type"] = (double) id;
    fConfigIsSet.insert("run.target.type");

    switch (id) {
    case 10:
    case 11:
        if (fConfigIsSet.count("run.target.z") == 0) {
            SetTarget(1, 1);
            SetTargetMass(1.00794 * kU);
        }

        if (fConfigIsSet.count("run.target.pf") == 0)
            SetTargetPF(0.55);

        break;

    case 20:
    case 21:
    case 22:
    case 23:
        if (fConfigIsSet.count("run.target.z") == 0) {
            SetTarget(6, 12);
            SetTargetMass(12.0107 * kU);
        }

        break;

    case 30:
    case 31:
    case 40:
        if (fConfigIsSet.count("run.target.z") == 0) {
            SetTarget(2, 4);
            SetTargetMass(4.002602 * kU);
        }

        break;
    }
}

void G2PRun::SetFieldType(const char *type)
{
    map<string, int> tempmap;
    tempmap["none"] = 0;
    tempmap["nofield"] = 0;
    tempmap["trans"] = 10;
    tempmap["transverse"] = 10;
    tempmap["long"] = 11;
    tempmap["longitudinal"] = 11;
    tempmap["gep"] = 20;

    if (tempmap.count(type) > 0)
        SetFieldType(tempmap[type]);
    else
        SetFieldType(0);
}

void G2PRun::SetFieldType(int id)
{
    fConfig["field.type"] = id;
    fConfigIsSet.insert("field.type");

    switch (id) {
    case 0:
        SetFieldRatio(0.0);
        break;

    case 10:
        if (fConfigIsSet.count("field.ratio") == 0)
            SetFieldRatio(2.5);

        if (fConfigIsSet.count("field.angle.alpha") == 0)
            SetFieldAngle(90 * kDEG, 90 * kDEG, -90 * kDEG);

        break;

    case 11:
        if (fConfigIsSet.count("field.ratio") == 0)
            SetFieldRatio(5.0);

        break;

    case 20:
        if (fConfigIsSet.count("field.ratio") == 0)
            SetFieldRatio(5.0);

        if (fConfigIsSet.count("field.angle.alpha") == 0)
            SetFieldAngle(90 * kDEG, 5.6 * kDEG, -90 * kDEG);

        break;
    }
}

void G2PRun::SetParticle(int id)
{
    static const char *const here = "SetParticle()";

    fConfig["run.particle.id"] = (double) id;

    switch (id) {
    case 11: // e-
        fConfig["run.particle.mass"] = 0.5109989281e-3;
        fConfig["run.particle.charge"] = -1 * e;
        break;

    case -11: // e+
        fConfig["run.particle.mass"] = 0.5109989281e-3;
        fConfig["run.particle.charge"] = e;
        break;

    case -211: // pi-
        fConfig["run.particle.mass"] = 139.5701835e-3;
        fConfig["run.particle.charge"] = -1 * e;
        break;

    case 211: // pi+
        fConfig["run.particle.mass"] = 139.5701835e-3;
        fConfig["run.particle.charge"] = e;
        break;

    default: // e-
        Warning(here, "Unknown pid, use electrons.");
        fConfig["run.particle.mass"] = 0.5109989281e-3;
        fConfig["run.particle.charge"] = -1 * e;
        break;
    }

    fConfigIsSet.insert("run.particle.id");
    fConfigIsSet.insert("run.particle.mass");
    fConfigIsSet.insert("run.particle.charge");
}

void G2PRun::SetBeamEnergy(double E)
{
    fConfig["run.e0"] = E;

    fConfigIsSet.insert("run.e0");
}

void G2PRun::SetHRSAngle(double angle)
{
    fConfig["run.hrs.angle"] = angle;

    fConfigIsSet.insert("run.hrs.angle");
}

void G2PRun::SetHRSMomentum(double P0)
{
    fConfig["run.hrs.p0"] = P0;

    fConfigIsSet.insert("run.hrs.p0");
}

void G2PRun::SetTarget(int Z, int A)
{
    fConfig["run.target.z"] = (double) Z;
    fConfig["run.target.a"] = (double) A;

    fConfigIsSet.insert("run.target.z");
    fConfigIsSet.insert("run.target.a");
}

void G2PRun::SetTargetMass(double M)
{
    fConfig["run.target.mass"] = M;

    fConfigIsSet.insert("run.target.mass");
}

void G2PRun::SetTargetOffset(double x, double y, double z)
{
    fConfig["run.target.offset.x"] = (double) x;
    fConfig["run.target.offset.y"] = (double) y;
    fConfig["run.target.offset.z"] = (double) z;

    fConfigIsSet.insert("run.target.offset.x");
    fConfigIsSet.insert("run.target.offset.y");
    fConfigIsSet.insert("run.target.offset.z");
}

void G2PRun::SetTargetPF(double pf)
{
    fConfig["run.target.pf"] = pf;

    fConfigIsSet.insert("run.target.pf");
}

void G2PRun::SetFieldRatio(double ratio)
{
    fConfig["field.ratio"] = ratio;

    fConfigIsSet.insert("field.ratio");
}

void G2PRun::SetFieldAngle(double alpha, double beta, double gamma)
{
    fConfig["field.angle.alpha"] = (double) alpha;
    fConfig["field.angle.beta"] = (double) beta;
    fConfig["field.angle.gamma"] = (double) gamma;

    fConfigIsSet.insert("field.angle.alpha");
    fConfigIsSet.insert("field.angle.beta");
    fConfigIsSet.insert("field.angle.gamma");
}

int G2PRun::ParseConfigFile()
{
    static const char *const here = "ParseConfigFile()";

    config_t cfg;
    config_init(&cfg);

    if (!config_read_file(&cfg, fConfigFile)) {
        Error(here, "Parse error at %s:%d - %s", config_error_file(&cfg), config_error_line(&cfg), config_error_text(&cfg));
        config_destroy(&cfg);
        return -1;
    }

    config_setting_t *root = config_root_setting(&cfg);

    if (ParseSetting("", root)) {
        config_destroy(&cfg);
        return -1;
    }

    config_destroy(&cfg);

    return 0;
}

int G2PRun::ParseSetting(const char *prefix, const config_setting_t *setting)
{
    static const char *const here = "ParseConfigFile()"; // NOTICE: this is not a typo

    int type = config_setting_type(setting);
    bool isgood = true;

    switch (type) {
    case CONFIG_TYPE_GROUP:
        for (int i = 0; i < config_setting_length(setting); i++) {
            config_setting_t *subsetting = config_setting_get_elem(setting, i);
            const char *newprefix;

            if (config_setting_name(setting) == NULL)
                newprefix = prefix;
            else
                newprefix = Form("%s%s.", prefix, config_setting_name(setting));

            if (ParseSetting(newprefix, subsetting)) {
                isgood = false;
                break;
            }
        }

        break;

    case CONFIG_TYPE_BOOL:
    case CONFIG_TYPE_INT:
    case CONFIG_TYPE_INT64:
    case CONFIG_TYPE_FLOAT: {
        string name = Form("%s%s", prefix, config_setting_name(setting));
        double value = GetValue(setting);

        if (fConfigIsSet.find(name) == fConfigIsSet.end())
            fConfig[name] = value;

        return 0;
    }

    default:
        Error(here, "Found illegal configuration, check the format.");
        return -1;
    }

    return (isgood ? 0 : -1);
}

double G2PRun::GetValue(const config_setting_t *setting)
{
    double value = 0;

    int type = config_setting_type(setting);

    switch (type) {
    case CONFIG_TYPE_INT:
        value = (double) config_setting_get_int(setting);
        break;

    case CONFIG_TYPE_INT64:
        value = (double) config_setting_get_int64(setting);
        break;

    case CONFIG_TYPE_FLOAT:
        value = config_setting_get_float(setting);
        break;

    case CONFIG_TYPE_BOOL:
        value = (double) config_setting_get_bool(setting);
        break;
    }

    return value;
}

void G2PRun::StaticSetSeed(unsigned n)
{
    pRand->SetSeed(n);
}

ClassImp(G2PRun)
