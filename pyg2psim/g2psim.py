#!/usr/bin/env python

from os.path import join, dirname, realpath
from time import time
from array import array

import ROOT
ROOT.gSystem.Load(join(dirname(realpath(__file__)), 'libG2PSim.so'))

from config_list import defaults
from config import expand

def run(**kwds):
    if not 'base' in kwds:
        kwds['base'] = defaults
    c = expand(kwds)

    kDEG = 3.14159265358979323846 / 180.0

    run = ROOT.G2PRun()
    n = c['n']
    run.SetNEvent(n)
    if 'conf' in c:
        run.SetConfigFile(c['conf'])
    if 'debug' in c:
        run.SetDebugLevel(c['debug'])
    if 'seed' in c:
        run.SetSeed(c['seed'])

    cr = c['run']
    run.SetBeamEnergy(cr['e0'])
    run.SetHRSAngle(cr['angle'] * kDEG)
    run.SetHRSMomentum(cr['p0'])
    hrs = cr['hrs']
    run.SetTargetType(cr['target'])
    run.SetFieldType(cr['field'])
    if 'particle' in cr:
        run.SetParticle(cr['particle'])

    if 'target' in c and isinstance(c['target'], dict):
        ct = c['target']
        if 'z' in ct and 'a' in ct:
            run.SetTarget(ct['z'], ct['a'])
        if 'mass' in ct:
            run.SetTargetMass(ct['mass'])
        if 'pf' in ct:
            run.SetTargetPF(ct['pf'])
        if 'offset' in ct:
            cto = ct['offset']
            run.SetTargetOffset(cto[0], cto[1], cto[2])

    if 'field' in c and isinstance(c['target'], dict):
        cf = c['field']
        if 'ratio' in cf:
            run.SetFieldRatio(cf['ratio'])
        if 'origin' in cf:
            cfo = cf['origin']
            run.SetFieldOrigin(cfo[0], cfo[1], cfo[2])
        if 'angle' in cf:
            cfa = cf['angle']
            run.SetFieldAngle(cfa[0] * kDEG, cfa[1] * kDEG, cfa[2] * kDEG)

    gun = ROOT.G2PEvGen()
    if 'generator' in c and isinstance(c['generator'], dict):
        cg = c['generator']
        if 'beam' in cg:
            cgb = cg['beam']
            if 'pos' in cgb and 'z' in cgb:
                cgbp = cgb['pos']
                cgbz = cgb['z']
                gun.SetBeamPos(cgbp[0], cgbp[1], cgbz)
            if 'angle' in cgb:
                cgba = cgb['angle']
                gun.SetTiltAngle(cgba[0], cgba[1])
        if 'raster' in cg:
            cgr = cg['raster']
            if 'fast' in cgr:
                cgrf = cgr['fast']
                gun.SetFastRasterSize(cgrf[0], cgrf[1])
            if 'slow' in cgr:
                cgrs = cgr['slow']
                gun.SetSlowRasterSize(cgrs[0], cgrs[1])
        if 'react' in cg:
            cgr = cg['react']
            if 't' in cgr:
                cgrt = cgr['t']
                gun.SetTargetTh(cgrt[0], cgrt[1])
            if 'p' in cgr:
                cgrp = cgr['p']
                gun.SetTargetPh(cgrp[0], cgrp[1])
            if 'z' in cgr:
                cgrz = cgr['z']
                gun.SetReactZ(cgrz[0], cgrz[1])
            if 'd' in cgr:
                cgrd = cgr['d']
                if isinstance(cgrd, list):
                    gun.SetDelta(cgrd[0], cgrd[1])
                elif isinstance(cgrd, str):
                    gun.SetDelta(cgrd)
    ROOT.gG2PApps.Add(gun)

    fwd = ROOT.G2PFwdHRS(hrs)
    if 'forward' in c and isinstance(c['forward'], dict):
        cf = c['forward']
        if 'sieve' in cf:
            fwd.SetSieve(cf['sieve'])
        if 'vdc' in cf:
            cfv = cf['vdc']
            fwd.SetVDCRes(cfv[0], cfv[1], cfv[2], cfv[3])
    ROOT.gG2PApps.Add(fwd)

    if 'bpm' in c and isinstance(c['bpm'], dict):
        bpm = ROOT.G2PBPM()
        cb = c['bpm']
        if 'res' in cb:
            cbr = cb['res']
            bpm.SetBPMRes(cbr[0], cbr[1])
        ROOT.gG2PApps.Add(bpm)

    if 'backward' in c and isinstance(c['backward'], dict):
        bwd = ROOT.G2PBwdHRS(hrs)
        cb = c['backward']
        if 'xpars' in cb:
            cbx = array('d', cb['xpars'])
            bwd.SetParsX(cbx)
        if 'ypars' in cb:
            cby = array('d', cb['ypars'])
            bwd.SetParsY(cby)
        if 'z' in cb:
            bwd.SetRecZ(cb['z'])
        ROOT.gG2PApps.Add(bwd)

    if 'phys' in c and isinstance(c['phys'], dict):
        cp = c['phys']
        phys = ROOT.G2PPhys(cp['model'])
        if 'pars' in cp:
            mpars = array('d', cp['pars'])
            phys.SetPars(mpars, len(cp['pars']))
        ROOT.gG2PApps.Add(phys)

    fname = 'test.root'
    if 'filename' in c:
        fname = c['filename']
    sim = ROOT.G2PSim(fname)

    begin = time()
    sim.Run()
    end = time()
    print('Average calculation time for one event: {0:.3f} ms'.format((end-begin) / n * 1000))

def optrun(**kwds):
    from config_list import defaults
    if not 'base' in kwds:
        kwds['base'] = defaults
    c = expand(kwds)

    kDEG = 3.14159265358979323846 / 180.0

    run = ROOT.G2PRun()
    n = c['n']
    run.SetNEvent(n)
    if 'conf' in c:
        run.SetConfigFile(c['conf'])
    else:
        run.SetConfigFile('optics.cfg')
    if 'debug' in c:
        run.SetDebugLevel(c['debug'])
    if 'seed' in c:
        run.SetSeed(c['seed'])

    cr = c['run']
    run.SetBeamEnergy(cr['e0'])
    run.SetHRSAngle(cr['angle'] * kDEG)
    run.SetHRSMomentum(cr['p0'])
    run.SetTargetType(cr['target'])
    run.SetFieldType(cr['field'])
    if 'particle' in cr:
        run.SetParticle(cr['particle'])

    co = c['optics']
    optics = ROOT.G2POptics(co['filename'])
    if 'foilz' in co:
        optics.SetFoilZ(co['foilz'])
    if 'bpmz' in co:
        optics.SetBPMZ(co['bpmz'])
    if 'eloss' in co:
        coe = array('d', co['eloss'])
        optics.SetEnergyLoss(coe, len(coe))
    ROOT.gG2PApps.Add(optics)

    fname = 'test.root'
    if 'filename' in c:
        fname = c['filename']
    sim = ROOT.G2PSim(fname)

    begin = time()
    sim.Run()
    end = time()
    print('Average calculation time for one event: {0:.3f} ms'.format((end-begin) / n * 1000))

def recrun(**kwds):
    pass
