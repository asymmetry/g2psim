#!/usr/bin/env python

import gzip, os, re
from os.path import exists, join
from math import tan
try: import cPickle as pkl
except: import pickle as pkl
import subprocess as sp
import numpy as np

import ROOT

from . import g2psim as sim
from . import config, config_list

def _zdump(value, filename):
    with gzip.open(filename, 'wb') as f:
        pkl.dump(value, f, -1)

def _zload(filename):
    with gzip.open(filename,'rb') as f:
        return pkl.load(f)

def gencut(run, conf, overwrite=False):
    c = _zload(conf)
    cutfile = join(c['cutpath'], 'g2p_{}.root.SieveCut.cut'.format(run))
    if exists(cutfile) and not overwrite:
        print('{} already exists.'.format(cutfile))
        return None

    print('Generating {} ......'.format(cutfile))
    with open(cutfile, 'w') as f:
        for i in range(7):
            subcutfile = join(c['cutpath'], 'g2p_{}.root.SieveCut.{}.{}.cut'.format(run // 10, run % 10, i))
            if exists(subcutfile):
                with open(subcutfile, 'r') as fi:
                    for j in fi:
                        f.write(j)
            else:
                for j in range(7):
                    f.write('fEvtHdr.fRun==0\n')
            f.write('\n')

def genf51(run, conf, overwrite=False):
    c = _zload(conf)
    f51file = join(c['f51path'], 'Sieve.{}.f51'.format(run))
    if exists(f51file) and not overwrite:
        print('{} already exists.'.format(f51file))
        return None
    tree2ascii = c['tree2ascii']
    vardef = c['vardef']
    cutfile = join(c['cutpath'], 'g2p_{}.root.SieveCut.cut'.format(run))
    gcutfile = join(c['cutpath'], 'g2p_{}.root.FullCut.root'.format(run // 10))
    rootfilelist = [join(c['rootpath'], 'g2p_{}.root'.format(x)) for x in c[run]['runlist']]
    if not all(exists(x) for x in [tree2ascii, vardef, cutfile, gcutfile] + rootfilelist):
        print('Some input files do not exist.')
        return None

    print('Generating {} ......'.format(f51file))
    command = [tree2ascii, '-pv']
    command += ['-d', vardef]
    command += ['-O', str(c['idoffset'] * c[run]['id'])]
    command += ['-S', str(c['eventlimit'])]
    command += ['-c', cutfile]
    command += ['-g', gcutfile]
    command += ['-o', f51file]
    command += rootfilelist
    sp.call(command)

def parsef51(run, conf, overwrite=False):
    c = _zload(conf)
    if 'f51pars' in c[run] and not overwrite:
        print('Found old f51 pars, set overwrite to reparse.')
        return None
    f51file = join(c['f51path'], 'Sieve.{}.f51'.format(run))
    if not exists(f51file):
        print('Require {} as input.'.format(f51file))
        return None

    print('Parsing {} ......'.format(f51file))
    data=[]
    with open(f51file, 'r') as f:
        for l in f:
            dl = [float(x) for x in l.strip().split()[5:10]]
            data.append(dl)
    data = np.array(data)
    c[run]['f51event'] = data.shape[0]

    if not 'f51pars' in c[run]:
        c[run]['f51pars'] = {}
    f51pars = c[run]['f51pars']
    dmean = data.mean(axis = 0)
    f51pars['e0'] = round(dmean[0] * 1e-3, 6)
    f51pars['xb'] = round(dmean[1] * 1e-3, 6)
    f51pars['tb'] = round(dmean[2], 6)
    f51pars['yb'] = round(dmean[3] * 1e-3, 6)
    f51pars['pb'] = round(dmean[4], 6)

    _zdump(c, conf)

def simf51(run, conf, overwrite=False):
    c = _zload(conf)
    simfile = join(c['simpath'], 'run_{}.root'.format(run))
    if exists(simfile) and not overwrite:
        print('{} already exists.'.format(simfile))
        return None
    if not 'f51pars' in c[run]:
        print('Require f51 pars.')
        return None

    f51pars = c[run]['f51pars']
    print('Simulating with e0 = {e0:.3e} xb = {xb:.3e} tb = {tb:.3e} yb = {yb:.3e} pb = {pb:.3e} ......'.format(**f51pars))
    pars = {
        'base': c['config'], 'n': 500000, 'filename': simfile,
        'run': {'e0': f51pars['e0'], 'p0': c[run]['p0']},
        'generator': {'beam': {
            'pos': [f51pars['xb'], f51pars['yb']],
            'angle': [f51pars['tb'], f51pars['pb']]}}
    }
    sim.run(**pars)

def calholeang(run, conf, overwrite=False):
    c = _zload(conf)
    if 'holeang' in c[run] and not overwrite:
        print('Found old hole angles, set overwrite to recalculate.')
        return None
    simfile = join(c['simpath'], 'run_{}.root'.format(run))
    if not exists(simfile):
        print('Require {} as input.'.format(simfile))
        return None

    print('Calculating angle of each sieve hole in {} ......'.format(run))
    if not 'holeang' in c[run]:
        c[run]['holeang'] = [None] * 49
    simconf = config.expand(getattr(config_list, c['config']))
    cgrt = simconf['generator']['react']['t']
    cgrp = simconf['generator']['react']['p']
    t = ROOT.TChain('T')
    t.Add(simfile)
    for i in range(49):
        hgrt = ROOT.TH1D('hgrt' + str(i), 'hgrt' + str(i), 200, cgrt[0] - 0.01, cgrt[1] + 0.01)
        hgrp = ROOT.TH1D('hgrp' + str(i), 'hgrp' + str(i), 100, cgrp[0] - 0.01, cgrp[1] + 0.01)
        t.Project('hgrt' + str(i), 'gen.react.t', 'fwd.id.hole==' + str(i))
        t.Project('hgrp' + str(i), 'gen.react.p', 'fwd.id.hole==' + str(i))
        ct = hgrt.GetBinCenter(hgrt.GetMaximumBin())
        cp = hgrp.GetBinCenter(hgrp.GetMaximumBin())
        c[run]['holeang'][i] = [round(ct, 5), round(cp, 5)]

    _zdump(c, conf)

def simf51hole(run, conf, holeid, overwrite=False):
    c = _zload(conf)
    simpath = join(c['simpath'], 'temp', str(run))
    if not exists(simpath):
        os.makedirs(simpath)
    simfile = join(simpath, 'run_{0}_{1:02d}.root'.format(run, holeid))
    if exists(simfile) and not overwrite:
        print('{} already exists.'.format(simfile))
        return None
    if not 'f51pars' in c[run] and not 'holeang' in c[run]:
        print('Require f51 pars and hole angles.')
        return None

    ct, cp = c[run]['holeang'][holeid]
    f51pars = c[run]['f51pars']
    print('Simulating around hole {0:d} = ({1:.2e}, {2:.2e}) ......'.format(holeid, ct, cp))
    pars = {
        'base': c['config'], 'n': 15000,
        'debug': 0, 'seed': 0, 'filename': simfile,
        'run': {'e0': f51pars['e0'], 'p0': c[run]['p0']},
        'generator': {
            'beam': {
            'pos': [f51pars['xb'], f51pars['yb']],
            'angle': [f51pars['tb'], f51pars['pb']]
            },
            'react': {
                't': [ct - 3.0e-3, ct + 3.0e-3],
                'p': [cp - 3.0e-3, cp + 3.0e-3]
            }
        },
        'bpm': 0, 'backward': 0, 'phys': 0
    }
    sim.run(**pars)

def caleloss(run, conf, overwrite=False):
    c = _zload(conf)
    simpath = join(c['simpath'], 'temp', str(run))
    if 'eloss' in c[run] and not overwrite:
        print('Found old hole energy loss, set overwrite to recalculate.')
        return None

    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    c1 = ROOT.TCanvas('c1','c1', 300 * 7, 200 * 7)
    c1.Divide(7, 7)
    if not 'eloss' in c[run]:
        c[run]['eloss'] = [0] * 49
    he = [None] * 49
    for i in range(49):
        c1.cd((7 - i % 7) * 7 - i // 7)
        print('Calculating energy loss of hole {} ......'.format(i))
        simfile = join(simpath, 'run_{0}_{1:02d}.root'.format(run, i))
        while True:
            simf51hole(run, conf, i, overwrite)
            te = ROOT.TChain('T')
            te.Add(simfile)
            he[i] = ROOT.TH1D('eloss' + str(i), 'eloss' + str(i), 77, 0.0003, 0.008)
            nhe = te.Draw('phys.eloss.b+phys.eloss.a>>eloss{}'.format(i), 'fwd.id.hole==' + str(i))
            if nhe <= 0:
                os.remove(simfile)
                os.mknod('check_{0}_{1:02d}'.format(run, i))
                break
            r = he[i].Fit('landau', 'QS')
            if int(r) == 0:
                eloss = r.Parameter(1)
                c[run]['eloss'][i] = round(eloss, 4)
                print('The energy loss of hole {0:d} is {1:.1f} MeV.'.format(i, eloss * 1000))
                he[i].SetTitle('Eloss = {0:.1f} MeV'.format(eloss * 1000))
                c1.Update()
                break
            else:
                print('The fit failed, redo simulation ......')
                os.remove(simfile)
    c1.Print(join(simpath, 'eloss_{}.png'.format(run)))
    ROOT.gROOT.SetBatch(ROOT.kFALSE)
    _zdump(c, conf)

def simoptics(run, conf, overwrite=False):
    c = _zload(conf)
    simfile = join(c['simpath'], 'ref_{}.root'.format(run))
    if exists(simfile) and not overwrite:
        print('{} already exists.'.format(simfile))
        return None
    f51file = join(c['f51path'], 'Sieve.{}.f51'.format(run))
    if not exists(f51file):
        print('Require {} as input.'.format(f51file))
        return None

    print('Simulating with {} ......'.format(f51file))
    pars = {
        'base': c['config'], 'n': c[run]['f51event'],
        'filename': simfile, 'run': {'p0': c[run]['p0']},
        'optics': {'data': f51file, 'eloss': c[run]['eloss']}
    }
    sim.optrun(**pars)

def getnewf51(run, conf, overwrite=False):
    c = _zload(conf)
    newf51file = join(c['f51path'], 'Sieve.{}.f51.new'.format(run))
    if exists(newf51file) and not overwrite:
        print('{} already exists.'.format(newf51file))
        return None
    simfile = join(c['simpath'], 'ref_{}.root'.format(run))
    if not exists(simfile):
        print('Require {} as input.'.format(simfile))
        return None

    scandef = c['scandef']
    print('Generating {} ......'.format(newf51file))
    t = ROOT.TChain('T')
    t.Add(simfile)
    N = t.GetEntries()
    with open(newf51file, 'w') as f:
        for i in range(N):
            t.GetEntry(i)
            for j in scandef:
                d = getattr(t, j[0])
                if j[2] == 1:
                    d = tan(d)
                if j[1] == 'I':
                    f.write('{0:<5d}  '.format(int(d)))
                else:
                    f.write('{0:13.6e}  '.format(d))
            f.write('\n')

def getfullf51(conf, overwrite=False):
    c = _zload(conf)
    fullf51file = join(c['f51path'], 'Sieve.full.f51.new')
    if exists(fullf51file) and not overwrite:
        print('{} already exists.'.format(fullf51file))
        return None
    f51list = [join(c['f51path'], 'Sieve.{}.f51.new'.format(x)) for x in c['fulllist']]
    if not all(exists(x) for x in f51list):
        print('Require ' + ',\n'.join(f51list) + ' as input.')
        return None

    print('Generating {} ......'.format(fullf51file))
    with open(fullf51file, 'w') as f:
        for i in f51list:
            with open(i, 'r') as fin:
                for j in fin:
                    f.write(j)

def gendata(run, conf, overwrite=False):
    c = _zload(conf)
    datafile = join(c['simpath'], 'g2p_{}.dat'.format(run))
    if exists(datafile) and not overwrite:
        print('{} already exists.'.format(datafile))
        return None
    tree2ascii = c['tree2ascii']
    vardef = c['vardef']
    cutfile = join(c['cutpath'], 'g2p_{}.root.SieveCut.cut'.format(run))
    dppat = re.compile(r'(\(L\.tr\.tg_dp>[-.0-9]* && L\.tr\.tg_dp<[-.0-9]*\) && )')
    hfppat = re.compile(r'(hfpcut_L_[0-9]_[0-9]_[0-9] && )')
    with open(cutfile, 'r') as f:
        with open('tempcut', 'w') as fout:
            for i in f:
                j = dppat.sub(r'', i)
                j = hfppat.sub(r'', j)
                fout.write(j)

    gcutfile = join(c['cutpath'], 'g2p_{}.root.FullCut.root'.format(run // 10))
    rootfilelist = [join(c['rootpath'], 'g2p_{}.root'.format(x)) for x in c[run]['runlist']]
    if not all(exists(x) for x in [tree2ascii, vardef, cutfile, gcutfile] + rootfilelist):
        print('Some input files do not exist.')
        return None

    print('Generating {} ......'.format(datafile))
    command = [tree2ascii, '-pv']
    command += ['-d', vardef]
    command += ['-c', 'tempcut']
    command += ['-g', gcutfile]
    command += ['-o', datafile]
    command += rootfilelist
    sp.call(command)

    if exists('tempcut'):
        os.remove('tempcut')

def recwithdb(run, conf, overwrite=False):
    c = _zload(conf)
    recfile = join(c['simpath'], 'g2p_{}.root'.format(run))
    if exists(recfile) and not overwrite:
        print('{} already exists.'.format(recfile))
        return None
    database = join(c['fitpath'], c['fitdir'], c['database'])
    if not exists(database):
        print('Require {} as input.'.format(databse))
        return None
    if not 'f51pars' in c[run]:
        print('Require f51 pars.')
        return None
    datafile = join(c['simpath'], 'g2p_{}.dat'.format(run))
    if not exists(datafile):
        print('Require {} as input.'.format(datafile))
        return None

    n = 0
    with open(datafile) as f:
        for l in f:
            n += 1
    print('Simulating with {} ......'.format(datafile))
    pars = {
        'base': c['config'], 'n': n, 'filename': recfile,
        'run': {'e0': c[run]['f51pars']['e0'], 'p0': c[run]['p0']},
        'dbrec': {'data': datafile, 'database': database}
    }
    sim.recrun(**pars)

def compare(run, conf, overwrite=False):
    c = _zload(conf)
    tplot = '{}_theta.png'.format(run)
    pplot = '{}_phi.png'.format(run)
    dplot = '{}_delta.png'.format(run)
    if all(exists(x) for x in [tplot, pplot, dplot]) and not overwrite:
        print('Plots already exists.')
        return None
    recfile = join(c['simpath'], 'g2p_{}.root'.format(run))
    simfile = join(c['simpath'], 'run_{}.root'.format(run))
    if not exists(recfile) or not exists(simfile):
        print('Require {} and {} as input.'.format(recfile, simfile))
        return None

    t1 = ROOT.TChain('T')
    t1.Add(recfile)
    t2 = ROOT.TChain('T')
    t2.Add(simfile)
    simconf = config.expand(getattr(config_list, c['config']))
    cgrt = simconf['generator']['react']['t']
    cgrp = simconf['generator']['react']['p']
    h11 = [None] * 49
    h12 = [None] * 49
    h13 = [None] * 49
    h21 = [None] * 49
    h22 = [None] * 49
    h23 = [None] * 49

    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    c1 = ROOT.TCanvas('c1', 'c1', 300 * 7, 200 * 7)
    c1.Divide(7, 7)
    c2 = ROOT.TCanvas('c2', 'c2', 300 * 7, 200 * 7)
    c2.Divide(7, 7)
    c3 = ROOT.TCanvas('c3', 'c3', 300 * 7, 200 * 7)
    c3.Divide(7, 7)
    for i in range(49):
        h11[i] = ROOT.TH1D('h11' + str(i), 'Theta', 200, cgrt[0] - 0.01, cgrt[1] + 0.01)
        h21[i] = h11[i].Clone('h21' + str(i))
        h12[i] = ROOT.TH1D('h12' + str(i), 'Phi', 200, cgrp[0] - 0.01, cgrp[1] + 0.01)
        h22[i] = h12[i].Clone('h22' + str(i))
        h13[i] = ROOT.TH1D('h13' + str(i), 'Delta', 200, c[run]['delta'] - 0.01, c[run]['delta'] + 0.005)
        h23[i] = h13[i].Clone('h23' + str(i))

        cut1 = 'data.id=={}'.format(i)
        cut2 = '(isgood&&fwd.id.hole=={})*phys.react.xs'.format(i)
        N1 = t1.Project('h11' + str(i), 'rec.t', cut1)
        N2 = t2.Project('h21' + str(i), 'rec.t', cut2)
        if N1 < 50 or N2 < 50:
            continue
        N1 = h11[i].GetMaximum()
        N2 = h21[i].GetMaximum()
        ratio = float(N1) / N2
        cut2 = '(isgood&&fwd.id.hole=={})*(phys.react.xs*{})'.format(i, ratio)

        c1.cd((7 - i % 7) * 7 - i // 7)
        t1.Draw('rec.t>>h11{}'.format(i), cut1)
        t2.Draw('rec.t>>h21{}'.format(i), cut2, 'same')
        h21[i].SetLineColor(2)
        c1.Update()
        c2.cd((7 - i % 7) * 7 - i // 7)
        t1.Draw('rec.p>>h12{}'.format(i), cut1)
        t2.Draw('rec.p>>h22{}'.format(i), cut2, 'same')
        h22[i].SetLineColor(2)
        c2.Update()
        c3.cd((7 - i % 7) * 7 - i // 7)
        t1.Draw('rec.d>>h13{}'.format(i), cut1)
        t2.Draw('rec.d>>h23{}'.format(i), cut2, 'same')
        h23[i].SetLineColor(2)
        c3.Update()

    c1.Print(tplot)
    c2.Print(pplot)
    c3.Print(dplot)
    ROOT.gROOT.SetBatch(ROOT.kFALSE)
