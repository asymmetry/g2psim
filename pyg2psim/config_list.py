# config_list.py

# defaults
defaults = {
    'n': 50000,
    'filename': 'test.root',
#    'conf': 'sim.cfg',
#    'debug': 1,
#    'seed': 1,
    'run': {
        'type': 'prod',
        'e0': 2.2535,
        'p0': 2.2495,
        'angle': 5.767,
#        'particle': 11,
        'hrs': '484816',
        'field': 'none'
    },
#    'target': {
#        'z': 1,
#        'a': 1,
#        'mass': 0.938890091,
#        'pf': 0.55,
#        'offset': [0, 0, 0]
#    },
#    'field': {
#        'ratio': 0,
#        'origin': [0, 0, 0],
#        'angle': [0, 0, 0]
#    },
    'generator': {
        'beam': {
            'pos': [0, 0],
            'z': 0,
            'angle': [0, 0]
        },
        'profile': [0, 0, 0],
        'raster': {
            'fast': [0, 0],
            'slow': [0.014, 0.014]
        },
        'react': {
            't': [-0.06, 0.06],
            'p': [-0.03, 0.03],
            'd': [-0.04, 0.04],
            'z': [-0.014135, 0.014135]
        }
    },
    'forward': {
        'sieve': 'out',
        'vdc': [1e-4, 5e-4, 1e-4, 5e-4]
    },
    'bpm': {
        'res': {
            'pos': [0.5e-3, 0.5e-3],
            'angle': [0.5e-3, 0.5e-3]
        }
    },
    'backward': {
        'xpars': [0, 0, 0],
        'ypars': [0, 0, 0],
        'z': 0
    },
    'phys': {
        'model': 'pbosted'
#        'pars': {1: 0}
    },
    'optics': {
        'data': 'Sieve.full.f51',
        'foilz': 0,
        'bpmz': 0,
        'eloss': [0] * 49
    }
}

optics_l = {
    'base': 'defaults',
    'run': {'type': 'optics'},
    'generator': {
        'beam': {'z': -13.6271e-3},
        'raster': {'slow': [0.0002, 0.0002]},
        'react': {
            'd': [-0.01, 0.005],
            'z': [-0.0141350, -0.0131191]
        }
    },
    'forward': {'sieve': 'in'},
    'backward': {'z': -13.6271e-3},
    'phys': {
        'model': 'elastic',
        'pars': {1: 1.0 / 3, 2: 2}
    },
    'optics': {
        'foilz': -13.6271e-3,
        'bpmz': -13.6271e-3
    }
}

# transverse 2.5T, 2.254GeV
run_l_22542590 = {
    'base': 'defaults',
    'run': {'field': 'trans'},
    'backward': {
        'xpars': [-0.00061, 3.13661, -4.87024],
        'ypars': [-0.08772, 0.21862, -17.5404]
    }
}

optics_l_22542590 = {
    'base': 'run_l_22542590',
    'run': optics_l['run'],
    'generator': {
        'beam': optics_l['generator']['beam'],
        'raster': optics_l['generator']['raster'],
        'react': {
            't': [-0.1, 0.02],
            'p': [-0.033, 0.032],
            'd': [-0.01, 0.005],
            'z': optics_l['generator']['react']['z']
        }
    },
    'forward': optics_l['forward'],
    'backward': optics_l['backward'],
    'phys': optics_l['phys'],
    'optics': optics_l['optics']
}

# transverse 2.5T, 1.706GeV
run_l_17062590 = {
    'base': 'run_l_22542590',
    'run': {
        'e0': 1.7105,
        'p0': 1.7044,
        'hrs': '400016'
    }
}

optics_l_17062590 = {
    'base': 'optics_l_22542590',
    'run': {
        'e0': 1.7105,
        'p0': 1.7044,
        'hrs': '400016'
    },
    'generator': {'react': {'t': [-0.12, 0]}}
}

# transverse 2.5T, 1.158GeV
run_l_11582590 = {
    'base': 'run_l_22542590',
    'run': {
        'e0': 1.1570,
        'p0': 1.1573,
        'hrs': '400016'
    }
}

optics_l_11582590 = {
    'base': 'optics_l_22542590',
    'run': {
        'e0': 1.1510, # NOTICE: not 1.1570
        'p0': 1.1573,
        'hrs': '400016'
    },
    'generator': {'react': {'t': [-0.16, -0.03]}}
}

# longitudinal
run_l_22545000 = {
    'base': 'defaults',
    'run': {
        'hrs': '400016',
        'field': 'long'
    },
    'backward': {
        'xpars': [0.00079, -0.75038, 178.902],
        'ypars': [0.15693, -0.39864, -176.049]
    }
}

optics_l_22545000 = {
    'base': 'run_l_22545000',
    'run': {'type': 'optics_C125'},
    'generator': {
        'beam': {'z': -12.5476e-3},
        'raster': {'slow': [0.0002, 0.0002]},
        'react': {
            'd': [-0.01, 0.005],
            'z': [-0.0141350, -0.0109601]
        }
    },
    'forward': {'sieve': 'in'},
    'backward': {'z': -12.5476e-3},
    'phys': {
        'model': 'elastic',
        'pars': {1: 1.0 / 3, 2: 2}
    },
    'optics': {
        'foilz': -12.5476e-3,
        'bpmz': -12.5476e-3
    }
}

# transverse 5.0T, 2.254GeV
run_l_22545090 = {
    'base': 'defaults',
    'run': {
        'hrs': '400016',
        'field': 'trans_5T'
    },
    'backward': {
        'xpars': [-0.03853, 6.34323, -4.97020],
        'ypars': [-0.35369, 0.86343, -29.1603]
    }
}

optics_l_22545090 = {
    'base': 'run_l_22545090',
    'run': optics_l_22545000['run'],
    'generator': {
        'beam': optics_l_22545000['generator']['beam'],
        'raster': optics_l_22545000['generator']['raster'],
        'react': {
            't': [-0.17, -0.05],
            'd': [-0.01, 0.005],
            'z': optics_l_22545000['generator']['react']['z']
        }
    },
    'forward': optics_l_22545000['forward'],
    'backward': optics_l_22545000['backward'],
    'phys': optics_l_22545000['phys'],
    'optics': optics_l_22545000['optics']
}

# transverse 5.0T, 3.352GeV
run_l_33525090 = {
    'base': 'run_l_22545090',
    'run': {
        'e0': 3.3514,
        'p0': 2.8200
    }
}

# right arm
optics_r = {
    'base': 'optics_l',
    'run': {'angle': -5.781}
}

run_r_22542590 = {
    'base': 'run_l_22542590',
    'run': {'angle': -5.781},
    'backward': {
        'xpars': [0.00317, 3.13296, 4.63844],
        'ypars': [0.08567, -0.21493, 16.5077]
    }
}

optics_r_22542590 = {
    'base': 'run_r_22542590',
    'run': optics_l['run'],
    'generator': {
        'beam': optics_l['generator']['beam'],
        'raster': optics_l['generator']['raster'],
        'react': {
            't': [-0.12, 0.02],
            'p': [-0.04, 0.04],
            'd': [-0.01, 0.005],
            'z': optics_l['generator']['react']['z']
        }
    },
    'forward': optics_l['forward'],
    'backward': optics_l['backward'],
    'phys': optics_l['phys'],
    'optics': optics_l['optics']
}

run_r_17062590 = {
    'base': 'run_r_22542590',
    'run': {
        'e0': 1.7105,
        'p0': 1.7044,
        'hrs': '400016'
    }
}

optics_r_17062590 = {
    'base': 'optics_r_22542590',
    'run': {
        'e0': 1.7105,
        'p0': 1.7044,
        'hrs': '400016'
    },
    'generator': {'react': {'t': [-0.12, 0]}}
}

run_r_11582590 = {
    'base': 'run_r_22542590',
    'run': {
        'e0': 1.1570,
        'p0': 1.1573,
        'hrs': '400016'
    }
}

optics_r_11582590 = {
    'base': 'optics_r_22542590',
    'run': {
        'e0': 1.1510, # NOTICE: not 1.1570
        'p0': 1.1573,
        'hrs': '400016'
    },
    'generator': {'react': {'t': [-0.16, -0.03]}}
}

run_r_22545000 = {
    'base': 'run_l_22545000',
    'run': {'angle': -5.781},
    'backward': {
        'xpars': [-0.00301, 0.76169, 177.993],
        'ypars': [-0.15435, 0.34665, -179.100]
    }
}

optics_r_22545000 = {
    'base': 'run_r_22545000',
    'run': optics_l_22545000['run'],
    'generator': optics_l_22545000['generator'],
    'forward': optics_l_22545000['forward'],
    'backward': optics_l_22545000['backward'],
    'phys': optics_l_22545000['phys'],
    'optics': optics_l_22545000['optics']
}

run_r_22545090 = {
    'base': 'run_l_22545090',
    'run': {'angle': -5.781},
    'backward': {
        'xpars': [-0.05071, 6.36930, 5.63490],
        'ypars': [0.34858, -0.85407, 30.4795]
    },
}

optics_r_22545090 = {
    'base': 'run_r_22545090',
    'run': optics_l_22545000['run'],
    'generator': {
        'beam': optics_l_22545000['generator']['beam'],
        'raster': optics_l_22545000['generator']['raster'],
        'react': {
            't': [-0.17, -0.05],
            'd': [-0.01, 0.005],
            'z': optics_l_22545000['generator']['react']['z']
        }
    },
    'forward': optics_l_22545000['forward'],
    'backward': optics_l_22545000['backward'],
    'phys': optics_l_22545000['phys'],
    'optics': optics_l_22545000['optics']
}

run_r_33525090 = {
    'base': 'run_r_22545090',
    'run': {
        'e0': 3.3514,
        'p0': 2.8200
    }
}
