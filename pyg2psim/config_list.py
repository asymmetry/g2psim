# config_list.py

# defaults
defaults = {
    'n': 50000,
    'filename': 'test.root',
#    'conf': 'sim.cfg',
#    'debug': 1,
#    'seed': 1,
    'run': {
        'e0': 2.2535,
        'p0': 2.2495,
        'angle': 5.767,
#        'particle': 11,
        'hrs': '484816',
        'target': 'prod',
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
        'res': [0.2e-3, 0.4e-3]
    },
    'backward': {
        'xpars': [0, 0, 0],
        'ypars': [0, 0, 0],
        'z': 0
    },
    'phys': {
#        'pars': [0, 0],
        'model': 'pbosted'
    },
    'optics': {
        'filename': 'Sieve.full.f51',
        'foilz': 0,
        'bpmz': 0,
        'eloss': [0] * 49
    }
}

# longitudinal setting
run_l_22545000 = {
    'base': defaults,
    'run': {
        'hrs': '400016',
        'field': 'long'
    },
    'backward': {
        'xpars': [0.00079, -0.75038, 178.902],
        'ypars': [0.15693, -0.39864, -176.049]
    },
}

optics_l_22545000 = {
    'base': run_l_22545000,
    'run': {'target': 'optics_C125'},
    'generator': {
        'beam': {'z': -12.5476e-3},
        'raster': {'slow': [0.0002, 0.0002]},
        'react': {
            'd': 'elastic',
            'z': [-0.0141350, -0.0109601]
        }
    },
    'forward': {'sieve': 'in'},
    'backward': {'z': -12.5476e-3},
    'phys': {
        'model': 'elastic',
        'pars': [2]
    },
    'optics': {
        'foilz': -12.5476e-3,
        'bpmz': -12.5476e-3
    }
}
