#!/usr/bin/env python

def expand(c):
    r = c.copy()
    while 'base' in r:
        base = r['base']
        del r['base']
        r = merge(base, r)
    return r

def merge(d, o):
    r = d.copy()
    for k, v in o.items():
        if k in r and isinstance(d[k], dict) and isinstance(v, dict):
            r[k] = merge(d[k], v)
        else:
            r[k] = v
    return r

def printf(c):
    r = expand(c)
    _printf(r, '')

def _printf(c, prefix):
    kl = c.keys()
    kl.sort()
    for k in kl:
        v = c[k]
        if isinstance(v, dict):
            nprefix = str(k) + '.'
            _printf(v, nprefix)
        else:
            print(prefix + str(k) + ': {0}'.format(v))
