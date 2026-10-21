#!/usr/bin/env python3
"""Pure exact loader for the frozen74-wall split-seed chamber."""
from scripts.fresh_h9_split_endpoint_walk_exact_20261013_base import *

INITIAL=PREFIX
PREFIX='data/enumerations/fresh_h9_split_second_frontier_20261013'

def current():
    G,I,h,old,scout,banks=source()
    family=read(INITIAL+'_family.json');geometry=read(INITIAL+'_family_geometry.json');case=geometry['cases'][0]
    for i in case['new_caps']:
        row=load_cap(family['caps'][i]);old[row['name']]=row
    assert list(old)==case['walls']and len(old)==74
    return G,I,h,old,scout,case,family

def field_record(f):
    return {'kind':'K'}if f is K else {'kind':'quadratic','base':field_record(f.base),'radicand':f.base.encode(f.d)}

def field_parse(r):
    if r['kind']=='K':return K
    assert r['kind']=='quadratic';base=field_parse(r['base']);return Quadratic(base,base.parse(r['radicand']))

def profile(record):
    F=field_parse(record['field'])
    return {'name':record['name'],'field':F,'c':tuple(map(F.parse,record['coefficients'])),
            'p':tuple(map(F.parse,record['covector'])),'q':scale(O,2)}

def selected():
    G,I,h,old,_,_,_=current();bank=read(PREFIX+'_caps.json');geometry=read(PREFIX+'_geometry.json')
    for record in bank['caps']:
        row=profile(record);old[row['name']]=row
    assert list(old)==geometry['walls']and len(old)==79
    return G,I,h,old,geometry
