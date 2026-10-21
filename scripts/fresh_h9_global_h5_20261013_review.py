#!/usr/bin/env python3
"""Separate small Fraction audit of the common-bank escape and provenance.

This does not repeat the 38,209 formula evaluations.  Those are replayed by
the separately reviewed primary verifier.  No project arithmetic is imported.
"""
import argparse
import hashlib
import json
from collections import Counter
from fractions import Fraction as Q
from math import isqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PRE = 'data/enumerations/fresh_h9_global_h5_20261013'
OLD = 'results/fresh_h9_split_second_frontier_independent_geometry_20261013_input.json'
COMPACT = 'data/enumerations/fresh_h9_compact_next_20261013_verify.json'
OUT = 'results/fresh_h9_global_h5_20261013_review.json'
ZERO = (Q(0), Q(0))
ONE = (Q(1), Q(0))

def read(p): return json.loads((ROOT / p).read_text())
def sha(p): return hashlib.sha256((ROOT / p).read_bytes()).hexdigest()
def f(x): return Q(x['a']), Q(x['b'])
def add(x,y): return x[0]+y[0], x[1]+y[1]
def neg(x): return -x[0], -x[1]
def mul(x,y): return x[0]*y[0]+x[1]*y[1], x[0]*y[1]+x[1]*y[0]+x[1]*y[1]
def scale(x,q): return x[0]*q,x[1]*q
def total(xs):
    ans=ZERO
    for x in xs: ans=add(ans,x)
    return ans
def dot(x,y): return total(mul(a,b) for a,b in zip(x,y,strict=True))
def mv(A,x): return [dot(r,x) for r in A]
def pair(x,A,y): return dot(x,mv(A,y))
def enc(x): return {'a':str(x[0]),'b':str(x[1])}
def sign(x):
    a,b=x[0]+x[1]/2,x[1]/2
    if not a: return (b>0)-(b<0)
    if not b or a*b>0: return (a>0)-(a<0)
    d=a*a-5*b*b
    assert d
    return ((a>0)-(a<0))*((d>0)-(d<0))

# Rational intervals use exact Fraction endpoints, with no shared rounding code.
def ia(x,y): return x[0]+y[0],x[1]+y[1]
def ni(x): return -x[1],-x[0]
def im(x,y):
    t=[a*b for a in x for b in y]
    return min(t),max(t)
def itotal(xs):
    ans=(Q(0),Q(0))
    for x in xs: ans=ia(ans,x)
    return ans
def idot(x,y): return itotal(im(a,b) for a,b in zip(x,y,strict=True))
def ipair(x,A,y): return idot(x,[idot(r,y) for r in A])
def ei(x): return list(map(str,x))
DEN=10**110
SQ5=isqrt(5*DEN*DEN)
PHI=(Q(DEN+SQ5,2*DEN),Q(DEN+SQ5+1,2*DEN))
def interval(x): return ia((x[0],x[0]),im((x[1],x[1]),PHI))
def frozen(x): return tuple(Q(int(t),10**80) for t in x)

def run(verify=False):
    d=read(PRE+'_verify.json'); recipes=read(PRE+'_recipes.json')
    assert d['status']=='PASS_EXACT_COMMON_ESCAPE_ALL_38209_LITERAL_FORMULA_CAPS'
    assert d['retained_walls']==125 and d['literal_caps']==38209 and d['distinct_incidence_groups']==121
    for p,s in d['sources'].items(): assert sha(p)==s,p
    for p,s in recipes['sources'].items(): assert sha(p)==s,p
    raw=read(OLD); G=[[f(x) for x in r] for r in raw['gram']]
    I=[[f(x) for x in r] for r in raw['inverse']]
    h=list(map(f,raw['center_covector'])); y=list(map(f,recipes['exact_escape_proposal']))
    for i,row in enumerate(G):
        for j in range(10): assert dot(row,[I[k][j] for k in range(10)])==(ONE if i==j else ZERO)
    q=pair(y,I,y); future=pair(y,I,h)
    assert q==f(d['strict_escape']['norm']) and future==f(d['strict_escape']['future_pairing'])
    assert sign(q)>0 and sign(future)<0 and sign(pair(h,I,h))<0
    assert len(d['strict_escape']['all125covector_pairings'])==125
    assert all(int(z[0])>0 for z in d['strict_escape']['all125covector_pairings'])

    compact=read(COMPACT)
    assert [z['name'] for z in compact['new_normals']]==['compact0','compact1','compact2','compact3']
    assert all(z['norm']==2 for z in compact['new_normals'])
    seed=next(r for r in raw['roots'] if r['name']=='seed2')
    assert seed['field']=={'kind':'golden'} and f(seed['norm'])==scale(ONE,2)
    fixed=[[interval(f(z)) for z in seed['coefficients']]]+[[frozen(z) for z in r] for r in compact['normal_coordinate_enclosures'][1:]]
    g=[[interval(z) for z in r] for r in G]
    # Their exact right/norm-two identities are inherited from the reviewed
    # compact-vertex proof. Here separately check the coordinate enclosures.
    for i,a in enumerate(fixed):
        for j,b in enumerate(fixed):
            p=ipair(a,g,b); expected=2 if i==j else 0
            assert p[0]<=expected<=p[1]
    x=[interval(neg(z)) for z in mv(I,y)]
    co=[im(ipair(x,g,e),(Q(1,2),Q(1,2))) for e in fixed]
    assert all(z[1]<0 for z in co)
    xf=[ia(a,ni(itotal(im(c,e[j]) for c,e in zip(co,fixed,strict=True)))) for j,a in enumerate(x)]
    fq=ipair(xf,g,xf); ff=idot(xf,[interval(z) for z in h])
    assert fq[0]>0 and ff[0]>0
    recorded_future=frozen(d['H5_escape']['future_covector_pairing'])
    assert recorded_future[0]>0 and max(recorded_future[0],ff[0])<=min(recorded_future[1],ff[1])
    recorded_norm=frozen(d['H5_escape']['norm'])
    assert max(recorded_norm[0],fq[0])<=min(recorded_norm[1],fq[1])
    assert all(int(z[1])<0 for z in d['H5_escape']['old_nonfixed_pairings'])

    count=Counter(); seen=set(); spans=set(); bysource={}
    for r in recipes['recipes']:
        key=(r['source'],r['row']); assert key not in seen;seen.add(key)
        count[r['source']]+=1; spans.add(tuple(r['incident'])); bysource[key]=r
        assert len(r['incident'])==len(set(r['incident']))==len(r['profile'])==9
        assert all(isinstance(m,int) and m>=2 for m in r['profile']) and r['branch'] in (-1,1)
    assert len(seen)==38209 and len(spans)==121 and len(count)==12 and dict(count)==d['bank_counts']
    assert Counter(str(r['branch']) for r in recipes['recipes'])==Counter(d['branch_counts'])
    # Alternate exact anchors are not uncensused new normals: each is the same
    # defining formula as one of the original-bank rows.
    for a in recipes['extra_anchors']:
        r=bysource[('fresh_h9_h5_closure_next_20261013_scout',a['source_row'])]
        for k in ('incident','profile','branch'): assert a[k]==r[k]
    assert int(d['maximum_target_pairing_upper'])<0
    assert all(int(z)<0 for z in d['bank_target_pairing_upper_bounds'].values())
    out={'status':'PASS_SEPARATE_ESCAPE_PROJECTION_PROVENANCE_AUDIT',
         'sources':{p:sha(p) for p in [PRE+'_verify.json',PRE+'_recipes.json',OLD,COMPACT,'scripts/'+Path(__file__).name]},
         'literal_rows':38209,'banks':12,'incidence_spans':121,
         'exact_H9_escape_norm':enc(q),'exact_H9_future_pairing':enc(future),
         'H5_escape_norm_interval':ei(fq),'H5_future_covector_pairing_interval':ei(ff),
         'negative_projection_coefficients':list(map(ei,co)),
         'all_three_alternate_anchor_formulas_found_in_bank':True,
         'scope':'Independent stdlib Fraction escape/projection and coverage-metadata checks, plus separately documented mathematical source review and full primary-verifier replay. Not a second implementation of all 38209 formula signs or inherited root certificates.'}
    if verify: assert out==read(OUT)
    else: (ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'])

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');run(ap.parse_args().verify)
