#!/usr/bin/env python3
"""Separate Fraction witness, consecutive-angle, and cover-logic replay.

The primary verifier separately replays all 153,812 actual normal signs and
1,198 actual cross-pairings. This checker does not duplicate that computation.
"""
import argparse
import json
import sys
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts import fresh_h9_global_h5_20261013_review as A

PRE=A.PRE
OUT='results/fresh_h9_global_h5_20261013_enlarged_review.json'
DEN=10**90
def rounded(z):
    return Q((z[0]*DEN)//1,DEN),Q(-((-z[1]*DEN)//1),DEN)
def atan(q,n):
    s=sum((-1)**j*q**(2*j+1)/Q(2*j+1) for j in range(n))
    t=q**(2*n+1)/Q(2*n+1)
    return (s-t,s) if n%2 else (s,s+t)
pa,pb=atan(Q(1,5),100);pc,pd=atan(Q(1,239),32)
PI=rounded((16*pa-4*pd,16*pb-4*pc))
@lru_cache(None)
def cos2(m):
    assert isinstance(m,int) and m>=2
    if m==2:return Q(0),Q(0)
    if m==3:return Q(1),Q(1)
    x=rounded((PI[0]/m,PI[1]/m));xx=rounded(A.im(x,x))
    term=(Q(1),Q(1));ans=term
    for j in range(1,71):
        z=A.ni(A.im(term,xx));c=(2*j)*(2*j-1)
        term=rounded((z[0]/c,z[1]/c));ans=rounded(A.ia(ans,term))
    err=Q(2)**142
    for j in range(1,143):err/=j
    return rounded((2*(ans[0]-err),2*(ans[1]+err)))

def run(verify=False):
    d=A.read(PRE+'_enlarged_verify.json');recipes=A.read(PRE+'_enlarged_recipes.json')
    assert d['status']=='PASS_EXACT_38209_PLUS119_BANK_EXCLUSION_WITH_OPTIONAL_OLD_WALL_DELETION'
    for p,s in d['sources'].items():assert A.sha(p)==s,p
    for p,s in recipes['sources'].items():assert A.sha(p)==s,p
    assert d['total_dictionary_entries']==38453 and d['complete_sign_checks']==153812
    assert d['old_optional_walls']==125 and d['original_formula_caps']==38209 and d['new_opposite_formula_caps']==119
    original=A.read(PRE+'_recipes.json');lookup={(r['source'],r['row']):r for r in original['recipes']}
    assert len(lookup)==38209 and len(recipes['new_recipes'])==119
    assert len({(r['source'],r['row']) for r in recipes['new_recipes']})==119
    for r in recipes['new_recipes']:
        old=lookup[(r['source'],r['row'])]
        assert r['incident']==old['incident'] and r['profile']==old['profile'] and r['branch']==-old['branch']
        assert {'seed2','compact1','compact2','compact3'}<=set(r['incident'])

    raw=A.read(A.OLD);G=[[A.f(z) for z in r] for r in raw['gram']]
    I=[[A.f(z) for z in r] for r in raw['inverse']];h=list(map(A.f,raw['center_covector']))
    g=[[A.interval(z) for z in r] for r in G]
    compact=A.read(A.COMPACT);seed=next(r for r in raw['roots'] if r['name']=='seed2')
    fixed=[[A.interval(A.f(z)) for z in seed['coefficients']]]+[[A.frozen(z) for z in r] for r in compact['normal_coordinate_enclosures'][1:]]
    witnesses=[]
    assert len(d['witnesses'])==4
    for w in d['witnesses']:
        y=list(map(A.f,w['covector']));q=A.pair(y,I,y);future=A.pair(y,I,h)
        assert q==A.f(w['norm']) and future==A.f(w['future_pairing']) and A.sign(q)>0 and A.sign(future)<0
        x=[A.interval(A.neg(z)) for z in A.mv(I,y)]
        co=[A.im(A.ipair(x,g,e),(Q(1,2),Q(1,2))) for e in fixed];assert all(z[1]<0 for z in co)
        xf=[A.ia(a,A.ni(A.itotal(A.im(c,e[j]) for c,e in zip(co,fixed,strict=True)))) for j,a in enumerate(x)]
        fq=A.ipair(xf,g,xf);ff=A.idot(xf,[A.interval(z) for z in h]);assert fq[0]>0 and ff[0]>0
        for actual,saved in [(fq,w['H5_projected_norm']),(ff,w['H5_future_covector_pairing'])]:
            saved=A.frozen(saved);assert max(actual[0],saved[0])<=min(actual[1],saved[1])
        assert len(w['cutters'])==len(set(w['cutters'])) and all(0<=j<38453 for j in w['cutters'])
        assert int(w['maximum_noncutter_pairing_upper'])<0<int(w['minimum_cutter_pairing_lower'])
        witnesses.append({'norm':A.enc(q),'future_pairing':A.enc(future),'H5_norm_interval':A.ei(fq),'H5_future_covector_interval':A.ei(ff),'cutters':w['cutters']})
    assert list(map(lambda w:len(w['cutters']),witnesses))==[15,54,1,2]

    kinds=Counter();gaps=Counter();branches=set();assignments=Counter();checks=0
    for row in d['coverage']:
        j=row['required_new_branch'];assert j not in branches and 0<=j<119;branches.add(j)
        wi=row['witness'];assert wi in range(4) and wi==recipes['assignment'][str(j)];assignments[wi]+=1
        cutters=set(witnesses[wi]['cutters']);assert 38334+j not in cutters
        assert len(row['cutter_conflicts'])==len(cutters)
        assert {p['cutter_entry'] for p in row['cutter_conflicts']}==cutters
        for p in row['cutter_conflicts']:
            interval=A.frozen(p['pairing']);assert interval[0]<=interval[1];kinds[p['kind']]+=1;checks+=1
            if p['kind']=='positive_pairing':assert interval[0]>0
            else:
                assert p['kind']=='between_consecutive_finite_angles'
                m,n=p['labels'];assert n==m+1
                # Recompute these cosine bounds with our own rational rounding
                # and Taylor routine; do not trust the saved cosine intervals.
                lo=A.ni(cos2(n));hi=A.ni(cos2(m))
                assert lo[1]<interval[0] and interval[1]<hi[0]
                assert interval[0]>-2
                gaps[(m,n)]+=1
    assert branches==set(range(119)) and checks==d['conflict_checks']==1198
    assert dict(assignments)=={0:68,1:2,2:28,3:21}
    out={'status':'PASS_SEPARATE_FOUR_WITNESS_COVER_AND_ALL_LABEL_GAP_AUDIT',
         'sources':{p:A.sha(p) for p in [PRE+'_enlarged_verify.json',PRE+'_enlarged_recipes.json',PRE+'_verify.json',A.OLD,A.COMPACT,'scripts/'+Path(A.__file__).name,'scripts/'+Path(__file__).name]},
         'dictionary_entries':38453,'opposite_branches':119,'cutters_per_witness':[15,54,1,2],
         'witnesses':witnesses,'assignments':{str(k):v for k,v in assignments.items()},
         'independent_conflict_checks':checks,'conflict_kinds':dict(kinds),'consecutive_label_gaps':{str(k):v for k,v in sorted(gaps.items())},
         'scope':'Separate standard-library Fraction witness/projection and cover replay, including independently recomputed cosine gap bounds. The complete actual normal signs and actual cross-pairings are checked by the separately source-reviewed and freshly replayed primary verifier; they are not implemented a second time here.'}
    if verify:assert out==A.read(OUT)
    else:(A.ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'])
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');run(p.parse_args().verify)
