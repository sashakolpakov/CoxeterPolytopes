#!/usr/bin/env python3
"""Independent saved incidence and numerical-value checks; no project imports."""
from collections import Counter
from hashlib import sha256
from itertools import combinations
from math import acos,cos,pi
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_b_region_split30_native_qr_20261018'


def native(x):
    tol=2e-7
    if x < -2-tol:return 'dotted'
    if abs(x+2)<=tol:return None
    if abs(x)<=tol:return 2
    if x>0 or x < -2:return None
    m=round(pi/acos(-x/2))
    if m<3:return None
    theta=lambda k:2*cos(pi/k)
    spacing=min(abs(theta(m)-theta(m-1)),abs(theta(m+1)-theta(m)))
    return m if abs(x+theta(m))<=min(tol,spacing/100)else None


def main():
    d=json.loads((ROOT/f'results/{TAG}.json').read_text())
    for path,h in d['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==h,path
    old=json.loads((ROOT/'results/fresh_b_region_lr_diagonal_20261017.json').read_text())
    assert d['deleted_walls']==[47]and d['retained_old_wall_count']==62
    assert d['exact_deletion_debit']==old['exact_deletion_debit']
    p=json.loads((ROOT/'results/fresh_nonarithmetic_release49_62_replacement_20261017.json').read_text())
    verts=p['retained_old_finite_supports']+[r['support']for r in p['new_finite_vertices']]
    assert len(verts)==199
    lost=[v for v in verts if 47 in v]
    assert lost==d['exact_deletion_debit']['lost_finite_supports']and len(lost)==2
    assert all(8 not in v for v in lost)
    fs=d['disk']['triangles'];bc=d['disk']['boundary'];es=Counter(tuple(sorted(e))for t in fs for e in combinations(t,2))
    assert(len(set(sum(fs,[]))),len(es),len(fs))==(9,20,12)
    assert (30,70)not in es and (68,69)not in es
    assert {e for e,n in es.items()if n==1}=={tuple(sorted(e))for e in zip(bc,bc[1:]+bc[:1])}
    assert set(es.values())=={1,2}
    for r in d['disk']['links']:
        v=r['vertex'];ee=[tuple(i for i in t if i!=v)for t in fs if v in t]
        assert sorted(ee)==sorted(map(tuple,r['edges']))
        degree=Counter(i for e in ee for i in e);seen={next(iter(degree))}
        while True:
            nxt=seen|{i for e in ee if set(e)&seen for i in e}
            if nxt==seen:break
            seen=nxt
        assert seen==set(degree)
        assert sorted(degree.values())==([1,1]+[2]*(len(degree)-2)if v in bc else[2]*len(degree))
    assert len(d['reconstructed_Q_candidates'])==440
    best={}
    for q in d['reconstructed_Q_candidates']:
        prev=old['near_native_rejections'][q['source_rejection']]
        assert prev['blockers'][0][0]==70 and len(prev['blockers'])==1
        assert q['event']==prev['event']and q['profile']==prev['profile']
        assert q['labels']==prev['labels']and q['branch']==prev['branch']
        assert q['obsolete_QR_product']>0 and abs(q['obsolete_QR_product']-prev['blockers'][0][1])<2e-5
        assert q['center_product']< -1e-8 and q['norm_residual']<2e-5
        k=q['event'],q['profile']
        if k not in best or(q['obsolete_QR_product'],q['index'])<(best[k]['obsolete_QR_product'],best[k]['index']):best[k]=q
    assert d['selected_recoveries']==[q['index']for k,q in sorted(best.items())]
    for l in d['boundary_ledgers']:
        assert not l['warnings']
        assert dict(Counter(r['kind']for r in l['rays']))==l['counts']
        assert {tuple(sorted(r['active']))for r in l['rays']if r['kind']=='finite'}=={tuple(sorted(x))for x in l['PD_finite_supports']}
        assert l['open_ends']==l['PD_open_ends']
        assert l['open_ends']==sum(1 if e['kind']=='one_open'else 2 if e['kind']=='two_open'else 0 for e in l['edges'])
    for r in d['Q_records']:
        l=d['boundary_ledgers'][r['score']]
        assert l['counts']=={'finite':29,'exterior':7}and l['open_ends']==11
    for r in d['J_records']+d['R_records']:
        cap=r.get('J',r.get('R'));l=d['boundary_ledgers'][cap['score']]
        assert cap['center_product']< -1e-8 and cap['norm_residual']<2e-5
        assert all(native(v)is not None for i,v in cap['native_products'])
        assert l['counts']=={'finite':31 if 'J'in r else 33,'exterior':7}and l['open_ends']==11
    stages=Counter()
    for r in d['R_records']:
        stages[r['bank_stage']]+=1
        vals=dict(r['R']['native_products'])
        missing={tuple(fs[i])for i,c in enumerate(r['finite_corners'])if c is None}
        if r['bank_stage']=='R_L_P_Q':
            assert native(vals[7])==native(vals[71])=='dotted'
            assert native(vals[67])==native(vals[69])==2
            assert missing=={(7,70,68),(7,70,71),(69,70,71)}
        else:
            assert r['bank_stage']=='R_Q_J_7'
            assert [native(vals[i])for i in(7,69,71)]==[5,3,2]
            assert native(vals[67])==native(vals[68])=='dotted'
            assert missing=={(7,70,68),(67,68,70),(67,69,70)}
        assert not r['completed_disk']
    assert stages=={'R_L_P_Q':284,'R_Q_J_7':1}
    assert d['counts']==dict(reconstructed_Q=440,selected_Q=66,scored_Q=66,native_J=96,native_R=285)
    assert d['bank_counts']==dict(norm_equations=3225,norm_branches=6450,native_caps=381)
    assert not d['improvements']and not d['unresolved']
    print('PASS_SAVED_SPLIT30_BOUNDARIES',len(d['sources']),'bindings',d['counts'],dict(stages),
          'all285_final_records_33finite_7exterior_11ends')


if __name__=='__main__':main()
