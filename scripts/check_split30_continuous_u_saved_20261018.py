#!/usr/bin/env python3
"""Saved boundary/value audit only: no project or numerical-library imports."""
from collections import Counter
from hashlib import sha256
from itertools import combinations,product
from math import acos,cos,pi
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_split30_continuous_u_20261018'

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
    for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    prior=json.loads((ROOT/'results/fresh_b_region_split30_native_qr_20261018.json').read_text())
    assert d['disk']==prior['disk']and d['exact_deletion_debit']==prior['exact_deletion_debit']
    primary=json.loads((ROOT/'results/fresh_nonarithmetic_release49_62_replacement_20261017.json').read_text())
    vertices=primary['retained_old_finite_supports']+[r['support']for r in primary['new_finite_vertices']]
    assert len(vertices)==199
    lost=[v for v in vertices if 47 in v]
    assert lost==d['exact_deletion_debit']['lost_finite_supports']and len(lost)==2
    assert not any(8 in v for v in lost)
    fs=d['disk']['triangles'];bc=d['disk']['boundary'];es=Counter(tuple(sorted(e))for t in fs for e in combinations(t,2))
    assert(len(set(sum(fs,[]))),len(es),len(fs))==(9,20,12)
    assert (30,70)not in es and (68,69)not in es
    assert {e for e,n in es.items()if n==1}=={tuple(sorted(e))for e in zip(bc,bc[1:]+bc[:1])}
    assert set(es.values())=={1,2}
    for v in set(sum(fs,[])):
        ee=[tuple(i for i in t if i!=v)for t in fs if v in t];ds=Counter(i for e in ee for i in e);seen={next(iter(ds))}
        while True:
            nxt=seen|{i for e in ee if set(e)&seen for i in e}
            if nxt==seen:break
            seen=nxt
        assert seen==set(ds)
        assert sorted(ds.values())==([1,1]+[2]*(len(ds)-2)if v in bc else[2]*len(ds))
    recipes=[(0,2,2),(0,2,3),(0,2,4),(0,3,2),(1,2,2),(2,2,2),(4,2,2)]
    labels=[2,3,4,5,6,8,12,24]
    expected={(k,qs,p,m7,mP,rs,3,2)for k,(p,m7,mP),qs,rs in product(labels,recipes,[-1,1],[-1,1])}
    expected|={(k,qs,0,3,2,rs,q30,q52)for k,qs,rs,(q30,q52)in product(labels,[-1,1],[-1,1],[(2,2),(2,3)])}
    assert len(expected)==288 and {tuple(r['key'])for r in d['branches']}==expected
    assert len(d['domain']['u_grid'])==122 and all(r['real_samples']==122 for r in d['branches'])
    assert len(d['brackets'])==43
    for b in d['brackets']:
        a,z,fa,fz=b['interval'];assert a<z and fa*fz<0
        assert tuple(b['key'])in expected and a in d['domain']['u_grid']and z in d['domain']['u_grid']
    represented=[e['bracket']for e in d['events']]+[e['bracket']for e in d['rejected_events']]
    assert sorted(represented)==list(range(43))
    for e in d['rejected_events']:assert e['blockers']and native(e['gap'])==2
    for led in d['boundary_ledgers']:
        assert not led['warnings']and dict(Counter(r['kind']for r in led['rays']))==led['counts']
        assert {tuple(sorted(r['active']))for r in led['rays']if r['kind']=='finite'}=={tuple(sorted(r))for r in led['PD_finite_supports']}
        assert led['open_ends']==led['PD_open_ends']==sum({'one_open':1,'two_open':2}.get(e['kind'],0)for e in led['edges'])
    assert d['boundary_ledgers'][d['control_score']]['counts']==dict(finite=26,exterior=6)
    assert d['boundary_ledgers'][d['deletion_score']]['counts']==dict(finite=24,exterior=6)
    jcount=0;patterns=Counter()
    oldext={(30,33,66),(30,59,66),(58,60,66),(59,60,66)}
    for e in d['events']:
        b=d['brackets'][e['bracket']];assert e['key']==b['key']and b['interval'][0]<=e['u']<=b['interval'][1]
        assert native(e['gap'])==2
        k,qs,p,m7,mP,rs,q30,q52=e['key'];assert(qs,p,m7,mP,rs)==(-1,0,3,2,-1)
        byid={r['node']:r for r in e['native']};assert set(byid)=={67,68,69,70}
        for r in e['native']:
            assert r['profile']==0 and not r['blockers']and r['center_product']< -1e-8 and r['norm_residual']<2e-5
            v=r['normal'];assert abs(-v[0]**2+sum(x*x for x in v[1:])-2)<2e-5
            assert len(r['native_products'])==59 and all(native(x)is not None for i,x in r['native_products'])
        q=dict(byid[69]['native_products']);r=dict(byid[70]['native_products'])
        assert[native(q[i])for i in(30,52,67)]==[q30,q52,k]
        assert[native(r[i])for i in(7,67,68,69)]==[3,2,2,2]
        ledger=d['boundary_ledgers'][e['score']]
        assert ledger['open_ends']==(13 if(q30,q52)==(2,3)else 10)
        assert e['J_bank']['symbolic_recipes']==[[0,[2,2,m],[]]for m in(2,3,4)]
        assert len(e['completed'])==len(e['J_bank']['caps'])
        for j,c in zip(e['J_bank']['caps'],e['completed']):
            assert c['cap']==j['index']and not c['completed_disk'];jcount+=1
            vals=dict(j['native_products']);assert all(native(x)is not None for x in vals.values())
            assert native(vals[69])=='dotted'and j['center_product']< -1e-8 and j['norm_residual']<2e-5
            missing={tuple(fs[i])for i,v in enumerate(c['finite_corners'])if v is None}
            assert missing=={(30,69,71),(69,70,71)}
            led=d['boundary_ledgers'][j['score']];assert led['counts']==dict(finite=34,exterior=6)and led['open_ends']==10
            ext={tuple(r['active'])for r in led['rays']if r['kind']=='exterior'}
            assert ext==oldext|{(30,69,71),(69,70,71)}
            patterns[tuple(sorted(ext))]+=1
    assert jcount==15 and len(patterns)==1
    assert d['counts']==dict(sample_real=35136,branches=288,brackets=43,native_events=15,native_J=15,completed_disks=0)
    assert d['J_bank_counts']==dict(norm_equations=45,norm_branches=90,native_caps=15)
    assert not d['improvements']and not d['unresolved']
    print('PASS_SAVED_CONTINUOUS_U',len(d['sources']),'bindings',d['counts'],'final_boundary_patterns',len(patterns))

if __name__=='__main__':main()
