#!/usr/bin/env python3
"""No-project-import replay of saved native pair values and boundary accounting."""
from collections import Counter
from hashlib import sha256
from itertools import combinations,product
import json
import math
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_wall59_remaining_low_seed_joins_20261018'
CYCLE=[30,67,58,66,33]
M=np.diag([-1.,1.,1.,1.]);J=np.array([[2.,1.],[1.,2.]])/3


def native(x):
    if x < -2-2e-7:return 'dotted'
    if abs(x+2)<=2e-7:return 'parallel_or_unresolved'
    if abs(x)<=2e-7:return 2
    if x>0 or x < -2:return None
    m=int(round(np.pi/np.arccos(-x/2)))
    if m<3:return None
    theta=lambda n:2*np.cos(np.pi/n)
    gap=min(abs(theta(m)-theta(m-1)),abs(theta(m+1)-theta(m)))
    return m if abs(x+theta(m))<=min(2e-7,gap/100)else None


def topology(t):
    fs=t['triangles'];bc=t['boundary'];es=Counter(tuple(sorted(e))for f in fs for e in combinations(f,2));vs=set(sum(fs,[]))
    n=len(t['caps']);assert n in(2,3)
    assert(len(vs),len(es),len(fs))==((7,13,7)if n==2 else(8,16,9))
    assert set(es.values())=={1,2}
    assert{e for e,k in es.items()if k==1}=={tuple(sorted(e))for e in zip(bc,bc[1:]+bc[:1])}
    for v in vs:
        links=[tuple(i for i in f if i!=v)for f in fs if v in f];deg=Counter(i for e in links for i in e);seen={next(iter(deg))}
        while True:
            new=seen|{b for a,b in links if a in seen}|{a for a,b in links if b in seen}
            if new==seen:break
            seen=new
        assert seen==set(deg)and sorted(deg.values())==([1,1]+[2]*(len(deg)-2)if v in bc else[2]*len(deg))


def ledger(l):
    assert not l['warnings']and l['PD_status']=='numerically_separated_acute_incidence'
    assert l['open_ends']==l['PD_open_ends']==sum(1 if e['kind']=='one_open'else 2 if e['kind']=='two_open'else 0 for e in l['edges'])
    assert Counter(r['kind']for r in l['rays'])==Counter(l['counts'])
    assert{tuple(r['active'])for r in l['rays']if r['kind']=='finite'}=={tuple(r)for r in l['PD_finite_supports']}


def run():
    path=ROOT/f'results/{TAG}.json';d=json.loads(path.read_text())
    for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h
    prior=json.loads((ROOT/'results/fresh_wall59_low_bond_companions_20261018.json').read_text())
    rep=json.loads((ROOT/'results/fresh_wall59_bond_replacement_20261018.json').read_text())
    profpath=ROOT/'results/fresh_nonarithmetic_b_release52_patch_20261017.json'
    profiles=json.loads(profpath.read_text())['profiles'];assert len(profiles)==7
    contacts=[]
    for pi,p in enumerate(profiles):
        expected=[0.,0.]
        if pi:expected[0 if pi<=3 else 1]=-2*math.cos(math.pi/(pi+2 if pi<=3 else pi-1))
        assert max(abs(a-b)for a,b in zip(expected,p['contacts']))<1e-14
        contacts.append(np.array(p['contacts']))
    oldlabels={tuple(r)for r in prior['seed_labels']}
    low=[r for r in rep['accepted']if r['release']==[59,60,61]and max(r['labels'])<=5]
    selected=[r for r in low if tuple(r['labels'])not in oldlabels]
    assert len(low)==14 and len(selected)==len(d['records'])==12
    assert d['selection']['previously_tested_labels']==[list(r)for r in sorted(oldlabels)]
    primary=json.loads((ROOT/'results/fresh_nonarithmetic_release49_62_replacement_20261017.json').read_text())
    vertices=primary['retained_old_finite_supports']+[r['support']for r in primary['new_finite_vertices']]
    assert len(vertices)==199
    debit=d['exact_deletion_debit'];lost=[v for v in vertices if set(v)&{59,60,61}]
    assert debit==prior['exact_deletion_debit']and lost==debit['lost_finite_supports']and len(lost)==7
    assert debit['compact_faces_by_dimension']=={'2':6,'3':1}and all(8 not in v for v in lost)
    for t in d['topologies']:topology(t)
    assert len(d['topologies'])==10
    for l in d['boundary_ledgers']:ledger(l)
    capcount=0;paircount=0;totals=Counter();bankdomains=0
    for seed,rec in zip(selected,d['records']):
        assert(rec['source_index'],rec['labels'],rec['normal'])==(seed['index'],seed['labels'],seed['normal'])
        l=d['boundary_ledgers'][rec['seed_score']];assert l['counts']==dict(finite=25,exterior=5)and l['open_ends']==9
        assert len(rec['banks'])==len(rec['pair_summaries'])==5
        for rot,b in enumerate(rec['banks']):
            assert b['frame']==[CYCLE[(rot+j)%5]for j in range(3)]and b['duties']==[[0,1],[1,2]]
            allowed=set()
            for pi,ms,high in b['symbolic_recipes']:
                for hs in product(d['high_label_samples'],repeat=len(high)):
                    labels=list(ms)
                    for i,m in zip(high,hs):labels[i]=m
                    allowed.add((pi,*labels))
            bankdomains+=1
            for i,c in enumerate(b['caps']):
                assert c['index']==i and c['frame']==b['frame']and(c['profile'],*c['labels'])in allowed
                assert c['center_product']<-1e-8 and c['norm_residual']<2e-5 and c['branch']in(-1,1)
                assert all(native(v)not in(None,'parallel_or_unresolved')for _,v in c['native_products'])
                capcount+=1
        for rot,summary in enumerate(rec['pair_summaries']):
            assert summary['rotation']==rot
            bx=rec['banks'][rot]['caps'];by=rec['banks'][(rot+2)%5]['caps'];hist=Counter();nearest=None;dg=sha256()
            assert summary['domain_size']==len(bx)*len(by)
            for ix,x in enumerate(bx):
                for iy,y in enumerate(by):
                    value=float(np.array(x['normal'])@M@np.array(y['normal']))+float(contacts[x['profile']]@J@contacts[y['profile']])
                    rel=native(value);hist[str(rel)]+=1;paircount+=1
                    dg.update(json.dumps([ix,iy,value,rel],separators=(',',':')).encode())
                    if rel not in(None,'parallel_or_unresolved')and(nearest is None or value>nearest['native_product']):
                        nearest=dict(x=ix,y=iy,native_product=value,relation=rel)
            assert dg.hexdigest()==summary['complete_pair_value_digest']
            assert dict(hist)==summary['relation_histogram']and summary['native_finite']==0
            assert not any(k not in('dotted','None')for k in hist)
            totals.update(hist)
            if nearest is None:assert summary['closest_native_pair']is None
            else:
                actual=summary['closest_native_pair'];assert{key:actual[key]for key in nearest}==nearest
                l=d['boundary_ledgers'][actual['score']];assert l['counts']==dict(finite=29,exterior=5)and l['open_ends']==9
                assert{(7,30,47),(30,47,52)}<={tuple(r['active'])for r in l['rays']if r['kind']=='exterior'}
    assert capcount==d['counts']['native_caps']==951 and paircount==d['complete_pair_tests']==18684
    assert not d['finite_junctions']and not d['dynamic_Z_banks']and not d['completed_disks']and not d['unresolved']
    return dict(status='PASS_SAVED_REMAINING_LOW_WALL59_JOINS',source_bindings=len(d['sources']),
        seeds=12,banks=bankdomains,native_caps=capcount,pair_tests=paircount,pair_histogram=dict(totals),pair_digest_checks=60,
        data_sha256=sha256(path.read_bytes()).hexdigest(),profile_input_sha256=sha256(profpath.read_bytes()).hexdigest(),
        scope='No project imports. Independent old199-support debit, complete remaining-seed selection, disk/link and cap-index/native-value checks; all18684 saved-normal pair products and all60 ordered pair digests recomputed with NumPy and explicit A2 offsets. Saved complete cone/PD ledger agreement only; no frame/norm-root, finite classifier, cone, interval-sign or ambient-volume replay.')


if __name__=='__main__':print(json.dumps(run(),sort_keys=True))
