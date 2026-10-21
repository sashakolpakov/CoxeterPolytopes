#!/usr/bin/env python3
"""No-project-import replay of saved native pair values and boundary accounting."""
from collections import Counter,defaultdict
from hashlib import sha256
from itertools import combinations,product
import json
import math
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
PREFIX='fresh_wall59_nonright55_prefix_20261018'
COMP='fresh_wall59_nonright55_companions_20261018'
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
    a=json.loads((ROOT/f'results/{PREFIX}.json').read_text())
    d=json.loads((ROOT/f'results/{COMP}.json').read_text())
    for doc in (a,d):
        for p,h in doc['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h
        assert not doc['unresolved']
        for l in doc['boundary_ledgers']:ledger(l)
    pri=json.loads((ROOT/'results/fresh_nonarithmetic_release49_62_replacement_20261017.json').read_text())
    vs=[tuple(v) for v in pri['retained_old_finite_supports']]+[tuple(v['support']) for v in pri['new_finite_vertices']]
    assert len(vs)==len(set(vs))==199
    lost={i for i,v in enumerate(vs) if set(v)&{59,60,61}}
    assert len(lost)==7 and all(8 not in vs[i] for i in lost)
    edges=Counter(e for v in vs for e in combinations(v,8));faces=defaultdict(list)
    for i,v in enumerate(vs):
        for k in range(8):
            for f in combinations(v,k):faces[f].append(i)
    touched=[]
    for f,indices in faces.items():
        if not(lost&set(indices)) or len(indices)<10-len(f):continue
        if all(edges[tuple(j for j in vs[i] if j!=k)]==2 for i in indices for k in vs[i] if k not in f):
            touched.append(dict(support=list(f),dimension=9-len(f),finite_supports=[list(vs[i]) for i in indices]))
    debit=a['exact_deletion_debit']
    assert debit==d['exact_deletion_debit']
    assert debit['lost_finite_supports']==[list(vs[i]) for i in sorted(lost)]
    assert touched==debit['compact_faces_touched']
    assert Counter(r['dimension'] for r in touched)=={2:6,3:1}
    contacts=[np.array([0.,0.])]
    for side in (0,1):
        for m in (3,4,5):
            c=np.zeros(2);c[side]=-2*np.cos(np.pi/m);contacts.append(c)
    squares=[2-float(c@J@c) for c in contacts]
    raw={(pi,a,b,c) for pi,a,b,c in product(range(7),(3,4,5),range(2,6),range(2,6))}
    keys={(r['profile'],*r['labels']) for r in a['recipes']}
    assert len(raw)==336 and len(keys)==len(a['recipes'])==len(a['equations'])==38 and keys<=raw
    assert [(r['profile'],r['labels'],r['corners']) for r in a['recipes']]==[(r['profile'],r['labels'],r['corners']) for r in a['equations']]
    accepted={(r['recipe'],r['branch']):r for r in a['accepted']};eligible=set();branches=0
    for i,r in enumerate(a['equations']):
        assert i==r['index'] and r['frame']==[55,30,58]
        assert sorted(x['branch'] for x in r['branches'])==[-1,1]
        for b in r['branches']:
            branches+=1;v=np.array(b['normal'])
            assert abs(float(v@M@v)-squares[r['profile']])<2e-5 and b['norm_residual']<2e-5
            # Six S6 contacts are prescribed by the seven native profiles;
            # the remaining54 retained old walls are checked explicitly.
            vals=dict(b['native_products']);assert len(vals)==54 and not({59,60,61}&set(vals))
            assert sorted(i for i,x in vals.items() if native(x)in(None,'parallel_or_unresolved'))==b['blockers']
            for old,m in zip(r['frame'],r['labels']):assert native(vals[old])==m
            if not b['blockers'] and b['center_product']<-1e-8:eligible.add((i,b['branch']))
    assert branches==76 and eligible==set(accepted) and len(eligible)==8
    for i,r in enumerate(a['accepted']):
        assert r['index']==i and r['profile']==0 and r['labels'][0]==3 and r['whole55_quadrilateral']
        b=a['equations'][r['recipe']];br=next(x for x in b['branches'] if x['branch']==r['branch'])
        for k,v in br.items():assert r[k]==v
        l=a['boundary_ledgers'][r['score']]
        assert l['counts']==dict(finite=25,exterior=5) and l['open_ends']==9
        rays=[x for x in l['rays'] if 55 in x['active']]
        assert all(x['kind']=='finite' for x in rays)
        assert {tuple(x['active']) for x in rays}=={(7,30,55),(7,55,58),(30,55,67),(55,58,67)}
    assert d['selection']['selected_labels']==[[3,3,4],[3,4,3],[3,5,5]] and len(d['records'])==3
    assert d['high_label_samples']==[] and d['low_labels']==[2,3,4,5]
    assert len(d['topologies'])==10
    for t in d['topologies']:topology(t)
    capcount=pairs=equations=branchcount=0;totals=Counter()
    for r in d['records']:
        prefix=a['accepted'][r['source_index']]
        assert r['labels']==prefix['labels'] and r['normal']==prefix['normal']
        assert len(r['banks'])==len(r['pair_summaries'])==5
        for rot,b in enumerate(r['banks']):
            assert b['frame']==[CYCLE[(rot+j)%5] for j in range(3)] and b['duties']==[[0,1],[1,2]]
            assert len(b['recipes'])==b['complete_low_recipe_count']
            equations+=len(b['recipes']);indices=[]
            for ri,rec in enumerate(b['recipes']):
                assert rec['index']==ri and max(rec['labels'])<=5 and min(rec['labels'])>=2
                for br in rec['branches']:
                    branchcount+=1
                    if 'cap'in br:
                        cap=b['caps'][br['cap']];indices.append(br['cap'])
                        assert cap['recipe']==ri and cap['branch']==br['branch']
            assert sorted(indices)==list(range(len(b['caps'])))
            for i,c in enumerate(b['caps']):
                capcount+=1;assert c['index']==i and c['center_product']<-1e-8 and c['norm_residual']<2e-5
                vals=dict(c['native_products']);assert len(vals)==55 and 67 in vals
                assert all(native(v)not in(None,'parallel_or_unresolved') for v in vals.values())
                assert abs(float(np.array(c['normal'])@M@np.array(c['normal']))-squares[c['profile']])<2e-5
                l=d['boundary_ledgers'][c['score']];assert l['counts']==dict(finite=27,exterior=5) and l['open_ends']==9
        for rot,ss in enumerate(r['pair_summaries']):
            bx=r['banks'][rot]['caps'];by=r['banks'][(rot+2)%5]['caps'];dg=sha256();hist=Counter();near=None
            assert ss['rotation']==rot and ss['domain_size']==len(bx)*len(by) and ss['native_finite']==0
            for ix,x in enumerate(bx):
                for iy,y in enumerate(by):
                    val=float(np.array(x['normal'])@M@np.array(y['normal']))+float(contacts[x['profile']]@J@contacts[y['profile']])
                    rel=native(val);assert rel=='dotted';hist[rel]+=1;pairs+=1
                    dg.update(json.dumps([ix,iy,val,rel],separators=(',',':')).encode())
                    if near is None or val>near['native_product']:near=dict(x=ix,y=iy,native_product=val,relation=rel)
            assert dg.hexdigest()==ss['complete_pair_value_digest'] and dict(hist)==ss['relation_histogram'];totals.update(hist)
            if near is None:assert ss['closest_native_pair']is None
            else:
                saved=ss['closest_native_pair'];assert {k:saved[k] for k in near}==near
                l=d['boundary_ledgers'][saved['score']];assert l['counts']==dict(finite=29,exterior=5) and l['open_ends']==9
    assert (equations,branchcount,capcount,pairs)==(513,1026,190,2968)
    assert not d['finite_junctions'] and not d['dynamic_Z_banks'] and not d['completed_disks']
    return dict(status='PASS_NONRIGHT55_SAVED_PREFIX_AND_COMPANIONS',prefix_recipes=38,prefix_branches=76,new_prefixes=8,
                selected_prefixes=3,companion_equations=equations,companion_branches=branchcount,caps=capcount,pairs=pairs,
                pair_histogram=dict(totals),pair_hashes=15,exact_debit=[7,6,1],
                scope='No project imports. Saved normalization/native values, accepted/index completeness, all2968 pair products and15hashes, independent199-support/compact-face debit, disk links and stored fullcone/PD equality. No normal/root/classifier/cone or exact physical interval replay.')


if __name__=='__main__':print(json.dumps(run(),sort_keys=True))
