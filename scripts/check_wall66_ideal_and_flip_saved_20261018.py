#!/usr/bin/env python3
"""Independent stdlib saved-data audit; no normals, cone or PSD recomputation."""
import json,math
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]


def native(x):
    if x < -2-2e-7:return 'dotted'
    if abs(x+2)<=2e-7:return 'parallel_or_unresolved'
    if abs(x)<=2e-7:return 2
    if x>0 or x < -2:return None
    m=round(math.pi/math.acos(-x/2))
    if m<3:return None
    th=lambda k:2*math.cos(math.pi/k)
    gap=min(abs(th(m)-th(m-1)),abs(th(m+1)-th(m)))
    return m if abs(x+th(m))<=min(2e-7,gap/100)else None


def topology(d):
    fs=d['triangles'];bc=d['boundary'];vs=set(sum(fs,[]));es=Counter(tuple(sorted(e))for f in fs for e in combinations(f,2))
    assert(len(vs),len(es),len(fs))==(9,19,11)and set(es.values())=={1,2}
    assert{e for e,n in es.items()if n==1}=={tuple(sorted(e))for e in zip(bc,bc[1:]+bc[:1])}
    for v in vs:
        links=[tuple(x for x in f if x!=v)for f in fs if v in f];deg=Counter(x for e in links for x in e);seen={next(iter(deg))}
        while True:
            new=seen|{b for a,b in links if a in seen}|{a for a,b in links if b in seen}
            if new==seen:break
            seen=new
        assert seen==set(deg)and sorted(deg.values())==([1,1]+[2]*(len(deg)-2)if v in bc else[2]*len(deg))
    return es


def load(stem):
    path=ROOT/f'results/{stem}.json';d=json.loads(path.read_text())
    for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h
    return d


def run():
    old=load('fresh_wall66_seam_subdivision_20261018')
    ideal=load('fresh_wall66_ideal_junction_pilot_20261018');flip=load('fresh_wall66_boundary_bridge_flip_20261018')
    p=json.loads((ROOT/'results/fresh_nonarithmetic_release49_62_replacement_20261017.json').read_text())
    vs=p['retained_old_finite_supports']+[r['support']for r in p['new_finite_vertices']];assert len(vs)==199
    lost=[v for v in vs if set(v)&{59,60,61}];assert len(lost)==7 and all(8 not in v for v in lost)
    for d in(ideal,flip):
        debit=d['exact_old_deletion_debit'];assert debit==old['exact_deletion_debit']and debit['lost_finite_supports']==lost
        assert debit['compact_faces_by_dimension']=={'2':6,'3':1};topology(d['topology'])
    edges=topology(flip['topology']);assert(68,71)not in edges and edges[67,70]==2
    assert ideal['profile']==0 and ideal['prefix_count']==flip['prefix_count']==16 and ideal['frozen_J_count']==95
    wanted={(r['prefix'],i)for r in ideal['domains']for i in range(len(r['recipes']))}
    tested=set();branches={};normbranches=0
    for r in ideal['tests']:
        k=(r['prefix'],r['recipe']);assert k not in tested;tested.add(k)
        for b in r['branches']:
            key=(*k,b['sign']);assert key not in branches;branches[key]=b;normbranches+=1
            assert b['norm_error']<2e-5
    assert wanted==tested and len(tested)==208 and normbranches==400
    expected=set();jtests=0
    for key,b in branches.items():
        if b['blockers']or b['center_product']>=-1e-8:assert not b['J_tests'];continue
        inds=[r['J']for r in b['J_tests']]
        assert inds==list(range(len(old['records'][key[0]]['J_bank']['caps'])))
        for r in b['J_tests']:
            jtests+=1;assert r['relation']==native(r['native_JZ'])
            if r['relation']not in(None,'parallel_or_unresolved'):expected.add((*key,r['J']))
    assert jtests==957
    actual=set()
    for r in ideal['native_candidates']:
        key=(r['prefix'],r['recipe'],r['branch'],r['J']);assert key not in actual;actual.add(key)
        ms=ideal['domains'][r['prefix']]['recipes'][r['recipe']]['labels'];parallel={i for i,m in zip([30,67,68],ms)if m==0}
        vals=dict(r['native_products']);assert len(vals)==57
        assert all(native(v)not in(None,'parallel_or_unresolved')or i in parallel and abs(v+2)<=2e-7 for i,v in vals.items())
        assert vals[69]==r['native_YZ'] and native(r['native_YZ'])==native(r['native_JZ'])=='dotted'
        assert r['mandatory_dotted_seams']==['Y/Z','J/Z']
        diag=ideal['complete_cone_diagnostics'][r['diagnostic']]
        assert diag['PD_incidence']['status']=='unresolved_possible_affine_supports'
        assert any(q['kind']=='unresolved'for q in diag['boundary']['rays'])
        assert r['status']=='UNPROMOTED_IDEAL_CANDIDATE_REQUIRES_AFFINE_INCIDENCE_CERTIFICATE'
    assert actual==expected and len(actual)==924 and not ideal['unresolved']
    nroots=0;ncaps=0;bankkeys=set()
    for i,r in enumerate(flip['records']):
        assert r['prefix']==i and r['selection']==old['records'][i]['selection']
        b=r['J_bank'];assert b['frame']==[66,67,68]and b['duties']==[[0,2],[1,2]]
        assert len(b['recipes'])==b['complete_low_recipe_count'];claims=[]
        for j,t in enumerate(b['recipes']):
            assert t['index']==j and all(m in[2,3,4,5]for m in t['labels'])
            assert sorted(q['branch']for q in t['branches'])==[-1,1];nroots+=2
            for q in t['branches']:
                if 'cap'in q:
                    ci=q['cap'];claims.append(ci);c=b['caps'][ci]
                    assert c['recipe']==j and c['branch']==q['branch']and c['labels']==t['labels']and c['profile']==t['profile']
        assert sorted(claims)==list(range(len(b['caps'])))
        for j,c in enumerate(b['caps']):
            assert c['index']==j and c['center_product']<-1e-8 and c['norm_residual']<2e-5
            assert all(native(v)not in(None,'parallel_or_unresolved')for _,v in c['native_products'])
            bankkeys.add((i,j));ncaps+=1
    assert nroots==1068 and ncaps==233
    keys=set()
    for r in flip['YJ_tests']:
        k=(r['prefix'],r['J']);assert k not in keys;keys.add(k)
        assert r['relation']==native(r['native_YJ'])=='dotted'and r['corner66_Y_J']is None
    assert keys==bankkeys and not flip['Z_banks']and not flip['completion_tests']and not flip['completions']and not flip['unresolved']
    for l in flip['boundary_ledgers']:
        assert not l['warnings']and l['PD_status']=='numerically_separated_acute_incidence'
        assert l['open_ends']==l['PD_open_ends']==sum(1 if e['kind']=='one_open'else 2 if e['kind']=='two_open'else 0 for e in l['edges'])
        assert{tuple(r['active'])for r in l['rays']if r['kind']=='finite'}=={tuple(r)for r in l['PD_finite_supports']}
    return dict(status='PASS_SAVED_WALL66_IDEAL_AND_BOUNDARY_FLIP',ideal_recipes=208,ideal_branches=400,
        ideal_J_tests=957,rejected_ideal_assemblies=924,flip_recipes=534,flip_branches=1068,rejected_flip_J=233,
        scope='Independent saved source/index/native-value/deletion/disk and separated-PD ledger checks. Affine diagnostics explicitly unresolved. No physical normals, norm roots, corner/PSD classifier, cones, interval signs or ambient-volume replay.')


if __name__=='__main__':print(json.dumps(run(),sort_keys=True))
