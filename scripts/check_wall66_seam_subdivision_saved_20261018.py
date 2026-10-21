#!/usr/bin/env python3
"""Stdlib saved-domain/ledger checks; no physical or corner-classifier imports."""
import json,math
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_wall66_seam_subdivision_20261018'

def native(x):
    tol=2e-7
    if x<-2-tol:return 'dotted'
    if abs(x+2)<=tol:return 'parallel_or_unresolved'
    if abs(x)<=tol:return 2
    if x>0 or x<-2:return None
    m=round(math.pi/math.acos(-x/2))
    if m<3:return None
    theta=lambda k:2*math.cos(math.pi/k)
    spacing=min(abs(theta(m)-theta(m-1)),abs(theta(m+1)-theta(m)))
    return m if abs(x+theta(m))<=min(tol,spacing/100)else None

def ledger(l):
    assert not l['warnings']and l['PD_status']=='numerically_separated_acute_incidence'
    assert l['open_ends']==l['PD_open_ends']
    assert Counter(r['kind']for r in l['rays'])==Counter(l['counts'])
    assert {tuple(sorted(r['active']))for r in l['rays']if r['kind']=='finite'}=={tuple(sorted(s))for s in l['PD_finite_supports']}
    assert l['open_ends']==sum(1 if e['kind']=='one_open'else 2 if e['kind']=='two_open'else 0 for e in l['edges'])

def topcheck(t):
    fs=t['triangles'];bc=t['boundary'];ee=Counter(tuple(sorted(e))for f in fs for e in combinations(f,2));vs=set(sum(fs,[]))
    assert(len(vs),len(ee),len(fs))==(9,19,11)and set(ee.values())=={1,2}and(68,69)not in ee
    assert {e for e,n in ee.items()if n==1}=={tuple(sorted(e))for e in zip(bc,bc[1:]+bc[:1])}
    for v in vs:
        links=[tuple(x for x in f if x!=v)for f in fs if v in f];dd=Counter(x for e in links for x in e);seen={next(iter(dd))}
        while True:
            ns=seen|{b for a,b in links if a in seen}|{a for a,b in links if b in seen}
            if ns==seen:break
            seen=ns
        assert seen==set(dd)and sorted(dd.values())==([1,1]+[2]*(len(dd)-2)if v in bc else[2]*len(dd))
    assert t['caps']==[68,69,70,71]

def run():
    path=ROOT/f'results/{TAG}.json';d=json.loads(path.read_text());old=json.loads((ROOT/'results/fresh_wall59_low_bond_companions_20261018.json').read_text());n=Counter()
    for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h;n['source_bindings']+=1
    assert d['exact_deletion_debit']==old['exact_deletion_debit']
    debit=d['exact_deletion_debit'];assert debit['deleted_walls']==[59,60,61]and debit['lost_finite_count']==7
    assert debit['compact_faces_by_dimension']=={'2':6,'3':1}and all(8 not in v for v in debit['lost_finite_supports'])
    primary=json.loads((ROOT/'results/fresh_nonarithmetic_release49_62_replacement_20261017.json').read_text())
    vs=primary['retained_old_finite_supports']+[r['support']for r in primary['new_finite_vertices']]
    assert len(vs)==199 and [v for v in vs if set(v)&{59,60,61}]==debit['lost_finite_supports']
    topcheck(d['topology']);groups={};eligible=0
    for r in old['compatible_pairs']:
        seed=old['seeds'][r['seed']];x=seed['banks'][1]['caps'][r['x']];y=seed['banks'][3]['caps'][r['y']]
        if max(x['labels']+y['labels'])>5:continue
        eligible+=1;key=(r['seed'],x['labels'][2],y['labels'][0])
        if key not in groups or r['native_product']>groups[key]['native_product']:groups[key]=r
    assert eligible==d['eligible_low_prefix_pairs']==3808
    assert d['selected_prefixes']==[dict(group=list(k),pair=r)for k,r in sorted(groups.items())]
    assert len(d['records'])==len(d['selected_prefixes'])==16
    for r in d['old66_contact_coverage']:
        b=old['seeds'][r['seed']]['banks'][r['bank']];j=b['frame'].index(66)
        assert r['native66_label_counts']=={str(k):v for k,v in Counter(c['labels'][j]for c in b['caps']).items()}
    for l in d['boundary_ledgers']:ledger(l);n['saved_ledgers']+=1
    def bankcheck(b,phase):
        seen=set();claimed=[]
        assert b['complete_low_recipe_count']==len(b['recipes'])
        for i,r in enumerate(b['recipes']):
            assert r['index']==i and 0<=r['profile']<7 and len(r['labels'])==3 and all(m in[2,3,4,5]for m in r['labels'])
            key=(r['profile'],*r['labels']);assert key not in seen;seen.add(key)
            assert r['norm']['status']=='two_real_norm_branches'
            assert sorted(x['branch']for x in r['branches'])==[-1,1]
            n['norm_equations']+=1
            for br in r['branches']:
                n['norm_branches']+=1;assert br['norm_residual']<=2e-5
                assert all(native(v)in(None,'parallel_or_unresolved')for j,v in br['blockers'])
                if 'cap'in br:
                    ci=br['cap'];claimed.append(ci);c=b['caps'][ci]
                    assert c['recipe']==i and c['profile']==r['profile']and c['labels']==r['labels']and c['branch']==br['branch']
                    assert c['center_product']==br['center_product']<0 and not br['blockers']
                    assert c['norm_residual']==br['norm_residual']
        assert sorted(claimed)==list(range(len(b['caps'])))
        for i,c in enumerate(b['caps']):
            assert c['index']==i and c['center_product']<0
            assert len(c['native_products'])==(57 if phase=='J' else 58)
            assert all(native(v)not in(None,'parallel_or_unresolved')for j,v in c['native_products'])
            l=d['boundary_ledgers'][c['score']];assert l['counts']==dict(finite=31 if phase=='J' else 33,exterior=5)and l['open_ends']==9
            assert {(7,30,47),(30,47,52)}<={tuple(sorted(r['active']))for r in l['rays']if r['kind']=='exterior'}
            n[phase+'_caps']+=1
        n[phase+'_banks']+=1
    for i,r in enumerate(d['records']):
        assert r['selection']==d['selected_prefixes'][i]
        assert r['J_bank']['frame']==[66,68,69]and r['J_bank']['duties']==[[0,1],[0,2]]
        bankcheck(r['J_bank'],'J');l=d['boundary_ledgers'][r['prefix_score']]
        assert l['counts']==dict(finite=29,exterior=5)and l['open_ends']==9
    zs=set();expectedz={(i,j)for i,r in enumerate(d['records'])for j in range(len(r['J_bank']['caps']))};allcaps=set()
    for i,r in enumerate(d['Z_banks']):
        key=(r['prefix'],r['J']);assert key not in zs;zs.add(key)
        assert r['bank']['frame']==[30,67,68]and r['bank']['duties']==[[0,1],[1,2]]
        bankcheck(r['bank'],'Z');allcaps|={(i,j)for j in range(len(r['bank']['caps']))}
    assert zs==expectedz and len(zs)==95
    tested=set()
    for r in d['completion_tests']:
        key=(r['Z_bank'],r['Z']);assert key not in tested;tested.add(key)
        c=d['Z_banks'][key[0]]['bank']['caps'][key[1]];vals=dict(c['native_products'])
        assert r['score']==c['score']and native(vals[69])==native(vals[70])=='dotted'
        assert r['missing_finite_edges']==[[69,71],[70,71]]and not r['whole_cap_ray_sets']
        expected=[f for f in d['topology']['triangles']if 71 in f and(69 in f or 70 in f)]
        assert r['nonfinite_required_corners']==expected
    assert tested==allcaps and len(tested)==2685 and not d['completions']and not d['unresolved']
    assert n['norm_equations']==d['counts']['norm_equations']==4049 and n['norm_branches']==d['counts']['norm_branches']==8098
    assert n['J_caps']==95 and n['Z_caps']==2685
    return dict(status='PASS_SAVED_WALL66_SEAM_SUBDIVISION',counts=dict(n),selected_prefixes=16,eligible_frozen_pairs=3808,complete_final_tests=len(tested),
        certificate_hash=sha256(path.read_bytes()).hexdigest(),source_hash=sha256(Path(__file__).read_bytes()).hexdigest(),
        scope='Independent stdlib frozen-input selection/oldsupport debit, disk/link, saved recipe/branch/cap index coverage, native value and final missing-edge checks, full saved cone/PD ledger agreement. No physical frame/norm-root/cone reconstruction, spherical classifier replay or exact newnormal/volume proof.')

if __name__=='__main__':print(json.dumps(run(),sort_keys=True))
