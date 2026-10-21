#!/usr/bin/env python3
"""Independent saved-record accounting; no project/numerical-package imports."""
import json,math
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_wall66_single_edge_conic_20261018'

def native(x):
    tol=2e-7
    if x<-2-tol:return 'dotted'
    if abs(x+2)<=tol:return 'parallel_or_unresolved'
    if abs(x)<=tol:return 2
    if x>0 or x<-2:return None
    m=round(math.pi/math.acos(-x/2))
    if m<3:return None
    th=lambda k:2*math.cos(math.pi/k)
    space=min(abs(th(m)-th(m-1)),abs(th(m+1)-th(m)))
    return m if abs(x+th(m))<=min(tol,space/100)else None

def pair(a,b):return-a[0]*b[0]+sum(x*y for x,y in zip(a[1:],b[1:]))

def ledger(l):
    assert not l['warnings']and l['open_ends']==l['PD_open_ends']
    assert l['PD_status']=='numerically_separated_acute_incidence'
    assert Counter(r['kind']for r in l['rays'])==Counter(l['counts'])
    assert{tuple(sorted(r['active']))for r in l['rays']if r['kind']=='finite'}=={tuple(sorted(r))for r in l['PD_finite_supports']}
    assert sum(1 if e['kind']=='one_open'else 2 if e['kind']=='two_open'else 0 for e in l['edges'])==l['open_ends']

def topology(t):
    fs=t['triangles'];edges=Counter(tuple(sorted(e))for f in fs for e in combinations(f,2));vs=set(sum(fs,[]));bc=t['boundary']
    assert(len(vs),len(edges),len(fs))==(9,18,10)
    assert len(vs)-len(edges)+len(fs)==1 and set(edges.values())=={1,2}
    assert{e for e,n in edges.items()if n==1}=={tuple(sorted(e))for e in zip(bc,bc[1:]+bc[:1])}
    for v in vs:
        le=[tuple(i for i in f if i!=v)for f in fs if v in f];degree=Counter(i for e in le for i in e);seen={next(iter(degree))}
        while True:
            new=seen|{b for a,b in le if a in seen}|{a for a,b in le if b in seen}
            if new==seen:break
            seen=new
        assert seen==set(degree)
        assert sorted(degree.values())==([1,1]+[2]*(len(degree)-2)if v in bc else[2]*len(degree))

def run():
    path=ROOT/f'results/{TAG}.json';d=json.loads(path.read_text());n=Counter()
    for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h;n['bindings']+=1
    assert d['retained_old_walls']==63 and d['deleted_walls']==[]
    assert d['exact_deletion_debit']==dict(finite_supports=0,compact_faces=0,wall8_supports=0)
    # Bind the existing exact label9 declaration; this does not reconstruct it.
    old=json.loads((ROOT/'results/fresh_nonarithmetic_pi8_quad_completion_20261017.json').read_text())
    old59=next(r for r in old['new_contact_checks']if r['wall']==59)
    c9=next(r for r in old59['contacts']if r['wall']==58)
    assert c9['label']==9 and 'equality'in c9
    assert d['symbolic_finite_corner_domain']==[[0,2,2],[0,3,2],[0,4,2],[0,5,2],[0,6,2],[1,2,2]]
    assert d['high_anchor_label_samples']==[6,8,12,24]
    assert len(d['sheets'])==18 and len(d['caps'])==78
    assert d['counts']=={'event_equations':7920,'positive_parameter_tests':7925,'native_or_center_rejected':7847}
    assert sum(x['count']for x in d['blockers'])==7847
    for l in d['boundary_ledgers']:ledger(l);n['ledgers']+=1
    base=d['boundary_ledgers'][d['control_score']]
    assert base['counts']==dict(finite=26,exterior=6)and base['open_ends']==10
    for i,s in enumerate(d['sheets']):
        assert s['index']==i and s['remaining_norm']>0
        assert s['sheet']==(-1 if i%2==0 else 1)
        assert len(s['coefficients'])==55
        assert s['acute_intervals']==[] if s['sheet']==-1 else s['acute_intervals']
        n['conic_sheets']+=1
    indices=[]
    for s in d['sheets']:
        indices.extend(e['cap']for e in s['events']if 'cap'in e)
    assert sorted(indices)==list(range(78))
    for i,c in enumerate(d['caps']):
        assert c['index']==i and c['parameter']>0 and c['center_product']<0 and c['norm_error']<2e-6
        s=d['sheets'][c['sheet_record']];t=c['parameter']
        v=[a+b*t+cc/t for a,b,cc in zip(s['part'],s['va'],s['vb'])]
        assert max(abs(a-b)for a,b in zip(v,c['normal']))<1e-8
        q=s['remaining_norm']+pair(s['part'],s['part'])
        assert abs(pair(v,v)-q)<2e-6
        values=dict(c['native_products']);assert len(values)==57
        assert all(native(x)not in(None,'parallel_or_unresolved')for x in values.values())
        finite=[[k,native(x)]for k,x in sorted(values.items())if isinstance(native(x),int)]
        assert finite==c['finite_old_contacts']
        assert all(k in[58,59,60]for k,m in finite)
        assert values[59]<=2e-7 and abs(values[60])<=2e-7
        if c['event']['kind']=='native_event':assert c['event']['wall']==58
        n[c['event']['kind']]+=1;n['caps']+=1
    assert(n['dotted_interval'],n['native_event'])==(9,69)
    q=d['companion_experiment'];topology(q['topology'])
    assert q['selected_seed_indices']==[1,9,10]and q['native58_59_label']==9 and native(q['native58_59_product'])==9
    for r in q['records']:
        seed=d['caps'][r['seed']];assert seed['finite_old_contacts']==r['contacts']and seed['profile']==0
        assert seed['score']==r['seed_score'];l=d['boundary_ledgers'][r['seed_score']]
        edges={tuple(sorted(e['sides'])):e['kind']for e in l['edges']};bc=q['topology']['boundary']
        ce={tuple(sorted(e))for e in zip(bc,bc[1:]+bc[:1])}
        assert all(edges[e]=='one_open'for e in ce)
        assert{e for e,k in edges.items()if k in['one_open','two_open']}-ce=={(7,30),(7,47),(30,52),(47,52)}
        rr=[z for z in l['rays']if 60 in z['active']]
        assert all(z['kind']=='finite'for z in rr)and{tuple(sorted(z['active']))for z in rr}=={(58,59,60),(58,60,67),(59,60,67)}
        for bi,b in enumerate(r['banks']):
            assert b['frame']==q['topology']['frames'][bi]and b['duties']==[[0,1],[1,2]]
            for ci,c in enumerate(b['caps']):
                assert c['index']==ci and c['center_product']<0 and c['norm_residual']<2e-5
                assert len(c['native_products'])==58 and all(native(x)not in(None,'parallel_or_unresolved')for j,x in c['native_products'])
                n['companion_caps']+=1
            n['banks']+=1
        expected_tests=sum(len(r['banks'][a]['caps'])*len(r['banks'][b]['caps'])for a,b in[(0,1),(1,2),(0,2)])
        assert sum(r['pair_histogram'].values())==expected_tests;n['pair_tests']+=expected_tests
        for p in r['pairs']:
            a,b=p['sides'];seen=set()
            for rec in p['records']:
                ia,ib=rec['indices'];assert(ia,ib)not in seen;seen.add((ia,ib))
                assert ia<len(r['banks'][a]['caps'])and ib<len(r['banks'][b]['caps'])
                assert native(rec['native_product'])==rec['relation']
                assert rec['relation']not in(None,'parallel_or_unresolved')
                n['compatible_pairs']+=1
        assert not r['triples']
    assert n['banks']==9 and n['pair_tests']==q['counts']['companion_pair_tests']
    assert not d['unresolved']and not d['improvements']and not q['improvements']
    assert all(l['open_ends']>=10 for l in d['boundary_ledgers'])
    return dict(status='PASS_SAVED_WALL66_SINGLE_EDGE_CONIC',counts=dict(n),source_hash=sha256(Path(__file__).read_bytes()).hexdigest(),certificate_hash=sha256(path.read_bytes()).hexdigest(),scope='Independent standard-library saved-record checks: bindings, exact inherited label9 declaration, zero deletion, conic coefficients/norm residuals and native-value gates, selected actual boundary and old60 triangle, disk links, banks/pair accounting and saved cone/PD equality. No physical normal reconstruction, root or cone rerun, interval proof, full finite classifier or ambient finite-volume certification.')

if __name__=='__main__':print(json.dumps(run(),sort_keys=True))
