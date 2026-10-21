#!/usr/bin/env python3
"""No-producer saved-data incidence, boundary and symbolic tangent check."""
import ast,json,re,hashlib
from pathlib import Path
from itertools import combinations
from collections import Counter,defaultdict
R=Path(__file__).resolve().parents[1];tag='fresh_release55_59_61_prism_replacement_20261017'
d=json.loads((R/f'results/{tag}.json').read_text())
for p,h in d['sources'].items():assert hashlib.sha256((R/p).read_bytes()).hexdigest()==h,p
ast.parse((R/f'scripts/{tag}.py').read_text())
links=re.findall(r'\]\(([^)]+)\)',(R/f'results/{tag}.md').read_text())
for p in links:assert (R/'results'/p).exists(),p
p=json.loads((R/'results/fresh_nonarithmetic_release49_62_replacement_20261017.json').read_text())
a=json.loads((R/'results/fresh_release58_66_anchor_completion_gate_20261017.json').read_text())
vs=[tuple(v)for v in p['retained_old_finite_supports']]+[tuple(v['support'])for v in p['new_finite_vertices']]
avs=[v for v in vs if not set(v)&{58,66}]+[tuple(v['support'])for v in a['new_finite_vertices']]
assert len(vs)==199 and len(avs)==192
# Rebuild face incidences independently from the saved complete finite supports.
def compact(vs):
    ed=Counter(e for v in vs for e in combinations(v,8));faces=defaultdict(set)
    for i,v in enumerate(vs):
        for k in range(8):
            for w in combinations(v,k):faces[w].add(i)
    return {w:ids for w,ids in faces.items()if len(ids)>=10-len(w)and all(ed[tuple(x for x in vs[i]if x!=j)]==2 for i in ids for j in vs[i]if j not in w)}
for v,cname in [(vs,'relative_control'),(avs,'relative_anchor')]:
    cs=compact(v)
    for row in d['exact_deletion_comparisons']:
        rec=row[cname];lost={i for i,s in enumerate(v)if set(s)&set(rec['deleted'])}
        assert sorted(v[i]for i in lost)==sorted(map(tuple,rec['lost_finite_supports']))
        touched={w:ids for w,ids in cs.items()if ids&lost}
        assert Counter(9-len(w)for w in touched)=={int(k):x for k,x in rec['compact_faces_by_dimension'].items()}
        assert len(lost)==rec['lost_finite_count']
        assert sum(8 in v[i]for i in lost)==rec['lost_supports_containing8']
        assert set(touched)=={tuple(r['support'])for r in rec['compact_faces_touched']}
        for rr in rec['compact_faces_touched']:assert sorted(v[i]for i in touched[tuple(rr['support'])])==sorted(map(tuple,rr['finite_supports']))
assert d['counts']=={'norm_equations':1778,'norm_branches':3012,'native_caps':298}
assert len(d['records'])==8 and sum(len(r['continuations'])for r in d['records'])==290
assert sum(len(r['coherent_pairs'])for r in d['records'])==80
assert sum(len(r['pair_pool'])*(len(r['pair_pool'])-1)//2 for r in d['records'])==1659
assert d['pair_relation_histogram']=={'None':1579,'dotted':80}
assert not d['improvements'] and len(d['boundary_ledgers'])==11
for b in d['boundary_ledgers']:
    assert not b['warnings'] and b['PD_status']=='numerically_separated_acute_incidence'
    assert {tuple(r['active'])for r in b['rays']if r['kind']=='finite'}==set(map(tuple,b['PD_finite_supports']))
    assert Counter(r['kind']for r in b['rays'])==b['counts']
    assert sum({'one_open':1,'two_open':2}.get(e['kind'],0)for e in b['edges'])==b['open_ends']==b['PD_open_ends']
for rec in d['records']:
    scores={rec['B']['score'],*[c['score']for c in rec['continuations']],*[p['score']for p in rec['coherent_pairs']]}
    for i in scores:
        b=d['boundary_ledgers'][i];assert b['open_ends']==11 and b['counts']['exterior']==7
        assert {(7,30,47),(30,47,52)}<={tuple(r['active'])for r in b['rays']if r['kind']=='exterior'}
    assert len(rec['derived_frames'])==6 and not rec['B']['numerical_deleted_wall_matches']
    for pair in rec['coherent_pairs']:assert pair['relation']=='dotted' and pair['native_product']<-2
for r in d['unresolved']:
    assert r['norm']['status']=='unresolved_tangent_norm'
    assert r['B_context']['index']in[1,3,4]
    assert r['projected_vector_diagnostic']['center_product']>50
    assert r['projected_vector_diagnostic']['blockers']
assert len(d['unresolved'])==5
print('PASS saved incidence debit /8 B/290 C/80 dotted pairs/11 ledgers/5 raw tangent flags; no numerical roots or cones replayed;',len(d['sources']),'hashes/AST',len(links),'links')
import sympy as sp
assert len(d['exact_tangent_gates'])==5 and not d['remaining_unresolved']
for t,u in zip(d['exact_tangent_gates'],d['unresolved']):
    assert t['B_context']==u['B_context']and t['frame']==u['frame']and t['labels']==u['labels']
    g=sp.Matrix([[sp.sympify(x)for x in r]for r in t['gram']]);b=sp.Matrix([sp.sympify(x)for x in t['rhs']]);c=sp.Matrix([sp.sympify(x)for x in t['projection_coefficients']])
    assert (g*c-b).applyfunc(sp.simplify)==sp.zeros(3,1)
    assert sp.simplify((c.T*g*c)[0])==2
    assert all(x<=0 for x in c)and any(x<0 for x in c)
    assert all(i!=60 or x==0 for i,x in zip(t['frame'],c))
    assert [sp.simplify(g[:k,:k].det())for k in(1,2,3)]==[sp.sympify(x)for x in t['positive_principal_minors']]
    assert all(sp.sympify(x)>0 for x in t['positive_principal_minors'])
print('PASS all5 symbolic tangent certificates; no producer imports; remaining_unresolved0')
