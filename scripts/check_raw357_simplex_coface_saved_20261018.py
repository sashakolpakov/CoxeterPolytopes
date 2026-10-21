#!/usr/bin/env python3
"""Independent saved combinatorial/arithmetic check; no producer imports.

Actual finite incidence and historical catalogue completeness remain inputs.
The geometric injectivity proof is in the report, not a numerical test here.
"""
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,permutations
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_raw357_simplex_coface_transfer_20261018'
d=json.loads((ROOT/f'results/{TAG}.json').read_text())
assert d['status']=='PASS_RAW357_COMPLETE_TYPED_SIMPLEX_COFACE_TRANSFER'
for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
def read(tag):return json.loads((ROOT/f'results/{tag}.json').read_text())
inc=read('fresh_raw357_extra_vertex_incidence_pilot_20261018')
base=read('fresh_sqrt5_prism_simplex_base_incidence_20261017')
forced=read('fresh_sqrt5_forced_base_simplex_domain_20261017')
model=read('fresh_all_marker_families_dual_continuation_20261018')
dual=next(c for c in model['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
alpha={f:F(x)for f,x in dual['face_multipliers']}
geo=d['geometry'];g=geo['native_gram'];right=[2,'0','0']
assert len(g)==10 and all(len(row)==10 for row in g)
assert all(g[i][j]==g[j][i]for i in range(10)for j in range(10))
assert g==inc['geometry']['native_gram']
assert geo['simplex']==inc['geometry']['compact_simplex_faces'][0]
assert geo['source_quad']==[[[1,3],[0,2]]]
assert (geo['deleted_cap'],geo['coface_facet'])==(2,0)
assert geo['coface_whole_A2']==[0,9]and geo['other_whole_A2']==[7,8]
assert geo['simplex']['fixed']==[0,7,8,9]and geo['simplex']['sides']==[1,3,4,5,6]
mp=dict(geo['core_isomorphism']);assert set(mp)=={1,3,4,5,6}and set(mp.values())==set(range(5))
core=base['cores'][1];assert geo['canonical_core']==core
edges={(i,j):e for i,j,e in core['edges']}
def edge(i,j):return edges.get(tuple(sorted((i,j))),right)
assert all(g[i][j]==edge(mp[i],mp[j])for i,j in combinations(mp,2))
auts=[list(p)for p in permutations(range(5))if all(edge(i,j)==edge(p[i],p[j])for i,j in combinations(range(5),2))]
assert geo['core_automorphisms']==auts==[list(range(5))]
for part in([0,9],[7,8]):
    assert g[part[0]][part[1]]==[3,'0','0']
    assert all(g[i][j]==right for i in part for j in set(range(10))-{2}-set(part))
assert sorted(mp[i]for i in(1,3))==geo['distinguished_edge_omissions']==[0,4]
assert geo['canonical_endpoint_omitted']==0 and geo['base_occurrence_multiplicity']==1
tpl=next(t for t in base['prism_templates']if t['canonical_core']==318)
bases=[x for x in tpl['bases']if x['deleted_normal']==2]
assert len(bases)==1 and bases[0]['base_core']==1
assert bases[0]['core_isomorphism']==geo['core_isomorphism']
assert bases[0]['opposite_base_facet']==0
family=inc['complete_family_inventory'];targets=base['base_types'];ctx={r[0]:r for r in base['contexts']}
oldmap=dict(forced['base_to_new_face'])
rows=d['family'];assert len(rows)==len(family)==47
rr={r['source_face']:r for r in rows};assert len(rr)==47 and set(rr)=={r['face']for r in family}
expected_columns={};expected_missing={};typed=set()
for source in family:
    f=source['face'];r=rr[f];old=ctx[source['source_context']]
    assert old[1:4]==[f,318,2]and source['required_quads']==1
    bb=[(t,m)for t,m in old[4]if targets[t]['core']==1]
    assert len(bb)==1 and bb[0][1]==1;t=bb[0][0];typ=targets[t]
    assert t not in typed;typed.add(t)
    sp=sorted(['A2','A2']+source['spectator_forest'])
    assert typ['spectator_forest']==sp==r['passive_forest']
    cap=2*Counter(sp)['A2'];assert cap in(4,6,8)
    fi=typ['old_face']if typ['old_face']is not None else oldmap.get(t)
    assert (r['source_context'],r['base_type'],r['capacity'],r['existing_simplex_face'])==(source['source_context'],t,cap,fi)
    assert r['base_occurrences_per_quad']==r['source_required_quads']==1
    assert r['simplex_required_triangle_payment']==2*10==20
    assert F(r['source_mass'])==F(source['current_mass'])==alpha.get(f,F())
    if fi is None:
        assert r['current_simplex_mass']is None and not alpha.get(f,F())
        expected_missing[t]=(f,sp)
    else:
        assert F(r['current_simplex_mass'])==alpha.get(fi,F())
        expected_columns[f]=(fi,t,cap,alpha.get(f,F())-cap*alpha.get(fi,F()))
cols=d['existing_target_transfer_columns'];missing=d['missing_simplex_rows']
assert len(cols)==len(expected_columns)==16 and len({c['source_face']for c in cols})==16
assert len({c['simplex_face']for c in cols})==16
for c in cols:
    f=c['source_face'];fi,t,cap,z=expected_columns[f]
    assert(c['simplex_face'],c['base_type'],c['capacity'])==(fi,t,cap)
    assert c['integer_face_coefficients']==[[f,1],[fi,-cap]]
    assert c['integer_vertex_capacities']==[]and c['vertex_cost']==0
    assert F(c['current_exact_dual_load'])==z and c['strict_current_violation']==(z>0)
assert len(missing)==len(expected_missing)==31 and len({r['base_type']for r in missing})==31
for r in missing:assert (r['source_face'],r['passive_forest'])==expected_missing[r['base_type']]and r['required']==20
violations={f:z for f,(fi,t,cap,z)in expected_columns.items()if z>0}
assert set(violations)=={70404,70406,70408,70414,70420,70422}
assert sum(violations.values(),F())==F(d['total_positive_load'])==F(37155891,10**10)
for r in d['inherited_capacity_comparisons']:
    old=[z for z in forced['new_shared_transfer_columns']if z['simplex_face']==r['simplex_face']]
    assert len(old)==1 and dict(old[0]['positive_quad_coefficients'])[r['source_face']]==2
    assert old[0]['negative_simplex_coefficient']==[r['simplex_face'],-45]
    assert(r['old_source_coefficient'],r['old_simplex_capacity'],r['new_source_coefficient'])==(2,45,1)
    assert r['new_simplex_capacity']==rr[r['source_face']]['capacity']
assert not d['model_changed']and not d['full_matrix_reconstructed']and d['optimizer_invocations']==0
print('PASS_RAW357_SAVED_COFACE_TRANSFER',len(rows),len(cols),len(missing),len(violations),flush=True)
