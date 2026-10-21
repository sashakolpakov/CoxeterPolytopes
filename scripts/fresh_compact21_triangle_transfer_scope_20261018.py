#!/usr/bin/env python3
"""Exact scope diagnostic for the new exceptional triangle constraint.

Literal weak duality only: no optimizer, geometry census, or new field bound.
"""
import argparse,json
from collections import defaultdict
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_compact21_triangle_transfer_scope_20261018'
TRI='fresh_compact29_three_triangle_star_constraint_20261018'
MODEL='fresh_sqrt5_dimension21_e7e8_context_20261016'
FIELD='fresh_compact21_single_field_20261016'

def run():
    triangle=json.loads((ROOT/f'results/{TRI}.json').read_text())
    d=json.loads((ROOT/f'results/{MODEL}_cap2.json').read_text())
    assert triangle['status']=='PASS_COMPACT29_UNIVERSAL_THREE_TRIANGLE_STAR_CONSTRAINT'
    assert triangle['triangle_factor_capacity']==1 and triangle['designated_triangle_upper']==3
    checked=0
    for doc in(triangle,d):
        for p,h in doc['sources'].items():
            assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p;checked+=1
    assert d['dimension']==21 and F(d['target'])==40
    assert d['certificate']['kind']=='exact_model_obstruction'
    z=d['variables'];coeff=defaultdict(F);value=F();forests=[];face_terms=0
    terms=d['certificate']['terms'];assert len(terms)==2490
    for i,t in enumerate(terms):
        y=F(t['multiplier']);assert y<0
        value+=y*F(t['rhs'])
        ids=[j for j,a in t['coefficients']];assert len(ids)==len(set(ids))
        for j,a in t['coefficients']:
            assert 0<=j<=z;coeff[j]+=y*F(a)
        m=t['metadata']
        if m['kind']=='forest':
            ts=sorted(m['types']);rank=lambda s:2 if s.startswith('I2(')else int(s[1:])
            assert sum(map(rank,ts))==21
            forests.append(dict(term=i,types=ts,positive_mass=str(-y)))
        else:assert m['kind']=='face';face_terms+=1
    assert all(v<=F(j==z)for j,v in coeff.items())
    assert value==F(d['certificate']['lower_bound'])==F(4119593261,102300000)>40
    target=[r for r in forests if r['types']==['E7']*3];assert not target
    solutions=[[a,b]for a in range(4)for b in range(3)if 7*a+8*b==21]
    assert solutions==[[3,0]]
    top=sorted(forests,key=lambda r:(-F(r['positive_mass']),r['term']))[:10]
    assert top[0]['types']==['A1','D5','E7','E8']
    assert any(r['types']==['D5','E8','E8']for r in top)
    paths=[f'scripts/{TAG}.py']
    for s,suffix in((TRI,''),(MODEL,'_cap2'),(FIELD,'')):
        paths += [f'scripts/{s}.py',f'results/{s}{suffix}.json',f'results/{s}.md']
    paths += ['results/fresh_sqrt5_dimension21_diagnostic_synthesis_20261016.md',
              'results/fresh_compact_global_bound_20261015.md']
    return dict(status='PASS_COMPACT21_TRIANGLE_ONLY_ALLOWANCE_SCOPE_DIAGNOSTIC',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        inherited_source_bindings_checked=checked,dimension=21,target='40',
        all_E7_E8_rank21_forests=[['E7']*3],all_E7_E8_rank_equation_solutions=solutions,
        exact_dual=dict(source_model=MODEL,lower_bound=str(value),excess=str(value-40),
            literal_terms=len(terms),face_terms=face_terms,forest_terms=len(forests),
            literal_columns_touched=len(coeff),nonzero_column_sums=sum(v!=0 for v in coeff.values()),
            normalization_mass=str(sum((F(r['positive_mass'])for r in forests),F())),
            target_forest_terms=target,target_forest_mass='0',largest_forest_multipliers=top,
            all_literal_weak_duality_inequalities=True),
        triangle_scope=dict(rank21_geometry_extension='The same proof applies to E7^3 by replacing positive index29 with21.',
            inherited_geometry_not_reexecuted=True,only_E7_E8_forests_transferred=True,
            no_bound_on_triangles_with_other_old_components=True),
        conclusion=dict(arbitrary_target_forest_allowance_only_change_cannot_break_this_dual=True,
            new_incidence_columns_tested=False,new_angle_payment_rule=False,new_dimension_bound=None),
        scope='Read-only literal dual arithmetic and forest-support diagnostic. No new optimization, complete face-coefficient reconstruction, geometric census, or nonexistence claim. Changing constraints beyond the E7^3 vertex row is outside this obstruction.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    p=ROOT/f'results/{TAG}.json';assert a.verify or not p.exists();d=run()
    if a.verify:assert json.loads(p.read_text())==d
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['exact_dual']['lower_bound'],d['exact_dual']['target_forest_mass'],flush=True)
