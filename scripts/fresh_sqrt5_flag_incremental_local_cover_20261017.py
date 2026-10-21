#!/usr/bin/env python3
"""Explicit singleton-packing repair after flag-driven simplex increases."""
import argparse
from collections import defaultdict
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_flag_incremental_local_cover_20261017'
FLAGS='fresh_sqrt5_flag_refined_base_capacity_20261017'
COVER='fresh_sqrt5_all_local_simplex_packing_cones_20261017'
REVIEW='review_sqrt5_all_local_packing_cover_20261017'
SUBSETS='fresh_sqrt5_all_subset_port_packing_20261017'


def digest(path):
    return sha256((ROOT/path).read_bytes()).hexdigest()


def read(tag):
    data=json.loads((ROOT/f'results/{tag}.json').read_text())
    for p,h in data['sources'].items():
        assert digest(p)==h,p
    return data


def run():
    flag,old,review,sub=map(read,[FLAGS,COVER,REVIEW,SUBSETS])
    assert old['complete_positive_simplex_forest_coverage'] and not old['open_forests']
    assert old['face_multipliers']==sub['face_multipliers']
    old_alpha={i:F(a) for i,a in old['face_multipliers']}
    cert=flag['full_optimized_test']['repaired_dual']
    alpha={i:F(a) for i,a in cert['face_multipliers']}
    assert set(alpha)==set(old_alpha)
    assert all(alpha[i]>=old_alpha[i]>0 for i in alpha)
    oldw,neww=defaultdict(F),defaultdict(F)
    for row in sub['positive_simplex_features']:
        i=row['face']
        for j,m in row['block_orbits']:
            assert m>0
            oldw[j]+=m*old_alpha[i]
            neww[j]+=m*alpha[i]
    columns={r['index']:r for r in sub['literal_block_orbits']}
    assert set(columns)==set(oldw)==set(neww)
    byvertex=defaultdict(list)
    for j,c in columns.items():
        assert c['literal_masks'] and len(c['literal_masks'])==len(set(c['literal_masks']))
        assert all(0<x<1<<13 and x.bit_count()==c['block_size'] for x in c['literal_masks'])
        byvertex[c['vertex']].append(j)
    local={r['vertex']:r for r in old['local_cones']}
    assert set(local)==set(byvertex) and len(local)==740
    original_mu={v:F(a) for v,a in cert['vertex_multipliers']}
    mu=dict(original_mu)
    rows=[]
    terms=added_terms=0
    for v,js in sorted(byvertex.items()):
        rec=local[v]
        assert rec['literal_block_orbits']==[columns[j] for j in sorted(js)]
        assert rec['block_orbit_marginals']==[[j,str(oldw[j])] for j in sorted(js)]
        covered={j:F(0) for j in js}
        mass=F(0)
        for item in rec['certificate']['packing_mixture']:
            a=F(item['mass']);assert a>0
            used=0
            for j,mask in item['blocks']:
                assert j in covered and mask in columns[j]['literal_masks']
                assert not used&mask
                used|=mask
                covered[j]+=a
            mass+=a;terms+=1
        assert all(covered[j]>=oldw[j] for j in js)
        assert mass==F(rec['certificate']['exact_cover_mass'])
        oldmass=mass
        increments=[]
        for j in sorted(js):
            deficit=max(F(0),neww[j]-covered[j])
            if deficit:
                # One literal block is a disjoint packing by itself.  Its
                # full automorphism average has orbit-count marginal one.
                mask=columns[j]['literal_masks'][0]
                increments.append(dict(mass=str(deficit),blocks=[[j,mask]]))
                covered[j]+=deficit;mass+=deficit;added_terms+=1
        assert all(covered[j]>=neww[j] for j in js)
        mu[v]=max(mu.get(v,F(0)),mass)
        rows.append(dict(vertex=v,old_cover_mass=str(oldmass),
            new_required_marginals=[[j,str(neww[j])] for j in sorted(js)],
            final_covered_marginals=[[j,str(covered[j])] for j in sorted(js)],
            singleton_packing_increments=increments,new_cover_mass=str(mass),
            old_flag_vertex_mass=str(original_mu.get(v,F(0))),repaired_vertex_mass=str(mu[v])))
    assert all(mu[v]>=a for v,a in original_mu.items())
    numerator=F(cert['objective_numerator']);denominator=sum(mu.values(),F(0))
    ratio=numerator/denominator
    assert ratio<=F(cert['bound'])
    scale=10**12
    floor=ratio.numerator*scale//ratio.denominator
    interval=[F(floor,scale),F(floor+1,scale)]
    assert interval[0]<=ratio<interval[1]
    # All flag signed inequalities depend only on alpha and hence survive.
    for col in flag['flag_transfer_columns']:
        load=sum((m*alpha.get(i,F(0)) for i,m in col['positive_face_coefficients']),F(0))
        assert load<=col['capacity']*alpha.get(col['simplex_face'],F(0))
    paths=[f'scripts/{TAG}.py',f'scripts/{FLAGS}.py',f'results/{FLAGS}.json',
        f'scripts/{COVER}.py',f'results/{COVER}.json',f'scripts/{REVIEW}.py',
        f'results/{REVIEW}.json',f'results/{SUBSETS}.json']
    return dict(status='PASS_EXACT_FLAG_ALPHA_INCREMENTAL_ALL_LOCAL_PACKING_COVER',
        sources={p:digest(p) for p in paths},optimizer_invocations=0,
        counts=dict(forests=len(rows),block_orbits=len(columns),inherited_packing_terms=terms,
            singleton_packing_increments=added_terms,
            increased_vertex_masses=sum(mu[v]>original_mu.get(v,F(0)) for v in mu)),
        full_flag_bound=cert['bound'],bound=str(ratio),objective_numerator=str(numerator),
        rational_bound_interval=list(map(str,interval)),
        normalizer=str(denominator),face_multipliers=cert['face_multipliers'],
        vertex_multipliers=[[v,str(a)] for v,a in sorted(mu.items())],forests=rows,
        scope='No numerical optimizer or count-vector enumeration. Revalidate every old explicit disjoint-packing mixture, then add one-block packings for each exact uncovered orbit marginal under the increased full flag alpha. The unchanged alpha preserves every flag signed inequality, and increasing vertex masses preserves the inherited full-model dual costs. Averaging these literal packings pays all nonnegative forest/block-orbit charges. The ratio is an explicit, potentially inefficient witness repair; failure to exceed12 is not feasibility or a dimension bound.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    result=run();path=ROOT/f'results/{TAG}.json'
    if a.verify:assert result==json.loads(path.read_text())
    else:path.write_text(json.dumps(result,sort_keys=True,separators=(',',':'))+'\n')
    print(result['status'],result['counts'],result['rational_bound_interval'],flush=True)
