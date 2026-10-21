#!/usr/bin/env python3
"""Extract the exact >12 dual supported only on determinant-flexible triangle stars."""
import argparse,hashlib,json
from fractions import Fraction as F
from pathlib import Path
import fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017 as B
import fresh_sqrt5_rank13_triangle_dual_core_20261017 as D
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_rank13_flexible_triangle_barrier_20261017'
INPUT='fresh_sqrt5_rank13_active_triangle_dual_20261017'
GATE='fresh_sqrt5_rank13_triangle_discriminant_gate_20261017'
def run():
    model=B.read(INPUT);gate=B.read(GATE);B.dual_bound=D.sparse_dual_bound
    restriction=next(r for r in gate['restricted_duals'] if r['discriminant_class'] is None)
    alpha={i:F(a) for i,a in restriction['face_multipliers']}
    records=[r for r in gate['records'] if r['face'] in alpha]
    assert len(records)==len(alpha)==2998
    assert all(not r['fixed_discriminant'] and any(n%2 for n in r['component_sizes']) for r in records)
    bound,beta,tight,num,den=D.sparse_dual_bound(model,alpha)
    assert str(bound)==restriction['bound'] and bound>12
    certificate=dict(kind='EXACT_RESTRICTED_MODEL_DUAL',bound=str(bound),objective_numerator=str(num),normalizer=str(den),
        face_multipliers=[[i,str(a)] for i,a in sorted(alpha.items())],vertex_multipliers=[[i,str(a)] for i,a in sorted(beta.items())],
        tight_orbits={str(i):v for i,v in tight.items()},excludes_target12=True)
    B.verify_certificate(model,certificate)
    paths=[f'scripts/{TAG}.py',f'scripts/{B.TAG}.py',f'scripts/{D.TAG}.py',f'results/{INPUT}.json',f'results/{GATE}.json']
    return dict(status='PASS_EXACT_DETERMINANT_FLEXIBLE_TRIANGLE_BARRIER',
        sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        counts=dict(triangles=len(alpha),positive_vertex_multipliers=len(beta)),certificate=certificate,
        witnesses=[dict(face=r['face'],provenance=r['provenance'],black=r['black'],edges=r['edges'],component_sizes=r['component_sizes']) for r in records],
        scope='Every retained triangle has an odd connected component, hence can match any admissible determinant squareclass by totally positive component scaling. This does not assert common Hasse invariants, one shared quadratic form, simultaneous incidence, or an ambient compact realization. The exact lower bound obstructs only nonnegative whole-forest corner payments on these witnessed local configurations.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=json.loads(json.dumps(run()))
    out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts'],d['certificate']['bound'],flush=True)
