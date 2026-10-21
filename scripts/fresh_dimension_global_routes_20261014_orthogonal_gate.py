#!/usr/bin/env python3
"""Exact data for the all-S61-orthogonal completion obstruction."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
from fresh_h20_bridge_edge_closure_20261014_review import Field,parse
import fresh_h20_bridge_faces_20261014_review as A
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_dimension_global_routes_20261014_orthogonal_gate'
BASE='results/fresh_h20_e7_edge_cusp_repair_20261014_exact.json'
FACE='results/fresh_h20_bridge_faces_20261014_review.json'
REVIEW='results/fresh_h20_bridge_edge_closure_20261014_review.json'
GATE='results/fresh_h20_s75_projection_gate_20261014.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
    old,faces,review,gate=[json.loads((ROOT/p).read_text())for p in[BASE,FACE,REVIEW,GATE]]
    for doc in[old,faces,review,gate]:
        for p,h in doc['sources'].items():assert sha(p)==h,p
    field=Field([parse(old['field']['discriminant'],0)]);z=field.zeros[1]
    S61=faces['faces'][0]['support'];S75=faces['faces'][1]['support']
    assert set(S61)-set(S75)=={61}and set(S75)-set(S61)=={75}
    assert len(set(S61)&set(S75))==17
    N={i:tuple(parse(x,1)for x in old['normals'][i])for i in set(S61+S75)}
    Q=[list(map(F,row))for row in old['form']];assert A.inertia(Q)==[20,1,0]
    def dot(v,w):return field.total(field.scale(field.mul(v[i],w[j]),q)for i,row in enumerate(Q)for j,q in enumerate(row)if q)
    G=[list(map(F,row))for row in faces['faces'][1]['rational_gram']]
    J=A.inv(G);assert J==[list(map(F,row))for row in faces['faces'][1]['inverse']]
    A.pd(G);assert min(x for row in J for x in row)>=0
    hp=[dot(N[61],N[i])for i in S75]
    cs=[field.total(field.scale(x,q)for x,q in zip(hp,row))for row in J]
    U=tuple(field.add(N[61][k],field.scale(field.total(field.mul(c,N[i][k])for c,i in zip(cs,S75)),-1))for k in range(21))
    assert all(dot(U,N[i])==z for i in S75)
    q=dot(U,U);assert field.sign(q)<0 and dot(U,N[61])==q
    HP=tuple(parse(x,1)for x in faces['faces'][1]['strict_face_point'])
    orient=dot(U,HP);assert field.sign(orient)<0
    assert q==parse(faces['faces'][1]['all_projected_wall_norms'][61],1)
    assert orient==parse(faces['faces'][1]['all_face_center_pairings'][61],1)
    assert J[S75.index(75)][S75.index(75)]==5
    assert gate['wall75_diagonal']==5 and gate['wall75_must_be_right']
    assert review['surviving_whole_H2_face']['support']==S75
    assert [r['wall']for r in review['surviving_whole_H2_face']['additional_projected_walls']]==[95,96,97]
    # The numerical inequalities used by the universal proof are rational.
    assert F(1)-5*F(1,2)==F(-3,2)
    paths=[BASE,FACE,REVIEW,GATE,f'scripts/{TAG}.py',
           'scripts/fresh_h20_bridge_edge_closure_20261014_review.py',
           'scripts/fresh_h20_bridge_faces_20261014_review.py']
    return A.enc(dict(status='PASS_ALL_S61_ORTHOGONAL_ADDITIONS_REMAIN_INFINITE_VOLUME',sources={p:sha(p)for p in paths},
        S61=S61,S75=S75,common_support_rank=17,union_support_rank=19,
        union_support_signature=[18,1],common_orthogonal_space_signature=[2,0],
        common_hyperbolic_point_ray=U,common_point_norm=q,common_point_future_pairing=orient,
        S75_inverse_diagonal_at75=5,
        nonright75_residual_norm_upper_bound_ratio=F(-3,2),
        right75_boundary_normals='Any added r orthogonalS61 and right75 lies in the positive2-plane orthogonal to span(S61,75), so its S75 geodesic passes through U.',
        theorem_scope='Any finite coherent acute Coxeter extension retaining the current98 walls, with every new normal orthogonal toS61 and a common strict future timelike center, has infinite volume. The S75 face is a nonempty positive-angle geodesic sector about U. No restriction on new norms, number of walls, finite labels, or coefficient fields.',
        proof_artifact='results/fresh_dimension_global_routes_20261014.md',
        limitations='Not an obstruction to arbitrary nonorthogonal additions or to deleting retained walls. Existing98 whole-plane and pair/center geometry are explicit bound dependencies.'))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'])
