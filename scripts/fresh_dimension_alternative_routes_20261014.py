#!/usr/bin/env python3
"""Small exact obstruction to a coupled H3-face lift of the known H6 source."""
import argparse,json
from hashlib import sha256
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_dimension_alternative_routes_20261014'
def run():
    phi=(1+s.sqrt(5))/2
    B=s.Matrix([[2,-phi,0],[-phi,2,-1],[0,-1,2]])
    J=B.inv().applyfunc(s.simplify)
    diagonal=[(9+3*s.sqrt(5))/4,(3+s.sqrt(5)),(5+s.sqrt(5))/4]
    assert list(J.diagonal())==diagonal
    assert all(x.is_positive for x in J)
    kappa=s.simplify(2-J[2,2]); assert kappa==(3-s.sqrt(5))/4
    assert (2-J[0,0]).is_negative and (2-J[1,1]).is_negative
    assert (2-2*J[2,2]).is_negative
    assert (2*kappa-1).is_negative
    assert (J[2,2]-kappa).is_positive
    # Literal finite nonright edges of Burcroff's Fig. 5 in the public
    # source's 1-based numbering. No dotted distances are needed here.
    edges=[(1,2,5),(2,3,3),(3,4,3),(4,5,3),(5,6,3),
           (9,10,5),(2,7,4),(2,8,4),(5,9,3),(6,10,5)]
    witnesses={str(i):next([a,b,m] for a,b,m in edges if i in (a,b)) for i in range(1,11)}
    assert len(witnesses)==10
    def enc(M):return [[str(s.simplify(x)) for x in row] for row in M.tolist()]
    p=f'scripts/{TAG}.py'
    return dict(status='PASS_H3_FACE_LIFT_ORTHOGONALITY_GATE',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()},
        H3_gram=enc(B),H3_inverse=enc(J),positive_residual=str(kappa),
        sole_nonorthogonal_face_profile=dict(pairings=[0,0,-1],label=3,residual_square=str(kappa)),
        all_label_proof='Positive inverse: two nonzero contacts or any contact at positions0,1 costs more than2. At position2 every label>=4 costs at least2*J22>2. Hence only terminal label3 can induce a face facet.',
        special_to_orthogonal_squared_cosine_lower=str(s.simplify(1/(2*kappa))),
        special_to_special_cosine_lower=str(s.simplify(J[2,2]/kappa)),
        source_nonright_edges=edges,source_every_facet_nonright_neighbor=witnesses,
        conclusion='Every facet of an H3-supported Coxeter face is orthogonal to H3 if its intrinsic chamber has no facet whose meeting neighbors are all right. In particular this holds for the displayed known compact H6 ten-facet source.',
        scope='Excludes changing H3 contacts of source-face facets while preserving this intrinsic source. Does not exclude outer caps disjoint from the face, changed intrinsic metrics, other sources, or a compact H9 completion.',
        primary_source='https://arxiv.org/pdf/2201.03437, Figure5 and Section8',
        machine_readable_diagram_source='https://raw.githubusercontent.com/johnamcleod/10-facets-6-d-polytopes/main/verify_polytope.py')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    out=ROOT/f'results/{TAG}.json';d=json.loads(json.dumps(run()))
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'])
