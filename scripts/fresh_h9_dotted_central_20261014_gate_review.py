#!/usr/bin/env python3
"""Independent literal-incidence and scalar-identity review of the dense gate.

No producer/project imports. The accompanying proof supplies the monotonicity
and finite-diagram implications, independently checked here at their thresholds.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_dotted_central_20261014_gate_review'
BASE='fresh_h9_dotted_central_20261014_gate'

def multiply(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out

def run():
    gate=json.loads((ROOT/f'results/{BASE}.json').read_text())
    for p,h in gate['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h
    catalogue='results/fresh_h9_tree_cube_wedges_20261014.json'
    g=json.loads((ROOT/catalogue).read_text())['representatives'][5]['geometry']
    missing={frozenset(x) for x in g['minimal_nonfaces']}
    vertices={frozenset(v) for v in combinations(range(18),9) if not any(s<=set(v) for s in missing)}
    archived=[frozenset(v['incident']) for v in g['vertices']]
    assert len(vertices)==444 and vertices==set(archived)
    checked=[]
    for i,row in enumerate(gate['maps']):
        A,B=row['dense_block'],row['other_block'];x,y=row['dense_leaf'],row['other_leaf']
        assert len(set(A+B+[x,y]))==6
        assert frozenset(A+[x]) in missing and frozenset(B+[y]) in missing
        expected=[('positive_four_block',A+B),('x_foreign_right',B+[x]),('y_foreign_right',A+[y])]
        expected += [('middle_label_upper',[x,y,a,b]) for a in A for b in B]
        assert len(expected)==len(row['actual_supports'])==7
        for (kind,support),r in zip(expected,row['actual_supports']):
            assert r['kind']==kind and r['support']==support
            assert frozenset(support)<=archived[r['actual_vertex']]
        conditions=sorted([sorted(A)+[6],sorted(B)+[7],sorted((x,A[0]))+[3],sorted((x,A[1]))+[3]])
        assert row['bounds']==conditions
        checked.append(dict(map=i,checked_faces=7,conditions=conditions))
    unique={tuple(map(tuple,row['bounds'])) for row in gate['maps']}
    assert len(checked)==44 and len(unique)==24
    assert {tuple(map(tuple,row['bounds'])) for row in gate['families']}==unique
    for row in gate['families']:
        assert not row['fixed'] and row['bounds']==gate['maps'][row['witness_map']]['bounds']
    # Direct Chebyshev coefficients, independent of the producer's recurrence.
    cubic=[F(-1),F(-4),F(4),F(8)]
    factor=multiply([-1,1],multiply(cubic,cubic))
    assert factor==list(map(F,[-1,-7,0,56,0,-112,0,64]))
    evaluate=lambda p,z:sum((v*z**i for i,v in enumerate(p)),F(0))
    assert evaluate(cubic,F(3,5))==-F(29,125) and evaluate(cubic,F(1))==7
    # Four-path determinant threshold and single-contact B projection threshold.
    assert (1-F(1,4))**2==F(9,16)
    assert F(1,4)/(1-F(4,5))==F(5,4)
    # Q(sqrt3) verifies the dense threshold and determinant lower bound.
    qmul=lambda a,b:(a[0]*b[0]+3*a[1]*b[1],a[0]*b[1]+a[1]*b[0])
    assert qmul((F(2),F(1)),(F(1),-F(1,2)))==(F(1,2),F(0))
    determinant=(F(1,4)-F(1,2),F(1,4))
    assert determinant==(-F(1,4),F(1,4)) and 3*determinant[1]**2>determinant[0]**2
    # Label5 exceeds the four-path threshold: (3+sqrt5)/8 >9/16.
    assert F(5)*2**2>3**2
    files=[f'scripts/{TAG}.py',f'scripts/{BASE}.py',f'results/{BASE}.json',catalogue]
    return dict(status='PASS_INDEPENDENT_LITERAL_AND_SCALAR_DENSE_HIGH_GATE_REVIEW',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in files},
        reconstructed_vertices=len(vertices),literal_maps=checked,actual_face_premises=44*7,
        distinct_condition_sets=len(unique),chebyshev_factor_coefficients=list(map(str,factor)),
        dense_projection_threshold=['2','1'],other_projection_strict_lower='5/4',
        residual_determinant_strict_lower=['-1/4','1/4'],
        scope='Independent reconstruction of all444 vertices,44 injective literal placements,308 face premises,24 condition sets and scalar thresholds. Analytic all-parameter reasoning is reviewed in the accompanying report. No full case5 exclusion or complete traversal is claimed.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['actual_face_premises'],flush=True)
