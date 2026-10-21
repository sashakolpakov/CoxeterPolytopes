#!/usr/bin/env python3
"""Complete compact tetrahedron Gram templates with a prescribed pi/5 pair."""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, permutations, product
from math import isqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_pi5_tetrahedron_templates_20261017'
ZERO = (F(0),)*4
ONE = (F(1), F(0), F(0), F(0))


def const(x):
    return (F(x), F(0), F(0), F(0))


def add(a, b):
    return tuple(x+y for x,y in zip(a,b))


def neg(a):
    return tuple(-x for x in a)


def mul(a, b):
    out = [F(0)]*4
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            if x and y:
                common = i & j
                out[i ^ j] += x*y*(2 if common & 1 else 1)*(5 if common & 2 else 1)
    return tuple(out)


def total(xs):
    result = ZERO
    for x in xs:
        result = add(result, x)
    return result


SCALE = 10**16
R2 = (F(isqrt(2*SCALE*SCALE), SCALE), F(isqrt(2*SCALE*SCALE)+1, SCALE))
R5 = (F(isqrt(5*SCALE*SCALE), SCALE), F(isqrt(5*SCALE*SCALE)+1, SCALE))
assert R2[0]**2 < 2 < R2[1]**2 and R5[0]**2 < 5 < R5[1]**2
BASIS_INTERVALS = [(F(1),F(1)),R2,R5,(R2[0]*R5[0],R2[1]*R5[1])]


def bounds(a):
    lo = hi = F(0)
    for c,(l,h) in zip(a,BASIS_INTERVALS):
        lo += c*(l if c>=0 else h)
        hi += c*(h if c>=0 else l)
    return lo,hi


def sign(a):
    if a == ZERO:
        return 0
    lo,hi = bounds(a)
    assert lo > 0 or hi < 0, ('Increase exact isolator precision',a,lo,hi)
    return 1 if lo>0 else -1


def determinant(g):
    out = ZERO
    for p in permutations(range(len(g))):
        term = ONE
        for i,j in enumerate(p):
            term = mul(term,g[i][j])
        if sum(p[i]>p[j] for i,j in combinations(range(len(g)),2))%2:
            term = neg(term)
        out = add(out,term)
    return out


def sub(g, rows, cols=None):
    if cols is None:
        cols=rows
    return [[g[i][j] for j in cols] for i in rows]


def encode(x):
    return list(map(str,x))


def run():
    pairs = list(combinations(range(4),2))
    values = {2:ZERO,3:const(-1),4:(F(0),F(-1),F(0),F(0)),
              5:(F(-1,2),F(0),F(-1,2),F(0))}
    accepted = []
    counts = Counter()
    representatives = []
    registry = {}
    marked_maps = [(a,b,c,d) for a,b in ((0,1),(1,0)) for c,d in ((2,3),(3,2))]
    for tail in product(range(2,6),repeat=5):
        labels = [5]+list(tail)
        gram = [[const(2) if i==j else ZERO for j in range(4)] for i in range(4)]
        lm = [[2]*4 for _ in range(4)]
        for (i,j),label in zip(pairs,labels):
            gram[i][j]=gram[j][i]=values[label]
            lm[i][j]=lm[j][i]=label
        minors = [determinant(sub(gram,[i for i in range(4) if i!=j])) for j in range(4)]
        if any(sign(v)<=0 for v in minors):
            counts['REJECT_NONSPHERICAL_VERTEX'] += 1
            continue
        det = determinant(gram)
        if sign(det)>=0:
            counts['REJECT_NONHYPERBOLIC_FULL_GRAM'] += 1
            continue
        # Independent finite-vertex cone verification: v_j=adj(G)e_j,
        # Gv_j=det(G)e_j and v_i.v_j=det(G)*adj(G)_ij<0.
        adj = []
        for i in range(4):
            row=[]
            for j in range(4):
                cofactor=determinant(sub(gram,[k for k in range(4) if k!=j],
                                         [k for k in range(4) if k!=i]))
                row.append(neg(cofactor) if (i+j)%2 else cofactor)
            adj.append(row)
        assert all(sign(v)>0 for row in adj for v in row)
        for i,j in product(range(4),repeat=2):
            assert total(mul(gram[i][k],adj[k][j]) for k in range(4)) == (det if i==j else ZERO)
        center_norm = mul(det,total(v for row in adj for v in row))
        assert sign(center_norm)<0
        vertex_norms = [mul(det,adj[i][i]) for i in range(4)]
        assert all(sign(v)<0 for v in vertex_norms)
        options = [(tuple(lm[p[i]][p[j]] for i,j in pairs),p) for p in marked_maps]
        canonical,mp = min(options)
        if canonical not in registry:
            registry[canonical]=len(representatives)
            representatives.append(dict(index=len(representatives),labels=list(canonical)))
        index = len(accepted)
        accepted.append(dict(index=index,labels=labels,
                             gram=[[encode(v) for v in row] for row in gram],
                             vertex_determinants=list(map(encode,minors)),determinant=encode(det),
                             adjugate=[[encode(v) for v in row] for row in adj],
                             finite_vertex_norms=list(map(encode,vertex_norms)),
                             strict_center_norm=encode(center_norm),
                             marked_representative=registry[canonical],
                             canonical_to_literal=list(mp)))
        counts['ACCEPT_COMPACT_COXETER_TETRAHEDRON'] += 1
    assert sum(counts.values()) == 4**5
    assert accepted
    # Every marked representative has its literal realization in the output.
    literal = {tuple(r['labels']) for r in accepted}
    assert all(tuple(r['labels']) in literal for r in representatives)
    # The fixed I2(5) pair is spherical; its complement in R^(3,1) is R^(1,1).
    fixed_determinant=add(const(4),neg(mul(values[5],values[5])))
    assert sign(fixed_determinant)>0
    paths=[f'scripts/{TAG}.py',f'results/{TAG}.md',
           'scripts/fresh_nonarithmetic_h3_edge730_domain_20261017.py',
           'results/fresh_nonarithmetic_h3_edge730_domain_20261017.json',
           'results/fresh_nonarithmetic_h3_edge730_domain_20261017.md']
    return dict(status='PASS_COMPLETE_PI5_COMPACT_TETRAHEDRON_TEMPLATES',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                exact_field_basis=['1','sqrt2','sqrt5','sqrt10'],
                radical_isolators=[[str(v) for v in pair] for pair in (R2,R5)],
                pair_order=[list(p) for p in pairs],old_root_labels=[7,30],
                full_search_size=4**5,counts=dict(counts),
                literal_templates=len(accepted),marked_classes=len(representatives),
                representatives=representatives,records=accepted,
                all_label_proof='A compact tetrahedron has every three-root principal block positive definite and full Gram Lorentzian. If a label is at least6, either endpoint attached nontrivially to a third root makes that triple determinant nonpositive, since cos(pi/m)^2+cos(pi/3)^2>=1. Thus both endpoints would be disconnected from the other two roots, making the full Gram positive definite. Therefore all labels are2..5.',
                compactness_proof='Every two-root block is positive definite and every three-root determinant is positive. The full negative determinant gives signature(3,1). The positive adjugate yields four future timelike dual rays with all other facet products zero and its own negative. These generate the whole inequality cone, so the resulting tetrahedron is compact and every wall is actual.',
                embedding_freedom='Fixing the two old norm-two roots with mutual product -phi leaves O(1,1) acting on their Lorentz-plane orthogonal complement. Each literal template therefore has a complete one-parameter Lorentz-boost family plus the discrete components; an ambient compatible embedding must still pass all native old gates.',
                limitations='Complete local tetrahedron shapes only. No ambient compatible embedding or closure of the original H9 face is asserted. Other compact polyhedra with more facets are not excluded by failure of these templates.')


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--verify',action='store_true')
    args=parser.parse_args()
    result=run()
    path=ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result==json.loads(path.read_text())
    else:
        path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['counts'],result['marked_classes'],flush=True)
