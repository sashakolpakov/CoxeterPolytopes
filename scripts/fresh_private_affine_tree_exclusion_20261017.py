#!/usr/bin/env python3
"""Exact local checks for the all-rank private affine-tree exclusion.

Finite probes corroborate the written resistance/extension proof; they are
not a cutoff substitute for that proof.  No predecessor arithmetic is imported.
"""
import argparse
from collections import deque
from fractions import Fraction as F
import hashlib
from itertools import combinations, product
import json
from pathlib import Path
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_private_affine_tree_exclusion_20261017'
SOURCES = [
    'results/fresh_private_affine_cycle_exclusion_20261017.md',
    'results/fresh_simple_finite_quad_vinberg_transfer_20261017.md',
    'results/fresh_simple_cusp_vinberg_weight_frontier_20261017.md',
    'scripts/fresh_simple_cusp_affine_tree_completion_20261017.py',
    'results/fresh_simple_cusp_affine_tree_completion_20261017.md',
]

def sha(p):
    return hashlib.sha256((ROOT / p).read_bytes()).hexdigest()

def graph(kind, n):
    """Edges carry finite labels, cusp entries are nu_i^2."""
    nu2 = [F(1)] * n
    if kind == 'C':
        edges = [(i, i+1, 4 if i in (0, n-2) else 3) for i in range(n-1)]
        nu2[0] = nu2[-1] = F(1, 2)
    elif kind == 'B':
        edges = [(i, i+1, 4 if i == n-4 else 3) for i in range(n-3)]
        edges += [(0, n-2, 3), (0, n-1, 3)]
        nu2[n-3] = F(1, 2)
        nu2[-2:] = [F(1, 4)] * 2
    elif kind == 'D':
        edges = [(i, i+1, 3) for i in range(n-5)]
        edges += [(0, n-4, 3), (0, n-3, 3),
                  (n-5, n-2, 3), (n-5, n-1, 3)]
        nu2[-4:] = [F(1, 4)] * 4
    elif kind == 'A':
        edges = [(i, (i+1) % n, 3) for i in range(n)]
    else:
        raise ValueError(kind)
    return edges, nu2

def adjacency(n, edges):
    out = [[] for _ in range(n)]
    for i, j, m in edges:
        out[i].append((j, m)); out[j].append((i, m))
    return out

def tree_key(n, edges):
    adj = adjacency(n, edges)
    assert len(edges) == n-1
    def rooted(i, parent):
        return tuple(sorted((m, rooted(j, i)) for j, m in adj[i] if j != parent))
    return min(rooted(i, -1) for i in range(n))

def gram(n, edges):
    G = 2*s.eye(n)
    for i, j, m in edges:
        assert m in (3, 4)
        G[i,j] = G[j,i] = -1 if m == 3 else -s.sqrt(2)
    return G

def arms(lengths):
    es = []; nxt = 1
    for length in lengths:
        p = 0
        for _ in range(length):
            es.append((p, nxt, 3)); p = nxt; nxt += 1
    return nxt, es

def affine_catalogue():
    rows = []
    for kind, lo in [('B',4), ('C',3), ('D',5)]:
        for n in range(lo, 10):
            edges, _ = graph(kind, n)
            rows.append((f'affine_{kind}_order{n}', n, edges))
    for name, aa in [('E6',(2,2,2)), ('E7',(1,3,3)), ('E8',(1,2,5))]:
        n, edges = arms(aa); rows.append((f'affine_{name}', n, edges))
    rows.append(('affine_F4', 5, [(0,1,3),(1,2,3),(2,3,4),(3,4,3)]))
    data = []; keys = {}
    for name, n, edges in rows:
        G = gram(n, edges)
        ker = G.nullspace(); assert len(ker) == 1
        v = ker[0]
        if v[0] < 0: v = -v
        assert all(x > 0 for x in v)
        # A positive null vector scales this acute Gram to a graph Laplacian.
        L = s.diag(*v)*G*s.diag(*v)
        assert all(s.simplify(sum(L.row(i))) == 0 for i in range(n))
        assert all(L[i,j] <= 0 for i in range(n) for j in range(n) if i != j)
        key = tree_key(n, edges); assert key not in keys
        keys[key] = name
        data.append(dict(name=name, order=n, edges=edges, positive_kernel=list(map(str,v))))
    return keys, data

def induced(edges, vertices):
    ix = {x:i for i,x in enumerate(sorted(vertices))}
    return [(ix[a],ix[b],m) for a,b,m in edges if a in ix and b in ix]

def leaf_witness(n, edges, at, label, keys):
    z = n; es = edges+[(at,z,label)]; adj = adjacency(n+1,es)
    level = {frozenset([z])}
    for size in range(1, 10):
        for U in sorted(level, key=lambda t:tuple(sorted(t))):
            ee = induced(es,U)
            if len(ee) == size-1:
                name = keys.get(tree_key(size,ee))
                if name:
                    return dict(at=at, label=label, vertices=sorted(U), affine_type=name)
        level = {U | {j} for U in level for i in U for j,_ in adj[i] if j not in U}
    raise AssertionError(('no small affine leaf witness',n,at,label))

def resistance_data(n, edges, nu2):
    adj = adjacency(n,edges)
    R = [[None]*n for _ in range(n)]
    L = [[F(0)]*n for _ in range(n)]
    for i,j,_ in edges:
        c = F(1,2) if nu2[i] == F(1,4) or nu2[j] == F(1,4) else F(1)
        L[i][i] += c; L[j][j] += c; L[i][j] -= c; L[j][i] -= c
    assert [L[i][i] for i in range(n)] == [2*x for x in nu2]
    for source in range(n):
        R[source][source] = F(0); todo = deque([source])
        while todo:
            i = todo.popleft()
            for j,_ in adj[i]:
                if R[source][j] is None:
                    R[source][j] = R[source][i]-1/L[i][j]; todo.append(j)
    checks = 0
    for i,j in combinations(range(n),2):
        # The grounded Green column, not merely its diagonal, is checked.
        v = [(R[k][j]+R[i][j]-R[k][i])/2 for k in range(n)]
        assert v[j] == 0 and v[i] == R[i][j]
        assert [sum(L[k][ell]*v[ell] for ell in range(n)) for k in range(n)] == [F(k==i)-F(k==j) for k in range(n)]
        checks += 1
    return R, checks

def terminal_type(nu2, i):
    return 'double' if nu2[i] == F(1,2) else 'fork' if nu2[i] == F(1,4) else 'main'

def two_contact_checks(n, edges, nu2, R):
    es = {frozenset((a,b)):m for a,b,m in edges}; rows=[]; trials=0
    for i,j in combinations(range(n),2):
        for ti,tj in product((1,2),repeat=2):
            trials += 1
            if ti*nu2[i]*R[i][j] >= 2 or tj*nu2[j]*R[i][j] >= 2: continue
            kinds=(terminal_type(nu2,i),terminal_type(nu2,j))
            if kinds == ('main','main'):
                assert R[i][j] == 1 and ti == tj == 1 and es[frozenset((i,j))] == 3
                outcome='forbidden_affine_A2'
            elif kinds == ('fork','fork'):
                assert R[i][j] == 4 and ti == tj == 1
                outcome='forbidden_affine_A3'
            else:
                assert set(kinds) == {'double','main'} and R[i][j] == 1
                a,b = (i,j) if kinds[0]=='double' else (j,i)
                ta,tb = (ti,tj) if a==i else (tj,ti)
                assert tb == 1 and ta in (1,2) and es[frozenset((a,b))] == 4
                outcome='terminal_triangle'
            rows.append(dict(pair=[i,j], contact_squares=[ti,tj], R=str(R[i][j]), outcome=outcome))
    return dict(trials=trials, surviving_rows=rows)

def glue_checks():
    rows=[]
    for kl,kr in product(('B','D'),repeat=2):
        for nl,nr in ((18,18),(18,23),(23,19)):
            # Start with affine C (for a B arm) or B (for a D arm),
            # remove its doubled terminal, then join the two long endpoints.
            arms_data=[]
            for kind,n in ((kl,nl),(kr,nr)):
                old='C' if kind=='B' else 'B'
                edges,nu2=graph(old,n)
                a=next(i for i,x in enumerate(nu2) if x==F(1,2))
                b=next(j if i==a else i for i,j,_ in edges if a in(i,j))
                keep=[i for i in range(n) if i!=a]; ix={x:i for i,x in enumerate(keep)}
                arms_data.append((len(keep),induced(edges,keep),ix[b]))
            n1,e1,b1=arms_data[0];n2,e2,b2=arms_data[1];z=n1+n2
            ee=e1+[(a+n1,b+n1,m)for a,b,m in e2]+[(b1,z,3),(b2+n1,z,3)]
            target='C' if kl==kr=='B' else 'D' if kl==kr=='D' else 'B'
            te,_=graph(target,z+1)
            assert tree_key(z+1,ee)==tree_key(z+1,te)
            rows.append(dict(left=kl,right=kr,left_order=n1,right_order=n2,
                             joined_order=z+1,affine_type=target))
    return rows

def schur_identity():
    # Grounding the affine block at its doubled terminal leaves residual
    # [[0,-h],[-h,c]], with h>0.  The inverse zz-entry is identically zero.
    h,c=s.symbols('h c',real=True)
    M=s.Matrix([[0,-h],[-h,c]])
    assert M.det()==-h*h and M.inv()[1,1]==0
    x,y=s.symbols('x y',real=True)
    assert (s.Matrix([[0,x]])*M.inv()*s.Matrix([0,y]))[0]==0
    return dict(residual_matrix=[['0','-h'],['-h','c']],determinant='-h**2',
                inverse_zz='0',outside_schur_correction='0')

def build():
    keys,catalogue=affine_catalogue(); probes=[]
    for n in (18,19,23):
        for kind in ('B','C','D'):
            edges,nu2=graph(kind,n);R,count=resistance_data(n,edges,nu2)
            leaves=[leaf_witness(n,edges,i,label,keys) for i in range(n) for label in (3,4)]
            twos=two_contact_checks(n,edges,nu2,R)
            if kind=='D': assert not any(r['outcome']=='terminal_triangle' for r in twos['surviving_rows'])
            probes.append(dict(kind=kind,order=n,green_columns=count,leaf_witnesses=leaves,two_contacts=twos))
            print('PRIVATE_TREE_LOCAL_PASS',kind,n,flush=True)
    # Also corroborate the simple-leaf obstruction for the second affine
    # component when it is a cycle; the private-cycle proof is inherited.
    cycle_leaves=[]
    for n in (18,19,23):
        edges,_=graph('A',n)
        cycle_leaves.append(dict(order=n,witness=leaf_witness(n,edges,0,3,keys)))
    return dict(status='PASS_PRIVATE_AFFINE_TREE_EXCLUSION_LOCAL_CHECKS',
                scope='All-rank theorem is in report; finite graph probes corroborate its identities only.',
                sources={p:sha(p)for p in SOURCES},source_sha256=sha('scripts/'+TAG+'.py'),
                affine_witness_catalogue=catalogue,probes=probes,cycle_leaf_checks=cycle_leaves,
                schur_identity=schur_identity(),terminal_gluings=glue_checks(),
                counts=dict(green_columns=sum(p['green_columns']for p in probes),
                            leaf_witnesses=sum(len(p['leaf_witnesses'])for p in probes)+len(cycle_leaves),
                            two_contact_trials=sum(p['two_contacts']['trials']for p in probes),
                            terminal_gluings=12))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=json.loads(json.dumps(build()));out=ROOT/'results'/f'{TAG}.json'
    if args.verify:
        assert result==json.loads(out.read_text())
    else:
        out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(result['status'],result['counts'],flush=True)

if __name__=='__main__': main()
