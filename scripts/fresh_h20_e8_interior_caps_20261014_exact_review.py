#!/usr/bin/env python3
"""Independent exact nested-quadratic replay of the frozen 87-wall extension.

No producer or project arithmetic imports.  Multiplication and sign are
recursive in Q(sqrt2)(sqrt79)(sqrt11), with sign determined by squaring.
The inherited 79-wall/source completeness remains a certified dependency.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_e8_interior_caps_20261014_exact_review'
PRIMARY = 'results/fresh_h20_e8_interior_caps_20261014_exact.json'
BASE = 'results/fresh_h20_e8_cusp_completion_20261014.json'
FACE = 'results/fresh_h20_contact_closure_20261014_h12_replacement.json'
SOURCE = 'results/fresh_h12_a12_source_20261014.json'
PRIMES = [2, 79, 11]
Z = (F(0),)*8


def scalar(x): return (F(x),)+(F(0),)*7
def add(x, y): return tuple(a+b for a, b in zip(x, y))
def neg(x): return tuple(-a for a in x)
def scale(x, q): return tuple(a*q for a in x)
def sub(x, y): return add(x, neg(y))


def multiply(x, y):
    """Nested quadratic multiplication; tuple length specifies the field."""
    if not any(x) or not any(y): return (F(0),)*len(x)
    if not any(x[1:]): return scale(y, x[0])
    if not any(y[1:]): return scale(x, y[0])
    n = len(x)//2
    if not n: return (x[0]*y[0],)
    a, b, c, d = x[:n], x[n:], y[:n], y[n:]
    p = PRIMES[n.bit_length()-1]
    ac, bd = multiply(a, c), multiply(b, d)
    # Karatsuba form gives the coefficient of the final square root.
    cross = sub(sub(multiply(add(a, b), add(c, d)), ac), bd)
    return add(ac, scale(bd, p))+cross


def sign(x):
    if not any(x): return 0
    if len(x) == 1: return 1 if x[0] > 0 else -1
    n = len(x)//2
    a, b = x[:n], x[n:]
    sa, sb = sign(a), sign(b)
    if not sb: return sa
    if not sa or sa == sb: return sb
    q = sub(multiply(a, a), scale(multiply(b, b), PRIMES[n.bit_length()-1]))
    # Independent square classes of 2,79,11 ensure a nonzero vector
    # cannot have zero relative norm.
    sq = sign(q)
    assert sq != 0
    return sa*sq


def field(x):
    assert len(x) in (4, 8)
    return tuple(map(F, x))+(F(0),)*(8-len(x))


def vector(x): return tuple(field(y) for y in x)
def vadd(x, y): return tuple(add(a, b) for a, b in zip(x, y))
def vscale(x, c): return tuple(multiply(a, c) for a in x)
def rational(x):
    assert not any(x[1:])
    return x[0]


def pd(a):
    """Rational LDL, returning false on a nonpositive pivot."""
    a = [list(row) for row in a]
    for i in range(len(a)):
        t = a[i][i]
        if t <= 0: return False
        for j in range(i+1, len(a)):
            for k in range(j, len(a)):
                a[k][j] = a[j][k] = a[j][k]-a[i][j]*a[i][k]/t
    return True


def inverse(a):
    n = len(a)
    b = [list(row)+[F(int(i == j)) for j in range(n)] for i, row in enumerate(a)]
    for i in range(n):
        j = next(j for j in range(i, n) if b[j][i])
        b[i], b[j] = b[j], b[i]
        q = b[i][i]; b[i] = [v/q for v in b[i]]
        for j in range(n):
            if j != i:
                q = b[j][i]
                b[j] = [v-q*w for v, w in zip(b[j], b[i])]
    return [row[n:] for row in b]


def digest(p): return sha256((ROOT/p).read_bytes()).hexdigest()


def run():
    d, base, face, source = [json.loads((ROOT/p).read_text()) for p in [PRIMARY, BASE, FACE, SOURCE]]
    for doc in (d, base, face, source):
        for p, h in doc['sources'].items(): assert digest(p) == h, p
    assert d['field_radicands'] == [1, 2, 79, 158, 11, 22, 869, 1738]
    Q = [[F(x) for x in row] for row in d['form']]
    assert Q == [[F(x) for x in row] for row in base['form']]
    sparse = [[(j, q) for j, q in enumerate(row) if q] for row in Q]
    form_cache = {}

    def dot(x, y):
        if y not in form_cache:
            form_cache[y] = tuple(tuple(sum(q*y[j][k] for j, q in row) for k in range(8)) for row in sparse)
        out = Z
        for a, b in zip(x, form_cache[y]):
            if any(a) and any(b): out = add(out, multiply(a, b))
        return out

    def linear(vectors, coefficients):
        out = (Z,)*21
        for v, c in zip(vectors, coefficients): out = vadd(out, vscale(v, c))
        return out

    N = list(map(vector, d['normals'])); assert len(N) == 87
    assert N[:79] == list(map(vector, base['normals']))
    H = vector(d['strict_center']); assert H == vector(base['strict_center'])
    Eold = vector(base['strict_escape']['point'])
    E = face['E8_support']; assert E == base['E8_support']
    B = [[rational(dot(N[i], N[j])) for j in E] for i in E]
    assert pd(B)
    I = inverse(B); assert I == [[F(x) for x in row] for row in face['E8_inverse']]
    weights = [linear([N[j] for j in E], list(map(scalar, row))) for row in I]
    assert all(dot(w, N[j]) == scalar(int(i == k)) for i, w in enumerate(weights) for k, j in enumerate(E))
    by_inc = {tuple(v['source_incident']):v for v in face['H12_face_finite_vertices']+face['H12_face_ideal_vertices']}
    V = [vector(by_inc[tuple(row['incident'])]['point']) for row in source['projective_vertices'][:2]]
    assert [dot(V[0], V[0]), dot(V[0], V[1]), dot(V[1], V[1])] == list(map(scalar, [-1, -1, F(-13, 121)]))
    literal_coefficients = [
        [(F(17,8),0,0,0,0,0,F(-1,24),0), (F(-2057,104),0,F(187,52),0,F(-99,52),0,F(121,312),0)],
        [(F(33,16),0,0,0,0,0,F(-1,48),0), (F(-3993,208),0,F(363,104),0,F(-99,104),0,F(121,624),0)]]
    prototypes = []
    for k, (profile, coefficients) in enumerate(zip([[3,6], [4,5]], literal_coefficients)):
        rec = d['new_recipes'][k]
        assert rec['E8_positions'] == profile and rec['E8_walls'] == [E[j] for j in profile]
        coeff = [field(row) for row in coefficients]
        assert coeff == [field(rec['a']), field(rec['b'])] and all(sign(c) > 0 for c in coeff)
        w = linear([weights[j] for j in profile], [scalar(1)]*2)
        y = linear(V, coeff)
        assert dot(w, y) == Z and dot(w, w) == field(rec['weight_norm'])
        assert sign(dot(y, y)) < 0 and sign(dot(y, H)) < 0
        r = vadd(vscale(w, scalar(-1)), y)
        assert dot(r, r) == scalar(2)
        assert [dot(r, N[j]) for j in E] == [scalar(-int(j in profile)) for j in range(8)]
        assert all(dot(r, N[j]) == Z for j in range(11))
        assert dot(r, N[19]) == scalar(-int(k == 0)) == field(rec['anchor_pairing'])
        assert dot(r, Eold) == field(rec['old_escape_pairing']) and sign(dot(r, Eold)) > 0
        prototypes.append(r)
    assert dot(*prototypes) == Z
    rebuilt = []
    for i, r in enumerate(prototypes):
        for s, t in [(1,1), (1,-1), (-1,1), (-1,-1)]:
            v = list(r); v[19] = scale(v[19], s); v[20] = scale(v[20], t)
            assert d['orbit'][len(rebuilt)] == dict(primitive=i, height_signs=[s,t])
            rebuilt.append(tuple(v))
    assert len(set(rebuilt)) == 8 and rebuilt == N[79:]

    # A complete Gram replay supplies every subsequent finite-neighbor test.
    G = [[Z]*87 for _ in range(87)]
    for i in range(87):
        for j in range(i+1): G[i][j] = G[j][i] = dot(N[i], N[j])
    norms = [rational(G[i][i]) for i in range(87)]
    assert norms == list(map(F, d['norms'])) and all(q > 0 for q in norms)
    assert all(sign(G[i][j]) <= 0 for i in range(87) for j in range(i))
    counts = Counter(); pairmap = {tuple(row['walls']):row for row in d['new_pairs']}
    assert len(pairmap) == 660
    for i in range(79,87):
        assert norms[i] == 2
        for j in range(i):
            g = G[i][j]; square = multiply(g, g); prod = norms[i]*norms[j]
            if g == Z: kind = 'right'
            elif square == scalar(prod/4): kind = 'pi3'
            else:
                assert sign(sub(square, scalar(prod))) > 0
                kind = 'dotted'
            assert pairmap[i,j]['relation'] == kind and field(pairmap[i,j]['pairing']) == g
            counts[kind] += 1
    assert dict(counts) == d['new_pair_counts'] == dict(right=144, pi3=20, dotted=496)

    Hq = dot(H,H); assert sign(Hq) < 0
    h = [dot(H, n) for n in N]; assert all(sign(x) < 0 for x in h)
    for i in range(87):
        p = vadd(H, vscale(N[i], scale(h[i], -1/norms[i])))
        assert p == vector(d['facet_points'][i])
        qp = sub(Hq, scale(multiply(h[i],h[i]), 1/norms[i]))
        assert sign(qp) < 0
        # Pairing formula, rather than another full coordinate multiplication.
        for j in range(87):
            g = sub(h[j], scale(multiply(h[i],G[i][j]), 1/norms[i]))
            assert g == Z if i == j else sign(g) < 0

    saved_vertices = {tuple(row['incident']):row for row in d['new_finite_vertices']}
    assert len(saved_vertices) == 8
    finite = [[i != j and sign(sub(multiply(G[i][j],G[i][j]), scalar(norms[i]*norms[j]))) < 0 for j in range(87)] for i in range(87)]
    accepted_all = set(); checks = 0
    for i, saved in zip(range(79,87), d['new_cap_finite_vertex_census']):
        neighbors = [j for j in range(87) if finite[i][j]]
        assert saved['wall'] == i and neighbors == saved['finite_neighbors'] and len(neighbors) == 21
        accepted = []
        for part in combinations(neighbors,19):
            J = tuple(sorted((i,)+part)); checks += 1
            if not all(finite[a][b] for a,b in combinations(J,2)): continue
            A = [[rational(G[a][b]) for b in J] for a in J]
            if pd(A): accepted.append(list(J)); accepted_all.add(J)
        assert saved['subsets_tested'] == 210 and accepted == saved['accepted']
        assert len(accepted) == (1 if i < 83 else 2)
    assert checks == 1680 and accepted_all == set(saved_vertices)
    types = Counter()
    for J, row in saved_vertices.items():
        A = [[rational(G[a][b]) for b in J] for a in J]
        assert A == [[F(x) for x in rr] for rr in row['gram']] and pd(A)
        # Reconstruct the center projection with our own rational elimination.
        AI = inverse(A)
        coef = [tuple(sum(AI[a][b]*h[J[b]][k] for b in range(20)) for k in range(8)) for a in range(20)]
        point = vadd(H, vscale(linear([N[j] for j in J],coef),scalar(-1)))
        assert point == vector(row['point']) and dot(point,point) == field(row['norm'])
        assert sign(dot(point,point)) < 0 and sign(dot(point,H)) < 0
        assert all(dot(point,N[j]) == Z if j in J else sign(dot(point,N[j])) < 0 for j in range(87))
        unseen = set(J); kinds=[]
        while unseen:
            todo=[min(unseen)]; block=set()
            while todo:
                a=todo.pop()
                if a in block: continue
                block.add(a); todo.extend(b for b in J if a != b and G[a][b] != Z)
            unseen-=block
            adj={a:[b for b in block if a!=b and G[a][b]!=Z] for a in block}
            assert all(norms[a]==2 for a in block)
            assert all(G[a][b]==scalar(-1) for a in block for b in adj[a])
            assert sum(map(len,adj.values()))==2*(len(block)-1)
            branch=[a for a in block if len(adj[a])==3]
            if not branch:
                assert max(map(len,adj.values()))<=2; kinds.append('A'+str(len(block)))
            else:
                assert len(branch)==1 and len(block)==7
                arms=[]
                for b in adj[branch[0]]:
                    previous,current,length=branch[0],b,1
                    while len(adj[current])==2:
                        nxt=next(a for a in adj[current] if a!=previous)
                        previous,current=current,nxt;length+=1
                    arms.append(length)
                assert sorted(arms)==[1,2,3]; kinds.append('E7')
        assert sorted(kinds)==row['spherical_type']; types[tuple(sorted(kinds))]+=1
    assert types == {('A11','A2','A7'):4, ('A11','A2','E7'):4}
    Jlist=[tuple(row['incident']) for row in d['new_finite_vertices']]
    edges=[]
    for a,b in combinations(range(8),2):
        S=sorted(set(Jlist[a]) & set(Jlist[b]))
        if len(S)==19:
            assert pd([[rational(G[i][j]) for j in S] for i in S])
            edges.append(dict(vertices=[a,b],support=S))
    assert len(edges)==4 and edges==d['new_finite_edges']

    for row, saved in zip(face['H12_face_finite_vertices'],d['preserved_H12_finite_vertices']):
        p=vector(row['point']); values=[dot(p,n) for n in N[79:]]
        assert row['incident']==saved['incident'] and values==[field(v) for v in saved['new_pairings']]
        assert all(sign(v)<0 for v in values)
    assert len(d['preserved_H12_finite_vertices'])==35
    for row,saved in zip(base['completed_cusp_links'],d['preserved_two_cusp_links']):
        p=vector(row['point']);values=[dot(p,n) for n in N[79:]]
        assert values==[field(v) for v in saved] and all(sign(v)<0 for v in values)
    assert len(d['preserved_two_cusp_links'])==2
    assert [dot(Eold,n) for n in N[79:]]==[field(v) for v in d['previous_escape']]
    escape=vector(d['strict_escape']['point']);eq=dot(escape,escape)
    assert eq==field(d['strict_escape']['norm'])==scalar(F(1336028658157091827,250000000000000000))
    assert sign(eq)>0 and sign(dot(escape,H))<0 and all(sign(dot(escape,n))<0 for n in N)
    assert dot(escape,H)==field(d['strict_escape']['center_pairing'])
    paths=[PRIMARY,BASE,FACE,SOURCE,'scripts/fresh_h20_e8_interior_caps_20261014_exact.py',f'scripts/{TAG}.py']
    return dict(status='PASS_INDEPENDENT_NESTED_QUADRATIC_87_WALL_REPLAY',sources={p:digest(p) for p in paths},arithmetic='stdlib Fraction; recursive quadratic multiplication and exact sign by relative norms, no project or producer imports',new_pairs=660,new_pair_counts=dict(counts),all_norms=87,actual_facets=87,new_cap_PD20_subsets=checks,new_finite_vertices=8,certified_finite_edges=4,preserved_finite_vertices=35,preserved_bounded_cusps=2,escape_norm=str(rational(eq)),scope='Reconstructs all8 new normals, the full87 Gram, all660 new Coxeter relations, strict center/facet projections, complete finite-corner census incident to new caps,8 physical vertices/4 connecting edges,35 inherited vertices/2 cusp links, and the new strict escape. Inherited79 geometry and original H12 finite-volume/vertex completeness are dependencies. No full87 ambient corner, edge or cusp census.')


if __name__ == '__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();path=ROOT/f'results/{TAG}.json'
    if args.verify: assert result==json.loads(path.read_text())
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['new_pair_counts'],result['escape_norm'])
