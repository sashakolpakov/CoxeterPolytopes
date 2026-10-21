#!/usr/bin/env python3
"""Exact small counterexample and contact-graph filter for cap coverings.

The H3 example and its complete lower-envelope cell rays are reconstructed
with standalone Q(sqrt(10)) arithmetic. The 860-branch graph inherits the
frozen cap norms/finite-contact identities, checks their source bindings,
recomputes the 301-parent graph, and checks every omitted graph edge with
the saved exact outward interval. It is not another 860-normal replay.
"""
import argparse
import hashlib
import itertools
import json
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h22_covering_cells_20261013'
OUT = ROOT / 'results' / (TAG + '.json')
PRIMARY = 'results/fresh_h22_multiquadratic_20261013.json'
RAW = 'data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
BANK = 'data/enumerations/fresh_h22_incidence_general_binary_20261013_bank.json'
CORE = 'results/fresh_h22_incidence_core_20261013.json'


@dataclass(frozen=True)
class K:
    a: F = F(0)
    b: F = F(0)

    def __post_init__(self):
        object.__setattr__(self, 'a', F(self.a))
        object.__setattr__(self, 'b', F(self.b))

    @staticmethod
    def cast(x):
        return x if isinstance(x, K) else K(x)

    def __add__(self, other):
        o = K.cast(other)
        return K(self.a + o.a, self.b + o.b)

    __radd__ = __add__

    def __neg__(self):
        return K(-self.a, -self.b)

    def __sub__(self, other):
        return self + -K.cast(other)

    def __rsub__(self, other):
        return K.cast(other) + -self

    def __mul__(self, other):
        o = K.cast(other)
        return K(self.a * o.a + 10 * self.b * o.b,
                 self.a * o.b + self.b * o.a)

    __rmul__ = __mul__

    def __truediv__(self, other):
        o = K.cast(other)
        d = o.a * o.a - 10 * o.b * o.b
        assert d != 0
        return self * K(o.a / d, -o.b / d)

    def sign(self):
        if not self.b:
            return (self.a > 0) - (self.a < 0)
        if not self.a:
            return (self.b > 0) - (self.b < 0)
        if self.a * self.b > 0:
            return (self.a > 0) - (self.a < 0)
        d = self.a * self.a - 10 * self.b * self.b
        assert d != 0  # sqrt(10) is irrational
        return ((self.a > 0) - (self.a < 0)) * ((d > 0) - (d < 0))

    def enc(self):
        return [str(self.a), str(self.b)]


def vec(v):
    return tuple(K.cast(x) for x in v)


def plus(u, v):
    return tuple(a + b for a, b in zip(u, v, strict=True))


def scale(v, s):
    return tuple(x * s for x in v)


def dot(u, v):
    return sum((a * b for a, b in zip(u, v, strict=True)), K())


def cov(u):
    return (-u[0],) + u[1:]


def pair(u, v):
    return dot(cov(u), v)


def rank(rows):
    a = [list(v) for v in rows]
    if not a:
        return 0
    r = 0
    for col in range(len(a[0])):
        p = next((j for j in range(r, len(a)) if a[j][col].sign()), None)
        if p is None:
            continue
        a[r], a[p] = a[p], a[r]
        z = a[r][col]
        a[r] = [x / z for x in a[r]]
        for j in range(len(a)):
            if j != r:
                z = a[j][col]
                a[j] = [x - z*y for x, y in zip(a[j], a[r], strict=True)]
        r += 1
        if r == len(a):
            break
    return r


def ev(v):
    return [x.enc() for x in v]


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2],
            a[0]*b[1]-a[1]*b[0])


def kernel4(rows):
    return tuple((-1)**j * dot([r[k] for k in range(4) if k != j],
                              cross([rows[1][k] for k in range(4) if k != j],
                                    [rows[2][k] for k in range(4) if k != j]))
                 for j in range(4) for r in [rows[0]])


def tiny_example():
    rt = K(0, 1)
    old = [vec(v) for v in [(0,1,0), (0,0,1), (1,-1,-1)]]
    caps = [vec(v) for v in [(3,0,0,rt), (3,0,0,-rt),
                             (3,-3,0,1), (3,-3,0,-1),
                             (3,0,-3,1), (3,0,-3,-1)]]
    walls = [v + (K(),) for v in old] + caps
    h = vec((3,-1,-1,0))
    assert pair(h,h) == K(-7)
    assert all(pair(h,r).sign() < 0 for r in walls)
    assert all(pair(r,r) == K(1) for r in walls)
    counts = Counter()
    gram = [[pair(r,s) for s in walls] for r in walls]
    for i in range(9):
        for j in range(i):
            p = gram[i][j]
            if not p.sign():
                counts['right'] += 1
            elif p == K(-1):
                counts['parallel'] += 1
            else:
                assert (p+1).sign() < 0
                counts['dotted'] += 1
    facets = []
    for i,r in enumerate(walls):
        p = plus(h, scale(r, -pair(h,r)))
        assert pair(p,p).sign() < 0 and pair(p,h).sign() < 0
        ps = [pair(p,s) for s in walls]
        assert [j for j,x in enumerate(ps) if not x.sign()] == [i]
        assert all(x.sign() < 0 for j,x in enumerate(ps) if j != i)
        facets.append(ev(p))
    finite = []
    for s in (1,-1):
        p = vec((1,0,0, s*3*rt/10))
        ps = [pair(p,r) for r in walls]
        inc = [0,1,3 if s == 1 else 4]
        assert pair(p,p) == K(F(-1,10))
        assert [j for j,x in enumerate(ps) if not x.sign()] == inc
        assert all(x.sign() < 0 for j,x in enumerate(ps) if j not in inc)
        finite.append({'point':ev(p), 'incident':inc, 'norm':pair(p,p).enc()})
    links = []
    for delta, blocks in [(vec((1,-1,0,0)), ((1,2),(5,6))),
                           (vec((1,0,-1,0)), ((0,2),(7,8)))]:
        ps = [pair(delta,r) for r in walls]
        inc = sorted(sum((list(b) for b in blocks), []))
        assert pair(delta,delta) == K() and pair(delta,h).sign() < 0
        assert [j for j,x in enumerate(ps) if not x.sign()] == inc
        assert all(x.sign() < 0 for j,x in enumerate(ps) if j not in inc)
        assert plus(walls[blocks[0][0]],walls[blocks[0][1]]) == delta
        assert plus(walls[blocks[1][0]],walls[blocks[1][1]]) == scale(delta,6)
        assert all(gram[i][j] == K() for i in blocks[0] for j in blocks[1])
        assert rank([walls[i] for i in inc]) == 3
        assert rank([[gram[i][j] for j in inc] for i in inc]) == 2
        links.append({'ray':ev(delta),'incident':inc,'affine_A1_blocks':[list(b) for b in blocks],
                      'positive_kernel_ray_multiples':[1,6], 'quotient_rank':2})
    x = vec((3,-1,-1,F(11,4)))
    xp = [pair(x,r) for r in walls]
    assert pair(x,x) == K(F(9,16)) and pair(x,h) == K(-7)
    assert all(p.sign() < 0 for p in xp)

    # The full projective cone has t>0 and a bounded t=1 section: the
    # old triangle bounds x,y and f+/f- bound z. Enumerating all independent
    # triples therefore certifies every vertex of this polyhedral slice.
    whole = set()
    for three in itertools.combinations(walls,3):
        z = kernel4([cov(r) for r in three])
        if not any(t.sign() for t in z):
            continue
        assert all(pair(z,r) == K() for r in three)
        for s in (1,-1):
            zz = scale(z,s)
            if all(pair(zz,r).sign() <= 0 for r in walls):
                assert zz[0].sign() > 0
                whole.add(tuple(t/zz[0] for t in zz))
    whole_counts = Counter(pair(z,z).sign() for z in whole)
    nulls = {z for z in whole if not pair(z,z).sign()}
    assert nulls == {vec((1,-1,0,0)),vec((1,0,-1,0))}
    assert whole_counts[-1] == 2

    # Complete 2D projective cell census: t=1 makes the old cone a compact
    # triangle, and every cell vertex has two independent active rows.
    upper = [caps[i] for i in (0,2,4)]
    bs = [tuple(z/r[-1] for z in r[:-1]) for r in upper]
    old_vertices = [vec(v) for v in [(1,0,0),(1,-1,0),(1,0,-1)]]
    for v in old_vertices:
        heights = [-pair(b,v) for b in bs]
        mm = heights[0]
        for z in heights[1:]:
            if (z-mm).sign() < 0:
                mm = z
        assert (mm*mm+pair(v,v)).sign() <= 0
    cells = []
    unique = set()
    fail = set()
    for i,b in enumerate(bs):
        rows = [cov(w) for w in old]
        rows += [cov(plus(c,scale(b,-1))) for j,c in enumerate(bs) if j != i]
        rays = set()
        for aa,bb in itertools.combinations(rows,2):
            z = cross(aa,bb)
            if not any(t.sign() for t in z):
                continue
            for s in (1,-1):
                zz = scale(z,s)
                if all(dot(a,zz).sign() <= 0 for a in rows):
                    assert zz[0].sign() > 0
                    rays.add(tuple(t/zz[0] for t in zz))
        assert rank(list(rays)) == 3
        center = tuple(sum((z[k] for z in rays),K()) for k in range(3))
        assert all(dot(a,center).sign() < 0 for a in rows)
        recs = []
        for z in sorted(rays,key=lambda v:json.dumps(ev(v))):
            q = pair(z,z)
            height = -pair(b,z)
            assert q.sign() <= 0 and height.sign() >= 0
            gap = height*height+q
            allh = [-pair(c,z) for c in bs]
            assert all((t-height).sign() >= 0 for t in allh)
            unique.add(z)
            if gap.sign() > 0:
                fail.add(z)
            recs.append({'ray':ev(z), 'norm':q.enc(), 'height':height.enc(),
                         'squared_covering_gap':gap.enc(),
                         'old_active':[j for j,w in enumerate(old) if not pair(w,z).sign()],
                         'tying_caps':[j for j,t in enumerate(allh) if t == height],
                         'covered':gap.sign() <= 0})
        cells.append({'cap':i,'rays':recs})
    assert len(unique) == 7 and len(fail) == 4
    assert sorted(len(c['rays']) for c in cells) == [3,3,5]
    return {'status':'PASS_EXACT_COXETER_COUNTEREXAMPLE',
            'field':'Q(sqrt(10)), positive square root; pairs [rational, sqrt(10) coefficient]',
            'form':[-1,1,1,1], 'norm_one_walls':[ev(r) for r in walls],
            'gram':[[x.enc() for x in row] for row in gram], 'pair_counts':dict(counts),
            'strict_center':ev(h),'actual_facet_points':facets,
            'closed_original_finite_edge_endpoints':finite,'bounded_original_cusps':links,
            'strict_escape':{'point':ev(x),'norm':pair(x,x).enc(),
                             'center_pair':pair(x,h).enc(),'wall_pairings':[p.enc() for p in xp]},
            'complete_projective_vertex_census':{
                'finite':whole_counts[-1],'ideal':whole_counts[0],
                'spacelike':whole_counts[1],
                'rays':[{'ray':ev(z),'norm':pair(z,z).enc()} for z in sorted(whole,key=lambda v:json.dumps(ev(v)))]},
            'lower_envelope_cells':cells, 'distinct_cell_rays':7,
            'uncovered_new_cell_rays':4,
            'scope':'All original finite vertices are covered and all original cusp links are bounded, but the nine-facet H3 Coxeter chamber has a strict spacelike escape. No finite-volume or novelty claim.'}


def read(path):
    return json.loads((ROOT/path).read_text())


def sha(path):
    return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def ipair(a,b,Q):
    return sum(x*sum(q*y for q,y in zip(row,b,strict=True)) for x,row in zip(a,Q,strict=True))


def graph_860():
    cert = read(PRIMARY)
    assert cert['status'].startswith('PASS_') and cert['facets'] == 860
    for path,digest in cert['sources'].items():
        assert sha(path) == digest
    old = read(RAW)['coxeter_input']
    Q = [[-int(x) for x in row] for row in old['basis_gram']]
    W = [[int(x) for x in row] for row in old['wall_coordinates_as_columns_transposed']]
    core = read(CORE)
    point_ids = core['incidence_parts'][0]
    vp = core['old_projections'][1]
    chosen = [c for c in read(BANK)['candidates'] if c['twice_height_square'] == 4]
    C = []
    for c in chosen:
        bits = [(c['mask']>>i)&1 for i in range(21)]
        assert sum(bits) == c['weight'] and c['A']**2-c['weight'] == 0
        twice = [c['A']*vp[j]-sum(bits[i]*W[point_ids[i]][j] for i in range(21)) for j in range(22)]
        assert twice == c['twice_old_projection'] and all(t%2 == 0 for t in twice)
        C.append([t//2 for t in twice])
    assert len(C) == len(set(map(tuple,C))) == 301
    covs = [[sum(Q[i][j]*v[j] for j in range(22)) for i in range(22)] for v in C]
    assert all(sum(x*y for x,y in zip(v,co,strict=True)) == 0 for v,co in zip(C,covs,strict=True))
    assert all(sum(x*y for x,y in zip(w,co,strict=True)) <= 0 for co in covs for w in W)
    ids = cert['selected_original_types']+[23,24]
    assert len(ids) == len(set(ids)) == 24
    adj = [set() for _ in range(325)]
    parent_counts = Counter()
    def edge(i,j):
        assert i != j
        adj[i].add(j)
        adj[j].add(i)
    for i in range(301):
        for j in range(i):
            p = sum(x*y for x,y in zip(C[i],covs[j],strict=True))+2
            assert p in (-2,-1,0)
            parent_counts[str(p)] += 1
            if p > -2:
                edge(i,j)
    S = int(cert['interval_denominator'])
    norms = {}
    inherited = 0
    for rec in cert['cap_records']:
        k = rec['type']
        norms[k] = F(rec['norm'])
        ix = 301+ids.index(k)
        finite = rec['finite_parent_labels']
        assert len(rec['all_parent_pairings']) == 812
        for j in range(210,511):
            lo,hi = map(int,rec['all_parent_pairings'][j])
            assert lo <= hi
            if str(j) in finite:
                assert finite[str(j)] in (2,3,4,6)
                edge(ix,j-210)
                inherited += 1
            else:
                assert hi < 0 and hi*hi > 2*norms[k]*S*S
    assert set(norms) == set(ids)
    checked = set()
    for row in cert['all_signed_cap_pairs']:
        i,j = row['i'],row['j']
        if i%2 or j%2:
            continue
        assert i > j
        a,b = i//2,j//2
        checked.add((a,b))
        lo,hi = map(int,row['pairing'])
        assert lo <= hi
        if row['kind'] == 'right_by_exact_frame_identity':
            assert lo <= 0 <= hi
            edge(301+a,301+b)
        else:
            assert row['kind'] == 'dotted'
            assert hi < 0 and hi*hi > norms[ids[a]]*norms[ids[b]]*S*S
    assert checked == {(i,j) for i in range(24) for j in range(i)}
    components = []
    seen = set()
    for i in range(325):
        if i in seen:
            continue
        stack = [i]
        seen.add(i)
        comp = []
        while stack:
            j = stack.pop()
            comp.append(j)
            for k in adj[j]-seen:
                seen.add(k)
                stack.append(k)
        comp.sort()
        components.append({'indices':comp,
                           'parent_null_caps':[j for j in comp if j < 301],
                           'additional_types':[ids[j-301] for j in comp if j >= 301]})
    components.sort(key=lambda c:(-len(c['indices']),c['indices']))
    edges = [(i,j) for i in range(325) for j in sorted(adj[i]) if j < i]
    assert len(edges) == 43900
    assert sorted(map(lambda c:len(c['indices']),components)) == [1]*17+[3,305]
    assert components[0]['additional_types'] == [0,2,6,22]
    assert components[1]['additional_types'] == [17,23,24]
    return {'status':'PASS_EXACT_860_UPPER_CAP_CONTACT_GRAPH',
            'source_sha256':{p:sha(p) for p in (PRIMARY,RAW,BANK,CORE)},
            'inherited_scope':'New-cap norms, heights, and exact finite-contact identities are inherited from the hash-bound primary 860 certificate. All parent-parent arithmetic and every absent new contact are checked here. This is not a separate normal reconstruction.',
            'upper_caps':325,'parent_pair_counts':dict(parent_counts),
            'finite_parent_new_contacts':inherited,
            'finite_contact_edges':len(edges),'components':components,
            'degrees':[len(a) for a in adj],
            'edge_list_sha256':hashlib.sha256(json.dumps(edges,separators=(',',':')).encode()).hexdigest(),
            'scope':'Nineteen finite-contact components obstruct full ball/horoball coverage. Any finite acute-compatible completion retaining these facets must add finite-angle bridges between these components. This does not exclude replacements or additions providing such bridges.'}


def run(verify):
    out = {'status':'PASS_COVERING_COUNTEREXAMPLE_AND_FINITE_CONTACT_FILTER',
           'counterexample':tiny_example(),'frozen_860_graph':graph_860()}
    if verify:
        assert out == json.loads(OUT.read_text())
    else:
        OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'])
    print('H3: 9 facets, 7 distinct cell rays, 4 uncovered new rays; strict escape 9/16.')
    print('H22: 325 upper caps, 43900 finite contacts, component sizes 305 + 3 + 17*1.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify',action='store_true')
    run(ap.parse_args().verify)
