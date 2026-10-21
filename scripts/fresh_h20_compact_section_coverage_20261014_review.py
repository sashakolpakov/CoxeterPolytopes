#!/usr/bin/env python3
"""Independent small witness replay; not a second census of all 162 recipes."""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import gcd, isqrt
from pathlib import Path
import json
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_compact_section_coverage_20261014'
RAW = 'data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
OUT = ROOT / f'results/{TAG}_review.json'
SCALE = 1 << 180
SAMPLES = [0, 27, 52, 55, 80, 107, 108, 135, 160, 60]


def clean(a):
    return {k:v for k,v in a.items() if v}


def load(a):
    return clean({int(k):F(v) for k,v in a.items()})


def rat(a):
    return {1:F(a)} if a else {}


def add(a,b):
    out = dict(a)
    for k,v in b.items():
        out[k] = out.get(k,F(0)) + v
    return clean(out)


def scale(a,c):
    return clean({k:v*c for k,v in a.items()})


def mul(a,b):
    out = {}
    for i,x in a.items():
        for j,y in b.items():
            g = gcd(i,j)
            k = (i//g)*(j//g)
            out[k] = out.get(k,F(0)) + g*x*y
    return clean(out)


def sign(a):
    if not a:
        return 0
    lo = hi = F(0)
    for n,c in a.items():
        b = isqrt(n*SCALE*SCALE)
        l = F(b,SCALE)
        h = l if b*b == n*SCALE*SCALE else F(b+1,SCALE)
        lo += c*(l if c>0 else h)
        hi += c*(h if c>0 else l)
    assert lo>0 or hi<0, 'Unresolved exact sign enclosure'
    return 1 if lo>0 else -1


def dot(a,b):
    assert len(a)==len(b)
    out = {}
    for x,y in zip(a,b):
        out = add(out,mul(x,y))
    return out


def mv(A,x):
    return [dot([rat(v) for v in row],x) for row in A]


def vs(x,c):
    return [mul(a,c) for a in x]


def va(x,y):
    assert len(x)==len(y)
    return [add(a,b) for a,b in zip(x,y)]


def scalar(a):
    assert not (set(a)-{1})
    return a.get(1,F(0))


def run():
    src = json.loads((ROOT/f'results/{TAG}.json').read_text())
    gate = json.loads((ROOT/f'results/{TAG}_gates.json').read_text())
    for data in [src,gate]:
        for p,h in data['sources'].items():
            assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    raw = json.loads((ROOT/RAW).read_text())['coxeter_input']
    Q = [[-int(v) for v in row] for row in raw['basis_gram']]
    W = [list(map(int,w)) for w in raw['wall_coordinates_as_columns_transposed']]
    assert len(W)==210 and all(len(w)==22 for w in W)
    C = [[sum(row[k]*w[k] for k in range(22)) for row in Q] for w in W]
    G = [[sum(w[k]*c[k] for k in range(22)) for c in C] for w in W]
    total = [sum(w[k] for w in W) for k in range(22)]
    g = gcd(*total)
    h = [F(v,g) for v in total]
    hc = [sum(row[k]*h[k] for k in range(22)) for row in Q]
    assert sum(h[k]*hc[k] for k in range(22))==-14
    assert max(sum(h[k]*c[k] for k in range(22)) for c in C)<0

    # Independent integer reconstruction of the modified-boundary input.
    table = Counter()
    short = [i for i in range(210) if G[i][i]==2]
    old_pairs = [(i,j) for i,j in combinations(short,2) if G[i][j]==0]
    for i,j in old_pairs:
        for k in range(210):
            if k in (i,j):
                continue
            d,a,b = G[k][k],G[i][k],G[j][k]
            if 2*d>a*a+b*b:
                table[(d,a,b)] += 1
    expected = {(r['norm'],r['pair_a'],r['pair_b']):r['count']
                for r in gate['spherical_boundary_table']}
    assert dict(table)==expected
    assert {p for p in table if p[1] and p[2]}=={(2,-1,-1)}
    for d,c in zip([G[i][i] for i in range(210)],C):
        assert all((2*v)%d==0 for v in c)

    selected = [src['cases'][i] for i in SAMPLES]
    assert {r['seed'] for r in selected}==set(range(9))
    assert {r['radical_sign'] for r in selected}=={-1,1}
    assert {r['parameter'] for r in selected}=={'1/5','1/2','4/5'}
    assert len({tuple(r['contact_squares']) for r in selected})==3
    assert len({r['defect'] for r in selected})==3
    radicands = set()
    for row in selected:
        for key in ['normal','positive_scaled_first_contact']:
            for v in row[key]:
                radicands.update(map(int,v))
    assert all(n>0 and all(v==1 for v in S.factorint(n).values()) for n in radicands)

    checked = []
    for row in selected:
        case = row['case']
        e = list(map(load,row['normal']))
        X = list(map(load,row['positive_scaled_first_contact']))
        P = [rat(F(v)) for v in row['initial_point']]
        assert all(len(v)==22 for v in [e,X,P])
        ec,xc = mv(Q,e),mv(Q,X)
        assert dot(e,ec)==rat(1)
        assert dot(X,ec)=={} and dot(P,ec)=={}
        qx = dot(X,xc)
        assert qx==load(row['contact_norm']) and sign(qx)<0
        future = dot(X,[rat(v) for v in hc])
        assert future==load(row['contact_center_pairing']) and sign(future)<0
        xp = mv(C,X)
        signs = list(map(sign,xp))
        assert signs==row['all210_wall_pairing_signs'] and max(signs)==0
        active = [i for i,v in enumerate(signs) if v==0]
        assert active==row['active_walls'] and len(active)==20

        profile = next(s for s in src['seeds']
                       if s['seed']==row['seed'] and s['contact_squares']==row['contact_squares'])
        E = profile['edge_walls']
        A = S.Matrix([[G[i][j] for j in E] for i in E])
        _,D = A.LDLdecomposition(hermitian=False)
        assert all(D[i,i]>0 for i in range(20))
        I = [[F(x) for x in rr] for rr in A.inv().tolist()]
        ep = mv(C,e)
        a = [ep[i] for i in E]
        v = mv(I,a)
        spatial = scalar(dot(a,v))
        gamma = 1-spatial
        assert 0<spatial<1 and gamma==F(profile['positive_residual_norm'])
        w = [dot([rat(W[j][d]) for j in E],v) for d in range(22)]
        z = va(e,vs(w,rat(-1)))
        assert dot(z,mv(Q,z))==rat(gamma)
        k = E.index(row['release'])
        wj = [rat(sum(I[k][l]*W[j][d] for l,j in enumerate(E))) for d in range(22)]
        dual = va(wj,vs(z,scale(v[k],-1/gamma)))
        dp = dot(dual,[rat(c) for c in C[row['new_wall']]])
        assert sign(dp)<0 and dot(dual,ec)=={}
        assert [dot(dual,[rat(c) for c in C[j]]) for j in E]==[rat(int(j==row['release'])) for j in E]
        pp = mv(C,P)
        assert [i for i,vv in enumerate(pp) if not vv]==E
        assert all(sign(pp[i])<0 for i in range(210) if i not in E)
        assert sign(dot(P,mv(Q,P)))<0
        rebuilt = va(vs(P,scale(dp,-1)),vs(dual,pp[row['new_wall']]))
        assert rebuilt==X
        assert active==sorted(set(E)-{row['release']}|{row['new_wall']})
        # H_E is PD by the positive rank-one Schur complement gamma.
        # The retained nineteen are independent, and dp != 0 supplies the
        # twentieth projected normal. Orthogonality to timelike X makes
        # this a PD20 actual vertex with every remaining wall strict.
        i,j = row['actual_ridge']
        assert i in active and j in active and i!=j
        pair = add(rat(G[i][j]),scale(mul(ep[i],ep[j]),-1))
        ni = add(rat(G[i][i]),scale(mul(ep[i],ep[i]),-1))
        nj = add(rat(G[j][j]),scale(mul(ep[j],ep[j]),-1))
        n,sq = mul(ni,nj),mul(pair,pair)
        assert pair==load(row['projected_pairing'])
        assert [ni,nj]==list(map(load,row['projected_norms']))
        assert sign(ni)>0 and sign(nj)>0 and sign(add(n,scale(sq,-1)))>0
        defect = row['defect']
        if defect=='obtuse':
            assert sign(pair)>0
        elif defect=='between_0_and_1/4':
            assert sign(pair)<0 and sign(sq)>0 and sign(add(n,scale(sq,-4)))>0
        else:
            assert defect=='between_1/4_and_1/2' and sign(pair)<0
            assert sign(add(scale(sq,4),scale(n,-1)))>0
            assert sign(add(n,scale(sq,-2)))>0
        checked.append(dict(case=case,seed=row['seed'],defect=defect,
                            actual_wall_count=len(active),wall_sign_count=len(signs)))
    sources = [RAW,f'scripts/{TAG}.py',f'results/{TAG}.json',
               f'scripts/{TAG}_gates.py',f'results/{TAG}_gates.json',
               f'scripts/{TAG}_review.py']
    return dict(status='PASS_SOURCE_PROOF_AND_TEN_INDEPENDENT_VERTEX_SAMPLES',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
                orthogonal_short_pairs=len(old_pairs),modified_boundary_type=[2,-1,-1],
                samples=checked,independent_wall_sign_count=210*len(checked),
                precision_bits=180,
                scope='Source/proof review of the compact orthogonal-short G2 obstruction and sparse-support scope. '
                      'Separate Fraction sparse-radical arithmetic replays ten saved A2 witnesses covering all nine source edges, '
                      'three ratios/positions/defects and both signs; SymPy supplies rational inverses only. '
                      'Not a second exhaustive 162-recipe replay or a compact-section classification.')


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--verify',action='store_true')
    args=ap.parse_args()
    result=run()
    if args.verify:
        assert result==json.loads(OUT.read_text())
    else:
        OUT.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'])
    print('Independent samples:',len(result['samples']),'wall signs:',result['independent_wall_sign_count'])
