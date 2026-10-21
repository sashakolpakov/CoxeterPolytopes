#!/usr/bin/env python3
"""Independent Fraction octic replay of the44/47-wall cap and cusp package.

No producer, helper, symbolic algebra, floating point or project imports.
Signs use rational square-root isolating intervals; inverses use Galois norms.
Finite20-supports use full matrix LDL, not the producer's forest algorithm.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
from math import isqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_new_projection_20261015_review'
PRIMARY = 'results/fresh_h20_new_projection_20261015_exact.json'
SOURCE = 'results/fresh_h20_joint_route_20261015_exact.json'
CUSPS = 'results/fresh_h20_new_projection_20261015_cusps.json'
GATE = 'results/fresh_h20_new_projection_20261015_null_gate.json'
ZERO = (F(0),) * 8
ONE = (F(1),) + (F(0),) * 7


def dec(x):
    if isinstance(x, (int, str, F)):
        return (F(x),) + (F(0),) * 7
    if isinstance(x, dict):
        assert set(x) == {'a', 'b'}
        return (F(x['a']), F(x['b'])) + (F(0),) * 6
    assert len(x) in (2, 4, 8)
    return tuple(map(F, x)) + (F(0),) * (8-len(x))


def add(x, y):
    return tuple(a+b for a, b in zip(x, y))


def neg(x):
    return tuple(-a for a in x)


def sub(x, y):
    return add(x, neg(y))


def mul(x, y):
    z = [F(0)] * 8
    for i, a in enumerate(x):
        if not a:
            continue
        for j, b in enumerate(y):
            if b:
                z[i ^ j] += a*b*(2 if i & j & 1 else 1)*(3 if i & j & 2 else 1)*(19 if i & j & 4 else 1)
    return tuple(z)


@lru_cache(None)
def inv(x):
    assert x != ZERO
    c = ONE
    for sigma in range(1, 8):
        c = mul(c, tuple(a*(-1 if (i & sigma).bit_count() % 2 else 1) for i, a in enumerate(x)))
    n = mul(x, c)
    assert n[1:] == (0,)*7 and n[0]
    return tuple(a/n[0] for a in c)


@lru_cache(None)
def radical_interval(d, bits):
    scale = 1 << bits
    q = isqrt(d * scale * scale)
    assert q*q <= d*scale*scale < (q+1)*(q+1)
    return F(q, scale), F(q+1, scale)


SIGN_PRECISIONS = Counter()


@lru_cache(None)
def sign(x):
    if x == ZERO:
        return 0
    for bits in (64, 128, 256, 512):
        lo = hi = x[0]
        for a, d in zip(x[1:], (2, 3, 6, 19, 38, 57, 114)):
            l, u = radical_interval(d, bits)
            lo += a * (l if a >= 0 else u)
            hi += a * (u if a >= 0 else l)
        if lo > 0 or hi < 0:
            SIGN_PRECISIONS[bits] += 1
            return 1 if lo > 0 else -1
    raise AssertionError(('unresolved exact sign', x))


def dot(x, y):
    z = ZERO
    for a, b in zip(x, y):
        z = add(z, mul(a, b))
    return z


def mv(M, x):
    return [dot(row, x) for row in M]


def inertia(M, require_positive=False):
    """Symmetric congruence; a zero diagonal is repaired by a literal basis sum."""
    a = [row[:] for row in M]
    counts = Counter()
    while a:
        n = len(a)
        j = next((i for i in range(n) if a[i][i] != ZERO), None)
        if j is None:
            pair = next(((i, j) for i, j in combinations(range(n), 2) if a[i][j] != ZERO), None)
            if pair is None:
                counts[0] += n
                break
            i, j = pair
            # Replace basis vector i by e_i+e_j; this is invertible.
            for k in range(n):
                a[i][k] = add(a[i][k], a[j][k])
            for k in range(n):
                a[k][i] = add(a[k][i], a[k][j])
            j = i
        if j:
            a[0], a[j] = a[j], a[0]
            for row in a:
                row[0], row[j] = row[j], row[0]
        s = sign(a[0][0])
        if require_positive and s <= 0:
            return False
        counts[s] += 1
        q = inv(a[0][0])
        a = [[sub(a[i][j], mul(mul(a[i][0], q), a[0][j])) for j in range(1, n)] for i in range(1, n)]
    if require_positive:
        return counts[0] == 0 and counts[-1] == 0
    return counts[1], counts[-1], counts[0]


def rank(M):
    a = [row[:] for row in M]
    k = 0
    for j in range(len(a[0])):
        i = next((i for i in range(k, len(a)) if a[i][j] != ZERO), None)
        if i is None:
            continue
        a[k], a[i] = a[i], a[k]
        q = inv(a[k][j])
        a[k] = [mul(v, q) for v in a[k]]
        for i in range(k+1, len(a)):
            q = a[i][j]
            if q != ZERO:
                a[i] = [sub(v, mul(q, w)) for v, w in zip(a[i], a[k])]
        k += 1
        if k == len(a):
            break
    return k


def run():
    def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
    seen=set()
    def read(p):
        d=json.loads((ROOT/p).read_text())
        if p not in seen:
            seen.add(p)
            for q,h in d.get('sources',{}).items():
                assert sha(q)==h,q
                if q.endswith('.json'):read(q)
        return d
    a,b,g,old=map(read,(PRIMARY,CUSPS,GATE,SOURCE))
    prior=read('results/fresh_h20_joint_route_20261015_exact_review.json')
    assert prior['status'].startswith('PASS')
    assert a['facets']==44 and b['facets']==47
    field=['1','sqrt2','sqrt3','sqrt6','sqrt19','sqrt38','sqrt57','sqrt114']
    assert a['field_basis']==b['field_basis']==field
    Q=[[dec(x)for x in row]for row in b['form']]
    N=[[dec(x)for x in row]for row in b['normals']]
    C=list(map(dec,b['strict_center']));E=list(map(dec,b['strict_escape']))
    assert len(Q)==len(C)==len(E)==21 and len(N)==47
    assert Q==[[dec(x)for x in row]for row in old['form']]
    assert N[:43]==[[dec(x)for x in row]for row in old['normals']]
    assert N[:44]==[[dec(x)for x in row]for row in a['normals']]
    assert C==list(map(dec,a['strict_center']))==list(map(dec,old['strict_center']))
    assert E==list(map(dec,a['strict_escape']))
    assert inertia(Q)==(20,1,0)
    cov=[mv(Q,n)for n in N];G=[[dot(n,c)for c in cov]for n in N]
    assert G==[[dec(x)for x in row]for row in b['gram']]
    assert [row[:44]for row in G[:44]]==[[dec(x)for x in row]for row in a['gram']]
    norms=[G[i][i]for i in range(47)]
    assert all(q==dec(2)for q in norms[42:])
    pairlabels={};counts=Counter();newcounts=Counter();capcounts=Counter()
    for i,j in combinations(range(47),2):
        p=G[i][j];assert sign(p)<=0
        t=mul(dec(4),mul(mul(p,p),inv(mul(norms[i],norms[j]))))
        if p==ZERO:label='right'
        elif t==ONE:label='pi3'
        elif t==dec(2):label='pi4'
        elif t==dec(3):label='pi6'
        elif t==dec(4):label='parallel'
        else:assert sign(sub(t,dec(4)))>0;label='dotted'
        pairlabels[i,j]=label;counts[label]+=1
        if j>=44:newcounts[label]+=1
        if j==43:capcounts[label]+=1
    assert dict(counts)==b['all_pair_counts'] and dict(newcounts)==b['new_pair_counts']
    assert dict(capcounts)==a['new_pair_counts']
    print('independent: all47 norms and1081 pairs pass',flush=True)
    # Reconstruct the new source projection by its defining full-rank system.
    recipe=a['cap_recipe'];Y=list(map(dec,recipe['source_projection']))
    assert Y==N[43][:19] and len(Y)==19
    ids=recipe['source_rights']+recipe['source_anchors']
    assert ids==[1,2,3,4,5,7,8,9,10,11,13,14,15,16,18,20,21,12,37]
    assert rank([cov[i][:19]for i in ids])==19
    assert [dot(Y,cov[i][:19])for i in ids]==[ZERO]*17+[dec(-1)]*2
    Q0=[r[:19]for r in Q[:19]]
    assert dot(Y,mv(Q0,Y))==dec('-15/2')
    assert N[43][19:]==[dec([0,0,0,0,'-1/2',0,0,0]),ZERO]
    assert sign(dot(Y,mv(Q0,C[:19])))<0
    cp=[dot(C,c)for c in cov];qC=dot(C,mv(Q,C))
    assert sign(qC)<0 and all(sign(p)<0 for p in cp)
    assert cp==list(map(dec,b['center_pairings'])) and qC==dec(b['center_norm'])
    # Every supplied facet witness is checked in the full47-wall chamber.
    for i,r in enumerate(b['facet_witnesses']):
        p=list(map(dec,r['point']))
        expected=[sub(x,mul(mul(cp[i],inv(norms[i])),y))for x,y in zip(C,N[i])]
        assert p==expected and dot(p,mv(Q,p))==dec(r['norm'])
        products=[dot(p,c)for c in cov]
        assert products==list(map(dec,r['pairings']))
        assert sign(dot(p,mv(Q,p)))<0
        assert [sign(x)for x in products]==[-1]*i+[0]+[-1]*(46-i)
    for i,row in enumerate(a['facet_points']):
        p=list(map(dec,row));products=[dot(p,c)for c in cov[:44]]
        assert sign(dot(p,mv(Q,p)))<0
        assert [sign(x)for x in products]==[-1]*i+[0]+[-1]*(43-i)
    print('independent: facet actuality passes',flush=True)
    # Full20x20 matrix congruence, independent of the producer's tree test.
    neighbors=[i for i in range(43)if pairlabels[i,43]in('right','pi3','pi4','pi6')]
    assert neighbors==a['complete_new_cap_vertex_census']['finite_neighbors'] and len(neighbors)==20
    positive=[];tested=0
    for subset in combinations(neighbors,19):
        S=list(subset)+[43];tested+=1
        if inertia([[G[i][j]for j in S]for i in S],True):positive.append(S)
    assert tested==20 and positive==a['complete_new_cap_vertex_census']['positive_supports']
    assert len(positive)==1
    supports=[]
    for row in a['finite_edge']['endpoints']+a['retained_previous_edge_vertices']:
        S=row['support'];v=list(map(dec,row['point']))
        assert inertia([[G[i][j]for j in S]for i in S],True)
        assert sign(dot(v,mv(Q,v)))<0 and dot(v,mv(Q,C))==dec(-1)
        assert [sign(dot(v,c))for c in cov]==[0 if i in S else -1 for i in range(47)]
        assert [dot(v,c)for c in cov[:44]]==list(map(dec,row['pairings']))
        supports.append(S)
    assert len({tuple(s)for s in supports})==4
    edge=sorted(set(supports[0])&set(supports[1]))
    assert edge==a['finite_edge']['shared_support'] and len(edge)==19
    assert inertia([[G[i][j]for j in edge]for i in edge],True)
    assert [r['support']for r in b['retained_finite_vertices']]==supports
    for row in b['retained_finite_vertices']:
        v=list(map(dec,row['point']))
        assert [dot(v,c)for c in cov]==list(map(dec,row['pairings']))
    oldE=list(map(dec,old['strict_escape']))
    assert dot(oldE,cov[43])==dec(a['old_escape_cut']) and sign(dot(oldE,cov[43]))>0
    ep=[dot(E,c)for c in cov];qE=dot(E,mv(Q,E))
    assert sign(qE)>0 and all(sign(p)<0 for p in ep)and sign(dot(E,mv(Q,C)))<0
    assert ep==list(map(dec,b['escape_pairings'])) and qE==dec(b['escape_norm'])
    assert ep[:44]==list(map(dec,a['escape_pairings']))
    face=list(map(dec,a['full_H18_source_face']['point']))
    assert [sign(dot(face,c))for c in cov[:44]]==[-1]*40+[0,0,-1,-1]
    cycle=mul(dec(4),mul(mul(G[0][36],G[0][36]),inv(mul(norms[0],norms[36]))))
    host_minor=sub(mul(norms[38],norms[39]),mul(G[38][39],G[38][39]))
    assert cycle==dec([12,8])==dec(a['conditional_nonquasi']['cartan_two_cycle'])
    assert host_minor==dec(-32)==dec(a['conditional_nonquasi']['determinant'])
    # Explicit three null-source caps and bounded Euclidean links.
    W=list(map(dec,b['cap_recipe']['common_transverse']))
    assert W[:19]==[ZERO]*19 and W[19:]==[dec([0,'-1/2',0,0])]*2
    assert dot(W,mv(Q,W))==dec(2)
    cusp_rays=[];component_count=0
    for k,row in enumerate(b['repaired_cusps']):
        i,j=row['source_pair'];c=[add(x,y)for x,y in zip(N[i],N[j])]
        cusp_rays.append(c)
        assert c==list(map(dec,row['ray'])) and N[44+k]==[add(x,y)for x,y in zip(c,W)]
        assert dot(c,mv(Q,c))==ZERO and sign(dot(c,mv(Q,C)))<0
        prods=[dot(c,v)for v in cov];inc=[i for i,p in enumerate(prods)if p==ZERO]
        assert inc==row['incident'] and len(inc)==23 and all(sign(x)<=0 for x in prods)
        assert prods==list(map(dec,row['pairings']))
        cg=[[G[i][j]for j in inc]for i in inc]
        assert inertia(cg)==(19,0,4)and rank([N[i]for i in inc])==20
        components=row['positive_affine_components'];parts=[]
        for comp in components:
            ids=comp['walls'];parts+=ids;ker=list(map(dec,comp['positive_kernel']))
            M=[[G[i][j]for j in ids]for i in ids]
            assert mv(M,ker)==[ZERO]*len(ids)and all(sign(x)>0 for x in ker)
            assert inertia([r[:-1]for r in M[:-1]],True)
            assert inertia(M)==(len(ids)-1,0,1)
            physical=[dot(ker,[N[i][j]for i in ids])for j in range(21)]
            mu=dec(comp['physical_kernel_multiple']);assert sign(mu)>0
            assert physical==[mul(mu,x)for x in c]
            component_count+=1
        assert sorted(parts)==inc and sorted(len(x['walls'])-1 for x in components)==[1,2,8,8]
        for x,y in combinations(components,2):assert all(G[i][j]==ZERO for i in x['walls']for j in y['walls'])
        assert row['minimum_scale_anchor']in row['pi3_source_anchors']
        anchor=row['minimum_scale_anchor'];assert norms[anchor]==dec(2)and prods[anchor]==dec(-1)
        assert [i for i in range(40)if pairlabels[i,44+k]=='pi3']==row['pi3_source_anchors']
        assert all(pairlabels[i,44+k]=='pi4'for i in(40,41))
        ce=dot(c,mv(Q,E));assert sign(ce)<0 and sign(sub(mul(ce,ce),dec(2)))>0
    assert component_count==12
    assert all(dot(c,mv(Q,e))==dec(-2)for c,e in combinations(cusp_rays,2))
    print('independent: one new finite vertex, finite edge and three rank19 cusp links pass',flush=True)
    # Complete literal parallel-pair ray list and every all-profile scale bound.
    NS=[n[:19]for n in N[:40]];ES=E[:19];ET=E[19:];HS=C[:19]
    assert sign(dot(ES,mv(Q0,ES)))<0 and sign(dot(ES,mv(Q0,HS)))<0
    assert mul(dec(2),dot(ET,ET))==ONE
    sqrt2=dec([0,1,0,0]);rays=[];rayrows=[]
    for i,j in combinations(range(40),2):
        if sign(G[i][j])>=0 or mul(G[i][j],G[i][j])!=mul(norms[i],norms[j]):continue
        assert norms[i]in(ONE,dec(2))and norms[j]in(ONE,dec(2))
        c=[add(mul(x,inv(ONE if norms[i]==ONE else sqrt2)),
               mul(y,inv(ONE if norms[j]==ONE else sqrt2)))for x,y in zip(NS[i],NS[j])]
        assert dot(c,mv(Q0,c))==ZERO and sign(dot(c,mv(Q0,HS)))<0
        products=[dot(c,v[:19])for v in cov[:40]];assert all(sign(x)<=0 for x in products)
        bounds=[]
        for k,p in enumerate(products):
            if p!=ZERO:
                first=ONE if norms[k]==dec(2)else mul(sqrt2,dec('1/2'))
                bounds.append((k,mul(first,inv(neg(p)))))
        anchor,L=bounds[0]
        for k,value in bounds[1:]:
            if sign(sub(value,L))>0:anchor,L=k,value
        ce=dot(c,mv(Q0,ES));margin=sub(mul(mul(L,L),mul(ce,ce)),dec(2))
        assert sign(L)>0 and sign(ce)<0 and sign(margin)>0
        row=g['ray_certificates'][len(rays)]
        assert row['source_parallel_pair']==[i,j]and row['mandatory_scale_anchor']==anchor
        assert list(map(dec,row['null_ray']))==c and list(map(dec,row['source_pairings']))==products
        assert dec(row['scale_lower_bound'])==L and dec(row['source_escape_pairing'])==ce
        assert dec(row['squared_escape_margin'])==margin
        rays.append(c);rayrows.append(dict(pair=[i,j],anchor=anchor,margin=[str(x)for x in margin]))
    assert len(rays)==31
    assert all(sign(dot(c,mv(Q0,e)))<0 for c,e in combinations(rays,2))
    margins=[dec(r['squared_escape_margin'])for r in g['ray_certificates']]
    least=g['minimum_margin_index'];assert all(sign(sub(x,margins[least]))>=0 for x in margins)
    assert dec(g['minimum_margin'])==margins[least]
    # At t=0, old support inequalities require p19,p20<0 for an interior;
    # acuteness requires w19,w20<=0, and w²=2 makes the new inequality impossible.
    assert N[40]==[ZERO]*19+[ONE,ZERO] and N[41]==[ZERO]*20+[ONE]
    paths=[PRIMARY,CUSPS,GATE,SOURCE,'results/fresh_h20_joint_route_20261015_exact_review.json',
           'scripts/fresh_h20_new_projection_20261015_exact.py',
           'scripts/fresh_h20_new_projection_20261015_cusps.py',
           'scripts/fresh_h20_new_projection_20261015_null_gate.py',f'scripts/{TAG}.py']
    return dict(status='PASS_INDEPENDENT_OCTIC47_GEOMETRY_AND31_NULL_RAY_PROFILE_GATE',
                sources={p:sha(p)for p in paths},arithmetic='Standalone Fraction bitmask field Q(sqrt2,sqrt3,sqrt19), Galois-norm inversion and integer-square-root interval signs; no producer/project imports.',
                signature=[20,1],walls=47,pairs=1081,pair_counts=dict(counts),
                actual_facets=47,new_cap_candidate_supports=tested,new_cap_finite_vertices=len(positive),
                verified_finite_vertices=4,verified_finite_edges=1,bounded_cusp_links=3,
                cusp_link_rank=19,positive_affine_components=12,parallel_source_rays=31,
                distinct_ray_pair_signs=465,ray_bounds=rayrows,
                strict_spacelike_escape=True,sign_precisions={str(k):v for k,v in sorted(SIGN_PRECISIONS.items())},
                scope='Independent full arithmetic for the47-wall package, the unique finite vertex on cap43, four retained finite vertices/one edge, three actual bounded cusp links and all-profile escape inequalities on the31 literal source parallel-pair rays. No full47-wall finite/ideal vertex census, no classification of all source null directions, no finite-volume claim, and no independent proof of the inherited external quasi-arithmetic criterion.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'])

