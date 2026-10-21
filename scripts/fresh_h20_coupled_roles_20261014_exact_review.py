#!/usr/bin/env python3
"""Independent old-field and normalized-frame replay of the seven-cap model."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from math import isqrt
from pathlib import Path
from fresh_h20_bridge_edge_closure_20261014_review import Field,parse
import fresh_h20_bridge_faces_20261014_review as B
ROOT=Path(__file__).resolve().parents[1]
BASE='fresh_h20_coupled_roles_20261014_exact';TAG=BASE+'_review'
FILES=['results/fresh_h20_e7_edge_cusp_repair_20261014_exact.json','results/fresh_h20_e7_bridge_face_20261014_exact.json','results/fresh_h20_bridge_edge_closure_20261014.json','results/fresh_h20_polar_pair_20261014.json','results/fresh_h20_e7_repaired_cusp_walk_20261014.json']
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
    old,bridge,edge,polar,walk=[json.loads((ROOT/p).read_text()) for p in FILES]
    src=f'results/{BASE}.json';d=json.loads((ROOT/src).read_text())
    for obj in [old,bridge,edge,polar,walk,d]:
        for p,h in obj['sources'].items():assert sha(p)==h,p
    ff=Field([parse(old['field']['discriminant'],0),parse(bridge['field']['new_positive_sqrt_radical'],1),parse(edge['field']['new_positive_sqrt_radical'],2)])
    add,mul,sc,total,sg=ff.add,ff.mul,ff.scale,ff.total,ff.sign;c=ff.const;zero=ff.zeros[3]
    for k,r in enumerate(ff.radicals):assert sg(r,k)>0
    Q=[list(map(F,row)) for row in old['form']];assert Q==[list(map(F,row)) for row in d['form']]
    rows=[[(j,q) for j,q in enumerate(row) if q] for row in Q];cache={}
    def vp(v,k):return tuple(ff.lift(parse(x,k),k) for x in v)
    def dot(v,w):
        if id(w) not in cache:cache[id(w)]=(w,tuple(total(sc(w[j],q) for j,q in row) for row in rows))
        return total(mul(x,y) for x,y in zip(v,cache[id(w)][1]))
    def lin(vs,cs):return tuple(total(mul(a,v[j]) for a,v in zip(cs,vs)) for j in range(21))
    def rat(x):
        for k in [3,2,1]:assert x[1]==ff.zeros[k-1];x=x[0]
        assert not any(x[1:]);return x[0]
    def eq(x,raw,k=3):assert x==ff.lift(parse(raw,k),k)
    N=[vp(v,1) for v in old['normals']]+[vp(bridge['new_normal'],2)]+[vp(v,3) for v in edge['new_normals']]+[tuple(c(x) for x in v) for v in polar['new_normals']]
    norms=[c(q) for q in old['norms']]+[c(2)]+[parse(edge['new_normal_norm'],3)]*2+[c(q) for q in polar['new_norms']]
    assert len(N)==len(norms)==100
    T=vp(polar['strict_center']['point'],1);assert sg(dot(T,T))<0
    root2=ff.lift((F(0),F(1))+(F(0),)*14,0)
    def projection(S,gramraw,inverseraw):
        roots=[lin([N[i]],[root2 if i==92 else c(1)]) for i in S]
        G=[[rat(dot(v,w)) for w in roots] for v in roots];B.pd(G);J=B.inv(G)
        assert G==[list(map(F,row)) for row in gramraw]
        assert J==[list(map(F,row)) for row in inverseraw]
        def project(v):
            pairs=[dot(v,n) for n in roots]
            return lin([v]+roots,[c(1)]+[sc(total(sc(g,q) for g,q in zip(pairs,row)),-1) for row in J])
        return project,G,J
    r=d['cap_R'];S=r['support'];assert S==walk['edge_germs'][0]['support']
    assert walk['edge_germs'][0]['type']=='open_opposite_end'
    project,G,J=projection(S,r['support_gram'],r['support_inverse'])
    H0=lin([N[i] for i in list(range(20))+[24]],[c(1)]*21);assert H0==vp(r['H0'],1)
    assert dot(H0,H0)==c(-2) and sg(dot(H0,T))<0
    nu=vp(r['nu'],1);assert nu==vp([[[str(x) for x in v],['0']*16] for v in old['completed_cusp']['point']],1)
    assert dot(nu,nu)==zero and all(dot(nu,N[i])==zero for i in S)
    W=project(H0);assert W==vp(r['W'],1) and all(dot(W,N[i])==zero for i in S)
    h=dot(W,nu);qW=dot(W,W);eq(h,r['h'],1);assert sg(h)<0 and sg(qW)<0
    Z=lin([W,nu],[sc(h,-2),qW]);assert Z==vp(r['Z'],1)
    assert dot(Z,Z)==zero and dot(Z,nu)==sc(mul(h,h),-2)
    R=lin([Z,nu],[c(1),c(-1)]);P=lin([Z,nu],[c(1),c(1)])
    assert R==vp(r['normal'],1) and P==vp(r['finite_vertex'],1)
    qR=sc(mul(h,h),4);assert dot(R,R)==qR and dot(P,P)==sc(qR,-1) and dot(P,R)==zero
    eq(qR,r['norm'],1);assert sg(qR)>0 and sg(dot(P,T))<0
    counts=Counter();rp=[];pp=[]
    for i,n in enumerate(N):
        g=dot(R,n);p=dot(P,n);rp.append(g);pp.append(p)
        if i in S:assert g==p==zero;counts['right']+=1
        else:
            assert sg(g)<0 and sg(p)<0 and sg(add(mul(g,g),sc(mul(qR,norms[i]),-1)))>0;counts['dotted']+=1
        eq(g,r['all100_pairings'][i]);eq(p,r['vertex_old100_pairings'][i])
    Eold=tuple(c(q) for q in polar['strict_escape']['rational_point']);assert sg(dot(R,Eold))>0
    eq(dot(R,Eold),r['old_escape_pairing'],1)
    N.append(R);norms.append(qR)
    b=d['hex_basis'];S75=b['support'];assert S75==[1,2,3,4,5,6,7,8,9,10,12,13,14,15,17,31,75,87]
    p75,G75,J75=projection(S75,b['support_gram'],b['support_inverse'])
    C=p75(lin([N[0],N[65],N[79]],[c(2),c(1),c(2)]));assert C==vp(b['C'],1)
    basis=[vp(b[x],1) for x in ['C','U','V']];squares=[ff.lift(parse(x,1),1) for x in b['squared_lengths']]
    dd,ss,tt=squares;assert all(sg(q)>0 for q in squares) and sg(dot(C,T))<0
    assert [[dot(v,w) for w in basis] for v in basis]==[[sc(dd,-1),zero,zero],[zero,ss,zero],[zero,zero,tt]]
    assert all(dot(v,N[i])==zero for v in basis for i in S75)
    # Independently reconstruct the denominator-free Gram-Schmidt recipe.
    v=p75(N[0]);U=lin([v,C],[dd,dot(v,C)]);assert U==basis[1]
    v=p75(N[b['third_basis_anchor']]);V=lin([v,C,U],[mul(dd,ss),mul(ss,dot(v,C)),sc(mul(dd,dot(v,U)),-1)])
    assert V==basis[2]
    coeff=[tuple(tuple(map(F,x)) for x in co) for co in d['hexagon']['frame_coefficients']]
    def frame(a,b):return B.ba(B.bs(B.bm(a[0],b[0]),-1),B.ba(B.bm(a[1],b[1]),B.bm(a[2],b[2])))
    gram=[[frame(a,b) for b in coeff] for a in coeff]
    expected=[[2 if i==j else 0 if (i-j)%6 in [1,5] else -6 if (i-j)%6==3 else -4 for j in range(6)] for i in range(6)]
    assert gram==[[B.c(x) for x in row] for row in expected] and expected==d['hexagon']['gram']
    @lru_cache(None)
    def bounds(x,bits):return ff.interval(x,3,bits)
    def sqrti(iv,bits):
        lo,hi=iv;assert lo>0;den=1<<bits
        a=isqrt(lo.numerator*den*den//lo.denominator);b=isqrt(hi.numerator*den*den//hi.denominator)+1
        assert F(a*a,den*den)<=lo and F(b*b,den*den)>hi
        return F(a,den),F(b,den)
    @lru_cache(None)
    def roots(bits):return tuple(sqrti(bounds(q,bits),bits) for q in squares)
    maxbits=0
    def normalized(co,pairs,dotted_norm=None):
        nonlocal maxbits
        if all(x==B.Z or y==zero for x,y in zip(co,pairs)):return 0
        for bits in [192,384,768,1536,3072]:
            lo=hi=F(0)
            for x,y,(a,b) in zip(co,pairs,roots(bits)):
                if x==B.Z or y==zero:continue
                l,u=bounds(mul(ff.lift(x,0),y),bits)
                vv=[l/a,l/b,u/a,u/b];lo+=min(vv);hi+=max(vv)
            if dotted_norm is not None:
                qlo,qhi=bounds(dotted_norm,bits)
                if hi<0 and hi*hi>2*qhi:maxbits=max(maxbits,bits);return -1
            elif lo>0 or hi<0:maxbits=max(maxbits,bits);return 1 if lo>0 else -1
        raise AssertionError('unresolved normalized-frame sign')
    bp=[[dot(v,n) for v in basis] for n in N]
    for k,co in enumerate(coeff):
        for i,(pairs,q) in enumerate(zip(bp,norms)):
            kind=d['hexagon']['new_to_old101_pair_checks'][k][i]['relation']
            if i in S75:assert pairs==[zero]*3 and kind=='right';counts['right']+=1
            else:assert normalized(co,tuple(pairs),q)<0 and kind=='dotted';counts['dotted']+=1
    counts['right']+=6;counts['dotted']+=9
    assert dict(counts)==d['new_pair_counts']=={'right':133,'dotted':588}
    assert all(normalized(co,tuple(dot(Eold,v) for v in basis))<0 for co in coeff)
    print('PASS_ALL721_NEW_PAIRS',flush=True)
    H=vp(d['strict_center']['point'],1);assert H==lin([C,T],[c(1),c(d['strict_center']['epsilon'])])
    assert sg(dot(H,H))<0 and sg(dot(H,T))<0 and all(sg(dot(H,n))<0 for n in N)
    hp=[dot(H,v) for v in basis];assert all(normalized(co,tuple(hp))<0 for co in coeff)
    vertices=[]
    for k,row in enumerate(d['hexagon']['finite_vertices']):
        co=tuple(tuple(map(F,x)) for x in row['frame_coefficients']);assert frame(co,co)==B.c(F(-1,3))
        assert co[0]==B.c(1)
        for i,pairs in enumerate(bp):assert normalized(co,tuple(pairs))==(0 if i in S75 else -1)
        pp=[frame(co,v) for v in coeff];ids=[i for i,x in enumerate(pp) if x==B.Z]
        assert ids==sorted([k,(k+1)%6]) and all(sg(ff.lift(x,0))<0 for i,x in enumerate(pp) if i not in ids)
        assert row['incident']==S75+[101+i for i in ids]
        vertices.append(row['incident'])
    assert len({tuple(x) for x in vertices})==6
    assert all(normalized(co,tuple(dot(P,v) for v in basis))<0 for co in coeff)
    # All actual facets follow from this strict common center, positive norms,
    # and the independently checked acute pair matrix by orthogonal projection.
    oldcusps=[old['completed_cusp']['point']]+[x['point'] for x in old['preserved_old_bounded_cusps']]
    for raw,row in zip(oldcusps,d['preserved_bounded_cusps']):
        nu0=tuple(ff.lift(tuple(map(F,x)),0) for x in raw)
        assert nu0==vp(row['point'],1) and dot(nu0,nu0)==zero and sg(dot(nu0,H))<0
        assert sg(dot(R,nu0))<0 and all(normalized(co,tuple(dot(nu0,v) for v in basis))<0 for co in coeff)
    neighbors=[S]+[S75+[101+(k-1)%6,101+(k+1)%6] for k in range(6)]
    assert neighbors==d['complete_new_finite_neighbors']
    assert d['complete_new_vertex_count']==7 and d['known_total_vertices_at_least']==141
    assert len(d['preserved_bounded_cusps'])==3
    print('PASS_SEVEN_ACTUAL_VERTICES_COMMON_CENTER_AND_THREE_CUSPS',flush=True)
    # Check the second face remains unaffected by these seven additions.
    S61=sorted(set(S75)-{75}|{61});roots61=[N[i] for i in S61]
    G61=[[rat(dot(v,w)) for w in roots61] for v in roots61];B.pd(G61);J61=B.inv(G61)
    assert all(x>=0 for row in J61 for x in row) and J61[S61.index(61)][S61.index(61)]==F(2,3)
    lim=d['S61_limitation'];assert lim['support']==S61 and J61==[list(map(F,row)) for row in lim['inverse_gram']]
    upper=add(qR,sc(mul(dot(R,N[61]),dot(R,N[61])),F(-2,3)))
    eq(upper,lim['R_projected_norm_upper_bound'],1);assert sg(add(upper,sc(qR,F(1,3))))<0
    assert F(lim['hex_norm_upper_bound'])==F(-2,3)
    pairs=[dot(H,n) for n in roots61]
    HC=lin([H]+roots61,[c(1)]+[sc(total(sc(g,q) for g,q in zip(pairs,row)),-1) for row in J61])
    assert sg(dot(HC,HC))<0 and sg(dot(HC,H))<0 and sg(dot(R,HC))<0
    assert all(normalized(co,tuple(dot(HC,v) for v in basis))<0 for co in coeff)
    assert all(x>=0 for row in J75 for x in row) and J75[S75.index(75)][S75.index(75)]==5
    R75upper=add(qR,sc(mul(dot(R,N[75]),dot(R,N[75])),-5))
    assert sg(add(R75upper,sc(qR,9)))<0 and sg(dot(R,p75(H)))<0
    # A literal rational witness independently proves ambient infinite volume.
    esc=d['strict_escape'];E=tuple(c(q) for q in esc['rational_point'])
    qE=rat(dot(E,E));assert qE==F(esc['norm'])>0
    eq(dot(E,H),esc['future_center_pairing'],1);assert sg(dot(E,H))<0
    for i,n in enumerate(N):
        g=dot(E,n);eq(g,esc['old101_pairings'][i]);assert sg(g)<0
    assert all(normalized(co,tuple(dot(E,v) for v in basis))<0 for co in coeff)
    print('PASS_RENEWED_STRICT_ESCAPE_AND_S61_REDUNDANCY',flush=True)
    paths=FILES+[src,f'scripts/{BASE}.py','scripts/fresh_h20_coupled_roles_20261014_field.py',f'scripts/{TAG}.py','scripts/fresh_h20_bridge_edge_closure_20261014_review.py','scripts/fresh_h20_bridge_faces_20261014_review.py']
    return dict(status='PASS_INDEPENDENT_107_WALL_COUPLED_GEOMETRY_AND_STRICT_ESCAPE',sources={p:sha(p) for p in paths},new_pair_counts=dict(counts),new_vertices=7,strict_old_cusps=3,max_normalization_bits=maxbits,escape_norm=str(qE),S61_all_seven_redundant=True,R_redundant_on_S75=True,
        implementation='No new producer arithmetic imports. Reuses previously independent sparse-radical tower; separately reconstructs projections/recipes and outward bounds for normalized Lorentz frame.',
        source_face='Exact compact right hexagon with all six actual finite endpoints; ambient chamber remains infinite volume by a strict rational spacelike escape.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],flush=True)
