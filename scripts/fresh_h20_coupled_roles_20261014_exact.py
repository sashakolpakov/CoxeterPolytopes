#!/usr/bin/env python3
"""A cusp-edge escape cutter and a compact right-angled S75 hexagon.

Exact vector recipes, formal Gram identities, and outward rational signs.
The normalized hexagon frame is kept unexpanded over its three square roots.
"""
import argparse,json
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
from collections import Counter
from fresh_h20_coupled_roles_20261014_field import Field,parse,enc,pd,B
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_coupled_roles_20261014_exact'
FILES=['results/fresh_h20_e7_edge_cusp_repair_20261014_exact.json','results/fresh_h20_e7_bridge_face_20261014_exact.json','results/fresh_h20_bridge_edge_closure_20261014.json','results/fresh_h20_polar_pair_20261014.json','results/fresh_h20_e7_repaired_cusp_walk_20261014.json']
S75=[1,2,3,4,5,6,7,8,9,10,12,13,14,15,17,31,75,87]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
    old,bridge,edge,polar,walk=[json.loads((ROOT/p).read_text())for p in FILES]
    for d in[old,bridge,edge,polar,walk]:
        for p,h in d['sources'].items():assert sha(p)==h,p
    A=Field(old,bridge,edge);z=A.z
    M=[tuple(parse(x,1)for x in v)for v in old['normals']]
    N=[A.vlift(v,1)for v in M]+[A.vlift(tuple(parse(x,2)for x in bridge['new_normal']),2)]
    N +=[tuple(parse(x,3)for x in v)for v in edge['new_normals']]
    N +=[tuple(A.c(x)for x in v)for v in polar['new_normals']]
    norms=[A.c(q)for q in old['norms']]+[A.c(2)]+[parse(edge['new_normal_norm'],3)]*2+[A.c(q)for q in polar['new_norms']]
    assert len(N)==len(norms)==100
    T=tuple(parse(x,1)for x in polar['strict_center']['point']);qT=A.dot(T,T,1)
    assert qT==parse(polar['strict_center']['norm'],1)and A.sign(qT,1)<0
    def projection(S):
        # The germ contains one B2 pair. Rescale its short-normal coordinate
        # by sqrt(2), preserving the subspace and making the Gram rational.
        scales=[(F(0),F(1))+(F(0),)*14 if i==92 else B.sc(1)for i in S]
        roots=[A.lin([M[i]],[A.lift(c,0,1)],1)for i,c in zip(S,scales)]
        G=[[A.rational(A.dot(v,w,1))for w in roots]for v in roots];assert pd(G)
        J=B.inv(G)
        def project(v):
            pairs=[A.dot(v,n,1)for n in roots]
            cs=[A.total((A.scale(g,q,1)for g,q in zip(pairs,row)),1)for row in J]
            ans=A.lin([v]+roots,[A.c(1,1)]+[A.scale(c,-1,1)for c in cs],1)
            assert all(A.dot(ans,M[i],1)==z[1]for i in S)
            return ans
        return project,G,J
    # A fresh literal support, not one of the four earlier excluded families.
    germ=walk['edge_germs'][0];S=germ['support'];assert germ['type']=='open_opposite_end'
    assert S==[2,3,4,5,6,7,13,18,31,65,79,83,84,88,89,90,91,92,94]
    project,G,J=projection(S)
    H0=tuple(A.total((M[i][j]for i in list(range(20))+[24]),1)for j in range(21))
    assert all(A.rational(x)==A.rational(x)for x in H0)
    assert A.dot(H0,H0,1)==A.c(-2,1)and A.sign(A.dot(H0,T,1),1)<0
    nu=tuple(A.lift(parse(x,0),0,1)for x in old['completed_cusp']['point'])
    assert A.dot(nu,nu,1)==z[1]and all(A.dot(nu,M[i],1)==z[1]for i in S)
    W=project(H0);h=A.dot(W,nu,1);w2=A.dot(W,W,1);assert A.sign(h,1)<0 and A.sign(w2,1)<0
    Z=A.lin([W,nu],[A.scale(h,-2,1),w2],1)
    h2=A.mul(h,h,1);assert A.dot(Z,Z,1)==z[1]and A.dot(Z,nu,1)==A.scale(h2,-2,1)
    assert A.sign(A.dot(Z,H0,1),1)<0
    R=A.lin([Z,nu],[A.c(1,1),A.c(-1,1)],1);P=A.lin([Z,nu],[A.c(1,1),A.c(1,1)],1)
    qR=A.scale(h2,4,1);assert A.dot(R,R,1)==qR and A.dot(P,P,1)==A.scale(qR,-1,1)
    assert A.dot(R,P,1)==z[1]and A.sign(A.dot(P,T,1),1)<0
    R3=A.vlift(R,1);P3=A.vlift(P,1);qR3=A.lift(qR,1)
    rp=[A.dot(R3,n)for n in N];pp=[A.dot(P3,n)for n in N]
    for i,(g,v,q)in enumerate(zip(rp,pp,norms)):
        if i in S:assert g==v==z[3]
        else:
            assert A.sign(g)<0 and A.sign(v)<0
            assert A.sign(A.add(A.mul(g,g),A.scale(A.mul(qR3,q),-1)))>0
    assert A.sign(A.dot(R,T,1),1)<0
    E=tuple(A.c(x,1)for x in polar['strict_escape']['rational_point'])
    cut=A.dot(R,E,1);assert A.sign(cut,1)>0
    # Regular right hexagon on the entire S75 plane.
    p75,G75,J75=projection(S75)
    C=p75(A.lin([M[0],M[65],M[79]],[A.c(2,1),A.c(1,1),A.c(2,1)],1))
    d=A.scale(A.dot(C,C,1),-1,1);assert A.sign(d,1)>0 and A.sign(A.dot(C,T,1),1)<0
    v=p75(M[0]);vc=A.dot(v,C,1);U=A.lin([v,C],[d,vc],1);s=A.dot(U,U,1);assert A.sign(s,1)>0
    basis_anchor=None
    for anchor in[11,28]:
        v=p75(M[anchor]);vc=A.dot(v,C,1);vu=A.dot(v,U,1)
        V=A.lin([v,C,U],[A.mul(d,s,1),A.mul(s,vc,1),A.scale(A.mul(d,vu,1),-1,1)],1)
        t=A.dot(V,V,1)
        if A.sign(t,1)>0:basis_anchor=anchor;break
    assert basis_anchor is not None
    basis=[C,U,V];squares=[d,s,t]
    gram=[[A.dot(v,w,1)for w in basis]for v in basis]
    assert gram==[[A.scale(d,-1,1),z[1],z[1]],[z[1],s,z[1]],[z[1],z[1],t]]
    assert all(A.dot(v,M[i],1)==z[1]for v in basis for i in S75)
    root2=(F(0),F(1))+(F(0),)*14;root3=(F(0),)*8+(F(1),)+(F(0),)*7
    root6=B.mul(root2,root3);one=B.sc(1);zero=B.Z
    coeff=[(root2,B.sc(b),B.scale(root3,c))for b,c in[(2,0),(1,1),(-1,1),(-2,0),(-1,-1),(1,-1)]]
    def frame_dot(a,b):return B.add(B.scale(B.mul(a[0],b[0]),-1),B.add(B.mul(a[1],b[1]),B.mul(a[2],b[2])))
    hexgram=[[frame_dot(a,b)for b in coeff]for a in coeff]
    expected=[[2 if i==j else 0 if (i-j)%6 in[1,5]else -6 if (i-j)%6==3 else -4 for j in range(6)]for i in range(6)]
    assert hexgram==[[B.sc(x)for x in row]for row in expected]
    N101=N+[R3];norm101=norms+[qR3]
    basis3=[A.vlift(v,1)for v in basis]
    bp=[[A.dot(v,n)for v in basis3]for n in N101]
    hexrows=[]
    for k,co in enumerate(coeff):
        row=[]
        for i,(pairs,q)in enumerate(zip(bp,norm101)):
            if i in S75:assert pairs==[z[3]]*3;row.append(dict(wall=i,relation='right',bounds=[F(0),F(0)]))
            else:row.append(dict(wall=i,relation='dotted',bounds=A.normalized_dotted(co,pairs,squares,q)))
        hexrows.append(row)
    old_E_hex=[]
    for co in coeff:
        sg,b=A.normalized_sign(co,[A.dot(E,v,1)for v in basis],squares,1)
        assert sg<0;old_E_hex.append(b)
    # An inward perturbation of the hexagon center is a common ambient center.
    eps=F(1,2**30);H=A.lin([C,T],[A.c(1,1),A.c(eps,1)],1);qH=A.dot(H,H,1)
    assert A.sign(qH,1)<0 and A.sign(A.dot(H,T,1),1)<0
    H3=A.vlift(H,1);hpair=[A.dot(H3,n)for n in N101];assert all(A.sign(g)<0 for g in hpair)
    hb=[A.dot(H,v,1)for v in basis];hexh=[]
    for co in coeff:
        sg,bound=A.normalized_sign(co,hb,squares,1);assert sg<0;hexh.append(bound)
    # Six actual finite vertices on S75, in the normalized orthogonal frame.
    pointcoeff=[(one,B.scale(root2,b),B.scale(root6,c))for b,c in[(F(1,2),F(1,6)),(0,F(1,3)),(F(-1,2),F(1,6)),(F(-1,2),F(-1,6)),(0,F(-1,3)),(F(1,2),F(-1,6))]]
    vertices=[]
    for k,co in enumerate(pointcoeff):
        assert frame_dot(co,co)==B.sc(F(-1,3))
        oldrows=[]
        for i,pairs in enumerate(bp):
            sg,bound=A.normalized_sign(co,pairs,squares)
            assert sg==0 if i in S75 else sg<0
            oldrows.append(bound)
        own=[frame_dot(co,x)for x in coeff];active=[i for i,x in enumerate(own)if x==zero]
        assert active==sorted([k,(k+1)%6])and all(B.sign(x)<0 for i,x in enumerate(own)if i not in active)
        vertices.append(dict(frame_coefficients=co,norm=F(-1,3),incident=S75+[101+i for i in active],old101_pair_bounds=oldrows,hex_pairings=own))
    phex=[]
    for co in coeff:
        sg,bound=A.normalized_sign(co,[A.dot(P,v,1)for v in basis],squares,1);assert sg<0;phex.append(bound)
    # All three old proper cusps remain strict at all seven additions.
    cuspdata=[]
    for raw in[old['completed_cusp']['point']]+[r['point']for r in old['preserved_old_bounded_cusps']]:
        nu1=tuple(A.lift(parse(x,0),0,1)for x in raw);g=A.dot(R,nu1,1);assert A.sign(g,1)<0
        pairs=[A.dot(nu1,v,1)for v in basis];bounds=[]
        for co in coeff:
            sg,b=A.normalized_sign(co,pairs,squares,1);assert sg<0;bounds.append(b)
        cuspdata.append(dict(point=nu1,R_pairing=g,hex_pair_bounds=bounds))
    # The new boundary roles do not close the other retained H2 face.
    S61=sorted(set(S75)-{75}|{61});_,G61,J61=projection(S61)
    assert all(q>=0 for row in J61 for q in row)
    assert J61[S61.index(61)][S61.index(61)]==F(2,3)
    R61=A.dot(R,M[61],1)
    R61_bound=A.add(qR,A.scale(A.mul(R61,R61,1),F(-2,3),1),1)
    assert A.sign(R61_bound,1)<0
    # All old support pairings are <=0. Entrywise nonnegative J61 makes
    # qproj<=q-(2/3)*g61^2<0 for all seven additions. The common strict
    # center, projected to S61, determines the future orientation.
    s61_limitation=dict(support=S61,inverse_gram=J61,R_projected_norm_upper_bound=R61_bound,hex_norm_upper_bound=F(-2,3),proof='Each hexagon normal is right toS75 and strictly dotted to61. All support pairings ofR are nonpositive. Entrywise nonnegative inverseGram and its61 diagonal2/3 give qproj<=q-(2/3)*g61²<-q/3. The projected common center fixes the future orientation. All seven additions are therefore redundant throughout the retainedS61 plane.')
    # Literal proposal from an independent numerical search. The scout is
    # deliberately not a dependency: every asserted property is rechecked here.
    escape_literal=['278120137/1000000','140149051/250000','849447027/1000000','742771047/1000000','637104417/1000000','665619/1250','214446999/500000','326293597/1000000','224700409/1000000','62062553/500000','4913027/200000','-73993233/1000000','-171529589/1000000','-218361597/1000000','-66046197/250000','-154503589/500000','-352768829/1000000','-79101397/200000','79305521/200000','-10590817/200000','19466647/500000']
    Enew=tuple(A.c(q,1)for q in escape_literal)
    qe=A.dot(Enew,Enew,1)
    assert qe==A.c(F(38017800175279,500000000000),1)
    future=A.dot(Enew,H,1);assert A.sign(future,1)<0
    E3=A.vlift(Enew,1);epairs=[A.dot(E3,n)for n in N101]
    assert all(A.sign(g)<0 for g in epairs)
    eframe=[A.dot(Enew,v,1)for v in basis];ehex=[]
    for co in coeff:
        sg,b=A.normalized_sign(co,eframe,squares,1);assert sg<0;ehex.append(b)
    newescape=dict(rational_point=escape_literal,norm=F(38017800175279,500000000000),future_center_pairing=future,old101_pairings=epairs,hex_pair_bounds=ehex,meaning='Strict future-oriented spacelike escape for all107 walls, hence infinite hyperbolic volume. The numerical proposal is not a certificate dependency.')
    counts=Counter(right=19+6*18+6,dotted=81+6*83+9);assert sum(counts.values())==721
    finite_neighbors=[S]+[S75+[101+(k-1)%6,101+(k+1)%6]for k in range(6)]
    paths=FILES+[f'scripts/{TAG}.py','scripts/fresh_h20_coupled_roles_20261014_field.py','scripts/fresh_h20_e7_causal_residuals_20261014_bridge_review.py','results/fresh_h20_polar_pair_20261014_census.json','results/acute_extension_preservation_20261014.md']
    return enc(dict(status='PASS_107_WALL_COUPLED_CUSP_EDGE_AND_COMPACT_S75_HEXAGON',sources={p:sha(p)for p in paths},form=A.Q,field=edge['field'],cap_R=dict(wall=100,support=S,support_gram=G,support_inverse=J,support_gram_rescaling='Root92 is multiplied bysqrt(2); all other support roots are unchanged. This rational Gram represents the same subspace.',H0=H0,nu=nu,W=W,h=h,Z=Z,normal=R,norm=qR,all100_pairings=rp,old_escape_pairing=cut,finite_vertex=P,finite_vertex_norm=A.scale(qR,-1,1),vertex_old100_pairings=pp,vertex_hex_pair_bounds=phex,formula='Z=-2hW+W²nu; R=Z-nu; P=Z+nu; R²=4h²=-P²'),hex_basis=dict(support=S75,support_gram=G75,support_inverse=J75,center_weights=[[0,2],[65,1],[79,2]],C=C,U=U,V=V,squared_lengths=squares,third_basis_anchor=basis_anchor,gram=gram,normalization='Frame C/sqrt(d), U/sqrt(s), V/sqrt(t); all three square roots positive'),hexagon=dict(walls=list(range(101,107)),frame_coefficients=coeff,coefficient_encoding='Each coefficient uses the16 squarefree basis; divide the corresponding C/U/V physical vector bysqrt(d/s/t). No expanded new compositum is assumed.',norms=[2]*6,gram=expected,new_to_old101_pair_checks=hexrows,old_escape_pair_bounds=old_E_hex,finite_vertices=vertices,compactness='The six normals are the regular right-angled hexagon at spatial directions k*pi/3 and inradius asinh(1). All six adjacent intersections are actual future timelike vertices; the Euclidean spatial halfplanes form a bounded hexagon inside the Klein unit disk.'),strict_center=dict(point=H,norm=qH,epsilon=eps,old101_pairings=hpair,hex_pair_bounds=hexh),new_pair_counts=dict(counts),complete_new_finite_neighbors=finite_neighbors,complete_new_vertex_count=7,known_total_vertices_at_least=141,local_census_proof='R has precisely19 finite neighbors, giving its unique vertex. Each hexagon wall has precisely18 support neighbors and its two adjacent hexagon walls. The latter pair is dotted, so only the two adjacent20-supports can be finite vertices. These give exactly6 distinct vertices, all strict atR. All134old vertices are strictly retained by acute extension.',finite_to_ideal_edge=dict(support=S,finite_wall=100,ideal_cusp='bound completed_cusp',proof='P=Z+nu andnu lie on the samePD19support; all omitted oldwalls are strict ontheopen segment, andnu has its retained boundedrank19link.'),preserved_bounded_cusps=cuspdata,S61_limitation=s61_limitation,strict_escape=newescape,maximum_interval_bits=A.maxbits,scope='Seven coherent additions with an actual cusp-edge endpoint and a compact right-angled S75 face. Exactly7 new finite vertices, at least141 overall; no complete old-vertex census. The explicit renewed strict spacelike escape proves the107-wall chamber still has infinite volume.'))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['new_pair_counts'],'bits',d['maximum_interval_bits'],flush=True)
