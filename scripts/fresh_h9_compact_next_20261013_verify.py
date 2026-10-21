#!/usr/bin/env python3
"""Exact four-cap insertion with an actual compact 3-5-3 tetrahedral face."""
import argparse,json
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts import fresh_h9_crossplane_repair_20261013_pair_verify as P
V=P.V;K=P.K
PREFIX='data/enumerations/fresh_h9_compact_next_20261013'
S=['crosspair0','seed1','seed2','seed6','wall17'];ANCHOR='nextH4'
NEW=['compact'+str(i)for i in range(4)]

def table(G):
    names,pairs=P.F.pair_table()
    for new,profile in zip(P.NEW,P.finite_profiles(G,names),strict=True):
        for n in names:pairs['|'.join(sorted((new,n)))]={'kind':'finite','label':profile[n]}if n in profile else{'kind':'dotted'}
    pairs['|'.join(P.NEW)]={'kind':'finite','label':2};names+=P.NEW
    return names,pairs

def gram(ss,pairs):
    vals={2:K.Z,3:K.neg(K.O),5:K.neg(K.P)}
    return[[K.scale(K.O,2)if a==b else vals[pairs['|'.join(sorted((a,b)))]['label']]for b in ss]for a in ss]

def reconstruct(d):
    assert d['source']==P.PREFIX+'_caps.json'and d['source_sha256']==V.sha(d['source'])
    previous=V.read(d['source']);prior=P.check(previous)
    G,I,h,g,old,names,extra,center,_=P.reconstruct(previous);old+=extra;names+=P.NEW
    assert d['orthogonal_old_walls']==S and d['anchor']==ANCHOR and d['anchor_label']==3
    assert d['frame_sources']==[0,3]and d['boost']=='13/4'
    expected,pairs=table(G);assert expected==names
    A=[old[names.index(n)]for n in S];M=gram(S,pairs);Mi=K.inverse(M);mii=[[V.golden(x)for x in row]for row in Mi]
    def ps(x):
        lam=V.mv(mii,[V.pair(x,g,r)for r in A])
        return[V.sub(x[i],V.total(V.mul(t,r[i])for t,r in zip(lam,A,strict=True)))for i in range(10)]
    p=[K.neg(K.P),K.Z,K.Z,K.Z,K.Z]
    assert [pairs['|'.join(sorted((ANCHOR,n)))]for n in S]==[{'kind':'finite','label':m}for m in(5,2,2,2,2)]
    qa=K.sub(K.scale(K.O,2),K.pair(p,Mi,p));assert qa==K.scale(K.sub(K.scale(K.O,4),K.scale(K.P,2)),Q(1,3))
    assert K.sign(qa)>0;a=ps(old[names.index(ANCHOR)]);qi=V.golden(qa)
    def project(x):
        z=ps(x);lam=V.div(V.pair(z,g,a),qi)
        return[V.sub(t,V.mul(lam,b))for t,b in zip(z,a,strict=True)]
    time=project(center);qt=V.pair(time,g,time);assert qt[1]<0
    time=[V.div(x,V.sqrt(V.neg(qt)))for x in time]
    target=tuple(map(K.parse,d['target']));assert d['target']==previous['escape']
    w=list(map(V.golden,[K.neg(x)for x in K.mv(I,target)]))
    direction=project(w);lam=V.pair(direction,g,time)
    direction=[V.add(x,V.mul(lam,t))for x,t in zip(direction,time,strict=True)]
    qd=V.pair(direction,g,direction);assert qd[0]>0
    direction=[V.div(x,V.sqrt(qd))for x in direction]
    frame=[time,direction];pivots=[V.enc(qt),V.enc(qd)]
    for j in d['frame_sources']:
        z=project([V.rat(int(i==j))for i in range(10)])
        for bi,b in enumerate(frame):
            lam=V.pair(z,g,b);lam=V.neg(lam)if bi==0 else lam
            z=[V.sub(x,V.mul(lam,y))for x,y in zip(z,b,strict=True)]
        q=V.pair(z,g,z);assert q[0]>0;pivots.append(V.enc(q));frame.append([V.div(x,V.sqrt(q))for x in z])
    nu=K.sub(K.div(K.O,qa),K.scale(K.O,2));assert K.sign(nu)>0
    k=Q(d['boost']);ct=(k+1/k)/2;st=(k-1/k)/2;assert ct*ct-st*st==1
    pt=[V.add(V.mul(V.rat(ct),t),V.mul(V.rat(st),e))for t,e in zip(time,direction,strict=True)]
    bt=[V.add(V.mul(V.rat(st),t),V.mul(V.rat(ct),e))for t,e in zip(time,direction,strict=True)]
    r=[V.add(V.neg(V.div(x,qi)),V.mul(V.sqrt(V.golden(nu)),y))for x,y in zip(a,pt,strict=True)]
    u,v=frame[2:];b=[V.mul(V.sqrt(V.rat(Q(2,3))),x)for x in bt]
    z0=[V.add(x,V.mul(V.sqrt(V.rat(Q(4,3))),y))for x,y in zip(b,u,strict=True)]
    z1=[V.add(V.sub(x,V.div(y,V.sqrt(V.rat(3)))),z)for x,y,z in zip(b,u,v,strict=True)]
    z2=[V.sub(V.sub(x,V.div(y,V.sqrt(V.rat(3)))),z)for x,y,z in zip(b,u,v,strict=True)]
    # Projection and orthonormal Gram-Schmidt identities give frame Gram
    # diag(-1,1,1,1). Then r²=1/qa-nu=2, r.anchor=-1;
    # each zj²=2, all six new cross pairings are exactly zero.
    inherited=prior['affine_and_vertex_census'];assert inherited['new_affine_circuits']==inherited['inherited_affine_circuits']==0
    return G,I,h,g,old,names,[r,z0,z1,z2],center,{'anchor_projection_norm':K.enc(qa),'nu':K.enc(nu),'frame_pivots':pivots,'frame_enclosures':[list(map(V.enc,row))for row in frame],'inherited_affine_circuits':0}

def check(d):
    G,I,h,g,old,names,caps,center,params=reconstruct(d);_,pairs=table(G);hi=list(map(V.golden,h));records=[]
    for i,r in enumerate(caps):
        cp=V.dot(r,hi);assert cp[0]>0
        q=V.pair(r,g,r);assert q[0]<=2*V.SCALE<=q[1]
        dotted=[]
        for a,n in zip(old,names,strict=True):
            if n in S+[ANCHOR]:
                p=V.pair(r,g,a);expected=-V.SCALE if i==0 and n==ANCHOR else 0
                assert p[0]<=expected<=p[1]
                continue
            p=V.pair(r,g,a);assert p[1]<-2*V.SCALE,(i,n,p);dotted.append({'old':n,'pairing':V.enc(p)})
        assert len(dotted)==84
        records.append({'name':NEW[i],'norm':2,'center_pairing_dual':V.enc(cp),'dotted_checks':dotted})
        for n in names:pairs['|'.join(sorted((NEW[i],n))) ]={'kind':'finite','label':3 if i==0 and n==ANCHOR else 2}if n in S+[ANCHOR]else{'kind':'dotted'}
        for j,n in enumerate(NEW[:i]):
            pp=V.pair(r,g,caps[j]);assert pp[0]<=0<=pp[1]
            pairs['|'.join(sorted((NEW[i],n))) ]={'kind':'finite','label':2}
    y=tuple(map(K.parse,d['target']));assert K.sign(K.pair(y,I,y))>0 and K.sign(K.pair(y,I,h))<0
    yi=list(map(V.golden,y));assert min(V.dot(r,yi)[0]for r in old)>0
    cuts=[V.dot(r,yi)for r in caps];assert cuts[0][1]<0
    allroots=old+caps;allnames=names+NEW
    positive=[];indefinite=0;checked=0
    for size in range(1,5):
        for ns in combinations(NEW,size):
            for mask in range(64):
                ss=list(ns)+[n for i,n in enumerate(S+[ANCHOR])if mask>>i&1];checked+=1
                typ=P.T.psd(gram(ss,pairs));assert typ!='positive_semidefinite',ss
                if typ=='positive_definite':positive.append(ss)
                else:indefinite+=1
    assert checked==960
    assert params['inherited_affine_circuits']==0
    vertices=[]
    for ss in positive:
        if len(ss)!=9:continue
        M=gram(ss,pairs);Mi=K.inverse(M);inc=[allroots[allnames.index(n)]for n in ss]
        lam=V.mv([[V.golden(x)for x in row]for row in Mi],[V.pair(center,g,r)for r in inc])
        assert all(p[1]<0 for p in lam)
        point=[V.sub(center[i],V.total(V.mul(t,r[i])for t,r in zip(lam,inc,strict=True)))for i in range(10)]
        q=V.pair(point,g,point);assert q[1]<0 and V.dot(point,hi)[0]>0
        pp=[]
        for r,n in zip(allroots,allnames,strict=True):
            if n in ss:continue
            p=V.pair(point,g,r);assert p[1]<0;pp.append({'wall':n,'pairing':V.enc(p)})
        vertices.append({'incident':ss,'gram':[[K.enc(x)for x in row]for row in M],'norm':V.enc(q),'point_enclosures':list(map(V.enc,point)),'strict_nonincident_pairings':pp})
    assert len(vertices)==4
    common=['seed1','seed2','wall17']+NEW[1:];tetra=['seed6','crosspair0',ANCHOR,NEW[0]];L=gram(tetra,pairs)
    assert all(pairs['|'.join(sorted((a,b)))]=={'kind':'finite','label':2}for a in common for b in tetra)
    assert all(P.T.psd(gram(list(s),pairs))=='positive_definite'for s in combinations(tetra,3))
    Li=K.inverse(L);assert all(K.sign(x)<0 for row in Li for x in row)
    assert {frozenset(v['incident'])for v in vertices}=={frozenset(common+[n for n in tetra if n!=omit])for omit in tetra}
    # Exhaust every possible ninth facet for the30 distinct outgoing edges.
    edge_sets=sorted({tuple(sorted(s))for v in vertices for s in combinations(v['incident'],8)})
    edge_rows=[]
    for ss in edge_sets:
        endings=[];dotted=0;other=0
        for n in allnames:
            if n in ss:continue
            if any(pairs['|'.join(sorted((n,a)))]['kind']=='dotted'for a in ss):dotted+=1;continue
            typ=P.T.psd(gram(list(ss)+[n],pairs));assert typ!='positive_semidefinite'
            if typ=='positive_definite':endings.append(n)
            else:other+=1
        assert len(endings)in(1,2)
        edge_rows.append({'common_facets':list(ss),'finite_end_walls':endings,'status':'finite'if len(endings)==2 else'open_to_light_cone','excluded_by_dotted_pair':dotted,'excluded_by_indefinite_gram':other})
    assert len(edge_rows)==30 and sum(x['status']=='finite'for x in edge_rows)==6
    out={'status':'PASS_EXACT_94_FACET_COMPACT_TETRAHEDRAL_FACE_WITH_VOLUME_GATE','facets':94,'retained_facets':90,
         'parameters':params,'new_normals':records,'normal_coordinate_enclosures':[list(map(V.enc,r))for r in caps],
         'new_pair_counts':{'right':29,'pi3':1,'dotted':336},'total_pair_count':len(pairs),'target_cut_pairings_dual':list(map(V.enc,cuts)),
         'affine_census':{'new_subsets':checked,'positive_definite':len(positive),'indefinite':indefinite,'new_affine_circuits':0,'inherited_affine_circuits':0},
         'actual_finite_vertices':vertices,'compact_H3_face':{'common_facets':common,'bounding_facets':tetra,'diagram':'3-5-3','gram':[[K.enc(x)for x in row]for row in L],
          'inverse':[[K.enc(x)for x in row]for row in Li],'vertices':4,'edges':6,
          'proof':'The four exact future timelike dual rays generate the complete four-wall halfspace cone in the common-face span; its projectivization is a compact tetrahedron. All four vertices satisfy every other inequality strictly, so no old wall truncates it.'},
         'new_vertex_edge_frontier':{'finite_edges':6,'open_rays':24,'ideal_vertices':0,'edges':edge_rows,
          'scope':'All outgoing rank8 spherical edges of the four new vertices, with every possible ninth wall checked. This is not the frontier of every inherited vertex.'},
         'scope':'Four coupled exact caps retain all90 walls, cut their saved gap, and create an actual compact Coxeter H3 face. The complete new-PSD census has no affine circuit. This does not establish finite volume or compactness of the H9 chamber.'}
    assert len(pairs)==94*93//2
    if 'escape'in d:
        y=tuple(map(K.parse,d['escape']));q=K.pair(y,I,y);future=K.pair(y,I,h);assert K.sign(q)>0 and K.sign(future)<0
        pp=[V.dot(r,list(map(V.golden,y)))for r in allroots];assert min(p[0]for p in pp)>0
        out['remaining_escape']={'norm':K.enc(q),'future_center_pairing':K.enc(future),'strict_pairings':list(map(V.enc,pp))}
        out['scope']+=' The saved strict future spacelike escape proves infinite volume.'
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=V.read(PREFIX+'_caps.json');out=check(d);out['source_sha256']=V.sha(PREFIX+'_caps.json')
    if a.verify:assert out==V.read(PREFIX+'_verify.json')
    else:(V.ROOT/(PREFIX+'_verify.json')).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'vertices',len(out['actual_finite_vertices']),'EXACT_ESCAPE'if'remaining_escape'in out else'VOLUME_UNTESTED')
