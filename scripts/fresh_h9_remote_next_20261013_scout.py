#!/usr/bin/env python3
"""Numerical finite-profile scout at the new H4 vertex, aimed at its escape."""
from itertools import combinations,product
import json
import numpy as np
from scipy.linalg import null_space
from scipy.spatial import HalfspaceIntersection
from scripts import fresh_h9_remote_next_20261013_verify as N
V=N.V
OUT=N.PREFIX+'_scout.json'

def run():
    source=N.PREFIX+'_caps.json';data=V.read(source);G,I,h,g,roots,names,cap,_,_,_=N.reconstruct(data)
    roots=roots+[cap];names=names+['nextH4'];C=np.array([[(x[0]+x[1])/(2*V.SCALE)for x in row]for row in roots])
    GG=np.array([[sum(V.golden(x))/(2*V.SCALE)for x in row]for row in G]);Cov=C@GG
    hh=np.array([sum(V.golden(x))/(2*V.SCALE)for x in h]);y=np.array([sum(V.golden(V.K.parse(x)))/(2*V.SCALE)for x in data['escape']])
    In=np.linalg.inv(GG);hc=hh/hh.sum();half=np.c_[C[:,9,None]-C[:,:9],-C[:,9]];half/=np.linalg.norm(half[:,:9],axis=1)[:,None]
    vv=HalfspaceIntersection(half,hc[:9],qhull_options='Qx').intersections;ys=np.c_[vv,1-vv.sum(axis=1)];qq=np.einsum('bi,ij,bj->b',ys,In,ys)
    vertices=[]
    for p in ys[qq<-1e-9]:
        active=tuple(int(i)for i in np.flatnonzero(np.abs(C@p)<2e-7));assert len(active)==9
        vertices.append(active)
    edges=sorted({tuple(j for j in ids if j!=i)for ids in vertices for i in ids})
    edge_data=[]
    for block in edges:
        M=C[list(block)]@GG@C[list(block)].T;assert np.linalg.eigvalsh(M).min()>1e-7
        edge_data.append((block,np.linalg.inv(M)))
    incident=[names.index(n)for n in N.S_NAMES+['wall17','coupled0','nextH4']]
    anchors=['seed'+str(i)for i in range(16)]+['wall'+str(i)for i in range(16,24)]+['coupled0','coupled1']+['remote'+str(i)for i in range(4)]
    alpha={2:0.,3:1.,4:2**.5,5:(1+5**.5)/2,6:3**.5,10:2*np.cos(np.pi/10)}
    angle_values=np.array([2*np.cos(np.pi/m)for m in(2,3,4,5,6,8,10,12,15,16,20,24,30)])
    anchor_values=[('finite',m,-v)for m,v in alpha.items()]+[('dotted',str(v),-v)for v in(2.5,3.,4.)]
    counters={'spherical_profiles':0,'quadratic_branches':0,'future_target_cutters':0,'fully_compatible_target_cutters':0};accepted=[];seen=set()
    for removed in incident:
        block=[j for j in incident if j!=removed];M=C[block]@GG@C[block].T;Mi=np.linalg.inv(M);profiles=[]
        for support_size in range(4):
            for support in combinations(range(8),support_size):
                for labels in product((3,4,5),repeat=support_size):
                    p=np.zeros(8)
                    for j,m in zip(support,labels,strict=True):p[j]=-alpha[m]
                    if 2-p@Mi@p>1e-7:profiles.append((support,labels,p))
        counters['spherical_profiles']+=len(profiles)
        for anchor in anchors:
            a=names.index(anchor)
            if a in block:continue
            A=Cov[block+[a]];ns=null_space(A)
            if ns.shape[1]!=1:continue
            v=ns[:,0];transform=np.linalg.pinv(A);aa=v@GG@v
            if abs(aa)<1e-10:continue
            for support,labels,p in profiles:
                for kind,m,ar in anchor_values:
                    rhs=np.r_[p,ar];u=transform@rhs
                    if np.linalg.norm(A@u-rhs)>1e-7:continue
                    bb=2*u@GG@v;cc=u@GG@u-2;D=bb*bb-4*aa*cc
                    if D<-1e-8:continue
                    for branch in(-1,1):
                        counters['quadratic_branches']+=1;t=(-bb+branch*np.sqrt(max(0,D)))/(2*aa);r=u+t*v
                        if r@hh<=1e-8 or -r@y<=1e-6:continue
                        counters['future_target_cutters']+=1;ps=Cov@r
                        if np.max(ps)>1e-6:continue
                        finite=ps>-2+1e-6
                        if np.any(np.min(np.abs(ps[finite,None]+angle_values[None,:]),axis=1)>1e-6):continue
                        counters['fully_compatible_target_cutters']+=1;key=tuple(np.round(r,7))
                        if key in seen:continue
                        seen.add(key);closed=[]
                        for edge,inv in edge_data:
                            q=ps[list(edge)]
                            if np.max(-q)>2-1e-7:continue
                            if 2-q@inv@q>1e-7:closed.append([names[i]for i in edge])
                        accepted.append({'removed_vertex_wall':names[removed],'edge':[names[i]for i in block],
                            'finite_profile':[[int(j),int(m)]for j,m in zip(support,labels,strict=True)],'anchor':anchor,
                            'anchor_kind':kind,'anchor_value':m,'branch':branch,'coefficients':list(map(float,r)),
                            'target_cut':float(-r@y),'closed_known_edges':closed,'closed_edge_count':len(closed)})
        print({'removed':names[removed],'profiles':len(profiles),'accepted_distinct':len(accepted)},flush=True)
    accepted.sort(key=lambda r:(-r['closed_edge_count'],-r['target_cut']))
    out={'status':'NUMERICAL_TARGETED_H4_FRONTIER_SCOUT','source':source,'source_sha256':V.sha(source),'counters':counters,
         'known_numeric_finite_vertices':len(vertices),'known_numeric_edges':len(edges),'distinct_accepted':len(accepted),'candidates':accepted,
         'scope':'Numerical proposals only. Sparse finite profiles, a stated anchor list and finite/dotted values are tested at the9 outgoing H4-vertex edges. No root, completeness, or volume claim until exact reconstruction.'}
    (V.ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n');print(counters,'distinct',len(accepted),'best',accepted[:1],flush=True)

if __name__=='__main__':run()
