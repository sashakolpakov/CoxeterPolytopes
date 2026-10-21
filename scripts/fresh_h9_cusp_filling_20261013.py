#!/usr/bin/env python3
"""Scout Gram-preserving cusp fillings after deleting walls of an H9 chamber.

All finite-angle mutations are declared. Schur diagonal quadratics propose
the changed old parallel pairings; cross Schur equations select subsets.
Numerical output is scouting only until a separate exact replay.
"""
import argparse,json,time
from collections import Counter
from itertools import combinations,product,chain
from pathlib import Path
import networkx as nx
import numpy as np
from scipy.spatial import HalfspaceIntersection,QhullError

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'data/enumerations/alexeev_h9_doubling_cubes.json'
OUT=ROOT/'data/enumerations/fresh_h9_cusp_filling_20261013.json'

def run(max_bases=212,broad=False):
    row=next(r for r in json.loads(SOURCE.read_text())['certified_types']if r['id']=='T1_D16_W4_7_12_14')
    G=-np.array(row['facet_gram'],dtype=float)/2
    parallel=[(i,j)for i in range(18)for j in range(i)if G[i,j]==-2]
    bases=[]
    for B in combinations(range(18),10):
        if any(i in B and j in B for i,j in parallel):continue
        vals=np.linalg.eigvalsh(G[np.ix_(B,B)])
        if vals[0]<-1e-8 and vals[1]>1e-8:bases.append(B)
    original_base_count=len(bases);orbit_sizes=[]
    if broad:
        import pynauty
        adjacency={i:[]for i in range(18)};colors=[set(),set()]
        for i,j in combinations(range(18),2):
            if G[i,j]==0:continue
            k=len(adjacency);adjacency[k]=[i,j];adjacency[i].append(k);adjacency[j].append(k);colors[int(G[i,j]==-2)].add(k)
        orbits={}
        for B in bases:
            graph=pynauty.Graph(len(adjacency),adjacency_dict=adjacency,vertex_coloring=[set(B),set(range(18))-set(B),*colors]);key=pynauty.certificate(graph)
            orbits.setdefault(key,[]).append(B)
        bases=[v[0]for v in orbits.values()];orbit_sizes=[len(v)for v in orbits.values()]
    weights={2:0.,3:1.,4:2**.5,5:(1+5**.5)/2,6:3**.5,8:2*np.cos(np.pi/8),10:2*np.cos(np.pi/10)}
    finite=[(m,2*np.cos(np.pi/m))for m in range(2,61)]
    stats=Counter();proposals=[];seen=set();start=time.monotonic()
    for bi,B in enumerate(bases[:max_bases]):
        ext=[i for i in range(18)if i not in B];partners=[next(k for k,j in enumerate(B)if G[i,j]==-2)for i in ext]
        original=G[np.ix_(B,B)];oldedges=[(i,j)for i in range(10)for j in range(i)if original[i,j]==-1]
        single=([(i,j,m)]for i,j in oldedges for m in(2,4,5,6,8,10))
        multiple=([(i,j,m)for(i,j),m in zip(edges,labels)]for k in range(2,5 if broad else 3)for edges in combinations(oldedges,k)for labels in product((2,4,5,6),repeat=k))
        # Add one new nonright angle as a distinct surgery on the zero pattern.
        added=([(i,j,m)]for i in range(10)for j in range(i)if original[i,j]==0 for m in(3,4,5,6))
        mutations=chain(single,multiple,added)
        for spec in mutations:
            stats['mutations']+=1;M=original.copy()
            for i,j,m in spec:M[i,j]=M[j,i]=-weights[m]
            vals=np.linalg.eigvalsh(M)
            if not(vals[0]<-1e-7 and vals[1]>1e-7):stats['wrong_basis_signature']+=1;continue
            I=np.linalg.inv(M);columns=G[np.ix_(B,ext)].copy();options=[]
            for ei,(original_id,p)in enumerate(zip(ext,partners)):
                b=columns[:,ei].copy();b[p]=0
                aa=I[p,p];bb=-2*(I@b)[p];cc=b@I@b-2
                if abs(aa)<1e-8:
                    if abs(bb)>1e-8:roots=[-cc/bb]
                    elif abs(cc)>1e-8:roots=[]
                    else:roots=[v for _,v in finite[:11]]+[2.1,2.5,3.,(1+5**.5)] ;stats['sampled_free_parameters']+=1
                else:
                    discr=bb*bb-4*aa*cc
                    roots=[]if discr< -1e-8 else[(-bb+ss*np.sqrt(max(0,discr)))/(2*aa)for ss in(-1,1)]
                for x in roots:
                    if x< -1e-7 or abs(x-2)<1e-7:continue
                    label='dotted'if x>2 else next((m for m,v in finite if abs(v-x)<1e-7),None)
                    if label is None:continue
                    c=b.copy();c[p]=-x
                    if abs(c@I@c-2)>1e-6:continue
                    options.append((ei,original_id,p,float(x),label,c))
            if len({r[0]for r in options})<4:stats['too_few_nonparallel_extensions']+=1;continue
            graph=nx.Graph();graph.add_nodes_from(range(len(options)))
            for i,j in combinations(range(len(options)),2):
                a,b=options[i],options[j]
                if a[0]!=b[0]and abs(a[5]@I@b[5]-G[a[1],b[1]])<1e-7:graph.add_edge(i,j)
            cliques=[sorted(c)for c in nx.find_cliques(graph)if len(c)>=4]
            if not cliques:stats['no_four_extension_clique']+=1;continue
            for clique in cliques:
                chosen=[options[i]for i in clique];ids=list(B)+[r[1]for r in chosen];H=G[np.ix_(ids,ids)].copy();H[:10,:10]=M
                for k,r in enumerate(chosen):H[:10,10+k]=r[5];H[10+k,:10]=r[5]
                key=tuple(np.round(H[np.triu_indices(len(H))],9))
                if key in seen:continue
                seen.add(key);stats['rank_ten_gram_proposals']+=1
                eig,V=np.linalg.eigh(H);v=V[:,0];v*=np.sign(sum(v))
                if np.min(v)<1e-7:stats['reducible_proposals']+=1;continue
                center=-eig[0]*v[:10];center/=sum(center);caps=np.array([I@r[5]for r in chosen])
                halves=np.vstack([-np.eye(10),-caps]);half=np.column_stack((halves[:,:9]-halves[:,9,None],halves[:,9]))
                try:
                    verts=HalfspaceIntersection(half,center[:9]).intersections;verts=np.column_stack((verts,1-np.sum(verts,axis=1)));norms=np.einsum('bi,ij,bj->b',verts,I,verts)
                    gate={'status':'NUMERICAL_COMPACT_CANDIDATE'if max(norms)<-1e-7 else'NUMERICAL_ESCAPE_OR_IDEAL','vertices':len(verts),'maximum_norm':float(max(norms)),'witness':verts[int(np.argmax(norms))].tolist()}
                except QhullError:gate={'status':'NUMERICAL_VOLUME_UNRESOLVED'}
                record={'basis':list(B),'retained_original_indices':ids,'mutated_finite_pairs':[[B[i],B[j],m]for i,j,m in spec],'parallel_replacements':[[r[1],B[r[2]],r[3],r[4]]for r in chosen],'gram':H.tolist(),'volume_gate':gate};proposals.append(record)
                print({'basis':bi,'facets':len(H),'gate':gate,'stats':dict(stats)},flush=True)
        if broad or bi%20==0:print({'bases_done':bi+1,'available_bases':len(bases),'seconds':round(time.monotonic()-start,2),'stats':dict(stats)},flush=True)
    result={'schema':'H9-cusp-filling-Gram-scout-v1','source':str(SOURCE.relative_to(ROOT)),'source_id':row['id'],'constant_Lorentz_bases':original_base_count,'bases_examined':min(len(bases),max_bases),'basis_orbit_sizes':orbit_sizes,'broad':broad,'counts':dict(stats),'proposals':proposals,'scope':'Numerical scout only; one/two (broad: up to four) old label-three mutations or one added nonright edge inside declared bases, then up to four facet deletions. Free diagonal parameters sampled. No exhaustive impossibility claim.'}
    output=OUT.with_name(OUT.stem+'_broad.json')if broad else OUT
    output.write_text(json.dumps(result,indent=2)+'\n');print({'done':True,'proposals':len(proposals),'counts':dict(stats)},flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--bases',type=int,default=212);p.add_argument('--broad',action='store_true');a=p.parse_args();run(a.bases,a.broad)
