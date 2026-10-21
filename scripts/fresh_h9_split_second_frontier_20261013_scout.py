#!/usr/bin/env python3
"""Finite-angle proposals on the outgoing edges of a newly reached vertex."""
import argparse
import itertools
import json
import time
import mpmath as mp
import numpy as np
from scripts.fresh_h9_split_second_frontier_20261013_base import *
from scripts.fresh_h9_split_edge_frontier_20261013 import frontier

def run(vertex='walk2720'):
    start=time.monotonic();mp.mp.dps=80
    G,I,h,old,scout,case,family=current();names=list(old);N=np.array(scout['source_snapshot']['N']);NB=N[list(BASIS)]
    W=np.array([[float(value(r['field'],x))for x in r['c']]for r in old.values()])@NB
    eta=np.diag([-1.]+[1.]*9);Gram=W@eta@W.T;phi=(1+5**.5)/2
    initial=list(scout['source_snapshot']['finite_vertices'])
    for proof in case['actual_finite_vertices']:
        F=old[proof['point_field_from_cap']]['field'];y=np.array([float(value(F,F.parse(x)))for x in proof['y']])
        In=np.array([[float(value(K,x))for x in r]for r in I]);q=-(In@y)@NB;q/=q[0]
        initial.append({'point':q[1:].tolist(),'active':[names.index(n)for n in proof['zero_walls']]})
    front=frontier(W,names,initial)
    open_edges=[e for e in front['edges']if e['classification']=='OPEN_TO_LIGHT_CONE_WITH_CURRENT_WALLS']
    endpoints=np.array([[1.]+r['light_cone_point']for r in open_edges])
    target=np.array([float(value(K,parse(x)))for x in case['geometry']['escapes'][0]['certificate']['y']])
    targetphysical=-(In@target)@NB
    snapshot={'walls':names,'physical_normals':W.tolist(),'N':N.tolist(),'finite_vertices':initial,
              'frontier':front,'strict_target':case['geometry']['escapes'][0]['certificate'],
              'source_sha256':{p:sha(p)for p in(INITIAL+'_family.json',INITIAL+'_family_geometry.json')}}
    path=ROOT/(PREFIX+'_source.json')
    if path.exists():assert json.loads(path.read_text())==snapshot
    else:path.write_text(json.dumps(snapshot,indent=2)+'\n')
    vi=next(i for i,v in enumerate(front['vertices'])if vertex in v['active_names'])
    selected_edges=[e for e in open_edges if e['start_vertex']==vi]
    profiles=[];labels=[]
    for count in range(4):
        for inds in itertools.combinations(range(8),count):
            for ms in itertools.product((3,5),repeat=count):
                p=np.zeros(8)
                for j,m in zip(inds,ms):p[j]=-1 if m==3 else-phi
                profiles.append(p);labels.append(list(zip(inds,ms)))
    profiles=np.array(profiles);fvals=-2*np.cos(np.pi/np.array([2,3,4,5,6,8,10,12,15,16,17,20,24,30,32,34,40,48,60]))
    def allowed(v):return(v<-2-1e-7)|np.any(abs(v[...,None]-fvals)<1e-7,axis=-1)
    E=old[vertex]['field']
    anchors=[]
    for i,n in enumerate(names):
        common,_=embed(E,old[n]['field'])
        if key(common)==key(E):anchors.append(i)
    out=[];seen=set();tested=0
    for edge in selected_edges:
        block=edge['edge_walls'];B=Gram[np.ix_(block,block)];ib=np.linalg.inv(B)
        qs=2-np.einsum('ij,jk,ik->i',profiles,ib,profiles);valid=np.flatnonzero(qs>1e-7)
        assert min(np.linalg.eigvalsh(B))>1e-7
        for anchor in anchors:
            if anchor in block:continue
            A=W[list(block)+[anchor]]@eta;_,ss,vh=np.linalg.svd(A,full_matrices=True)
            if min(ss)<1e-8:continue
            direction=vh[-1];inv=np.linalg.pinv(A)
            targets=np.column_stack([np.repeat(profiles[valid],3,axis=0),np.tile([0.,-1.,-phi],len(valid))])
            line0=targets@inv.T
            qa=direction@eta@direction;qb=2*line0@eta@direction;qc=np.einsum('ij,jk,ik->i',line0,eta,line0)-2
            disc=qb*qb-4*qa*qc
            for branch in(-1,1):
                inds=np.flatnonzero(disc>=-1e-9)
                if abs(qa)<1e-10:
                    if branch==1:continue
                    inds=inds[abs(qb[inds])>1e-8];t=-qc[inds]/qb[inds]
                else:t=(-qb[inds]+branch*np.sqrt(np.maximum(0,disc[inds])))/(2*qa)
                roots=line0[inds]+t[:,None]*direction;tested+=len(roots)
                vals=roots@eta@W.T;good=np.all(allowed(vals),axis=1)&(roots[:,0]>1e-7)
                for local in np.flatnonzero(good):
                    row=roots[local];v=vals[local];kk=tuple(np.round(v[:16],7))
                    if kk in seen:continue
                    seen.add(kk)
                    # A cap must actually cut at least one currently open ideal endpoint.
                    cuts=(endpoints@eta@row)>1e-7
                    if not cuts.any():continue
                    affine=False
                    for begin in(0,8):
                        sums={0.}
                        for z in v[begin:begin+8]:
                            if z>=-2-1e-7:sums|={s+z*z for s in list(sums)}
                        if any(abs(s-4)<1e-7 for s in sums):affine=True;break
                    if affine:continue
                    pi=int(valid[inds[local]//3]);angle=(2,3,5)[int(inds[local]%3)]
                    out.append({'edge':block,'edge_names':[names[i]for i in block],
                                'released_wall':edge['released_wall_name'],'anchor':anchor,'anchor_name':names[anchor],
                                'finite_profile':labels[pi],'anchor_label':angle,'quadratic_branch':branch,
                                'normal_physical':row.tolist(),'basis_covector':(NB@eta@row).tolist(),
                                'seed_covector':v[:16].tolist(),'open_edge_cut_indices':np.flatnonzero(cuts).tolist(),
                                'cuts_strict_target':float(targetphysical@eta@row)>1e-7,
                                'strict_target_pairing':float(targetphysical@eta@row)})
        print({'released':edge['released_wall_name'],'candidates':len(out),'tested':tested},flush=True)
    out.sort(key=lambda r:(-len(r['open_edge_cut_indices']),-r['strict_target_pairing']))
    target_path='data/enumerations/fresh_h9_split_multiple_edge_caps_20261013_escape_scout.json'
    if(ROOT/target_path).exists():
        target_data=read(target_path);tx=np.r_[-1.,target_data['selected_point']]
        for r in out:r['remote_escape_pairing']=float(np.array(r['normal_physical'])@tx)
    result={'status':'NUMERICAL_SECOND_VERTEX_EDGE_CAPS','source':PREFIX+'_source.json','source_sha256':sha(PREFIX+'_source.json'),
            'selected_vertex':vertex,'edges_considered':len(selected_edges),'anchor_indices':anchors,
            'tested_roots':tested,'candidates':out,'seconds':round(time.monotonic()-start,2),
            'scope':'No-removal proposals on open outgoing edges of one actual finite vertex. Profiles have at most3 attachments with labels3/5; anchors share its exact field. Numerical proposals only, not an exhaustive edge or metric search.'}
    suffix='_scout.json'if vertex=='walk2720'else '_opposite_scout.json'
    (ROOT/(PREFIX+suffix)).write_text(json.dumps(result,indent=2)+'\n')
    print({'status':result['status'],'candidates':len(out),'best_coverage':len(out[0]['open_edge_cut_indices'])if out else 0,
           'strict_target_cutters':sum(r['cuts_strict_target']for r in out),'seconds':result['seconds']},flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--vertex',default='walk2720');a=p.parse_args();run(a.vertex)
