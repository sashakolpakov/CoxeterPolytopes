#!/usr/bin/env python3
"""Producer-free exact linear-metric cap, affine census, and escape replay."""
import argparse
from itertools import combinations
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts import fresh_h9_remote_next_20261013_bridge_verify as B
from scripts import fresh_h9_remote_next_20261013_bridge_frontier as F
from scripts.fresh_h9_remote_next_20261013_projection_obstruction import inverse
from scripts.fresh_h9_remote_caps_20261013_frontier import psd
V=B.V;K=V.K
PREFIX='data/enumerations/fresh_h9_crossplane_repair_20261013'
CONSTRAINTS=[['seed1',2],['seed2',2],['seed3',5],['seed4',3],['seed5',3],['seed6',3],['seed8',5],['wall17',2],['nextH4',5]]

def reconstruct(data):
    assert data['source']==B.PREFIX+'_caps.json'and data['source_sha256']==V.sha(data['source'])
    previous=V.read(data['source']);B.check(previous)
    G,I,h,g,old,names,caps,center,_=B.reconstruct(previous);old+=caps;names+=['bridge0','bridge1']
    assert data['constraints']==CONSTRAINTS and data['branch']==-1
    A=[old[names.index(n)]for n,m in CONSTRAINTS]
    M=[[V.rat(2)if i==j else V.pair(a,g,b)for j,b in enumerate(A)]for i,a in enumerate(A)]
    Mi,pivots=inverse(M);p=[V.rat(0)if m==2 else V.rat(-1)if m==3 else V.neg(V.PHI)for n,m in CONSTRAINTS]
    lam=V.mv(Mi,p);u=[V.total(V.mul(t,a[i])for t,a in zip(lam,A,strict=True))for i in range(10)]
    hc=[V.pair(center,g,a)for a in A];lc=V.mv(Mi,hc)
    v=[V.sub(center[i],V.total(V.mul(t,a[i])for t,a in zip(lc,A,strict=True)))for i in range(10)]
    qu=V.dot(p,lam);qv=V.sub(V.pair(center,g,center),V.dot(hc,lc));assert qv[0]>0
    rem=V.sub(V.rat(2),qu);assert rem[0]>0
    t=V.neg(V.sqrt(V.div(rem,qv)));r=[V.add(x,V.mul(t,y))for x,y in zip(u,v,strict=True)]
    # M inverse identities give (u,A)=p and v perpendicular to A, hence u.v=0.
    # Thus r²=qu+t²qv=2 exactly and all nine prescribed pairings hold exactly.
    params={'inverse_nonzero_pivots':pivots,'u_squared':V.enc(qu),'v_squared':V.enc(qv),'norm_remainder':V.enc(rem),'coefficient':V.enc(t)}
    return G,I,h,g,old,names,r,center,params

def census(data,roots,names,r,center,g):
    oldnames,pairs=F.pair_table();assert oldnames==names
    prescribed=dict(CONSTRAINTS);new='crossplane0'
    for n in names:pairs['|'.join(sorted((n,new))) ]={'kind':'finite','label':prescribed[n]}if n in prescribed else{'kind':'dotted'}
    allnames=names+[new];allroots=roots+[r]
    def pair(a,b):return pairs['|'.join(sorted((a,b)))]
    def gram(ss):
        vals={2:K.Z,3:K.neg(K.O),5:K.neg(K.P)}
        return[[K.scale(K.O,2)if a==b else vals[pair(a,b)['label']]for b in ss]for a in ss]
    neighbors=[n for n,m in CONSTRAINTS]+[new];positive=[];indefinite=0;dotted=0
    for mask in range(1<<len(neighbors)):
        ss=[n for i,n in enumerate(neighbors)if mask>>i&1]
        if new not in ss:continue
        if any(pair(a,b)['kind']=='dotted'for a,b in combinations(ss,2)):dotted+=1;continue
        kind=psd(gram(ss));assert kind!='positive_semidefinite',('affine subset',ss)
        if kind=='positive_definite':positive.append(ss)
        else:indefinite+=1
    maximal=[ss for ss in positive if not any(set(ss)<set(other)for other in positive)]
    vertices=[]
    for ss in positive:
        if len(ss)!=9:continue
        Mi=K.inverse(gram(ss));inc=[allroots[allnames.index(n)]for n in ss]
        lam=V.mv([[V.golden(x)for x in row]for row in Mi],[V.pair(center,g,a)for a in inc])
        point=[V.sub(center[i],V.total(V.mul(x,a[i])for x,a in zip(lam,inc,strict=True)))for i in range(10)]
        assert V.pair(point,g,point)[1]<0
        records=[]
        for a,n in zip(allroots,allnames,strict=True):
            if n in ss:continue
            p=V.pair(a,g,point);assert p[1]<0;records.append({'wall':n,'pairing':V.enc(p)})
        vertices.append({'incident':ss,'point_enclosures':list(map(V.enc,point)),'nonincident_checks':records})
    inherited=V.read(B.PREFIX+'_frontier.json');assert inherited['inherited_affine_circuits']==inherited['new_affine_circuits']==0
    return{'new_psd_subsets_checked':512,'positive_definite_subsets':len(positive),'dotted_subsets':dotted,'indefinite_subsets':indefinite,
           'new_affine_circuits':0,'inherited_affine_circuits':0,'maximal_positive_definite_blocks_containing_new_facet':maximal,
           'new_finite_vertices':vertices,'all_pair_count':len(pairs)}

def check(data):
    G,I,h,g,roots,names,r,center,params=reconstruct(data);hi=list(map(V.golden,h));cp=V.dot(r,hi);assert cp[0]>0
    nr=V.pair(r,g,r);assert nr[0]<=2*V.SCALE<=nr[1]
    prescribed=dict(CONSTRAINTS);pairs=[]
    for old,n in zip(roots,names,strict=True):
        if n in prescribed:continue
        p=V.pair(r,g,old);assert p[1]<-2*V.SCALE,(n,p);pairs.append({'old':n,'pairing':V.enc(p)})
    assert len(pairs)==79
    assert data['target_point']==V.read(data['source'])['escape']
    target=tuple(map(K.parse,data['target_point']));assert K.sign(K.pair(target,I,target))>0 and K.sign(K.pair(target,I,h))<0
    cut=V.dot(r,list(map(V.golden,target)));assert cut[1]<0
    out={'status':'PASS_EXACT_CROSSPLANE_METRIC_REPAIR','facets':89,'norm':2,'parameters':params,
         'normal_coordinate_enclosures':list(map(V.enc,r)),'center_pairing_dual':V.enc(cp),'target_pairing_dual':V.enc(cut),
         'new_pair_counts':{'right':3,'pi3':3,'pi5':3,'dotted':79},'dotted_checks':pairs,
         'affine_and_vertex_census':census(data,roots,names,r,center,g),
         'scope':'One exact metric-coupled cap preserves all88 old walls and cuts their saved escape. Its prescribed Gram relations use only labels2,3,5. The complete new-subset census finds no positive affine circuit. Absence of affine circuits does not imply compactness or finite volume.'}
    if'escape'in data:
        y=tuple(map(K.parse,data['escape']));q=K.pair(y,I,y);future=K.pair(y,I,h);assert K.sign(q)>0 and K.sign(future)<0
        values=[V.dot(a,list(map(V.golden,y)))for a in roots+[r]];assert min(p[0]for p in values)>0
        out['remaining_escape']={'norm':K.enc(q),'future_center_pairing':K.enc(future),'strict_pairings':list(map(V.enc,values))}
        out['scope']+=' The saved strict future spacelike escape certifies infinite volume for this89-facet chamber.'
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=V.read(PREFIX+'_caps.json');out=check(d);out['source_sha256']=V.sha(PREFIX+'_caps.json')
    if a.verify:assert V.read(PREFIX+'_verify.json')==out
    else:(V.ROOT/(PREFIX+'_verify.json')).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'vertices',len(out['affine_and_vertex_census']['new_finite_vertices']),'EXACT_ESCAPE'if'remaining_escape'in out else'VOLUME_UNTESTED')
