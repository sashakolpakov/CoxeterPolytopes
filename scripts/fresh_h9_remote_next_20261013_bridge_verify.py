#!/usr/bin/env python3
"""Exact two-cap bridge across an indefinite remote face."""
import argparse
import json
from scripts import fresh_h9_remote_next_20261013_verify as N
V=N.V;K=V.K
PREFIX=N.PREFIX+'_bridge'

def reconstruct(data):
    assert data['source']==N.PREFIX+'_caps.json'and data['source_sha256']==V.sha(data['source'])
    previous=V.read(data['source']);N.check(previous)
    G,I,h,g,old,names,r,_,_,_=N.reconstruct(previous);roots=old+[r];names=names+['nextH4']
    assert data['orthogonal_old_walls']==N.S_NAMES+['wall22','nextH4']and data['finite_anchor']=='wall38'and data['finite_label']==3
    assert data['branch']==-1
    S=[roots[names.index(n)]for n in N.S_NAMES];A=roots[names.index('wall22')];B=r
    old_source=V.read(V.OLD);old_pairs=V.read(old_source['source'])['pair_checks']
    for n in N.S_NAMES:
        assert old_pairs['|'.join(sorted((n,'wall22')))]=={'kind':'finite','label':2}
    p=V.pair(A,g,B);assert p[1]<-2*V.SCALE;det=V.sub(V.rat(4),V.mul(p,p));assert det[1]<0
    def project(x):
        w=list(x)
        for s in S:
            a=V.div(V.pair(w,g,s),V.rat(2));w=[V.sub(t,V.mul(a,z))for t,z in zip(w,s,strict=True)]
        wa=V.pair(w,g,A);wb=V.pair(w,g,B)
        a=V.div(V.sub(V.mul(V.rat(2),wa),V.mul(p,wb)),det)
        b=V.div(V.sub(V.mul(V.rat(2),wb),V.mul(p,wa)),det)
        return[V.sub(t,V.add(V.mul(a,x),V.mul(b,y)))for t,x,y in zip(w,A,B,strict=True)]
    center=[V.neg(x)for x in V.mv([[V.golden(x)for x in row]for row in I],list(map(V.golden,h)))]
    u=project(roots[names.index('wall38')]);d=V.pair(u,g,u);assert d[0]>0
    v=project(center);z=V.div(V.pair(v,g,u),d);v=[V.sub(x,V.mul(z,y))for x,y in zip(v,u,strict=True)]
    e=V.pair(v,g,v);assert e[0]>0
    q=V.sub(V.rat(2),V.div(V.rat(1),d));assert q[0]>0
    t=V.sqrt(V.div(q,e));s=V.sqrt(V.div(q,d));z=V.div(V.rat(1),V.sqrt(V.mul(d,e)))
    r0=[V.sub(V.neg(V.div(x,d)),V.mul(t,y))for x,y in zip(u,v,strict=True)]
    r1=[V.add(V.neg(V.mul(s,x)),V.mul(z,y))for x,y in zip(u,v,strict=True)]
    # Exact projection identities give u.v=0,u²=d,v²=e. Hence both roots
    # have square2, their pairing is0, and r0.wall38=-1. No interval
    # residual is used to decide an equality.
    params={'old_dotted_pairing':V.enc(p),'lorentz_block_determinant':V.enc(det),'d':V.enc(d),'e':V.enc(e),'q':V.enc(q)}
    return G,I,h,g,roots,names,[r0,r1],center,params

def check(data):
    G,I,h,g,roots,names,caps,center,params=reconstruct(data);hi=list(map(V.golden,h));checks=[]
    for i,r in enumerate(caps):
        assert V.dot(r,hi)[0]>0
        for old,n in zip(roots,names,strict=True):
            if n in data['orthogonal_old_walls']or(i==0 and n=='wall38'):continue
            p=V.pair(r,g,old);assert p[1]<-2*V.SCALE,('not dotted',i,n,p)
            checks.append({'new':i,'old':n,'pairing':V.enc(p)})
    target=tuple(map(K.parse,data['target_point']));assert data['target_point']==V.read(data['source'])['escape']
    assert K.sign(K.pair(target,I,target))>0 and K.sign(K.pair(target,I,h))<0
    assert all(V.dot(r,list(map(V.golden,target)))[0]>0 for r in roots)
    cuts=[V.dot(r,list(map(V.golden,target)))for r in caps];assert all(p[1]<0 for p in cuts)
    allroots=roots+caps;allnames=names+['bridge0','bridge1'];vertices=[]
    for endpoint in('wall22','nextH4'):
        incident=N.S_NAMES+[endpoint,'bridge0','bridge1'];inc=[allroots[allnames.index(n)]for n in incident]
        lam=[V.div(V.pair(center,g,r),V.rat(2))for r in inc];assert all(p[1]<0 for p in lam)
        point=[V.sub(center[i],V.total(V.mul(t,r[i])for t,r in zip(lam,inc,strict=True)))for i in range(10)]
        assert V.pair(point,g,point)[1]<0 and V.dot(point,hi)[0]>0
        ps=[]
        for r,n in zip(allroots,allnames,strict=True):
            if n in incident:continue
            p=V.pair(point,g,r);assert p[1]<0;ps.append({'wall':n,'pairing':V.enc(p)})
        vertices.append({'endpoint_wall':endpoint,'incident':incident,'gram':'2I9','point_enclosures':list(map(V.enc,point)),
                         'nonincident_checks':ps})
    out={'status':'PASS_EXACT_INDEFINITE_FACE_BRIDGE','facets':88,'new_caps':2,'norms':[2,2],'new_mutual_label':2,
         'parameters':params,'normal_coordinate_enclosures':[list(map(V.enc,r))for r in caps],
         'new_old_pair_counts':{'right':16,'pi3':1,'dotted':155},'dotted_checks':checks,'target_cut_pairings':list(map(V.enc,cuts)),
         'actual_vertices':vertices,'finite_edge':{'common_walls':N.S_NAMES+['bridge0','bridge1'],'end_walls':['wall22','nextH4']},
         'scope':'Two exact caps in the positive normal plane of an eight-wall Lorentzian span retain all86 old facets, cut the saved escape and create two actual9A1 vertices joined by a finite edge. This is a specific finite-angle choice in a circle family, not an exhaustive cap census.'}
    assert len(checks)==155
    if 'escape'in data:
        y=tuple(map(K.parse,data['escape']));q=K.pair(y,I,y);future=K.pair(y,I,h)
        assert K.sign(q)>0 and K.sign(future)<0
        vals=[V.dot(r,list(map(V.golden,y)))for r in allroots];assert min(p[0]for p in vals)>0
        out['remaining_escape']={'norm':K.enc(q),'future_center_pairing':K.enc(future),'strict_pairings':list(map(V.enc,vals))}
        out['scope']+=' The saved strict future spacelike escape proves that this88-facet chamber still has infinite volume.'
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();data=V.read(PREFIX+'_caps.json');out=check(data)
    out['source_sha256']=V.sha(PREFIX+'_caps.json')
    if a.verify:assert V.read(PREFIX+'_verify.json')==out
    else:(V.ROOT/(PREFIX+'_verify.json')).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['facets'],'facets','EXACT_ESCAPE'if'remaining_escape'in out else'volume untested')
