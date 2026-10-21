#!/usr/bin/env python3
"""Exact metric-coupled H4 endpoint cap extending the remote85 chamber."""
import argparse
import json
from fractions import Fraction as Q
from scripts import fresh_h9_remote_caps_20261013_verify as V
from scripts import fresh_h9_remote_caps_20261013_coupled_verify as W
K=V.K
PREFIX='data/enumerations/fresh_h9_remote_next_20261013'
S_NAMES=['seed'+str(i)for i in(1,2,3,5,6,7)]

def reconstruct(data):
    assert data['source']==W.PREFIX+'_caps.json'and data['source_sha256']==V.sha(data['source'])
    previous=V.read(data['source']);previous_certificate=W.check(previous)
    G,I,h,g,old,names,pair,P,poly=W.reconstruct(previous)
    roots=old+pair;names=names+['coupled0','coupled1'];assert len(roots)==85
    assert data['orthogonal_walls']==S_NAMES+['wall17']and data['pi3_wall']=='coupled0'and data['pi10_wall']=='coupled1'
    assert data['negative_unit_seed']=='seed0'and data['norm']==2
    # The previous exact frame gives P²=-1,R²=C²=1 and pairwise orthogonality.
    # All three are orthogonal to the norm-two seed0 and to S_NAMES.
    R=[V.div(V.add(x,y),V.rat(2))for x,y in zip(pair[0],pair[1],strict=True)]
    C=[V.div(V.sub(x,y),V.rat(2))for x,y in zip(pair[0],pair[1],strict=True)]
    U=roots[names.index('seed0')];anchor=roots[names.index('wall17')]
    old_input=V.read(V.OLD);old_pairs=V.read(old_input['source'])['pair_checks']
    assert old_pairs['seed0|wall17']=={'kind':'finite','label':2}
    assert 'seed0'in W.S_NAMES and set(S_NAMES)|{'seed0'}==set(W.S_NAMES)
    L=V.sqrt(V.golden(K.add(K.scale(K.O,2),K.P)))
    b=V.div(V.neg(V.add(V.rat(1),L)),V.rat(2));c=V.div(V.sub(L,V.rat(1)),V.rat(2))
    pp,pr,pc=[V.pair(x,g,anchor)for x in(P,R,C)];assert pp[1]<0
    a=V.div(V.neg(V.add(V.mul(b,pr),V.mul(c,pc))),pp)
    D=V.div(V.sub(V.add(V.rat(2),V.mul(a,a)),V.add(V.mul(b,b),V.mul(c,c))),V.rat(2));assert D[0]>0
    u=V.neg(V.sqrt(D))
    cap=[V.total(V.mul(z,v[i])for z,v in((a,P),(b,R),(c,C),(u,U)))for i in range(10)]
    # Norm two is exact by 2u²=2+a²-b²-c². Prescribed pairings are:
    # cap.(R+C)=b+c=-1, cap.(R-C)=b-c=-L, and cap.wall17=0.
    # L>0 and L²=2+phi give precisely 2cos(pi/10), not its conjugate.
    return G,I,h,g,roots,names,cap,previous_certificate,{'a':V.enc(a),'b':V.enc(b),'c':V.enc(c),'u':V.enc(u),
        'L':V.enc(L),'u_squared':V.enc(D),'anchor_P_pairing':V.enc(pp)},[V.enc(x)for x in cap]

def check(data):
    G,I,h,g,roots,names,cap,previous_certificate,parameters,coordinates=reconstruct(data)
    hi=list(map(V.golden,h));assert V.dot(cap,hi)[0]>0
    records=[]
    for r,n in zip(roots,names,strict=True):
        if n in S_NAMES+['wall17','coupled0','coupled1']:continue
        p=V.pair(cap,g,r);assert p[1]<-2*V.SCALE,('non-dotted unprescribed pair',n,p)
        records.append({'old_wall':n,'pairing':V.enc(p)})
    assert len(records)==76
    target=tuple(map(K.parse,data['target_point']));assert data['target_point']==V.read(data['source'])['escape']
    assert K.sign(K.pair(target,I,target))>0 and K.sign(K.pair(target,I,h))<0
    assert all(V.dot(r,list(map(V.golden,target)))[0]>0 for r in roots)
    cut=V.dot(cap,list(map(V.golden,target)));assert cut[1]<0
    # Actual endpoint: seed3--5--wall17--3--coupled0--3--newcap,
    # together with the other five orthogonal seed normals.
    incident=S_NAMES+['wall17','coupled0','nextH4']
    M=[[K.scale(K.O,2)if i==j else K.neg(K.P)if {i,j}=={2,6}
        else K.neg(K.O)if {i,j}in({6,7},{7,8})else K.Z for j in range(9)]for i in range(9)]
    reduced=[row[:]for row in M]
    for j in range(9):
        assert K.sign(reduced[j][j])>0
        for k in range(j+1,9):
            for l in range(j+1,9):reduced[k][l]=K.sub(reduced[k][l],K.div(K.mul(reduced[k][j],reduced[j][l]),reduced[j][j]))
    inv=K.inverse(M);assert all(K.sign(x)>=0 for row in inv for x in row)
    center=[V.neg(x)for x in V.mv([[V.golden(x)for x in row]for row in I],hi)]
    allroots=roots+[cap];allnames=names+['nextH4'];inc=[allroots[allnames.index(n)]for n in incident]
    lam=V.mv([[V.golden(x)for x in row]for row in inv],[V.pair(center,g,r)for r in inc]);assert all(x[1]<0 for x in lam)
    vertex=[V.sub(center[i],V.total(V.mul(t,r[i])for t,r in zip(lam,inc,strict=True)))for i in range(10)]
    assert V.pair(vertex,g,vertex)[1]<0 and V.dot(vertex,hi)[0]>0
    vertex_checks=[]
    for r,n in zip(allroots,allnames,strict=True):
        if n in incident:continue
        p=V.pair(vertex,g,r);assert p[1]<0;vertex_checks.append({'wall':n,'pairing':V.enc(p)})
    previous_point=[tuple(map(int,row))for row in previous_certificate['attached_vertex']['point_enclosures']]
    previous_pair=V.pair(cap,g,previous_point);assert previous_pair[1]<0
    out={'status':'PASS_EXACT_REMOTE_H4_ENDPOINT','facets':86,'source_sha256':V.sha(data['source']),
         'parameters':parameters,'normal_coordinate_enclosures':coordinates,'norm':2,
         'new_pair_counts':{'right':7,'pi3':1,'pi10':1,'dotted':76},'dotted_checks':records,
         'strict_center_pairing_in_dual_convention':V.enc(V.dot(cap,hi)),'target_cut_pairing_in_dual_convention':V.enc(cut),
         'actual_vertex':{'type':'H4+5A1','incident':incident,'gram':[[K.enc(x)for x in row]for row in M],
                          'point_enclosures':list(map(V.enc,vertex)),'nonincident_checks':vertex_checks},
         'old_H3_vertex_new_pairing':V.enc(previous_pair),
         'closed_old_edge':{'walls':S_NAMES+['wall17','coupled0'],'old_end_wall':'seed0','new_end_wall':'nextH4'},
         'scope':'One exact cap retains all85 previous walls, cuts their saved strict escape, and closes an actual edge with an H4+5A1 endpoint. Its norm and three metric constraints follow from the explicit orthogonal-frame formula. No exhaustive cap-family or finite-volume claim.'}
    if 'escape'in data:
        y=tuple(map(K.parse,data['escape']));q=K.pair(y,I,y);future=K.pair(y,I,h)
        assert K.sign(q)>0 and K.sign(future)<0
        vals=[V.dot(r,list(map(V.golden,y)))for r in allroots];assert min(x[0]for x in vals)>0
        out['remaining_escape']={'norm':K.enc(q),'future_center_pairing':K.enc(future),'strict_pairings':list(map(V.enc,vals))}
        out['scope']+=' A new strict future spacelike escape proves that the86-facet chamber still has infinite volume.'
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();data=V.read(PREFIX+'_caps.json');out=check(data)
    out['artifact_sha256']=V.sha(PREFIX+'_caps.json')
    if a.verify:assert V.read(PREFIX+'_verify.json')==out
    else:(V.ROOT/(PREFIX+'_verify.json')).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['facets'],'facets','EXACT_ESCAPE'if'remaining_escape'in out else'volume untested')
