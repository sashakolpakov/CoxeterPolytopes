#!/usr/bin/env python3
"""Producer-free exact coupled two-cap extension of the remote83 chamber."""
import argparse
import json
from fractions import Fraction as Q
from scripts import fresh_h9_remote_caps_20261013_verify as V
K=V.K
PREFIX=V.PREFIX+'_coupled'
S_NAMES=['seed0','seed1','seed2','seed3','seed5','seed6','seed7']

def reconstruct(data):
    assert data['source']==V.PREFIX+'_caps.json'and data['source_sha256']==V.sha(data['source'])
    previous=V.read(data['source']);V.check(previous)
    G,I,h,old,S0,caps0,point0,old_roots=V.reconstruct(previous)
    center=tuple(K.neg(x)for x in K.mv(I,h));g=[[V.golden(x)for x in row]for row in G]
    roots=old_roots+caps0;names=[r['name']for r in old['roots']]+['remote'+str(i)for i in range(4)]
    B=K.matrix(data['frame']);norms=tuple(map(K.parse,data['frame_norms']))
    assert len(B)==len(norms)==3 and all(len(r)==10 for r in B)
    assert K.sign(norms[0])<0 and all(K.sign(x)>0 for x in norms[1:])
    assert all(K.pair(B[i],G,B[j])==(norms[i]if i==j else K.Z)for i in range(3)for j in range(3))
    S=data['orthogonal_old_walls'];assert S==S_NAMES
    sr={r['name']:r for r in old['roots']};seeds=[tuple(map(K.parse,sr[n]['coefficients']))for n in S]
    assert all(sr[n]['field']=={'kind':'golden'}for n in S)
    assert all(K.pair(a,G,b)==(K.scale(K.O,2)if i==j else K.Z)for i,a in enumerate(seeds)for j,b in enumerate(seeds))
    assert all(K.pair(a,G,b)==K.Z for a in seeds for b in B)
    assert K.sign(K.pair(B[0],G,center))<0 and all(K.pair(b,G,center)==K.Z for b in B[1:])
    ds=[K.div(K.O,K.neg(norms[0])if i==0 else q)for i,q in enumerate(norms)]
    assert data['radicands']==list(map(K.enc,ds))
    vectors=[[V.mul(V.sqrt(V.golden(d)),V.golden(x))for x in row]for d,row in zip(ds,B,strict=True)]
    anchor=data['finite_anchor'];assert anchor=='wall17'and data['finite_cap_sign']==1 and data['finite_label']==3 and data['boost_branch']==-1
    ar=roots[names.index(anchor)]
    pairings=[V.pair(r,g,ar)for r in vectors]
    A=V.div(V.add(pairings[0],pairings[1]),V.rat(2))
    Bcoef=V.div(V.sub(pairings[1],pairings[0]),V.rat(2))
    C=V.add(pairings[2],V.rat(1))
    D=V.sub(V.mul(C,C),V.mul(V.rat(4),V.mul(A,Bcoef)))
    assert A[1]<0 and D[0]>0
    # The minus radical branch is selected exactly. It is the unique root >1.
    k=V.div(V.sub(V.neg(C),V.sqrt(D)),V.mul(V.rat(2),A))
    other=V.div(V.add(V.neg(C),V.sqrt(D)),V.mul(V.rat(2),A))
    assert k[0]>V.SCALE and other[1]<V.SCALE
    ct=V.div(V.add(V.mul(k,k),V.rat(1)),V.mul(V.rat(2),k));st=V.div(V.sub(V.mul(k,k),V.rat(1)),V.mul(V.rat(2),k))
    common=[V.add(V.mul(st,x),V.mul(ct,y))for x,y in zip(vectors[0],vectors[1],strict=True)]
    caps=[[V.add(x,V.mul(V.rat(s),z))for x,z in zip(common,vectors[2],strict=True)]for s in(1,-1)]
    point=[V.add(V.mul(ct,x),V.mul(st,y))for x,y in zip(vectors[0],vectors[1],strict=True)]
    # Exact identities: scaled frame squares(-1,1,1), ct²-st²=1.
    # Therefore two caps have norm2 and pairing0, P²=-1 and P is orthogonal
    # to both. The defining quadratic is equivalent to (cap0,anchor)=-1.
    # All these are identities, not approximate interval equality tests.
    return G,I,h,g,roots,names,caps,point,{'A':V.enc(A),'B':V.enc(Bcoef),'C':V.enc(C),'D':V.enc(D),'positive_boost':V.enc(k),'other_boost':V.enc(other)}

def check(data):
    G,I,h,g,old,names,caps,point,poly=reconstruct(data)
    hi=list(map(V.golden,h));center=[V.neg(x)for x in V.mv([[V.golden(x)for x in r]for r in I],hi)]
    checks=[]
    for i,r in enumerate(caps):
        assert V.dot(r,hi)[0]>0
        for oldr,n in zip(old,names,strict=True):
            if n in S_NAMES:continue
            if i==0 and n=='wall17':continue
            p=V.pair(r,g,oldr);assert p[1]<-2*V.SCALE,('not dotted',i,n,p)
            checks.append({'new':i,'old':n,'pairing':V.enc(p)})
    vertex_checks=[]
    for r,n in zip(old,names,strict=True):
        if n in S_NAMES:continue
        p=V.pair(point,g,r);assert p[1]<0;vertex_checks.append({'old':n,'pairing':V.enc(p)})
    assert V.dot(point,hi)[0]>0
    # Second actual vertex: H3 on seed3--5--wall17--3--cap0, six A1's.
    incident=S_NAMES+['wall17','coupled0']
    M=[[K.scale(K.O,2)if i==j else K.neg(K.P)if {i,j}=={3,7} else K.neg(K.O)if {i,j}=={7,8} else K.Z
        for j in range(9)]for i in range(9)]
    # Match every old anchor pairing to source, independently of named type.
    orig=V.read(V.OLD);source={r['name']:r for r in orig['roots']};anchor=source['wall17']
    # Source values may have a quadratic encoding; enclosure isolates the
    # asserted old Coxeter pair, while its exact old pair proof is hash anchored.
    oldpairs=V.read(orig['source'])['pair_checks']
    for j,n in enumerate(S_NAMES):
        assert oldpairs['|'.join(sorted((n,'wall17')))]=={'kind':'finite','label':5 if j==3 else 2}
    m=[r[:]for r in M]
    for j in range(9):
        assert K.sign(m[j][j])>0
        for a in range(j+1,9):
            for b in range(j+1,9):m[a][b]=K.sub(m[a][b],K.div(K.mul(m[a][j],m[j][b]),m[j][j]))
    inv=K.inverse(M);assert all(K.sign(x)>=0 for r in inv for x in r)
    allroots=old+caps;allnames=names+['coupled0','coupled1'];inc=[allroots[allnames.index(n)]for n in incident]
    products=[V.pair(center,g,r)for r in inc];lam=V.mv([[V.golden(x)for x in row]for row in inv],products)
    assert all(x[1]<0 for x in lam)
    second=[V.sub(center[i],V.total(V.mul(t,r[i])for t,r in zip(lam,inc,strict=True)))for i in range(10)]
    assert V.pair(second,g,second)[1]<0 and V.dot(second,hi)[0]>0
    second_checks=[]
    for r,n in zip(allroots,allnames,strict=True):
        if n in incident:continue
        p=V.pair(second,g,r);assert p[1]<0,(n,p);second_checks.append({'wall':n,'pairing':V.enc(p)})
    original_incident=['seed'+str(i)for i in range(8)]+['wall17']
    for j in range(8):
        assert oldpairs['|'.join(sorted(('seed'+str(j),'wall17')))]=={'kind':'finite','label':5 if j==3 else 2}
    original_gram=[[K.scale(K.O,2)if i==j else K.neg(K.P)if {i,j}=={3,8} else K.Z for j in range(9)]for i in range(9)]
    oi=K.inverse(original_gram);assert all(K.sign(x)>=0 for row in oi for x in row)
    assert K.sign(K.sub(K.scale(K.O,4),K.mul(K.P,K.P)))>0 # Only nontrivial 2x2 block.
    original_roots=[allroots[allnames.index(n)]for n in original_incident]
    olam=V.mv([[V.golden(x)for x in row]for row in oi],[V.pair(center,g,r)for r in original_roots])
    original_point=[V.sub(center[i],V.total(V.mul(t,r[i])for t,r in zip(olam,original_roots,strict=True)))for i in range(10)]
    assert V.pair(original_point,g,original_point)[1]<0 and V.dot(original_point,hi)[0]>0
    original_checks=[]
    for r,n in zip(allroots,allnames,strict=True):
        if n in original_incident:continue
        p=V.pair(original_point,g,r);assert p[1]<0;original_checks.append({'wall':n,'pairing':V.enc(p)})
    target=tuple(map(K.parse,data['target_point']));assert list(map(K.enc,target))==V.read(data['source'])['escape']
    assert K.sign(K.pair(target,I,target))>0 and K.sign(K.pair(target,I,h))<0
    assert all(V.dot(r,list(map(V.golden,target)))[0]>0 for r in old)
    cuts=[V.dot(r,list(map(V.golden,target)))for r in caps];assert all(x[1]<0 for x in cuts)
    out={'status':'PASS_EXACT_FINITE_ANGLE_REMOTE_PAIR','total_facets':85,'new_facets':2,
         'norms':[2,2],'mutual_label':2,'new_old_right_pairings':14,'new_old_finite_extra':{'new':0,'old':'wall17','label':3},
         'new_old_dotted_pairings':len(checks),'dotted_checks':checks,'boost_equation':poly,
         'target_cut_pairings':list(map(V.enc,cuts)),
         'joint_vertex':{'type':'9A1','norm':'-1','nonincident_checks':vertex_checks},
         'attached_vertex':{'type':'H3+6A1','incident':incident,'gram':[[K.enc(x)for x in row]for row in M],
                            'point_enclosures':list(map(V.enc,second)),'nonincident_checks':second_checks},
         'original_vertex':{'type':'H2+7A1','incident':original_incident,'point_enclosures':list(map(V.enc,original_point)),
                            'nonincident_checks':original_checks},
         'finite_edges':[{'common_facets':S_NAMES+['coupled0'],'ends':['joint_vertex','attached_vertex']},
                         {'common_facets':S_NAMES+['wall17'],'ends':['original_vertex','attached_vertex']}],
         'scope':'Exact two-cap extension preserves all83 walls, repairs the saved second gap, and has two actual finite vertices joined by a finite edge. It attaches to the existing rank-eight spherical face S+wall17. The remote83 certificate and frozen79 angle proof are dependencies.'}
    if'escape'in data:
        y=tuple(map(K.parse,data['escape']));q=K.pair(y,I,y);future=K.pair(y,I,h)
        assert K.sign(q)>0 and K.sign(future)<0
        vals=[V.dot(r,list(map(V.golden,y)))for r in allroots];assert min(x[0]for x in vals)>0
        out['remaining_escape']={'norm':K.enc(q),'future_center_pairing':K.enc(future),'strict_pairings':list(map(V.enc,vals))}
        out['scope']+=' A saved future strict spacelike point proves this85-facet chamber still has infinite volume.'
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();data=V.read(PREFIX+'_caps.json')
    out=check(data);out['source_sha256']=V.sha(PREFIX+'_caps.json')
    if a.verify:assert V.read(PREFIX+'_verify.json')==out
    else:(V.ROOT/(PREFIX+'_verify.json')).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['total_facets'],'facets','EXACT_ESCAPE'if'remaining_escape'in out else'volume untested')
