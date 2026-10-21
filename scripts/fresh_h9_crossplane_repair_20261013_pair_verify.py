#!/usr/bin/env python3
"""Producer-free exact orthogonal-pair replacement and full new-affine replay."""
import argparse
from itertools import combinations
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts import fresh_h9_crossplane_repair_20261013_verify as T
V=T.V;K=V.K;B=T.B;F=T.F
PREFIX=T.PREFIX+'_pair'
CONSTRAINTS=[['seed1',2],['seed2',2],['seed3',5],['seed4',3],['seed6',3],['seed8',5],['wall17',2],['nextH4',5]]
NEW=['crosspair0','crosspair1']

def reconstruct(data):
    assert data['source']==B.PREFIX+'_caps.json' and data['source_sha256']==V.sha(data['source'])
    previous=V.read(data['source']);B.check(previous)
    G,I,h,g,old,names,caps,center,_=B.reconstruct(previous);old+=caps;names+=['bridge0','bridge1']
    assert names[6:8]==['seed6','seed7']
    for j in(6,7):assert old[j]==[V.rat(int(i==j))for i in range(10)]
    assert data['constraints']==CONSTRAINTS and data['difference_pairing']==2 and data['branch'] in(-1,1)
    D=[V.sub(a,b)for a,b in zip(old[6],old[7],strict=True)]
    # The difference is an exact Q(phi) vector of norm4.
    dq=[K.sub(a,b)for a,b in zip([K.O if i==6 else K.Z for i in range(10)],
                                [K.O if i==7 else K.Z for i in range(10)],strict=True)]
    assert K.pair(dq,G,dq)==K.scale(K.O,4)
    A=[old[names.index(n)]for n,m in CONSTRAINTS]+[D]
    M=[[V.rat(4 if i==8 else 2)if i==j else V.pair(a,g,b)for j,b in enumerate(A)]for i,a in enumerate(A)]
    Mi,pivots=T.inverse(M)
    p=[V.rat(0)if m==2 else V.rat(-1)if m==3 else V.neg(V.PHI)for n,m in CONSTRAINTS]+[V.rat(2)]
    lam=V.mv(Mi,p);u=[V.total(V.mul(t,a[i])for t,a in zip(lam,A,strict=True))for i in range(10)]
    hc=[V.pair(center,g,a)for a in A];lc=V.mv(Mi,hc)
    v=[V.sub(center[i],V.total(V.mul(t,a[i])for t,a in zip(lc,A,strict=True)))for i in range(10)]
    qu=V.dot(p,lam);qv=V.sub(V.pair(center,g,center),V.dot(hc,lc));assert qv[0]>0
    rem=V.sub(V.rat(2),qu);assert rem[0]>0
    t=V.sqrt(V.div(rem,qv));t=t if data['branch']==1 else V.neg(t)
    r=[V.add(x,V.mul(t,y))for x,y in zip(u,v,strict=True)]
    # Reflection in D gives s=r-D: s²=2 and r.s=0 exactly.
    s=[V.sub(a,b)for a,b in zip(r,D,strict=True)]
    params={'inverse_nonzero_pivots':pivots,'u_squared':V.enc(qu),'v_squared':V.enc(qv),
            'norm_remainder':V.enc(rem),'coefficient':V.enc(t)}
    return G,I,h,g,old,names,[r,s],center,params

def finite_profiles(G,names):
    _,pairs=F.pair_table();first=dict(CONSTRAINTS)
    second={('seed7'if n=='seed6'else n):m for n,m in CONSTRAINTS}
    # Every unchanged prescribed wall is perpendicular to D. Exact seed
    # pairings handle seed8 even though its two separate relations are dotted.
    for n,m in CONSTRAINTS:
        if n=='seed6':continue
        if n.startswith('seed'):
            j=int(n[4:]);assert G[j][6]==G[j][7]
        else:
            a=pairs['|'.join(sorted((n,'seed6')))];b=pairs['|'.join(sorted((n,'seed7')))]
            assert a==b and a['kind']=='finite'
    return [first,second]

def census(names,profiles):
    oldnames,pairs=F.pair_table();assert oldnames==names
    for new,profile in zip(NEW,profiles,strict=True):
        for n in names:pairs['|'.join(sorted((new,n))) ]={'kind':'finite','label':profile[n]}if n in profile else{'kind':'dotted'}
    pairs['|'.join(NEW)]={'kind':'finite','label':2}
    def pair(a,b):return pairs['|'.join(sorted((a,b)))]
    vals={2:K.Z,3:K.neg(K.O),5:K.neg(K.P)}
    checked=0;dotted=0;indefinite=0;positive=[];affine=[]
    for count in(1,2):
        for new in combinations(NEW,count):
            common=sorted(set.intersection(*(set(profiles[NEW.index(n)])for n in new)))
            for mask in range(1<<len(common)):
                ss=list(new)+[n for i,n in enumerate(common)if mask>>i&1];checked+=1
                if any(pair(a,b)['kind']=='dotted'for a,b in combinations(ss,2)):dotted+=1;continue
                M=[[K.scale(K.O,2)if a==b else vals[pair(a,b)['label']]for b in ss]for a in ss]
                kind=T.psd(M)
                if kind=='positive_semidefinite':affine.append(ss)
                elif kind=='positive_definite':positive.append(ss)
                else:indefinite+=1
    inherited=V.read(B.PREFIX+'_frontier.json');assert inherited['inherited_affine_circuits']==inherited['new_affine_circuits']==0
    assert not affine,affine
    vertices=[s for s in positive if len(s)==9];assert not vertices
    maximal=[s for s in positive if not any(set(s)<set(o)for o in positive)]
    return {'all_pair_count':len(pairs),'new_subset_count':checked,'dotted_subsets':dotted,
            'indefinite_subsets':indefinite,'positive_definite_subsets':len(positive),
            'new_affine_circuits':0,'inherited_affine_circuits':0,
            'new_finite_vertices':vertices,'maximal_positive_definite_blocks_meeting_pair':maximal,
            'scope':'Every possible new positive-semidefinite subset contains one or both new walls; all its old walls must occur in their common finite-neighbor list. All such subsets are checked exactly over Q(phi).'}

def check(data):
    G,I,h,g,old,names,caps,center,params=reconstruct(data);profiles=finite_profiles(G,names)
    hi=list(map(V.golden,h));records=[]
    for r,profile,name in zip(caps,profiles,NEW,strict=True):
        cp=V.dot(r,hi);assert cp[0]>0
        norm=V.pair(r,g,r);assert norm[0]<=2*V.SCALE<=norm[1]
        dotted=[]
        for a,n in zip(old,names,strict=True):
            if n in profile:continue
            p=V.pair(a,g,r);assert p[1]<-2*V.SCALE,(name,n,p);dotted.append({'wall':n,'pairing':V.enc(p)})
        assert len(dotted)==80
        records.append({'name':name,'normal_enclosures':list(map(V.enc,r)),
                        'norm':2,'finite_relations':profile,'center_pairing_dual':V.enc(cp),'dotted_checks':dotted})
    cross=V.pair(caps[0],g,caps[1]);assert cross[0]<=0<=cross[1]
    targets=[V.read(data['source'])['escape'],V.read(T.PREFIX+'_caps.json')['escape']]
    assert data['targets']==targets
    cuts=[];target_feasibility=[]
    for point in targets:
        y=tuple(map(K.parse,point));assert K.sign(K.pair(y,I,y))>0 and K.sign(K.pair(y,I,h))<0
        oldpp=[V.dot(r,list(map(V.golden,y)))for r in old];assert all(p[0]>0 for p in oldpp)
        target_feasibility.append(list(map(V.enc,oldpp)))
        pp=[V.dot(r,list(map(V.golden,y)))for r in caps];assert all(p[1]<0 for p in pp)
        cuts.append(list(map(V.enc,pp)))
    out={'status':'PASS_EXACT_ORTHOGONAL_PAIR_REPAIR','facets':90,'old_facets_retained':88,
         'replaces':'The separate crossplane0 cap is omitted; this is an alternative extension of the original88.',
         'parameters':params,'new_normals':records,'new_pair':{'pairing':0,'label':2},
         'target_pairings_dual':cuts,'targets_old88_strict_pairings':target_feasibility,
         'affine_and_vertex_census':census(names,profiles),
         'scope':'The exact coupled pair cuts both saved symmetric gaps. All90 inequalities are actual facets by strict-center projection and nonpositive pairings. No positive affine circuits occur, but this does not establish compactness or finite volume.'}
    if 'escape'in data:
        y=tuple(map(K.parse,data['escape']));q=K.pair(y,I,y);future=K.pair(y,I,h)
        assert K.sign(q)>0 and K.sign(future)<0
        pp=[V.dot(a,list(map(V.golden,y)))for a in old+caps];assert min(p[0]for p in pp)>0
        out['remaining_escape']={'norm':K.enc(q),'future_center_pairing':K.enc(future),'strict_pairings':list(map(V.enc,pp))}
        out['scope']+=' A strict future spacelike escape proves that this90-facet chamber has infinite volume.'
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=V.read(PREFIX+'_caps.json');out=check(d);out['source_sha256']=V.sha(PREFIX+'_caps.json')
    if a.verify:assert out==V.read(PREFIX+'_verify.json')
    else:(V.ROOT/(PREFIX+'_verify.json')).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'EXACT_ESCAPE'if'remaining_escape'in out else'VOLUME_UNTESTED')
