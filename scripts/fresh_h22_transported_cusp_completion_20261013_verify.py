#!/usr/bin/env python3
"""Exact complete symmetry-orbit repair of 483840 specified H22 cusp links.

Checks the equivariant frame arrays, all changed pairs by orbit reduction,
and a full-support positive relation at the representative cusp. No optimizer
or floating-point predicate is used.
"""
import argparse,hashlib,json
from fractions import Fraction as F
from collections import Counter
import numpy as np
from fresh_h22_new_cusp_completion_20261013_verify import setup,P,read,sha,RAW
from fresh_h22_new_cusp_completion_20261013_tight import shifted,DEN as OLD_DEN
from fresh_h22_incidence_field_obstruction_20261013 import rank
PRE='fresh_h22_transported_cusp_completion_20261013'
OUT='results/'+PRE+'_verify.json'
FRAME='data/enumerations/'+PRE+'_frames.npz'
META='data/enumerations/'+PRE+'_frames.json'
PROPOSAL='data/enumerations/'+PRE+'_local_proposal.json'
OLD='results/fresh_h22_incidence_binary_20261013.json'
SHIFTS=[3,4,3,4,4,4,4,4,5,8,7,9,6]
SPEC=[(i,1)for i in range(5)]+[(i,s)for i in range(5,9)for s in(1,-1)]
def sign(a,b=0,d=1):
    a,b=F(a),F(b)
    if not b:return (a>0)-(a<0)
    if not a:return (b>0)-(b<0)
    if (a>0)==(b>0):return 1 if a>0 else-1
    c=a*a-d*b*b
    if c==0:return 0
    return (1 if a>0 else-1) if c>0 else (1 if b>0 else-1)
def dotted(a,b,d,den,q,n):
    a,b=int(a),int(b);assert den>0 and q>0 and n>0
    assert sign(a,b,d)<0,(a,b,d,den,q,n)
    assert sign(F(a*a+d*b*b)-q*n*den*den,2*a*b,d)>0,(a,b,d,den,q,n)
def label(p,q,n):
    p=F(p);assert p<=0
    if not p:return 'right'
    ratio=p*p/(q*n)
    if ratio in(F(1,4),F(1,2),F(3,4)):
        return {F(1,4):'pi/3',F(1,2):'pi/4',F(3,4):'pi/6'}[ratio]
    assert ratio>=1,(p,q,n,ratio)
    return 'parallel'if ratio==1 else'dotted'
def load():
    lower=read('results/fresh_h22_new_cusp_completion_20261013_lower.json')
    assert sha(lower['source'])==lower['source_sha256']
    assert all(sha(p)==h for p,h in read(lower['source'])['sources'].items())
    old=read(OLD);assert all(sha(p)==h for p,h in old['sources'].items())
    _,Q,W,h,meta,gens,D,V,B,A,ids,sizes=setup()
    fm=read(META);assert sha(FRAME)==fm['frames_sha256']
    z=np.load(P/FRAME);bases=z['old_base_numerators'];height=z['height_numerators'];E=z['six_delta_old'];den=int(z['denominator']);seed=fm['seed']
    q0=[F(x)for x in fm['base_norms']];q=[q0[j]for j,s in SPEC]
    C=np.array([s*bases[seed,j]+t*(den//6)*E[seed]for (j,s),t in zip(SPEC,SHIFTS)])
    Z=np.array([s*height[j]+2*t*den for (j,s),t in zip(SPEC,SHIFTS)])
    return Q,W,h,meta,gens,D,V,B,A,ids,sizes,z,bases,height,E,den,seed,q,C,Z
def families(bases,height,E,den):
    for k,((j,s),t)in enumerate(zip(SPEC,SHIFTS)):
        yield k,s*bases[:,j]+t*(den//6)*E,int(s*height[j]+2*t*den)
def run():
    Q,W,h,meta,gens,D,V,B,A,ids,sizes,z,bases,height,E,den,seed,q,C,Z=load();N=len(E)
    assert N==241920 and bases.shape==(N,9,22) and den==32760
    assert int(abs(E).max())*22*max(int(abs(g).max())for g in gens)<2**63
    assert 22*int(abs(bases).max())*max(int(abs(g).max())for g in gens)<2**63
    # Independent direct closure, not just a spanning tree of proposals.
    permutations=z['permutations']
    for j,(g,m)in enumerate(zip(gens,meta)):
        assert np.array_equal(g.T@Q@g,Q)and np.array_equal(g@h,h)
        assert np.array_equal(W@g.T,W[m['wall_permutation']])
        perm=permutations[j];assert len(np.unique(perm))==N
        assert np.array_equal(E@g.T,E[perm])
        assert np.array_equal(bases@g.T,bases[perm])
    parents=z['parents'];edge=z['generator_indices'];depth=z['depth']
    assert np.count_nonzero(depth==0)==1 and depth[seed]==0 and np.all(depth>=0)
    for i in np.flatnonzero(depth>0):
        assert depth[parents[i]]==depth[i]-1 and permutations[edge[i],parents[i]]==i
    prop=read(PROPOSAL);assert prop['orbit_index']==seed
    oldbases=[[F(x)for x in r]for r in prop['bases']]
    assert all(F(int(C[k,j]),den)==oldbases[k][j]+SHIFTS[k]*F(int(E[seed,j]),6)for k in range(13)for j in range(22))
    assert all(F(int(Z[k]),den)==oldbases[k][22]+2*SHIFTS[k]for k in range(13))
    assert np.all(Z>0)
    cv=C@Q;assert 22*int(abs(C).max())*int(abs(cv).max())+11*int(Z.max())**2<2**63
    G=C@cv.T+11*Z[:,None]*Z
    assert [F(int(G[i,i]),den*den)for i in range(13)]==q
    assert np.all(C@(Q@h)<0)
    local=Counter(label(F(int(G[i,j]),den*den),q[i],q[j])for i in range(13)for j in range(i))
    cov=cv.T;pair_bounds={}
    def products(rows):
        assert 22*int(abs(rows).max())*int(abs(cov).max())<2**63
        return rows@cov
    # Old210 has the only nontrivial finite old-root contacts.
    wn=np.einsum('ij,ij->i',W,W@Q);wp=products(W);oldcounts=[]
    incident=prop['target']['incident'];extra=prop['target']['extra']
    for k in range(13):
        cc=Counter()
        for i in range(210):
            kind=label(F(int(wp[i,k]),den),q[k],F(int(wn[i])));cc[kind]+=1
            if i not in incident:assert kind=='dotted'
        oldcounts.append(dict(cc))
    # The new four and560 rational-height walls are all nonincident.
    for name,path,key in [('core4','results/fresh_h22_incidence_core_20261013.json','old_projections'),('binary560',OLD,'accepted_old_projections')]:
        rows=np.array(read(path)[key],dtype=np.int64);mx=products(rows).max(axis=0)
        for k in range(13):dotted(mx[k],2*Z[k],11,den,q[k],F(2))
        pair_bounds[name]=mx.tolist()
    # U: p=(a+Z sqrt110)/(2den); K: p=(4a+four Z sqrt11)/(4den).
    uf='data/enumerations/fresh_h22_symmetry_orbit_20261013.npz';U=np.load(P/uf)['twice_old_projections'];mx=products(U).max(axis=0)
    for k in range(13):dotted(mx[k],Z[k],110,2*den,q[k],F(1))
    pair_bounds['U']=mx.tolist()
    kf='data/enumerations/fresh_h22_cusp_orbits_20261013.npz';ks=np.load(P/kf);kp=products(ks['old_null_projections']);km={}
    for four in(1,4):
        mx=kp[ks['four_heights']==four].max(axis=0);km[str(four)]=mx.tolist()
        for k in range(13):dotted(4*int(mx[k]),four*Z[k],11,4*den,q[k],F(four*four,16))
    pair_bounds['K']=km
    rf='data/enumerations/fresh_h22_gap_face_orbit_20261013.npz';R=np.load(P/rf)['six_old_projections'];mx=products(R).max(axis=0)
    for k in range(13):dotted(int(mx[k])+66*int(Z[k]),0,1,6*den,q[k],F(2,3))
    pair_bounds['R']=mx.tolist()
    # Exactly one earlier repair wall is incident to this representative cusp.
    high_index=int(np.load(P/('data/enumerations/'+PRE+'_cusps.npz'))['highest_indices'][seed])
    oldcap=None;repair={};an=np.einsum('ij,ij->i',A,A@Q)
    for name,rows,t,ix,index in shifted(D,V,B,A,ids,sizes,Q,True):
        pp=products(rows);records={}
        nn=an if index==0 else np.zeros(len(rows),dtype=np.int8)
        for height_value in np.unique(t):
            for normclass in np.unique(nn[t==height_value]):
                mask=(t==height_value)&(nn==normclass)
                normold=F(int(normclass))if index==0 else[None,F(279,14),F(279,14),F(55,62496),F(55,62496)][index]
                if index==0 and mask[high_index]:
                    oldcap=[F(int(x),OLD_DEN)for x in rows[high_index]]+[F(int(t[high_index]),10)]
                    for k in range(13):
                        value=10*int(pp[high_index,k])+11*int(height_value)*int(Z[k])*OLD_DEN
                        assert value==0
                        dotted(10*int(pp[high_index,k])-11*int(height_value)*int(Z[k])*OLD_DEN,0,1,10*OLD_DEN*den,q[k],normold)
                    mask[high_index]=False
                mx=pp[mask].max(axis=0)
                for k in range(13):dotted(10*int(mx[k])+11*int(height_value)*int(Z[k])*OLD_DEN,0,1,10*OLD_DEN*den,q[k],normold)
                records[f't={int(height_value)},q={normold}']=mx.tolist()
        repair[name]=records
    assert oldcap is not None
    print('retained bank pair gates PASS',flush=True)
    # Every nonlocal new pair is strictly dotted. Same-cusp pairs are local G.
    cross={}
    for k,rows,zz in families(bases,height,E,den):
        pp=products(rows);bound=22*int(abs(rows).max())*int(abs(cov).max())+11*abs(zz)*int(Z.max());assert bound<2**63
        positive=pp+11*zz*Z;negative=pp-11*zz*Z
        assert np.array_equal(positive[seed],G[k])
        positive[seed]=np.iinfo(np.int64).min
        mp=positive.max(axis=0);mm=negative.max(axis=0)
        for j in range(13):
            dotted(mp[j],0,1,den*den,q[k],q[j]);dotted(mm[j],0,1,den*den,q[k],q[j])
        cross[str(k)]={'same_sign_nonlocal_maxima':mp.tolist(),'opposite_sign_maxima':mm.tolist()}
    print('complete new orbit pair gates PASS',flush=True)
    # Full-support positive relation and exact rank at the representative.
    Q11=[list(map(int,r))+[0]for r in Q]+[[0]*22+[11]]
    delta=[F(int(x),6)for x in E[seed]]+[F(2)]
    assert [oldcap[j]+(int(W[extra,j])if j<22 else 0)for j in range(23)]==delta
    normals=[list(map(F,W[i]))+[F(0)]for i in incident]+[oldcap]
    normals.extend([[F(int(x),den)for x in C[k]]+[F(int(Z[k]),den)]for k in range(13)])
    def dot(a,b):return sum(x*y for x,y in zip(a,b))
    def mv(a):return [dot(r,a)for r in Q11]
    assert len(normals)==31 and rank(normals)==22
    gram=[[dot(a,mv(b))for b in normals]for a in normals];assert rank(gram)==21
    assert all(dot(delta,mv(r))==0 for r in normals)
    assert [sum(r[j]for r in normals)for j in range(23)]==[(1+sum(SHIFTS))*x for x in delta]
    target=read('results/fresh_h22_incidence_binary_20261013_escape.json');X=target['numerator'];xc=Q.astype(object)@np.array(X[:22],dtype=object)
    cuts=[]
    for k in range(13):
        a=sum(int(u)*int(v)for u,v in zip(C[k],xc));b=int(Z[k])*X[22]
        if sign(a,b,11)>0:cuts.append({'profile':k,'rational_numerator':str(a),'sqrt11_numerator':str(b),'denominator':den*target['denominator']})
    assert cuts
    sources=[RAW,OLD,FRAME,META,PROPOSAL,'data/enumerations/'+PRE+'_cusps.npz','results/fresh_h22_new_cusp_completion_20261013_lower.json','results/fresh_h22_new_cusp_completion_20261013_verify.json','results/fresh_h22_incidence_core_20261013_full.json']
    return {'status':'PASS_EXACT_COMPLETE_TRANSPORTED_H22_CUSP_REPAIR','sources':{p:sha(p)for p in sources},'positive_cusps':N,'completed_signed_cusps':2*N,'new_caps':26*N,'facets':read(OLD)['facets']+26*N,'denominator':den,'shifts':SHIFTS,'representative_old_numerators':C.tolist(),'representative_height_numerators':Z.tolist(),'new_norms':list(map(str,q)),'old210_pair_counts':oldcounts,'other_retained_pair_bounds':pair_bounds,'old_repair_pair_bounds':repair,'local_new_pair_counts':dict(local),'nonlocal_pair_bounds':cross,'bounded_link':{'incident_old_indices':incident,'incident_old_repair_entry':high_index,'incident_normals':31,'normal_span_rank':22,'quotient_rank':21,'horosphere_dimension':21,'all_positive_weights':1,'sum_of_normals_delta_multiple':1+sum(SHIFTS)},'previous_escape_cutters':cuts,'scope':'Coherent13-wall repair at every one of241920 positive cusp rays and their negative-height copies. The complete same-generator orbit is checked. All nonincident retained and nonlocal new pairs are dotted; all specified links are bounded. Previous escape is cut. A global volume certificate is still required.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
    if a.verify:assert out==read(OUT)
    else:(P/OUT).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['facets'],flush=True)
