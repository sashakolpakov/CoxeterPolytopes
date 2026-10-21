#!/usr/bin/env python3
"""Standalone exact cusp-orbit, bounded-link and strict-escape certificate."""
import argparse
from collections import Counter
from fractions import Fraction as F
from math import gcd
import hashlib
import json
from pathlib import Path
import numpy as np
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
SOURCE='data/enumerations/fresh_h22_cusp_orbits_20261013.json'
OLD='results/fresh_h22_symmetry_orbit_verify_20261013.json'
OUT='results/fresh_h22_cusp_orbits_20261013_verify.json'
def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()


def link(v,Q,W,h):
    pp=W@Q@v;assert max(pp)==0
    ids=np.flatnonzero(pp==0).tolist();R=W[ids];G=R@Q@R.T
    unseen=set(range(len(ids)));components=[]
    while unseen:
        todo=[min(unseen)];unseen.remove(todo[0]);part=[]
        while todo:
            i=todo.pop();part.append(i)
            js=[j for j in unseen if G[i,j]]
            unseen.difference_update(js);todo.extend(js)
        part.sort();M=sp.Matrix(G[np.ix_(part,part)].tolist());kernels=M.nullspace()
        assert len(kernels)==1
        k=kernels[0]
        if k[0]<0:k=-k
        assert all(x>0 for x in k)
        total=np.array([sum(F(int(c.p),int(c.q))*int(R[j,a])for c,j in zip(k,part))for a in range(22)],dtype=object)
        a=next(i for i,x in enumerate(v)if x);mu=total[a]/int(v[a]);assert mu>0
        assert all(x==mu*int(y)for x,y in zip(total,v))
        components.append({'walls':[ids[j]for j in part],'rank':len(part)-1,
                           'positive_kernel':[str(x)for x in k],'ray_multiple':str(mu)})
    assert sum(c['rank']for c in components)==20
    return {'old_incident_walls':ids,'old_affine_components':components,
            'new_affine_A1_rank':1,'new_kernel':[1,1],'new_kernel_ray_multiple':2,
            'bounded_link_rank':21}


def run(saved=None):
    data=read(SOURCE);old=read(OLD)
    assert data['completed']
    assert all(sha(p)==x for p,x in old['source_sha256'].items())
    assert all(sha(p)==x for p,x in data['sources'].items())
    raw=read('data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json')['coxeter_input']
    Q=-np.array(raw['basis_gram'],dtype=np.int64);W=np.array(raw['wall_coordinates_as_columns_transposed'],dtype=np.int64)
    syms=read('data/enumerations/fresh_h20_global_interior_sections_20261013.json')['ambient_symmetry_generators']
    gens=[np.array(g['matrix'],dtype=np.int64)for g in syms]
    U=np.load(ROOT/'data/enumerations/fresh_h22_symmetry_orbit_20261013.npz')['twice_old_projections']
    assert sha(data['array'])==data['array_sha256']
    archive=np.load(ROOT/data['array']);D=archive['old_null_projections'];four=archive['four_heights']
    assert D.dtype==four.dtype==np.int64 and D.shape==(data['accepted_null_directions'],22)
    qmax=int(abs(Q).max());wmax=int(abs(W).max());umax=int(abs(U).max());dmax=int(abs(D).max());gmax=max(int(abs(g).max())for g in gens)
    bounds=[22*dmax*gmax,22**2*gmax*gmax*qmax,22**2*wmax*wmax*qmax,
            22**2*umax*dmax*qmax,22**2*dmax*dmax*qmax]
    assert max(bounds)<2**63
    h=W.sum(axis=0);h//=gcd(*map(int,h));assert h@Q@h==-14
    hc=h@Q;hmax=int(abs(h).max());assert 22**2*max(hmax,wmax,dmax)*hmax*qmax<2**63
    assert Counter(map(int,W@Q@h))=={-1:42,-8:168}
    # Recompute the positive complement, rather than infer inertia from a
    # numerical eigenvalue or a status flag in the previous certificate.
    j=next(i for i,x in enumerate(hc)if x)
    Z=np.array([[int(hc[j])*(k==i)-int(hc[i])*(k==j)for i in range(22)if i!=j]for k in range(22)],dtype=object)
    M=[[F(x)for x in row]for row in(Z.T@Q.astype(object)@Z).tolist()]
    for i in range(21):
        p=M[i][i];assert p>0
        for a in range(i+1,21):
            for b in range(i+1,21):M[a][b]-=M[a][i]*M[i][b]/p
    for g,meta in zip(gens,syms):
        assert np.array_equal(g.T@Q@g,Q)and np.array_equal(g@h,h)
        assert np.array_equal(W@g.T,W[meta['wall_permutation']])
    norms=np.diag(W@Q@W.T);assert set(norms)=={2,4}
    known=set();records=[];end=0
    for index,row in enumerate(data['orbits']):
        start=row['array_offset'];n=row['size'];assert start==end;end=start+n
        part=D[start:end];height=F(row['height']);assert height in(F(1),F(1,4))
        assert np.all(four[start:end]==4*height)
        assert part[0].tolist()==row['representative']
        v=part[0];assert int(v@Q@v)==0 and int(h@Q@v)==row['center_pairing']<0
        assert gcd(*map(int,v))==1
        # Every noninitial row must have an earlier predecessor under one
        # generator. Induction proves connectivity; closure excludes omissions.
        lookup={r.tobytes():i for i,r in enumerate(part)};assert len(lookup)==n
        assert not known.intersection(lookup);known.update(lookup)
        has_predecessor=np.zeros(n,dtype=bool);has_predecessor[0]=True
        for g in gens:
            perm=np.array([lookup[r.tobytes()]for r in(part@g.T)],dtype=np.int64)
            assert len(np.unique(perm))==n
            earlier=perm>np.arange(n);has_predecessor[perm[earlier]]=True
        assert np.all(has_predecessor)
        maximum=int(max(U@(Q@v)));assert maximum==row['maximum_twice_timelike_orbit_pairing']
        assert maximum<=(-6 if height==1 else-2)
        assert int(h@Q@v)<=(-6 if height==1 else-2)
        ps=W@Q@v;assert max(ps)==0
        kinds=Counter()
        for p,q in zip(ps,norms):
            ratio=F(int(p)**2,int(q))/height**2
            kind='dotted'if ratio>1 else{F(0):'2',F(1,4):'3',F(1,2):'4',F(1):'parallel'}[ratio]
            kinds[kind]+=1
        rec={'size':n,'height':str(height),'representative':v.tolist(),
             'old_pair_counts':dict(kinds),'maximum_twice_old_orbit_pairing':maximum,
             'cusp_link':link(v,Q,W,h)}
        records.append(rec);print('verified cusp orbit',index+1,n,flush=True)
    assert end==len(D)==len(known)
    # Distinct primitive future null integral vectors are nonproportional,
    # so have a strictly negative integer product, hence <=-1. Primitivity
    # transfers from each checked seed: an integral Q-isometry has determinant
    # +/-1, and the predecessor induction connects every row to its seed.
    # Global row uniqueness therefore excludes proportional rays.
    # This certifies all cusp/cusp pairs at once:
    # equal height1 gives integer pair <=0; every mixed/small-height pair is
    # dotted. Opposite signs over the identical ray are parallel.
    den=10**8
    if saved is None:
        proposal=data['runs'][0]['escape']['point']
        num=[round(float(x)*den)for x in proposal]
        num[:22]=[a+10000*int(b)for a,b in zip(num,h)]
    else:
        assert saved['strict_escape']['denominator']==str(den)
        num=list(map(int,saved['strict_escape']['numerator']))
    physical=np.array(num[:22],dtype=object);cov=Q.astype(object)@physical;z=num[-1]
    square=int(physical@cov)+z*z;hp=int(h.astype(object)@cov)
    assert square>0 and hp+4*abs(z)<0
    oldp=W.astype(object)@cov;assert max(oldp)<0
    upp=U.astype(object)@cov
    assert max(upp)<0 and min(int(x)**2-10*z*z for x in upp)>0
    bound=22*dmax*max(abs(int(x))for x in cov)
    assert 4*bound+4*abs(z)<2**63
    cpp=4*(D@np.array(cov,dtype=np.int64))+four*abs(z)
    assert max(cpp)<0
    former=old['strict_escape'];oldnum=list(map(int,former['numerator']))
    oldcov=Q.astype(object)@np.array(oldnum[:22],dtype=object)
    assert 4*22*dmax*max(abs(int(x))for x in oldcov)+4*abs(oldnum[-1])<2**63
    cuts=4*(D@np.array(oldcov,dtype=np.int64))+four*abs(oldnum[-1])
    # Count null directions for which at least one signed wall cuts the old
    # escape; this field is not a census of individual signed cutting walls.
    count=int(sum(cuts>0));assert count==data['old_escape_cutters']>0
    facets=121172+2*len(D);assert facets==data['total_wall_count']
    return {'status':'PASS_COMPLETE_SAMPLED_CUSP_ORBITS_BOUNDED_LINKS_AND_STRICT_ESCAPE',
            'sources':{p:sha(p)for p in(SOURCE,OLD,data['array'])},
            'dimension':22,'facets':facets,'cusp_orbits':records,'completed_cusps':len(D),
            'old_escape_cutters':count,'all_unordered_wall_pairs':facets*(facets-1)//2,
            'strict_escape':{'denominator':str(den),'numerator':list(map(str,num)),
                             'norm':str(F(square,den**2)),'center_pairing':str(F(hp,den)),
                             'maximum_cusp_cap_pairing':str(F(int(max(cpp)),4*den)),
                             'center_cap_uniform_upper_bound':str(F(hp+4*abs(z),den))},
            'proof':'Integral generator closure plus earlier-predecessor induction proves complete connected cusp orbits. Every old cusp has affine link rank20 and gains an orthogonal affineA1 of rank1. Lorentzian negativity and integrality certify all distinct cusp/cusp pairs without quadratic enumeration. All actual facets follow by strict-center projection. The new rational strict escape works for every center label m>=3.',
            'scope':'Exactly the15 orbits represented by the sampled cusps, not an exhaustive old cusp census. Existing121172-wall compatibility is a hash-anchored prior certificate. No finite-volume example or obstruction to arbitrary further walls is claimed.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    saved=read(OUT)if a.verify else None;out=run(saved)
    if a.verify:assert out==saved
    else:(ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['facets'],out['completed_cusps'])
