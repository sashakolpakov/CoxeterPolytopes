#!/usr/bin/env python3
"""Transport nine cusp base directions through the exact 241920-element orbit.

This constructs proposed frames. Global new-new Coxeter checks are separate.
"""
from collections import deque
from fractions import Fraction as F
import hashlib,json
from math import lcm
import numpy as np
from fresh_h22_new_cusp_completion_20261013_verify import setup,P,read
PRE='fresh_h22_transported_cusp_completion_20261013'
def run():
    _,Q,W,h,meta,gens,*_=setup()
    source='data/enumerations/'+PRE+'_cusps.npz';z=np.load(P/source);D=z['six_old'];N=len(D)
    proposal=read('data/enumerations/'+PRE+'_local_proposal.json');seed=proposal['orbit_index']
    assert N==241920 and np.all(np.einsum('ij,ij->i',D,D@Q)==-1584)
    lookup={r.tobytes():i for i,r in enumerate(D)};assert len(lookup)==N
    permutations=[]
    for g in gens:
        assert np.array_equal(g.T@Q@g,Q) and np.array_equal(g@h,h)
        rows=D@g.T;perm=np.array([lookup[r.tobytes()]for r in rows],dtype=np.int32)
        assert len(np.unique(perm))==N;permutations.append(perm)
    parents=np.full(N,-1,dtype=np.int32);edges=np.full(N,-1,dtype=np.int8);depth=np.full(N,-1,dtype=np.int16)
    depth[seed]=0;queue=deque([seed])
    while queue:
        i=queue.popleft()
        for j,perm in enumerate(permutations):
            k=int(perm[i])
            if depth[k]<0:depth[k]=depth[i]+1;parents[k]=i;edges[k]=j;queue.append(k)
    assert np.all(depth>=0)
    selected=[0,1,2,3,4,5,7,9,11]
    bases=[[F(x)for x in proposal['bases'][i]]for i in selected]
    den=lcm(*(x.denominator for r in bases for x in r));assert den==32760
    B=np.empty((N,9,22),dtype=np.int64);B[seed]=[[int(den*x)for x in r[:22]]for r in bases]
    heights=np.array([int(den*r[22])for r in bases],dtype=np.int64)
    for d in range(1,int(max(depth))+1):
        for j,g in enumerate(gens):
            ids=np.flatnonzero((depth==d)&(edges==j));B[ids]=B[parents[ids]]@g.T
        print('depth',d,'reached',int(sum(depth<=d)),flush=True)
    assert 22*int(abs(B).max())*max(int(abs(g).max())for g in gens)<2**63
    for j,g in enumerate(gens):
        assert np.array_equal(B@g.T,B[permutations[j]])
    q=[F(proposal['base_norms'][i])for i in selected]
    for j in range(9):
        cov=B[:,j]@Q
        assert 22*int(abs(B[:,j]).max())*int(abs(cov).max())+11*int(heights[j])**2<2**63
        norm=np.einsum('ij,ij->i',B[:,j],cov)+11*heights[j]**2
        assert np.all(norm==int(q[j]*den*den))
        assert np.all(np.einsum('ij,ij->i',D,cov)+132*heights[j]==0)
    target=P/('data/enumerations/'+PRE+'_frames.npz')
    np.savez_compressed(target,old_base_numerators=B,height_numerators=heights,six_delta_old=D,parents=parents,generator_indices=edges,depth=depth,permutations=np.array(permutations),selected_base_indices=np.array(selected),denominator=den)
    out={'status':'PASS_EXACT_EQUIVARIANT_THIRTEEN_WALL_CUSP_FRAMES','positive_cusps':N,'signed_cusps':2*N,'nine_stored_bases':selected,'thirteen_bases':'first5 plus both signs of last4','base_denominator':den,'base_norms':list(map(str,q)),'height_numerators':heights.tolist(),'seed':seed,'maximum_tree_depth':int(max(depth)),'frames_sha256':hashlib.sha256(target.read_bytes()).hexdigest(),'scope':'Exact equivariant local frames; no global compatibility or volume conclusion yet.'}
    (P/('data/enumerations/'+PRE+'_frames.json')).write_text(json.dumps(out,indent=2)+'\n');print(out['status'],flush=True)
if __name__=='__main__':run()
