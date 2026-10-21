#!/usr/bin/env python3
"""Reproduce three bounded numerical diagnostics of the remote split-seed gap."""
from collections import Counter
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.spatial import HalfspaceIntersection

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'data/enumerations/fresh_h9_split_endpoint_walk_optional_anchor_20261013_scout.json'
FRONTIER=ROOT/'results/fresh_h9_split_edge_frontier_20261013.json'
OUT=ROOT/'data/enumerations/fresh_h9_split_remote_diagnostics_20261013_scout.json'


def bits(mask):
    while mask:
        b=mask&-mask;yield b.bit_length()-1;mask-=b


def run():
    source=json.loads(SOURCE.read_text());snapshot=source['source_snapshot']
    W=np.array(snapshot['physical_normals']);metric=np.diag([-1]+[1]*9);G=W@metric@W.T
    hull=HalfspaceIntersection(np.c_[W[:,1:],-W[:,0]],np.zeros(9))
    V=hull.intersections;q=np.sum(V*V,axis=1);outside=np.flatnonzero(q>1+1e-7)
    ids=[i for i,c in enumerate(source['candidates'])if not c['conflicts']]
    R=np.array([source['candidates'][i]['normal_physical']for i in ids]);residual=np.full(len(V),-np.inf)
    for i in range(0,len(R),100):
        residual=np.maximum(residual,np.max(R[i:i+100,1:]@V.T-R[i:i+100,0,None],axis=0))
    uncut=[int(i)for i in outside if residual[i]<1e-7]
    assert uncut
    chosen=max(uncut,key=lambda i:q[i]);active=hull.dual_facets[chosen]
    escaped={'old_outside_vertices':len(outside),'surviving_vertices':len(uncut),
             'point':V[chosen].tolist(),'squared_radius':float(q[chosen]),
             'maximum_cap_residual':float(residual[chosen]),'active':active}

    adjacency=[sum(1<<int(i)for i in np.flatnonzero((row>-2+1e-7)&(row<1e-7)))for row in G]
    counts=Counter();blocks=[]
    def grow(inds,inverse,mask):
        counts[len(inds)]+=1
        if len(inds)==8:blocks.append(inds);return
        for j in bits(mask):
            p=G[inds,j];ip=inverse@p;schur=2-p@ip
            if schur<=1e-7:continue
            expanded=np.block([[inverse+np.outer(ip,ip)/schur,-ip[:,None]/schur],
                               [-ip[None,:]/schur,np.array([[1/schur]])]])
            grow(inds+[j],expanded,mask&adjacency[j]&~((1<<(j+1))-1))
    grow([],np.zeros((0,0)),(1<<len(W))-1)
    old_frontier=json.loads(FRONTIER.read_text())['cases'][0]
    reached={tuple(e['edge_walls'])for e in old_frontier['edges']}
    unreached=[b for b in blocks if tuple(b)not in reached]

    # Nine old walls at this affine spacelike vertex need not form a
    # spherical block. Test finite coupling profiles directly instead.
    A=W[active]@metric;_,singular,vh=np.linalg.svd(A,full_matrices=True)
    assert min(singular)>1e-8
    n=vh[-1];pinv=np.linalg.pinv(A);profiles=[]
    for count in range(5):
        for support in itertools.combinations(range(9),count):
            for labels in itertools.product([3,4,5],repeat=count):
                row=np.zeros(9)
                for i,m in zip(support,labels):row[i]=-2*np.cos(np.pi/m)
                profiles.append(row)
    initial=np.array(profiles)@pinv.T;aa=n@metric@n;bb=2*initial@metric@n
    cc=np.einsum('ij,jk,ik->i',initial,metric,initial)-2;dd=bb*bb-4*aa*cc
    assert abs(aa)>1e-8
    finite=-2*np.cos(np.pi/np.arange(2,121));seed_pass=[]
    for sign in(-1,1):
        rows=np.flatnonzero(dd>=-1e-8)
        t=(-bb[rows]+sign*np.sqrt(np.maximum(0,dd[rows])))/(2*aa)
        normals=initial[rows]+t[:,None]*n;pair=normals@metric@W[:16].T
        ok=np.all((pair<-2-1e-7)|np.any(abs(pair[:,:,None]-finite)<1e-7,axis=2),axis=1)&(normals[:,0]>0)
        seed_pass.append(int(sum(ok)))
    result={'status':'NUMERICAL_REMOTE_GAP_DIAGNOSTICS',
            'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(SOURCE,FRONTIER)},
            'escape':escaped,'elliptic_subset_counts':dict(counts),'rank_eight_blocks':blocks,
            'unreached_rank_eight_blocks':unreached,
            'remote_nine_wall_profiles':{'profile_count':len(profiles),'seed_pass_by_branch':seed_pass},
            'scope':'Floating geometry and bounded profiles only. The exact common-escape certificate is separate. Absence of candidates does not exclude other angle profiles, parameters, supports, or new generations.'}
    OUT.write_text(json.dumps(result,indent=2)+'\n')
    print({'outside_vertices':len(outside),'survive_all_caps':len(uncut),'rank8':len(blocks),
           'unreached_rank8':len(unreached),'remote_profiles':len(profiles),'seed_pass':seed_pass},flush=True)


if __name__=='__main__':run()
