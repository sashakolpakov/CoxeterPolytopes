#!/usr/bin/env python3
"""Numerical two-angle H20 repairs at the surviving one-angle bank gaps."""
import argparse
import itertools
import json
from pathlib import Path
from collections import Counter
import numpy as np
import fresh_dimension_bridge_20261013 as K

ROOT=Path(__file__).resolve().parents[1]
SOURCE='data/enumerations/fresh_dimension_bridge_angle_repairs_20261013_scout.json'
OUT=ROOT/'data/enumerations/fresh_dimension_bridge_two_angle_20261013_scout.json'


def run(labels):
    old=json.loads((ROOT/SOURCE).read_text());data=json.loads((ROOT/old['source']).read_text())
    base=json.loads((ROOT/data['source']).read_text());ev=lambda r:np.array([K.value(K.dec(x))for x in r])
    Q=np.array([ev(r)for r in base['form']]);R=np.array([ev(r)for r in base['normals']]+[ev(r['vector'])for r in data['cuts']])
    h=ev(base['strict_center']);norms=np.einsum('ij,jk,ik->i',R,Q,R);cov=R@Q
    V=np.array([ev(r['vector'])for r in data['rejected_polar_proposals']])
    caps=np.array([r['normal']for r in old['candidates']]);survive=np.flatnonzero(np.max(caps@Q@V.T,axis=0)<1e-7)
    pairlabels=list(itertools.product(labels,repeat=2));cs=np.cos(np.pi/np.array(pairlabels))
    finite=np.cos(np.pi/np.arange(2,121));found=[];counts=Counter();tested=[]
    for ti in survive:
        target=data['rejected_polar_proposals'][ti]
        if len(target['forbidden_indices'])!=1:continue
        bad=target['forbidden_indices'][0];v=V[ti];tested.append(int(ti));startcount=len(found)
        for released in itertools.combinations(target['active_indices'],2):
            kept=[i for i in target['active_indices']if i not in released]
            for anchor in released:
                ids=kept+[bad,anchor];A=cov[ids];_,s,vh=np.linalg.svd(A,full_matrices=True)
                if min(s)<1e-8:counts['rank_deficient']+=1;continue
                direction=vh[-1];pinv=np.linalg.pinv(A)[:,-2:]
                rhs=-cs*np.sqrt(2*norms[[bad,anchor]])
                initial=rhs@pinv.T;aa=direction@Q@direction;bb=2*initial@Q@direction
                cc=np.einsum('ij,jk,ik->i',initial,Q,initial)-2;disc=bb*bb-4*aa*cc
                counts['quadratics']+=len(pairlabels)
                if abs(aa)<1e-10:counts['null_direction_skipped']+=1;continue
                for sign in(-1,1):
                    rows=np.flatnonzero(disc>=-1e-8)
                    t=(-bb[rows]+sign*np.sqrt(np.maximum(0,disc[rows])))/(2*aa)
                    normals=initial[rows]+t[:,None]*direction
                    positive=(normals@Q@h< -1e-7)&(normals@Q@v>1e-7)
                    rows=rows[positive];normals=normals[positive]
                    if not len(rows):continue
                    counts['oriented_cutters']+=len(rows)
                    ratio=-(normals@cov.T)/np.sqrt(2*norms)
                    ix=np.clip(np.searchsorted(finite,ratio),1,len(finite)-1)
                    err=np.minimum(abs(ratio-finite[ix]),abs(ratio-finite[ix-1]))
                    passing=np.all((ratio>=1-1e-7)|(err<1e-7),axis=1)
                    for ri in np.flatnonzero(passing):
                        r=normals[ri];assert abs(r@Q@r-2)<1e-6
                        found.append({'target_index':int(ti),'released':list(released),'anchors':[bad,anchor],
                            'labels':list(pairlabels[rows[ri]]),'normal':r.tolist(),'target_pairing':float(r@Q@v),
                            'ratios':ratio[ri].tolist()})
        print({'target':int(ti),'new_candidates':len(found)-startcount,'total':len(found),'counts':dict(counts)},flush=True)
        OUT.write_text(json.dumps({'status':'NUMERICAL_TWO_ANGLE_H20_CAP_SCOUT','source':SOURCE,
             'base_source':old['source'],'labels':labels,'tested_targets':tested,'counts':dict(counts),'candidates':found,
             'scope':'For single-forbidden-pair polar gaps surviving all one-angle caps, retain18 active orthogonalities, impose finite pairings with the forbidden wall and one released old wall, and solve norm2. Listed labels, floating geometry, all82 old pairing checks. Null quadratic directions skipped. No exhaustive cap or volume claim.'},indent=2)+'\n')
    print({'done':True,'candidates':len(found)},flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--labels',nargs='+',type=int,default=[3,4,6]);a=p.parse_args();run(a.labels)
