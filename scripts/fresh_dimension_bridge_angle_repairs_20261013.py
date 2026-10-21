#!/usr/bin/env python3
"""Numerical non-right repairs of rejected H20 polar truncations."""
import json
import hashlib
from collections import Counter
from pathlib import Path
import numpy as np
import fresh_dimension_bridge_20261013 as K

ROOT=Path(__file__).resolve().parents[1]
SOURCE='data/enumerations/fresh_dimension_bridge_polar_20261013.json'
OUT=ROOT/'data/enumerations/fresh_dimension_bridge_angle_repairs_20261013_scout.json'


def run():
    data=json.loads((ROOT/SOURCE).read_text());base=json.loads((ROOT/data['source']).read_text())
    vec=lambda r:np.array([K.value(K.dec(x))for x in r])
    Q=np.array([vec(r)for r in base['form']]);R=np.array([vec(r)for r in base['normals']]+[vec(r['vector'])for r in data['cuts']])
    norms=np.einsum('ij,jk,ik->i',R,Q,R);h=vec(base['strict_center']);cov=R@Q
    cosines=np.cos(np.pi/np.arange(2,121));counts=Counter();candidates=[];tested=[]
    for index,target in enumerate(data['rejected_polar_proposals']):
        v=vec(target['vector'])
        if max(cov@v)>1e-6:continue
        active=target['active_indices'];tested.append(index)
        for released in active:
            kept=[i for i in active if i!=released]
            for anchor in target['forbidden_indices']:
                ids=kept+[anchor];A=cov[ids];_,s,vh=np.linalg.svd(A,full_matrices=True)
                if min(s)<1e-8:counts['rank_deficient']+=1;continue
                d=vh[-1];e=np.linalg.pinv(A)[:,-1];qa=d@Q@d
                for m in (2,3,4,5,6,8,10,12):
                    r0=e*(-np.sqrt(2*norms[anchor])*np.cos(np.pi/m))
                    qb=2*r0@Q@d;qc=r0@Q@r0-2;disc=qb*qb-4*qa*qc
                    counts['quadratics']+=1
                    if abs(qa)<1e-10:ts=[]if abs(qb)<1e-10 else[-qc/qb]
                    elif disc< -1e-9:continue
                    else:ts=[(-qb+sgn*np.sqrt(max(0,disc)))/(2*qa)for sgn in(-1,1)]
                    for t in ts:
                        r=r0+t*d
                        if r@Q@h>=-1e-7 or r@Q@v<=1e-7:continue
                        counts['oriented_cutters']+=1
                        pp=cov@r;ratio=-pp/np.sqrt(2*norms)
                        ok=(ratio>=1-1e-7)|(np.min(abs(ratio[:,None]-cosines),axis=1)<1e-7)
                        if not all(ok):counts['bad_pair_rejected']+=1;continue
                        assert abs(r@Q@r-2)<1e-6
                        record={'target_index':index,'released':released,'anchor':anchor,'label':m,
                                'normal':r.tolist(),'target_pairing':float(r@Q@v),
                                'ratios':ratio.tolist(),'max_norm_error':abs(float(r@Q@r)-2)}
                        candidates.append(record)
        print({'target':index,'candidates':len(candidates),'counts':dict(counts)},flush=True)
    out={'status':'NUMERICAL_H20_NONRIGHT_POLAR_REPAIR_SCOUT','source':SOURCE,
         'source_sha256':hashlib.sha256((ROOT/SOURCE).read_bytes()).hexdigest(),
         'tested_targets':tested,'counts':dict(counts),'candidates':candidates,
         'scope':'For each surviving sampled rejected polar ray, retain19 of its20 active wall orthogonalities and prescribe one finite angle with one forbidden wall, labels2,3,4,5,6,8,10,12; solve the norm2 quadratic and test all82 old relations using labels2..120 or dotted pairs. Numerical only; this is not all caps or all affine vertices.'}
    OUT.write_text(json.dumps(out,indent=2)+'\n');print({'done':True,'candidates':len(candidates)},flush=True)


if __name__=='__main__':run()
