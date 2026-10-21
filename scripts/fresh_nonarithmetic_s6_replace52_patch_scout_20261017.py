#!/usr/bin/env python3
"""Change the old30--52 angle while solving a three-cap S6 patch.

All other62 frozen normals are retained. The new52 is right to its former
eight-wall support, with a changed finite label to30. Candidate caps C,A,B
are reconstructed from small literal Grams; only discovery arithmetic is used.
"""
import argparse
import json
import math
from collections import Counter
from hashlib import sha256
from pathlib import Path

import numpy as np
import fresh_nonarithmetic_s6_split_bridge_scout_20261017 as F

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_nonarithmetic_s6_replace52_patch_scout_20261017'
S=F.S;M=F.METRIC
OLD_LABELS=list(range(4,33))+[48,64,128]
C44_LABELS=list(range(2,65))+[96,128,192,256]
C52_LABELS=list(range(2,17))+[20,24,32,48,64]


def theta(m):return 0. if m==2 else 2*math.cos(math.pi/m)


def run():
    d=F.read(F.TAG);original=np.array(d['old_normals']);h=np.array([1.]+[0.]*9)
    ids=S+[8,34,30];frame=original[ids];gram=frame@M@frame.T
    inverse=np.linalg.inv(gram);kernel=h-(h@M@frame.T)@inverse@frame
    kernel/=math.sqrt(kernel@M@kernel)
    contacts=np.array([0.]*8+[-theta(4)]);part=contacts@inverse@frame
    sign=min([-1,1],key=lambda s:np.linalg.norm(part+s*math.sqrt(2-part@M@part)*kernel-original[52]))
    counts=Counter();stages=[];near=[];hits=[]
    grid=np.array([theta(m)for m in range(2,1025)])
    other52=[j for j in range(63)if j not in ids+[52]]
    otherC=[j for j in range(63)if j not in S+[8,44,52]]
    otherA=[j for j in range(63)if j not in S+[30,44]]
    otherB=[j for j in range(63)if j not in S+[30,52]]
    for oldlabel in OLD_LABELS:
        N=original.copy();contacts=np.array([0.]*8+[-theta(oldlabel)])
        part=contacts@inverse@frame;rho=2-part@M@part
        N[52]=part+sign*math.sqrt(rho)*kernel
        if max(N[52]@M@N[other52].T)>-2+1e-7:
            counts['replacement52_old_contact_rejected']+=1
            continue
        f=N[[8,44,52]]
        x=float(N[44]@M@N[52]);b=float(N[8]@M@N[30]);t=theta(oldlabel)
        g3=np.array([[2.,0,0],[0,2,x],[0,x,2]])
        j3=np.linalg.inv(g3);g30=np.array([b,0,-t]);w30=N[30]-g30@j3@f
        w2=2-g30@j3@g30
        assert w2>0
        local=Counter()
        for a_label in C44_LABELS:
            for b_label in C52_LABELS:
                ac,bc=theta(a_label),theta(b_label)
                cc=np.array([0.,-ac,-bc]);p=cc@j3@f;residual=2-cc@j3@cc
                if residual<=0:continue
                for branch in [-1,1]:
                    c=p+branch*math.sqrt(residual/w2)*w30
                    c30=float(cc@j3@g30+branch*math.sqrt(residual*w2))
                    if c30>-2+1e-7 or max(c@M@N[otherC].T)>-2+1e-7:
                        local['C_unprescribed_contact_rejected']+=1
                        continue
                    g4=np.array([[2.,0,-t,c30],[0,2.,x,-ac],[-t,x,2.,-bc],[c30,-ac,-bc,2.]])
                    j4=np.linalg.inv(g4)
                    if j4[1,1]<=1e-10 or j4[2,2]<=1e-10:
                        local['polar_normal_not_spacelike']+=1
                        continue
                    full=np.vstack([N[[30,44,52]],c])
                    ca=-math.sqrt(2/j4[2,2])*(j4[2]@full)
                    cb=-math.sqrt(2/j4[1,1])*(j4[1]@full)
                    if max(ca@M@N[otherA].T)>-2+2e-6 or max(cb@M@N[otherB].T)>-2+2e-6:
                        local['polar_old_contact_rejected']+=1
                        continue
                    value=2*j4[1,2]/math.sqrt(j4[1,1]*j4[2,2])
                    if not -2+1e-8<value<1e-8:
                        local['AB_not_finite']+=1
                        continue
                    local['finite_patch_candidates']+=1
                    m=int(np.argmin(abs(grid+value)))+2;error=float(abs(value+theta(m)))
                    record=dict(replacement52_label=oldlabel,C44_label=a_label,C52_label=b_label,
                                C_branch=branch,AB_product=float(value),nearest_AB_label=m,angle_residual=error,
                                new52=N[52].tolist(),normals=[ca.tolist(),cb.tolist(),c.tolist()])
                    near.append(record);near.sort(key=lambda r:r['angle_residual']);del near[50:]
                    if error<1e-9:
                        local['possible_angle_hits']+=1;hits.append(record)
        counts.update(local);stages.append(dict(replacement52_label=oldlabel,counts=dict(local)))
        print('REPLACED52_PATCH',oldlabel,dict(local),'best residual',near[0]['angle_residual']if near else None,flush=True)
    paths=[f'scripts/{TAG}.py',f'scripts/{F.TAG}.py',f'results/{F.TAG}.json']
    return dict(status='COMPLETED_NUMERICAL_REPLACED52_THREE_CAP_PATCH_SCOUT',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
                retained_frozen_wall_count=62,replaced_wall=52,
                replacement_labels=OLD_LABELS,C44_labels=C44_LABELS,C52_labels=C52_LABELS,
                mutual_AB_recognition_through=1024,counts=dict(counts),stages=stages,hits=hits,closest_candidates=near,
                scope='Bounded numerical exploration that changes one earlier wall. All unprescribed contacts with the62 other fixed walls are restricted to be dotted. A/B are right polars of C plus their separate old endpoint frames. Numerical angle coincidences need exact certificates; no finite-volume or all-label exclusion follows.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['counts'],'hits',len(result['hits']),flush=True)
