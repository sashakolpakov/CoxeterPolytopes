#!/usr/bin/env python3
"""Bounded mutual-angle solving between the continuous one-parameter families."""
import argparse,hashlib,json,itertools
from collections import defaultdict,Counter
from pathlib import Path
import numpy as np
from fresh_h20_two_anchor_families_20261016 import point,real_roots
from fresh_h20_quad_adjacent_h3_closure_20261016 import geometry158,face_data,poly,gate_mask
from fresh_h20_release144_cap_20261016_residual_scout import close
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_two_anchor_families_20261016_pairs'
BASE='results/fresh_h20_two_anchor_families_20261016.json'
def mask(ids):return sum(1<<int(j)for j in ids)
def inrange(t,f):
    return any((a is None or t>a+1e-8*max(1,abs(a)))and(b is None or t<b-1e-8*max(1,abs(b)))for a,b in f['intervals'])
def run(limit):
    d=json.loads((ROOT/BASE).read_text());G=geometry158();M=G['metric'];X=G['X'];q=G['norms']
    H=(G['H']@G['g']@G['F'].T)@M;base,V,A,b,ids=face_data(G)
    Z=np.array([p['z']for p in d['joint_samples']]);families=d['families'];samples=d['candidates']
    groups=defaultdict(list)
    for p in samples:groups[mask(p['joint_sample_cuts'])].append(p)
    # Two geometrically separated representatives per coverage class, then
    # a declared finite rank by sample coverage. These choices are discovery.
    anchors=[]
    for m,pp in groups.items():
        pp.sort(key=lambda p:p['orthonormal_chart'][0]);anchors.extend(pp[:1])
        if len(pp)>1:anchors.extend(pp[-1:])
    anchors.sort(key=lambda p:(len(p['joint_sample_cuts']),-p['orthonormal_chart'][0]),reverse=True)
    anchors=anchors[:limit];best={};stats=Counter();num=np.array([f['numerator']for f in families]);den=np.array([f['denominator']for f in families])
    for no,a in enumerate(anchors):
        r=np.array(a['orthonormal_chart']);am=mask(a['joint_sample_cuts']);products=np.einsum('fij,j->fi',num,M@r)
        for label in [2,3,4,5,6,8,10]:
            c=0 if label==2 else np.cos(np.pi/label);P=products+2*c*den
            for f,coef in zip(families,P):
                for t in real_roots(coef):
                    stats['quadratic_roots']+=1
                    if not inrange(t,f):continue
                    rr=point(f,t)
                    if abs(rr@M@rr-2)>1e-4 or rr[0]<=0 or rr@M@H>=0:continue
                    aa=-rr[1:]/rr[0];cuts=np.where(Z@aa>1+1e-7)[0];bm=mask(cuts);union=am|bm
                    if union.bit_count()<15:continue
                    if not gate_mask((-(rr@M@X.T)/np.sqrt(2*q))[None,:])[0]:continue
                    if abs(r@M@rr+2*c)>1e-5:continue
                    stats['accepted_strong_pairs']+=1
                    # Keep a bounded list of best geometric representatives
                    # for each literal union of old diagnostic points.
                    margin=float(np.min(np.maximum(Z@np.array(a['residual_constraint'])-1,Z@aa-1)))
                    key=(union,label)
                    item=dict(anchor=a['id'],family=f['id'],parameter=float(t),mutual_label=label,
                      partner_orthonormal_chart=rr.tolist(),partner_constraint=aa.tolist(),
                      joint_sample_cuts=cuts.tolist(),union_cuts=[j for j in range(len(Z))if union>>j&1],
                      minimum_cover_margin=margin)
                    if key not in best or margin>best[key]['minimum_cover_margin']:best[key]=item
        if no%20==0:print('ANCHOR',no,len(anchors),'best',len(best),dict(stats),flush=True)
    pairs=sorted(best.values(),key=lambda p:(len(p['union_cuts']),p['minimum_cover_margin']),reverse=True)
    results=[]
    for p in pairs[:64]:
        ar=np.array(samples[p['anchor']]['residual_constraint']);br=np.array(p['partner_constraint'])
        face=poly(base,V,np.r_[A,[ar@V,br@V]],np.r_[b,1-ar@base,1-br@base],ids+[158,159])
        results.append(dict(pair=p,face=face))
    source=[f'scripts/{TAG}.py','scripts/fresh_h20_two_anchor_families_20261016.py',BASE]
    return dict(status='PASS_BOUNDED_CONTINUOUS_FAMILY_MUTUAL_ANGLE_SOLVER',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in source},
       anchor_limit=limit,anchor_ids=[p['id']for p in anchors],mutual_labels=[2,3,4,5,6,8,10],stats=dict(stats),
       strong_pairs=pairs,full_face_checks=results,scope='Numerical bounded choice of sampled first cap, all saved second-cap parameter intervals, specified mutual-angle labels. All158 old gates and solved mutual gates rechecked. No exact cap or complete cluster/nonexistence claim.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--limit',type=int,default=160);ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run(a.limit);p=ROOT/f'results/{TAG}.json'
    if a.verify:assert close(d,json.loads(p.read_text()))
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],len(d['strong_pairs']),min((p['face'].get('maximum_squared_radius',float('inf'))for p in d['full_face_checks']),default=None),flush=True)
