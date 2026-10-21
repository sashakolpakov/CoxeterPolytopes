#!/usr/bin/env python3
"""Bounded third-cap search with one solved finite new-to-new contact."""
import argparse,hashlib,json
from collections import Counter
from pathlib import Path
import numpy as np
from fresh_h20_two_anchor_families_20261016 import point,real_roots
from fresh_h20_two_anchor_families_20261016_pairs import inrange
from fresh_h20_quad_adjacent_h3_closure_20261016 import geometry158,face_data,poly,gate_mask
from fresh_h20_release144_cap_20261016_residual_scout import close
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_two_anchor_families_20261016_finite_third'
BASE='results/fresh_h20_two_anchor_families_20261016.json';PAIRS='results/fresh_h20_two_anchor_families_20261016_pairs.json'
def run():
    D=json.loads((ROOT/BASE).read_text());B=json.loads((ROOT/PAIRS).read_text());G=geometry158();M=G['metric'];X=G['X'];q=G['norms']
    base,V,A,b,ids=face_data(G);Z=np.array([p['z']for p in D['joint_samples']]);out=[];stats=Counter()
    pairs=sorted(B['full_face_checks'],key=lambda p:p['face']['maximum_squared_radius'])[:32]
    num=np.array([f['numerator']for f in D['families']]);den=np.array([f['denominator']for f in D['families']])
    for pi,item in enumerate(pairs):
        p=item['pair'];first=D['candidates'][p['anchor']];R=np.array([first['orthonormal_chart'],p['partner_orthonormal_chart']]);AR=-R[:,1:]/R[:,0,None]
        missing=np.where(np.max(Z@AR.T-1,axis=1)<=1e-7)[0];seen=set();local=[]
        for partner in [0,1]:
            products=np.einsum('fij,j->fi',num,M@R[partner])
            for label in [2,3,4,5,6,8,10]:
                c=0 if label==2 else np.cos(np.pi/label)
                for f,pol in zip(D['families'],products+2*c*den):
                    for t in real_roots(pol):
                        stats['quadratic_roots']+=1
                        if not inrange(t,f):continue
                        r=point(f,t)
                        if abs(r@M@r-2)>1e-4 or r[0]<=0:continue
                        co=-R@M@r/2
                        if not gate_mask(co[None,:])[0]:continue
                        stats['new_pair_gates_pass']+=1
                        ar=-r[1:]/r[0]
                        if min(Z[missing]@ar-1)<=1e-7:continue
                        if not gate_mask((-(r@M@X.T)/np.sqrt(2*q))[None,:])[0]:continue
                        key=tuple(np.round(r,7))
                        if key in seen:continue
                        seen.add(key);margin=float(min(np.max(np.c_[Z@AR.T,Z@ar]-1,axis=1)))
                        local.append(dict(family=f['id'],parameter=float(t),finite_partner=partner,finite_label=label,
                            mutual_cosines=co.tolist(),orthonormal_chart=r.tolist(),constraint=ar.tolist(),joint_margin=margin))
        local.sort(key=lambda p:p['joint_margin'],reverse=True)
        for third in local[:24]:
            aa=np.r_[AR,[third['constraint']]];face=poly(base,V,np.r_[A,aa@V],np.r_[b,np.ones(3)-aa@base],ids+[158,159,160])
            out.append(dict(pair=p,third=third,face=face))
        print('PAIR',pi,'covering_thirds',len(local),'total',len(out),flush=True)
    src=[f'scripts/{TAG}.py','scripts/fresh_h20_two_anchor_families_20261016.py','scripts/fresh_h20_two_anchor_families_20261016_pairs.py',BASE,PAIRS]
    return dict(status='PASS_BOUNDED_FINITE_THIRD_CAP_QUADRATIC_SOLVER',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in src},stats=dict(stats),selected_pairs=[p['pair']for p in pairs],full_face_checks=out,
      scope='Numerical32 selected pairs and all2,730 saved families, prescribed mutual labels2..6,8,10 to either earlier cap; other new contact and all158 old gates tested. Only triples covering21 saved points retained. No exact cap or universal nonexistence claim.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
    if a.verify:assert close(d,json.loads(p.read_text()))
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],len(d['full_face_checks']),flush=True)
