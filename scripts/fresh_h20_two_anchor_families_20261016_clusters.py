#!/usr/bin/env python3
"""Continuous third-cap interval intersections for a bounded compatible pair bank."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from fresh_h20_two_anchor_families_20261016 import point,real_roots,middle,samples
from fresh_h20_quad_adjacent_h3_closure_20261016 import geometry158,face_data,poly,gate_mask
from fresh_h20_release144_cap_20261016_residual_scout import close
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_two_anchor_families_20261016_clusters'
BASE='results/fresh_h20_two_anchor_families_20261016.json'
PAIRS='results/fresh_h20_two_anchor_families_20261016_pairs.json'
def intersection(old,P):
    roots=sorted(set(r for c in P for r in real_roots(c)if np.isfinite(r)))
    out=[];scale=np.maximum(np.max(abs(P),axis=1),1e-30)
    for lo,hi in old:
        cuts=[lo]+[t for t in roots if(lo is None or t>lo)and(hi is None or t<hi)]+[hi]
        for a,b in zip(cuts,cuts[1:]):
            t=middle(a,b)
            if min((P@np.array([1.,t,t*t]))/scale)>-1e-10:out.append([a,b])
    return out
def run(limit):
    D=json.loads((ROOT/BASE).read_text());B=json.loads((ROOT/PAIRS).read_text());G=geometry158();M=G['metric'];X=G['X'];q=G['norms']
    base,V,A,b,ids=face_data(G);Z=np.array([p['z']for p in D['joint_samples']]);out=[];intervals=[]
    pairs=sorted(B['full_face_checks'],key=lambda p:p['face']['maximum_squared_radius'])[:limit]
    for pi,item in enumerate(pairs):
        p=item['pair'];first=D['candidates'][p['anchor']];R=np.array([first['orthonormal_chart'],p['partner_orthonormal_chart']]);AR=-R[:,1:]/R[:,0,None]
        missing=np.where(np.max(Z@AR.T-1,axis=1)<=1e-7)[0];local=[]
        for f in D['families']:
            N=np.array(f['numerator']);den=np.array(f['denominator'])
            extra=np.r_[(-(N@M@R.T)).T-2*den[None,:],-(N@np.c_[np.ones(len(missing)),Z[missing]].T).T-1e-7*N[:,0]]
            allowed=intersection(f['intervals'],extra)
            if not allowed:continue
            intervals.append(dict(pair_index=pi,family=f['id'],intervals=allowed))
            # Parameter intersection is continuous; the selected representative
            # and all geometric scoring remain bounded numerical discovery.
            for ab in allowed:
                for t in samples(ab)[:3]:
                    r=point(f,t)
                    if abs(r@M@r-2)>1e-4 or r[0]<=0:continue
                    if not gate_mask((-(r@M@X.T)/np.sqrt(2*q))[None,:])[0]:continue
                    if min(-R@M@r/2)<1-1e-7:continue
                    ar=-r[1:]/r[0];margin=float(min(np.max(np.c_[Z@AR.T,Z@ar]-1,axis=1)))
                    if margin<=1e-8:continue
                    local.append(dict(family=f['id'],parameter=float(t),orthonormal_chart=r.tolist(),constraint=ar.tolist(),joint_margin=margin))
        local.sort(key=lambda z:z['joint_margin'],reverse=True)
        for c in local[:24]:
            AR3=np.r_[AR,[c['constraint']]]
            face=poly(base,V,np.r_[A,AR3@V],np.r_[b,np.ones(3)-AR3@base],ids+[158,159,160])
            out.append(dict(pair=p,third=c,face=face))
        print('PAIR',pi,'families',sum(x['pair_index']==pi for x in intervals),'samples',len(local),'census',min((x['face']['maximum_squared_radius']for x in out if x['pair']==p),default=None),flush=True)
    sources=[f'scripts/{TAG}.py','scripts/fresh_h20_two_anchor_families_20261016.py',BASE,PAIRS]
    return dict(status='PASS_BOUNDED_CONTINUOUS_THIRD_CAP_COVER_INTERSECTIONS',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in sources},pair_limit=limit,
       selected_pairs=[p['pair']for p in pairs],feasible_third_family_intervals=intervals,full_face_checks=out,
       scope='Bounded selected first pairs; every saved third family interval intersected numerically with dotted contacts to both new caps and cuts of all remaining21 saved samples. Representative and complete-face scoring bounded; no all-cap exclusion or exact finite-volume claim.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--limit',type=int,default=32);ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run(a.limit);p=ROOT/f'results/{TAG}.json'
    if a.verify:assert close(d,json.loads(p.read_text()))
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],len(d['feasible_third_family_intervals']),len(d['full_face_checks']),min((r['face']['maximum_squared_radius']for r in d['full_face_checks']),default=None),flush=True)
