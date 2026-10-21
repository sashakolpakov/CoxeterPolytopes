#!/usr/bin/env python3
"""Bounded simultaneous parameter recovery for all-acute five-plus-cap covers."""
import argparse,hashlib,itertools,json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from fresh_h20_two_anchor_families_20261016 import point
from fresh_h20_quad_adjacent_h3_closure_20261016 import geometry158,face_data,poly,gate_mask
from fresh_h20_release144_cap_20261016_residual_scout import close
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_two_anchor_families_20261016_coordinated_angles'
BASE='results/fresh_h20_two_anchor_families_20261016.json';COV='results/fresh_h20_two_anchor_families_20261016_acute_cover.json'
def run():
    D=json.loads((ROOT/BASE).read_text());C=json.loads((ROOT/COV).read_text());G=geometry158();M=G['metric'];base,V,A,b,ids=face_data(G)
    Z=np.array([x['z']for x in D['joint_samples']]);records=[];success=[]
    lab=[2,3,4,5,6,7,8,9,10,12,16,0];cs=np.array([0 if m==2 else 1 if m==0 else np.cos(np.pi/m)for m in lab])
    for row in C['steps'][:4]:
        if 'roots'not in row:continue
        ff=[D['families'][r['family']]for r in row['roots']];t0=np.array([r['parameter']for r in row['roots']]);lo=[];hi=[]
        for f,t in zip(ff,t0):
            a,z=next((a,z)for a,z in f['intervals']if(a is None or t>a)and(z is None or t<z))
            lo.append(-np.inf if a is None else a+1e-9*max(1,abs(a)));hi.append(np.inf if z is None else z-1e-9*max(1,abs(z)))
        def roots(t):return np.array([point(f,v)for f,v in zip(ff,t)])
        R0=roots(t0);co0=-R0@M@R0.T/2
        bad=[tuple(p['positions'])for p in row['noncoxeter_new_pairs']]
        options=[list(np.argsort(abs(cs-co0[i,j]))[:2])for i,j in bad]
        finite=[(i,j)for i in range(len(ff))for j in range(i+1,len(ff))if co0[i,j]<1-1e-7]
        dots=[(i,j)for i in range(len(ff))for j in range(i+1,len(ff))if co0[i,j]>=1-1e-7]
        for choices in itertools.product(*options):
            target={pair:cs[c]for pair,c in zip(bad,choices)}
            for pair in finite:
                if pair not in target:target[pair]=cs[np.argmin(abs(cs-co0[pair]))]
            def residual(t):
                R=roots(t);co=-R@M@R.T/2;ar=-R[:,1:]/R[:,0,None]
                return np.r_[[co[p]-v for p,v in target.items()],
                    [max(1+1e-7-co[p],0)for p in dots],
                    100*np.maximum(1e-7-np.max(Z@ar.T-1,axis=1),0)]
            fit=least_squares(residual,t0,bounds=(lo,hi),max_nfev=1200,ftol=1e-12,xtol=1e-12,gtol=1e-12,x_scale='jac')
            err=float(max(abs(residual(fit.x))));rec=dict(cover_iteration=row['iteration'],bad_pairs=[list(x)for x in bad],target_labels=[lab[c]for c in choices],maximum_residual=err,nfev=fit.nfev,parameters=fit.x.tolist())
            records.append(rec)
            if err>1e-7:continue
            R=roots(fit.x);co=-R@M@R.T/2;np.fill_diagonal(co,0)
            if not gate_mask(co).all():continue
            if not gate_mask(-(R@M@G['X'].T)/np.sqrt(2*G['norms'])).all():continue
            if max(abs(np.einsum('ij,jk,ik->i',R,M,R)-2))>1e-4:continue
            ar=-R[:,1:]/R[:,0,None];face=poly(base,V,np.r_[A,ar@V],np.r_[b,1-ar@base],ids+list(range(158,158+len(R))))
            success.append(dict(**rec,family_ids=[f['id']for f in ff],roots=R.tolist(),constraints=ar.tolist(),mutual_cosines=co.tolist(),face=face))
        print('COVER',row['iteration'],'patterns',len(list(itertools.product(*options))),'min_error',min(r['maximum_residual']for r in records if r['cover_iteration']==row['iteration']),'success',len(success),flush=True)
    src=[f'scripts/{TAG}.py','scripts/fresh_h20_two_anchor_families_20261016.py',BASE,COV]
    return dict(status='PASS_BOUNDED_SIMULTANEOUS_CLUSTER_ANGLE_RECOVERY',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in src},mutual_label_domain=lab,records=records,successes=success,
       scope='Numerical least-squares from4 selected all-acute covers, two nearest target angles per initially non-Coxeter contact, max1200evaluations each. Failure is local optimizer evidence only; success still requires exact normalization/angle/corner and remaining-face certification.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
    if a.verify:assert close(d,json.loads(p.read_text()))
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],len(d['records']),len(d['successes']),flush=True)
