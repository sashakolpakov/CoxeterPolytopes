#!/usr/bin/env python3
"""Bounded dotted-contact slices of two-parameter H20 normal families.

Only 151 and one genuine old H3 side have prescribed finite angles. A
third contact is a continuous dotted cosine, sliced at five stated values.
Each conic slice is partitioned at every saved ideal-sample pairing zero.
All arithmetic here is numerical; this is discovery, not an exclusion.
"""
import argparse,hashlib,itertools,json
from collections import Counter
from pathlib import Path
import numpy as np
from scipy.linalg import null_space
from fresh_h20_two_anchor_families_20261016 import charts,allowed_intervals,constant_gate,point,real_roots,middle
from fresh_h20_quad_adjacent_h3_closure_20261016 import geometry158,face_data,poly,gate_mask
from fresh_h20_release144_cap_20261016_residual_scout import close

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_one_anchor_dotted_families_20261016'
POINTS='results/fresh_h20_continuous_sample_joint_cover_20261016.json'
SIDES=[95,142,145,150,152,156,157]
OTHERS=[31,95,132,142,145,150,152,156,157]
LABELS=[2,3,4,5,6,8,10]
SLICES=[1.01,1.1,1.5,2.,3.]

def inside(t,a,b):return (a is None or a<t) and (b is None or t<b)

def cells(f,Z):
    num=np.array(f['numerator']);den=np.array(f['denominator'])
    # Positive means that this cap cuts the chart point. The denominator
    # is positive on every chart domain used by the inherited constructor.
    P=-(num[:,0,None]+num[:,1:]@Z.T).T
    out=[]
    for interval in f['intervals']:
        a,b=interval;rr=sorted(r for c in P for r in real_roots(c) if np.isfinite(r) and inside(r,a,b))
        distinct=[]
        for r in rr:
            if not distinct or abs(r-distinct[-1])>1e-11*max(1,abs(r)):distinct.append(r)
        ends=[a,*distinct,b]
        for lo,hi in zip(ends,ends[1:]):
            t=middle(lo,hi)
            if not inside(t,lo,hi):continue
            values=P@np.array([1.,t,t*t])
            mask=np.where(values>0)[0].tolist()
            if not mask:continue
            out.append(dict(interval=[lo,hi],parameter=t,positive_sample_indices=mask))
    return out,P

def dynamic_point(f,c,t):
    """Stable two-variable reconstruction within a non-null saved chart."""
    p=np.array(f['projection_base'])+c*np.array(f['projection_slope'])
    model=next(x for x in charts(p,np.array(f['kernel']),np.diag([-1,1,1,1,1]))
        if x['kind']==f['kind'] and x.get('branch')==f.get('branch'))
    return point(model,t)

def run():
    G=geometry158();X,M,q=G['X'],G['metric'],G['norms'];H=(G['H']@G['g']@G['F'].T)@M
    J=json.loads((ROOT/POINTS).read_text());Z=np.array([s['z'] for s in J['samples']])
    base,V,A,b,ids=face_data(G);old=poly(base,V,A,b,ids)
    families=[];pool=[];frames=[];seen=set();stats=Counter()
    for anchor in SIDES:
      for other in OTHERS:
        if anchor==other:continue
        frame=[151,anchor,other];U=X[frame];B=U@M@U.T
        if min(abs(np.linalg.eigvalsh(B)))<1e-8:stats['singular_frames']+=1;continue
        inv=np.linalg.inv(B);K=null_space(U@M);start=len(families);ps=np.array([0.,0.,-np.sqrt(2*q[other])])@inv@U
        for m151,m in itertools.product([2,3],LABELS):
          cv=np.array([0 if m151==2 else np.cos(np.pi/m151),0 if m==2 else np.cos(np.pi/m),0.])
          p0=(-np.sqrt(2*q[frame])*cv)@inv@U
          for c in SLICES:
            stats['patterns']+=1;p=p0+c*ps
            for model in charts(p,K,M):
                num=np.array(model['numerator']);den=np.array(model['denominator']);di=np.argmax(abs(den))
                co=-(num@M@X.T).T/np.sqrt(2*q)[:,None];fixed=[];P=[];who=[];bad=False
                for j,cc in enumerate(co):
                    value=cc[di]/den[di]
                    if j in frame or max(abs(cc-value*den))<1e-7:
                        kind=constant_gate(value)
                        if kind is None:bad=True;break
                        fixed.append(dict(wall=j,cosine=float(value),kind=kind))
                    else:P.append(cc-den);who.append(j)
                if bad:stats['bad_fixed_contact']+=1;continue
                P.extend([-num@M@H,num[:,0]]);who+=['center','positive_time'];P=np.array(P)
                ranges=allowed_intervals(P,model['positive'])
                if not ranges:stats['empty_slice']+=1;continue
                model.update(id=len(families),frame=frame,finite_labels=[m151,m],dotted_cosine=c,
                    projection_base=p0.tolist(),projection_slope=ps.tolist(),kernel=K.tolist(),
                    intervals=ranges,fixed_contacts=fixed,inequality_rows=who,inequality_coefficients=P.tolist())
                split,cutP=cells(model,Z);model['coverage_cells']=split;model['point_pairing_coefficients']=cutP.tolist()
                families.append(model)
                for k,cell in enumerate(split):
                    t=cell['parameter'];r=point(model,t)
                    if abs(r@M@r-2)>2e-4:stats['sample_norm_reject']+=1;continue
                    if not gate_mask((-(r@M@X.T)/np.sqrt(2*q))[None,:])[0] or r@M@H>=0 or r[0]<=0:
                        stats['sample_gate_reject']+=1;continue
                    ar=-r[1:]/r[0];cut=np.where(Z@ar>1+1e-9)[0]
                    if not len(cut):stats['empty_strict_coverage']+=1;continue
                    key=tuple(np.round(r,7))
                    if key in seen:stats['duplicate_sample']+=1;continue
                    seen.add(key)
                    if model['kind']!='null_line':assert np.linalg.norm(dynamic_point(model,c,t)-r)<1e-6
                    pool.append(dict(id=len(pool),family=model['id'],cell=k,dotted_cosine=c,parameter=t,
                        orthonormal_chart=r.tolist(),residual_constraint=ar.tolist(),joint_sample_cuts=cut.tolist()))
        frames.append(dict(frame=frame,families=len(families)-start))
        print('FRAME',frame,'families',len(families)-start,'total',len(families),'samples',len(pool),flush=True)
    src=[f'scripts/{TAG}.py',POINTS,'scripts/fresh_h20_two_anchor_families_20261016.py',
        'scripts/fresh_h20_quad_adjacent_h3_closure_20261016.py','scripts/fresh_h20_release144_cap_20261016_residual_scout.py']
    return dict(status='PASS_BOUNDED_ONE_ANCHOR_DOTTED_SLICE_DIAGNOSTIC',
        sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in src},
        anchor_domain=SIDES,dotted_anchor_domain=OTHERS,label151=[2,3],anchor_labels=LABELS,dotted_slices=SLICES,
        stats=dict(stats),frames=frames,families=families,candidates=pool,joint_samples=J['samples'],old_face=old,
        scope='Two free parameters in principle; only five dotted-cosine slices sampled. Each slice is numerically split at all saved sample pairing zeros. Unprescribed variable old contacts are dotted. No exact interval certificate, complete continuous search, or new Coxeter chamber claim.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
    if a.verify:assert close(d,json.loads(p.read_text()))
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],len(d['families']),len(d['candidates']),flush=True)
