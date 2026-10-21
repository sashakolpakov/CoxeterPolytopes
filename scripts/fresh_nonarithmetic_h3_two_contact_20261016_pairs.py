#!/usr/bin/env python3
"""Solve mutual finite-angle equations inside numerical conic intervals."""
import argparse,json,itertools,math
from collections import Counter
from hashlib import sha256
import numpy as np
import fresh_nonarithmetic_h3_two_contact_20261016 as C
ROOT=C.ROOT;TAG=C.TAG+'_pairs';ETA=C.ETA

def inside(t,intervals,positive):
    if not math.isfinite(t)or(positive and t<=0):return False
    return any((a is None or t>a+1e-8*max(1,abs(a)))and(b is None or t<b-1e-8*max(1,abs(b)))for a,b in intervals)

def clip(current,v,wall):
    vals=np.array([r['v']@ETA@v for r in current]);sg=np.sign(np.where(abs(vals)<2e-8,0,vals));neg=np.flatnonzero(sg<0);pos=np.flatnonzero(sg>0);zero=np.flatnonzero(sg==0)
    out={r['mask']:r for i,r in enumerate(current)if sg[i]<0}
    for i in zero:
        r=dict(current[i]);r['mask']|=1<<wall;out[r['mask']]=r
    for i,j in itertools.product(neg,pos):
        mask=current[i]['mask']&current[j]['mask']
        if mask.bit_count()<8:continue
        face=[k for k,r in enumerate(current)if r['mask']&mask==mask]
        if face!=sorted([int(i),int(j)]):continue
        z=vals[j]*current[i]['v']-vals[i]*current[j]['v']
        if abs(z[0])<1e-12:raise ArithmeticError('zero chart')
        z/=z[0];q=z@ETA@z;kind='finite_vertices'if q<-1e-7 else'hyperideal_rays'if q>1e-7 else'unresolved_null'
        mask|=1<<wall
        if mask in out:raise ArithmeticError('mask collision')
        out[mask]=dict(mask=mask,v=z,kind=kind)
    return sorted(out.values(),key=lambda r:r['mask'])

def run(maxseeds,maxlabel):
    path=ROOT/f'results/{C.TAG}.json';d=json.loads(path.read_text());NC=np.array(d['old_projected_normals'])
    blocked=json.loads((ROOT/f'results/{C.TAG}_affine_gate.json').read_text())['excluded_family_indices']
    families=[f for f in d['families']if f['index']not in blocked]
    hyp=[r for r in d['actual_face_rays']if r['kind']=='hyperideal_rays'];points=np.array([np.r_[1.,np.array(r['point'][1:])/np.linalg.norm(r['point'][1:])]for r in hyp])
    groups=[{0,1,3},{2,4,5,6}]
    initial=[dict(mask=sum(1<<i for i in r['active']),v=np.array(r['point']),kind=r['kind'])for r in d['actual_face_rays']]
    seeds=[]
    for f in families:
        if max(f['labels'])>maxlabel:continue
        for rec in f['intervals']:
            for s in rec['samples']:
                v=C.evaluate(f,s['parameter']);pv=points@ETA@v;cover=set(np.flatnonzero(pv>1e-7).tolist())
                if not cover:continue
                seeds.append(dict(family=f['index'],parameter=s['parameter'],v=v,coverage=cover,score=float(sum(np.maximum(pv,0)))))
    seeds.sort(key=lambda s:(-len(s['coverage']),-s['score'],s['family'],s['parameter']))
    buckets={}
    for s in seeds:buckets.setdefault(tuple(sorted(s['coverage'])),[]).append(s)
    chosen=[];per=Counter();positions={k:0 for k in buckets}
    while len(chosen)<maxseeds:
        changed=False
        for coverage in sorted(buckets):
            rows=buckets[coverage]
            while positions[coverage]<len(rows):
                s=rows[positions[coverage]];positions[coverage]+=1;key=(s['family'],coverage)
                if per[key]>=2:continue
                per[key]+=1;chosen.append(s);changed=True;break
            if len(chosen)>=maxseeds:break
        if not changed:break
    counts=Counter();best={};bridges={};union_counts=Counter()
    for si,s in enumerate(chosen):
        u=s['v'];cov=ETA@u
        for f in families:
            if max(f['labels'])>maxlabel:continue
            intervals=[r['interval']for r in f['intervals']]
            for m in range(2,11):
                goal=0. if m==2 else-2*math.cos(math.pi/m);qq=C.poly(f,cov,-goal);counts['mutual_equations']+=1
                if max(abs(qq))<1e-10:counts['identically_fixed_mutual']+=1;continue
                for t in C.roots(qq):
                    if not inside(t,intervals,f['kind']=='hyperbola'):continue
                    v=C.evaluate(f,t);prod=u@ETA@v
                    if abs(prod-goal)>1e-6:continue
                    g=NC@ETA@v
                    if v[0]<1e-8 or abs(v@ETA@v-2)>2e-5 or not all(C.compatible(x,2e-5)for x in g):continue
                    coverage=set(np.flatnonzero(points@ETA@v>1e-7).tolist());union=coverage|s['coverage'];counts['compatible_finite_pairs']+=1
                    union_counts[len(union)]+=1
                    if coverage-s['coverage'] and s['coverage']-coverage:
                        key=tuple(sorted(union));record=dict(first=dict(family=s['family'],parameter=s['parameter'],vector=u,coverage=sorted(s['coverage'])),
                            second=dict(family=f['index'],parameter=t,vector=v,coverage=sorted(coverage)),mutual_label=m)
                        score=float(sum(np.maximum(points@ETA@u,points@ETA@v)))
                        bridges.setdefault(key,[]).append((score,record))
                        bridges[key].sort(key=lambda z:-z[0]);bridges[key]=bridges[key][:12]
                    targets=[i for i,group in enumerate(groups)if group<=union]
                    if not targets:continue
                    counts['initial_region_covered_pairs']+=1
                    try:rr=clip(clip(initial,u,55),v,56)
                    except ArithmeticError:counts['unresolved_clip']+=1;continue
                    kinds=dict(Counter(r['kind']for r in rr));exterior=[r for r in rr if r['kind']!='finite_vertices']
                    rec=dict(first=dict(family=s['family'],parameter=s['parameter'],vector=u,coverage=sorted(s['coverage'])),
                             second=dict(family=f['index'],parameter=t,vector=v,coverage=sorted(coverage)),mutual_label=m,
                             counts=kinds,covered_initial_regions=targets,
                             remaining_exterior=[dict(active=[i for i in range(r['mask'].bit_length())if r['mask']>>i&1],point=r['v'],kind=r['kind'])for r in exterior])
                    score=(len(exterior),sum(max(0,r['v']@ETA@r['v'])for r in exterior))
                    signature=(s['family'],f['index'],m,round(t,7),round(s['parameter'],7))
                    best[signature]=(score,rec)
        if (si+1)%10==0:print('seed',si+1,'counts',dict(counts),'best',min((r[0]for r in best.values()),default=None),flush=True)
    ordered=sorted(best.values(),key=lambda z:z[0]);records=[r for _,r in ordered[:100]]
    paths=[f'scripts/{TAG}.py',f'scripts/{C.TAG}.py',f'results/{C.TAG}.json',f'results/{C.TAG}_affine_gate.json']
    return C.clean(dict(status='BOUNDED_NUMERICAL_TWO_CONTACT_MUTUAL_ANGLE_SEARCH',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
                        max_seeds=maxseeds,max_anchor_label=maxlabel,mutual_labels=list(range(2,11)),counts=dict(counts),
                        seed_count=len(chosen),seed_coverage_counts=dict(Counter(str(sorted(s['coverage']))for s in chosen)),eligible_sample_count=len(seeds),candidate_pair_count=len(best),best_pairs=records,
                        union_coverage_counts={str(k):v for k,v in union_counts.items()},bridge_pairs=[r for key in sorted(bridges)for _,r in bridges[key]],
                        scope='For bounded chosen first-cap samples, solves every mutual-label2..10 quadratic on each bounded admissible second family. Whole actual H3 cone clipping tests each initial-region covering pair. Discovery only: no exact contacts, degeneracy resolution, persistent-cusp exclusion, complete family packing or compactness proof.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--maxseeds',type=int,default=160);ap.add_argument('--maxlabel',type=int,default=6);ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    d=run(a.maxseeds,a.maxlabel);p=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(p.read_text());print('PASS saved pair search')
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['candidate_pair_count'],flush=True)
