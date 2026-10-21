#!/usr/bin/env python3
"""Bounded numerical symmetric short-arc cap bank and compatibility scout.

No finite-label or all-parameter exclusion is claimed. All accepted rows
require exactification before promotion. Normalization is norm two.
"""
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
from itertools import combinations
from collections import Counter
import numpy as np

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_dotted_arc_bank_20261016_scout'
INPUT='results/fresh_h20_dotted_arc_cap_20261016_scout.json'

def classify(g,tol=2e-7):
    if abs(g)<tol:return 2
    if g<-2-tol:return 'dotted'
    if abs(g+2)<tol:return 'parallel'
    if g>tol:return 'obtuse'
    m=np.pi/np.arccos(-g/2);n=int(round(m))
    if n>=3 and abs(g+2*np.cos(np.pi/n))<tol:return n
    return 'forbidden'

def run():
    inp=json.loads((ROOT/INPUT).read_text())['inputs']
    base=json.loads((ROOT/'results/fresh_h20_s61_multicap_20261014_exact.json').read_text())
    P=np.array(inp['P'],float);W=np.array(inp['W'],float);J=np.array(inp['J'],float)
    S=inp['support'];ep=np.array(inp['ep'],float);ae=np.array(inp['AE'],float)
    def mid(q):return float((F(q[0])+F(q[1]))/2)
    co=np.array([[mid(x)for x in v]for v in base['chain']['frame_coefficient_bounds']]);aa=mid(base['parameter']['isolating_interval'])
    oldep=np.array([mid(q)for q in base['strict_escape']['frame_pair_bounds']])
    # The old escape's support pairings are reconstructed in the high-precision
    # discovery load separately if needed; the exact chosen cap already records
    # its sign. This bank targets the stronger saved escape and new simple one.
    ez=np.array([1,-.01,0]);eb=np.full(18,-1e-5);eb[S.index(5)]=-1
    Je=J@eb
    labels=[2]+list(range(3,21))+[24,32,64,128,0]
    values={m:0. if m==2 else 2. if m==0 else 2*np.cos(np.pi/m)for m in labels}
    rows=[];counts=Counter();trials=0
    for center in range(108,132):
        for ma,A in values.items():
            for mb,B in values.items():
                if B<=A:continue
                y=np.r_[((B-A)*aa*aa+2*B)/(2*aa),(B-A)*co[center-107,1:]/2]
                q=(A*A-B*B)/2-B*B/(aa*aa)
                if q>=0:continue
                vp=P@y
                for j in range(18):
                    trials+=1;c=np.sqrt((2-q)/J[j,j]);g=vp-c*W[j]
                    if g.max()>2e-7:counts['obtuse']+=1;continue
                    rel=[classify(x)for x in g]
                    if 'forbidden'in rel or 'obtuse'in rel:counts['forbidden']+=1;continue
                    counts['accepted']+=1
                    cuts=dict(bank847=float(ep@y+c*ae[j]),simple_new=float(-y[0]*ez[0]+y[1:]@ez[1:]-c*Je[j]))
                    rows.append(dict(center=center,middle_label=ma or 'parallel',side_label=mb or 'parallel',
                                     support_wall=S[j],dual_index=j,frame=y.tolist(),dual_coefficient=float(c),
                                     residual_norm=float(q),relations=rel,pairs=g.tolist(),escape_pairings=cuts))
    seed=next(i for i,r in enumerate(rows)if r['center']==129 and r['middle_label']==3 and r['side_label']==4 and r['support_wall']==31)
    compatibility=Counter();seed_compatible=[];cutter_pairs=[]
    Y=np.array([r['frame']for r in rows]);cs=np.array([r['dual_coefficient']for r in rows])
    js=np.array([r['dual_index']for r in rows]);cutters={i for i,r in enumerate(rows)if r['escape_pairings']['bank847']>1e-7}
    for i,r in enumerate(rows):
        gs=-Y[i,0]*Y[:i,0]+Y[:i,1:]@Y[i,1:]+cs[i]*cs[:i]*J[js[i],js[:i]]
        dotted=gs<-2-2e-7;right=np.abs(gs)<2e-7;parallel=np.abs(gs+2)<2e-7
        ok=dotted|right|parallel
        compatibility['dotted']+=int(dotted.sum());compatibility['right']+=int(right.sum());compatibility['parallel']+=int(parallel.sum())
        for k in np.flatnonzero((gs<0)&(gs>-2)&~ok):
            rel=classify(float(gs[k]))
            if rel!='forbidden':ok[k]=True;compatibility[f'pi/{rel}']+=1
        if i==seed:seed_compatible +=list(map(int,np.flatnonzero(ok)))
        elif i>seed and ok[seed]:seed_compatible.append(i)
        if i in cutters:
            for k in np.flatnonzero(ok):
                if int(k)in cutters:cutter_pairs.append(dict(rows=[int(k),i],relation=classify(float(gs[k])),pair=float(gs[k])))
    sources=[INPUT,'results/fresh_h20_s61_multicap_20261014_exact.json',
             'results/fresh_h20_dotted_arc_cap_20261016.json',f'scripts/{TAG}.py']
    return dict(status='BOUNDED_NUMERICAL_SYMMETRIC_ARC_CAP_AND_COMPATIBILITY_BANK',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in sources},labels=labels,
                centers=list(range(108,132)),support=S,trials=trials,counts=dict(counts),rows=rows,
                exact_cap_row=seed,compatible_pair_counts=dict(compatibility),cutter_compatible_pairs=cutter_pairs,
                compatible_with_exact_cap=seed_compatible,
                summary=dict(accepted=len(rows),support_counts=dict(Counter(r['support_wall']for r in rows)),
                    bank_escape_cutters=sum(r['escape_pairings']['bank847']>1e-7 for r in rows),
                    simple_escape_cutters=sum(r['escape_pairings']['simple_new']>1e-7 for r in rows),
                    compatible_pairs=sum(compatibility.values()),compatible_with_exact_cap=len(seed_compatible)),
                scope='Finite label and symmetric three-contact recipe bank only. Float comparisons are discovery diagnostics, never exact Coxeter/compatibility claims or all-parameter exclusions.')

if __name__=='__main__':
    d=run();(ROOT/f'results/{TAG}.json').write_text(json.dumps(d,indent=2)+'\n');print(d['status'],d['summary'],flush=True)
