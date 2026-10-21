#!/usr/bin/env python3
"""Bounded numerical compatible-cover search on a complete H3 face."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from fresh_h20_nonright151_h3_closure_20261016 import geometry158,face_data,poly,gate_mask
from fresh_h20_release144_cap_20261016_residual_scout import close
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_nonright151_h3_joint_cover_20261016'
INPUT='results/fresh_h20_nonright151_h3_closure_20261016.json'

def bits(row):
    return int.from_bytes(np.packbits(row,bitorder='little').tobytes(),'little')

def clique_cover(cov,adjacency):
    """Complete finite backtracking, unless the explicit node guard is reached."""
    sample_masks=[bits(row) for row in cov]
    coverage=[bits(row) for row in cov.T]
    full=(1<<len(cov))-1;failed=set();nodes=0
    def dfs(allowed,need):
        nonlocal nodes
        nodes+=1
        if nodes>2000000:raise RuntimeError('node_limit')
        if not need:return []
        key=(allowed,need)
        if key in failed:return None
        q=need;choices=None
        while q:
            bit=q&-q;j=bit.bit_length()-1;q-=bit
            available=allowed&sample_masks[j]
            if not available:failed.add(key);return None
            if choices is None or available.bit_count()<choices.bit_count():choices=available
        order=[]
        while choices:
            bit=choices&-choices;i=bit.bit_length()-1;choices-=bit
            order.append(i)
        order.sort(key=lambda i:(-(coverage[i]&need).bit_count(),i))
        for i in order:
            tail=dfs(allowed&adjacency[i],need&~coverage[i])
            if tail is not None:return [i]+tail
        failed.add(key);return None
    try:return dfs((1<<len(coverage))-1,full),dict(nodes=nodes,failed_states=len(failed),complete=True)
    except RuntimeError:return None,dict(nodes=nodes,failed_states=len(failed),complete=False)

def run():
    p=json.loads((ROOT/INPUT).read_text());G=geometry158();base,V,A,b,ids=face_data(G)
    pool=p['candidates'];n=len(pool);assert n==8009
    roots=np.array([x['orthonormal_chart'] for x in pool]);ar=np.array([x['residual_constraint'] for x in pool])
    comp=np.empty((n,n),dtype=bool)
    for start in range(0,n,128):
        co=-(roots[start:start+128]@G['metric']@roots.T)/2
        for j,row in enumerate(co):comp[start+j]=gate_mask(row[:,None])
    assert np.array_equal(comp,comp.T)
    old=p['old_face']['vertices'];finite=np.array([x['z'] for x in old if x['radius2']<1-1e-6])
    center=finite.mean(axis=0);assert center@center<1
    samples=[np.array(x['z']) for x in old if x['radius2']>1+1e-6]
    samplekind=['old_exterior_vertex']*len(samples)
    adjacency=[bits(row) for row in comp]
    steps=[];final=None;reason=None
    for step in range(1,101):
        cov=(np.array(samples)@ar.T>1+1e-7)
        if not np.all(cov.any(axis=1)):
            reason='uncovered_feasible_ideal_sample';break
        selected,search=clique_cover(cov,adjacency)
        if selected is None:
            reason='infeasible_bounded_compatible_cover' if search['complete'] else 'node_limit'
            steps.append(dict(iteration=step,search=search,samples=len(samples)));break
        assert all(comp[a,c] for i,a in enumerate(selected) for c in selected[i+1:])
        assert np.all(cov[:,selected].any(axis=1))
        rec=dict(iteration=step,selected=selected,search=search,samples=len(samples))
        aa=ar[selected]@V;bb=1-ar[selected]@base
        face=poly(base,V,np.r_[A,aa],np.r_[b,bb],ids+list(range(158,158+len(selected))))
        rec['face_counts']={k:face.get(k) for k in ['finite','near_sphere','exterior','maximum_squared_radius']};steps.append(rec)
        print('STEP',step,'compatible',selected,rec['face_counts'],flush=True)
        assert not face.get('degenerate')
        if face['exterior']==0 and face['near_sphere']==0:
            final=dict(selected=selected,face=face);reason='compact_face_candidate';break
        worst=max(face['vertices'],key=lambda x:x['radius2']);z=np.array(worst['z'])
        if worst['radius2']<=1+1e-6:
            # Near-sphere vertices are not silently classified as compact.
            reason='near_sphere_requires_exact_classification';final=dict(selected=selected,face=face);break
        delta=z-center;coeff=[delta@delta,2*center@delta,center@center-1]
        ts=np.roots(coeff);t=next(float(x) for x in ts if np.isreal(x) and 0<float(x)<1)
        ideal=center+t*delta;assert abs(ideal@ideal-1)<1e-7
        yy=(ideal-base)@V
        assert max(A@yy-b)<1e-6 and max(ar[selected]@ideal-1)<1e-6
        assert not any(np.linalg.norm(ideal-x)<1e-7 for x in samples)
        samples.append(ideal);samplekind.append('feasible_ideal_boundary_sample')
    else:reason='iteration_limit'
    src=[f'scripts/{TAG}.py',INPUT,'scripts/fresh_h20_nonright151_h3_closure_20261016.py','results/fresh_h20_compact_quad_cap_20261016.json']
    return dict(status='PASS_BOUNDED_H3_COMPATIBLE_COVER_DIAGNOSTIC',
       sources={x:hashlib.sha256((ROOT/x).read_bytes()).hexdigest() for x in src},
       candidate_count=n,termination=reason,steps=steps,compatible_pairs=int(comp.sum()//2),
       samples=[dict(z=x.tolist(),kind=k) for x,k in zip(samples,samplekind)],final=final,
       scope='Numerical finite-clique covers with actual all158 angle gates, complete mutual compatibility, and full3D face intersection after each compatible cover. Explicit2million-node search guard and100sample iterations; numerical tolerances, no exact completion or universal impossibility theorem.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert close(d,json.loads(path.read_text()))
    else:path.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['termination'],len(d['steps']),flush=True)
