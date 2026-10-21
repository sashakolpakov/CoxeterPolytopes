#!/usr/bin/env python3
"""Numerical finite-clique closure using all continuous-family samples.

Adjacency is computed on demand, so no quadratic dense matrix or arbitrary
coverage-ranking truncation of the candidate bank is needed.
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from fresh_h20_quad_adjacent_h3_closure_20261016 import geometry158,face_data,poly,gate_mask
from fresh_h20_release144_cap_20261016_residual_scout import close

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_continuous_sample_joint_cover_20261016'
BASE='results/fresh_h20_two_anchor_families_20261016.json'
PAIRS='results/fresh_h20_two_anchor_families_20261016_pairs.json'

def bits(row):return int.from_bytes(np.packbits(row,bitorder='little').tobytes(),'little')

class Cover:
    def __init__(self,roots,metric):
        self.roots=roots;self.metric=metric;self.cache={};self.normals_seen=[]
    def neighbors(self,i):
        if i not in self.cache:
            co=-(self.roots@(self.metric@self.roots[i]))/2
            good=gate_mask(co[:,None]);good[i]=False
            self.cache[i]=bits(good);self.normals_seen.append(i)
        return self.cache[i]
    def search(self,cov):
        samples=[bits(x)for x in cov];coverage=[bits(x)for x in cov.T]
        allneed=(1<<len(samples))-1;failed=set();nodes=0
        def dfs(allowed,need):
            nonlocal nodes
            nodes+=1
            if nodes>2000000:raise RuntimeError('node_guard')
            if not need:return []
            key=(allowed,need)
            if key in failed:return None
            todo=need;choices=None
            while todo:
                v=todo&-todo;j=v.bit_length()-1;todo-=v
                possibilities=allowed&samples[j]
                if not possibilities:failed.add(key);return None
                if choices is None or possibilities.bit_count()<choices.bit_count():choices=possibilities
            order=[]
            while choices:
                v=choices&-choices;j=v.bit_length()-1;choices-=v;order.append(j)
            order.sort(key=lambda j:(-(coverage[j]&need).bit_count(),j))
            for j in order:
                tail=dfs(allowed&self.neighbors(j),need&~coverage[j])
                if tail is not None:return [j]+tail
            failed.add(key);return None
        try:
            answer=dfs((1<<len(coverage))-1,allneed);complete=True
        except RuntimeError:answer=None;complete=False
        return answer,dict(nodes=nodes,failed_states=len(failed),complete=complete,adjacency_rows_computed=len(self.cache))

def run():
    base_doc=json.loads((ROOT/BASE).read_text());pair_doc=json.loads((ROOT/PAIRS).read_text())
    G=geometry158();base,V,A,b,ids=face_data(G)
    pool=[dict(origin='family_sample',origin_id=p['id'],family=p['family'],parameter=p['parameter'],
               orthonormal_chart=p['orthonormal_chart'],residual_constraint=p['residual_constraint']) for p in base_doc['candidates']]
    for i,p in enumerate(pair_doc['strong_pairs']):
        pool.append(dict(origin='mutual_angle_solution',origin_id=i,family=p['family'],parameter=p['parameter'],
                   orthonormal_chart=p['partner_orthonormal_chart'],residual_constraint=p['partner_constraint']))
    roots=np.array([p['orthonormal_chart']for p in pool]);ar=np.array([p['residual_constraint']for p in pool])
    assert len(pool)==42203
    solver=Cover(roots,G['metric'])
    samples=[np.array(p['z'])for p in base_doc['joint_samples']]
    samplekind=[p['kind']for p in base_doc['joint_samples']]
    old=base_doc['old_face'];finite=np.array([p['z']for p in old['vertices']if p['radius2']<1-1e-6])
    center=finite.mean(axis=0);assert center@center<1
    steps=[];final=None;reason=None
    for iteration in range(1,101):
        cov=np.array(samples)@ar.T>1+1e-7
        if not np.all(cov.any(axis=1)):
            reason='sample_uncovered_by_entire_bank';break
        selected,search=solver.search(cov)
        if selected is None:
            reason='infeasible_bounded_cover'if search['complete']else'node_guard'
            steps.append(dict(iteration=iteration,search=search,samples=len(samples)));break
        assert np.all(cov[:,selected].any(axis=1))
        # Check both orientations independently of cached row membership.
        mutual=-(roots[selected]@G['metric']@roots[selected].T)/2
        np.fill_diagonal(mutual,0)
        assert gate_mask(mutual).all()
        aa=ar[selected]@V;bb=1-ar[selected]@base
        face=poly(base,V,np.r_[A,aa],np.r_[b,bb],ids+list(range(158,158+len(selected))))
        assert not face.get('degenerate')
        counts={k:face[k]for k in ['finite','near_sphere','exterior','maximum_squared_radius']}
        rec=dict(iteration=iteration,selected=selected,search=search,samples=len(samples),face_counts=counts)
        steps.append(rec);print('STEP',iteration,selected,counts,search,flush=True)
        if not face['exterior']and not face['near_sphere']:
            reason='compact_face_candidate';final=dict(selected=selected,roots=[pool[i]for i in selected],face=face);break
        worst=max(face['vertices'],key=lambda p:p['radius2'])
        if worst['radius2']<=1+1e-6:
            reason='near_sphere_requires_exact_classification';final=dict(selected=selected,roots=[pool[i]for i in selected],face=face);break
        z=np.array(worst['z']);direction=z-center
        ts=np.roots([direction@direction,2*center@direction,center@center-1])
        t=next(float(x)for x in ts if np.isreal(x)and 0<float(x)<1)
        ideal=center+t*direction;yy=(ideal-base)@V
        assert abs(ideal@ideal-1)<1e-7 and max(A@yy-b)<1e-6
        assert max(ar[selected]@ideal-1)<1e-6
        assert not any(np.linalg.norm(ideal-p)<1e-7 for p in samples)
        samples.append(ideal);samplekind.append('new_feasible_ideal_boundary_sample')
    else:reason='iteration_guard'
    sources=[f'scripts/{TAG}.py',BASE,PAIRS,'scripts/fresh_h20_quad_adjacent_h3_closure_20261016.py',
             'results/fresh_h20_compact_quad_cap_20261016.json']
    return dict(status='PASS_BOUNDED_ALL_CONTINUOUS_SAMPLES_JOINT_COVER',
       sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in sources},candidate_count=len(pool),
       termination=reason,steps=steps,samples=[dict(z=p.tolist(),kind=k)for p,k in zip(samples,samplekind)],
       tested_adjacency_rows=solver.normals_seen,final=final,
       scope='Numerical finite candidate bank, all41984saved samples plus219solved partners. No ranking truncation; on-demand complete-clique search subject to explicit2million-node/100sample guards. Every compatible cover gets a fullH3intersection and replacement ideal samples. No exact or continuous-domain nonexistence claim.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run()
    p=ROOT/f'results/{TAG}.json'
    if args.verify:assert close(d,json.loads(p.read_text()))
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['termination'],len(d['steps']),flush=True)
