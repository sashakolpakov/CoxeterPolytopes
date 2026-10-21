#!/usr/bin/env python3
"""Geometric relaxation: old gates Coxeter, new-to-new contacts merely acute."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from fresh_h20_continuous_sample_joint_cover_20261016 import Cover,bits
from fresh_h20_quad_adjacent_h3_closure_20261016 import geometry158,face_data,poly,gate_mask
from fresh_h20_release144_cap_20261016_residual_scout import close
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_two_anchor_families_20261016_acute_cover'
BASE='results/fresh_h20_two_anchor_families_20261016.json';PAIRS='results/fresh_h20_two_anchor_families_20261016_pairs.json'
class Acute(Cover):
    def neighbors(self,i):
        if i not in self.cache:
            co=-(self.roots@(self.metric@self.roots[i]))/2;good=co>=-1e-8;good[i]=False
            self.cache[i]=bits(good);self.normals_seen.append(i)
        return self.cache[i]
def run(iterations):
    D=json.loads((ROOT/BASE).read_text());P=json.loads((ROOT/PAIRS).read_text());G=geometry158();M=G['metric']
    base,V,A,b,ids=face_data(G)
    pool=[dict(origin='sample',origin_id=p['id'],family=p['family'],parameter=p['parameter'],orthonormal_chart=p['orthonormal_chart'],residual_constraint=p['residual_constraint'])for p in D['candidates']]
    pool += [dict(origin='solved_pair_partner',origin_id=i,family=p['family'],parameter=p['parameter'],orthonormal_chart=p['partner_orthonormal_chart'],residual_constraint=p['partner_constraint'])for i,p in enumerate(P['strong_pairs'])]
    R=np.array([p['orthonormal_chart']for p in pool]);AR=-R[:,1:]/R[:,0,None];solver=Acute(R,M)
    Z=[np.array(s['z'])for s in D['joint_samples']];finite=np.array([v['z']for v in D['old_face']['vertices']if v['radius2']<1-1e-6]);center=finite.mean(axis=0)
    steps=[];reason='iteration_guard';final=None
    for i in range(iterations):
        cov=np.array(Z)@AR.T>1+1e-7
        if not cov.any(axis=1).all():reason='uncovered_sample';break
        selected,meta=solver.search(cov)
        if selected is None:reason='no_sampled_all_acute_cover'if meta['complete']else'node_guard';steps.append(dict(iteration=i,search=meta));break
        co=-R[selected]@M@R[selected].T/2;np.fill_diagonal(co,0);assert co.min()>-1e-6
        aa=AR[selected];face=poly(base,V,np.r_[A,aa@V],np.r_[b,1-aa@base],ids+list(range(158,158+len(selected))))
        pairs=[]
        for j in range(len(selected)):
            for k in range(j+1,len(selected)):
                v=co[j,k]
                if not gate_mask(np.array([[v]]))[0]:pairs.append(dict(positions=[j,k],cosine=float(v),approx_label=float(np.pi/np.arccos(v))))
        row=dict(iteration=i,selected=selected,roots=[pool[j]for j in selected],search=meta,face=face,noncoxeter_new_pairs=pairs)
        steps.append(row);print('STEP',i,'caps',len(selected),'badpairs',len(pairs),'max',face['maximum_squared_radius'],flush=True)
        if not face['exterior']and not face['near_sphere']:reason='compact_geometric_face_candidate';final=row;break
        worst=max(face['vertices'],key=lambda x:x['radius2'])
        if worst['radius2']<=1+1e-6:reason='near_sphere';final=row;break
        v=np.array(worst['z'])-center;t=next(float(t)for t in np.roots([v@v,2*center@v,center@center-1])if np.isreal(t)and 0<float(t)<1)
        z=center+t*v;assert abs(z@z-1)<1e-7 and np.max(aa@z-1)<1e-6;Z.append(z)
    src=[f'scripts/{TAG}.py','scripts/fresh_h20_continuous_sample_joint_cover_20261016.py',BASE,PAIRS]
    return dict(status='PASS_BOUNDED_ALL_ACUTE_CLUSTER_RELAXATION',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in src},candidate_count=len(pool),iteration_limit=iterations,termination=reason,steps=steps,final=final,samples=[z.tolist()for z in Z],
      scope='Numerical finite bank relaxation: every old contact passes source Coxeter gates, while new-new only acute. Even a compact geometric face is not a Coxeter construction until every flagged mutual angle is solved and all geometry exactified. No continuous-domain exclusion claim.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--iterations',type=int,default=100);ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run(a.iterations);p=ROOT/f'results/{TAG}.json'
    if a.verify:assert close(d,json.loads(p.read_text()))
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['termination'],len(d['steps']),flush=True)
