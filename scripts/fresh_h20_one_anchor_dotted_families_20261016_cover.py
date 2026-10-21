#!/usr/bin/env python3
"""Bounded all-acute joint covers including the new two-variable slices."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from fresh_h20_continuous_sample_joint_cover_20261016 import Cover,bits
from fresh_h20_quad_adjacent_h3_closure_20261016 import geometry158,face_data,poly,gate_mask
from fresh_h20_near_closed_cluster_recovery_20261016 import adjacency
from fresh_h20_release144_cap_20261016_residual_scout import close
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_one_anchor_dotted_families_20261016_cover'
BASE='results/fresh_h20_two_anchor_families_20261016.json'
NEW='results/fresh_h20_one_anchor_dotted_families_20261016.json'
ACUTE='results/fresh_h20_two_anchor_families_20261016_acute_cover.json'
PAIRS='results/fresh_h20_two_anchor_families_20261016_pairs.json'
class Acute(Cover):
    def neighbors(self,i):
        if i not in self.cache:
            co=-(self.roots@(self.metric@self.roots[i]))/2;ok=co>=-1e-8;ok[i]=False
            self.cache[i]=bits(ok);self.normals_seen.append(i)
        return self.cache[i]
def run(limit):
    D=json.loads((ROOT/BASE).read_text());N=json.loads((ROOT/NEW).read_text());P=json.loads((ROOT/PAIRS).read_text());K=json.loads((ROOT/ACUTE).read_text())
    G=geometry158();M=G['metric'];base,V,A,b,ids=face_data(G)
    pool=[dict(origin='dotted_slice',**r)for r in N['candidates']]
    pool += [dict(origin='old_sample',**r)for r in D['candidates']]
    pool += [dict(origin='old_pair',id=i,family=r['family'],parameter=r['parameter'],orthonormal_chart=r['partner_orthonormal_chart'],residual_constraint=r['partner_constraint'])for i,r in enumerate(P['strong_pairs'])]
    R=np.array([x['orthonormal_chart']for x in pool]);AR=-R[:,1:]/R[:,0,None]
    Z=[np.array(z)for z in K['samples']];solver=Acute(R,M)
    center=np.mean([v['z']for v in D['old_face']['vertices']if v['radius2']<1-1e-6],axis=0)
    steps=[];final=None;reason='iteration_guard'
    for k in range(limit):
        cov=np.array(Z)@AR.T>1+1e-7
        if not cov.any(axis=1).all():reason='uncovered_sample';break
        chosen,meta=solver.search(cov)
        if chosen is None:reason='no_all_acute_sampled_cover'if meta['complete']else'node_guard';break
        rr=R[chosen];co=-rr@M@rr.T/2;np.fill_diagonal(co,0);assert co.min()>-1e-6
        aa=AR[chosen];face=poly(base,V,np.r_[A,aa@V],np.r_[b,1-aa@base],ids+list(range(158,158+len(chosen))))
        adj,oa,old,q=adjacency(rr,G)
        row=dict(iteration=k,selected=chosen,roots=[pool[i]for i in chosen],search=meta,face=face,
            potential_degrees=(adj.sum(1)+oa.sum(1)).tolist(),
            intrinsic_new_pairs=[[i,j]for i in range(len(chosen))for j in range(i+1,len(chosen))if adj[i,j]],
            noncoxeter_new_pairs=[dict(positions=[i,j],cosine=float(co[i,j]))for i in range(len(chosen))for j in range(i+1,len(chosen))if not gate_mask(np.array([[co[i,j]]]))[0]])
        steps.append(row)
        print('STEP',k,'caps',len(chosen),'two_variable',sum(pool[i]['origin']=='dotted_slice'for i in chosen),'degree',min(row['potential_degrees']),'max',face['maximum_squared_radius'],flush=True)
        if not face['exterior']and not face['near_sphere']:reason='compact_acute_geometric_candidate';final=row;break
        worst=max(face['vertices'],key=lambda v:v['radius2'])
        if worst['radius2']<=1+1e-6:reason='near_sphere';final=row;break
        direction=np.array(worst['z'])-center
        t=next(float(t)for t in np.roots([direction@direction,2*center@direction,center@center-1])if np.isreal(t)and 0<float(t)<1)
        z=center+t*direction;assert abs(z@z-1)<1e-7 and max(aa@z-1)<1e-6;Z.append(z)
    src=[f'scripts/{TAG}.py',BASE,NEW,ACUTE,PAIRS,'scripts/fresh_h20_continuous_sample_joint_cover_20261016.py','scripts/fresh_h20_near_closed_cluster_recovery_20261016.py']
    return dict(status='PASS_BOUNDED_TWO_VARIABLE_SAMPLE_ACUTE_COVERS',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in src},candidate_count=len(pool),iteration_limit=limit,termination=reason,steps=steps,final=final,samples=[z.tolist()for z in Z],scope='Finite sample-bank all-acute relaxation, old Coxeter gates retained. New-new contacts need later exact angle recovery. Explicit iteration and inherited node guards; no continuous exclusion or finite-volume claim.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--iterations',type=int,default=30);a=ap.parse_args();d=run(a.iterations);p=ROOT/f'results/{TAG}.json'
    if a.verify:assert close(d,json.loads(p.read_text()))
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['termination'],len(d['steps']),flush=True)
