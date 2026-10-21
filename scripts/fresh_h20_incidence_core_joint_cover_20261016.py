#!/usr/bin/env python3
"""Numerical clique cover with simultaneous selected-facet degree obligations."""
import argparse,hashlib,json
from collections import Counter
from pathlib import Path
import numpy as np
from fresh_h20_intrinsic_facet_incidence_20261016 import bank
from fresh_h20_continuous_sample_joint_cover_20261016 import bits
from fresh_h20_quad_adjacent_h3_closure_20261016 import geometry158,face_data,poly,gate_mask
from fresh_h20_release144_cap_20261016_residual_scout import close

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_incidence_core_joint_cover_20261016'
CORE='results/fresh_h20_intrinsic_facet_incidence_20261016.json'
OLD='results/fresh_h20_continuous_sample_joint_cover_20261016.json'

def members(mask):
    while mask:
        v=mask&-mask;mask-=v;yield v.bit_length()-1

class Search:
    def __init__(self,compat,finite,olddegree):
        self.compat=[bits(x)for x in compat]
        self.finite=[bits(x)for x in finite]
        self.olddegree=list(map(int,olddegree));self.n=len(olddegree)
    def solve(self,cov):
        point=[bits(x)for x in cov];covers=[bits(x)for x in cov.T]
        nodes=0;failed=set();prunes=Counter()
        def dfs(allowed,need,selected):
            nonlocal nodes
            nodes+=1
            if nodes>2000000:raise RuntimeError('node_guard')
            # Selected facets must be able to reach degree3. Unselected
            # candidates incapable of doing so may also be deleted safely.
            changed=True
            while changed:
                changed=False;allpossible=allowed|selected
                for v in members(selected):
                    if self.olddegree[v]+(self.finite[v]&allpossible).bit_count()<3:
                        prunes['selected_degree_impossible']+=1;return None
                remove=0
                for v in members(allowed):
                    if self.olddegree[v]+(self.finite[v]&allpossible).bit_count()<3:remove|=1<<v
                if remove:
                    allowed&=~remove;changed=True;prunes['unselected_degree_removed']+=remove.bit_count()
            key=(allowed,need,selected)
            if key in failed:return None
            choices=None
            for p in members(need):
                opts=allowed&point[p]
                if not opts:prunes['uncovered_point']+=1;failed.add(key);return None
                if choices is None or opts.bit_count()<choices.bit_count():choices=opts
            for v in members(selected):
                deficit=3-self.olddegree[v]-(self.finite[v]&selected).bit_count()
                if deficit>0:
                    opts=self.finite[v]&allowed
                    if opts.bit_count()<deficit:raise AssertionError('missed degree pruning')
                    if choices is None or opts.bit_count()<choices.bit_count():choices=opts
            if choices is None:return list(members(selected))
            order=sorted(members(choices),key=lambda j:(-(covers[j]&need).bit_count(),
                         -(self.finite[j]&selected).bit_count(),j))
            for j in order:
                result=dfs(allowed&self.compat[j],need&~covers[j],selected|(1<<j))
                if result is not None:return result
            failed.add(key);return None
        try:
            ans=dfs((1<<self.n)-1,(1<<len(point))-1,0);complete=True
        except RuntimeError:ans=None;complete=False
        return ans,dict(nodes=nodes,failed_states=len(failed),complete=complete,prunes=dict(prunes))

def run():
    c=json.loads((ROOT/CORE).read_text());old=json.loads((ROOT/OLD).read_text())
    pool,_=bank();keep=c['surviving_degree_core'];rows=[pool[i]for i in keep]
    R=np.array([r['root']for r in rows]);labs=np.array([r['label151']for r in rows])
    h=-2*np.cos(np.pi/labs);h[labs==2]=0;qr=2-h*h/2
    G=geometry158();M=G['metric'];X=G['X'].astype(np.longdouble);u=X[151];qu=u@M@u
    hu=X@M@u;oq=np.einsum('ij,jk,ik->i',X,M,X)-hu*hu/qu
    oi=np.where(oq>1e-7)[0]
    op=R.astype(np.longdouble)@M@X[oi].T-h[:,None]*hu[oi]/qu
    olddegree=(qr[:,None]*oq[oi]-op*op>1e-6).sum(axis=1)
    co=-R@M@R.T/2;np.fill_diagonal(co,0)
    valid=(abs(co)<2e-6)|(co>=1-2e-6)
    inside=(co>0)&(co<1-2e-6);angle=np.zeros_like(co);angle[inside]=np.pi/np.arccos(co[inside])
    valid|=inside&(angle>=2.9999)&(abs(angle-np.rint(angle))<2e-4)
    valid&=valid.T;np.fill_diagonal(valid,False)
    finite=valid&(qr[:,None]*qr[None,:]-(2*co+h[:,None]*h[None,:]/2)**2>1e-6)
    search=Search(valid,finite,olddegree)
    base,V,A,b,ids=face_data(G);ar=-R[:,1:]/R[:,0,None]
    samples=[np.array(x['z'])for x in old['samples']]
    oldface=poly(base,V,A,b,ids)
    center=np.mean([x['z']for x in oldface['vertices']if x['radius2']<1-1e-6],axis=0)
    steps=[];final=None;reason=None
    for it in range(1,101):
        cov=np.array(samples)@ar.T>1+1e-7
        chosen,stats=search.solve(cov)
        if chosen is None:
            reason='finite_bank_compact_degree_cover_infeasible'if stats['complete']else'node_guard'
            steps.append(dict(iteration=it,samples=len(samples),search=stats));break
        assert np.all(cov[:,chosen].any(axis=1))
        assert min(olddegree[chosen]+finite[np.ix_(chosen,chosen)].sum(axis=1))>=3
        mutual=-R[chosen]@M@R[chosen].T/2;np.fill_diagonal(mutual,0)
        assert gate_mask(mutual).all()
        face=poly(base,V,np.r_[A,ar[chosen]@V],np.r_[b,1-ar[chosen]@base],ids+list(range(158,158+len(chosen))))
        assert not face.get('degenerate')
        rec=dict(iteration=it,samples=len(samples),selected=[keep[i]for i in chosen],search=stats,
                 face_counts={k:face[k]for k in ['finite','near_sphere','exterior','maximum_squared_radius']})
        steps.append(rec);print('STEP',it,rec['face_counts'],stats,flush=True)
        if not face['exterior']:
            reason='compact_face_candidate'if not face['near_sphere']else'near_sphere_requires_exact_classification'
            final=dict(roots=[rows[i]for i in chosen],face=face);break
        worst=max(face['vertices'],key=lambda p:p['radius2']);z=np.array(worst['z']);v=z-center
        t=next(float(x)for x in np.roots([v@v,2*center@v,center@center-1])if np.isreal(x)and 0<float(x)<1)
        ideal=center+t*v
        assert abs(ideal@ideal-1)<1e-7 and max(A@((ideal-base)@V)-b)<1e-6
        assert max(ar[chosen]@ideal-1)<1e-6
        assert not any(np.linalg.norm(ideal-x)<1e-7 for x in samples)
        samples.append(ideal)
    else:reason='iteration_guard'
    src=[f'scripts/{TAG}.py','scripts/fresh_h20_intrinsic_facet_incidence_20261016.py',CORE,OLD,
         'scripts/fresh_h20_continuous_sample_joint_cover_20261016.py',
         'scripts/fresh_h20_quad_adjacent_h3_closure_20261016.py']
    return dict(status='PASS_NUMERICAL_INCIDENCE_CORE_JOINT_COVER',
       sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in src},
       full_bank=c['candidate_count'],core_count=len(keep),termination=reason,steps=steps,
       samples=[x.tolist()for x in samples],final=final,
       scope='Complete finite-graph clique and sampled-point cover search only when complete=true; selected new facets must simultaneously have at least3 finite neighbors, using degree-core pruning. Floating-point Coxeter and incidence gates; compact target H3 only. No continuous-family or finite-volume H20 nonexistence claim.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    d=run();p=ROOT/f'results/{TAG}.json'
    if a.verify:assert close(d,json.loads(p.read_text()))
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['termination'],len(d['steps']),flush=True)
