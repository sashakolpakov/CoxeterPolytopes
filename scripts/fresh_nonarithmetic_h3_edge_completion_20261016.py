#!/usr/bin/env python3
"""Bounded numerical edge-endpoint continuation of the actual H9 H3 face.

All arithmetic here is discovery arithmetic. No exact geometry or complete
continuous/nonexistence assertion is exported by this scout.
"""
import argparse,json,itertools,math,time
from collections import Counter
from hashlib import sha256
import numpy as np
from scipy.linalg import null_space
import fresh_nonarithmetic_h3_two_contact_20261016 as C
import fresh_nonarithmetic_h3_two_contact_20261016_pairs as B
ROOT=C.ROOT;ETA=C.ETA;TAG='fresh_nonarithmetic_h3_edge_completion_20261016'
LABELS=(2,3,4,5,6,8,10)
def product(m):return 0. if m==2 else -2*math.cos(math.pi/m)

def ledger(rays,N):
    sides=sorted({j for r in rays for j in range(r['mask'].bit_length())if r['mask']>>j&1 and j not in C.S})
    edges=[];ends=[]
    for i,j in itertools.combinations(sides,2):
        ids=[k for k,r in enumerate(rays)if r['mask']>>i&1 and r['mask']>>j&1]
        if len(ids)!=2:continue
        G=N[[i,j]]@ETA@N[[i,j]].T
        if min(np.linalg.eigvalsh(G))<=1e-7:continue
        x,y=[rays[k]['v']for k in ids];v=y-x
        roots=C.roots([v@ETA@v,2*x@ETA@v,x@ETA@x])
        nulls=[x+t*v for t in sorted(set(roots))if 1e-8<t<1-1e-8]
        e=dict(sides=[i,j],ray_indices=ids,endpoint_kinds=[rays[k]['kind']for k in ids],ideal_points=nulls)
        edges.append(e)
        for p in nulls:ends.append(dict(sides=[i,j],point=p))
    return edges,ends,sides

def family_frames(N,pair,labels):
    A=N[list(pair)];G=A@ETA@A.T
    if min(np.linalg.eigvalsh(G))<=1e-7:return[]
    c=np.array([product(m)for m in labels]);p=c@np.linalg.inv(G)@A;rho=2-p@ETA@p
    # Positive triple Gram is mandatory for a finite edge endpoint.
    if rho<=1e-7:return[]
    K=null_space(A@ETA);ev,U=np.linalg.eigh(K.T@ETA@K)
    if ev[0]>=-1e-8 or ev[1]<=1e-8:return[]
    K=K@U@np.diag(1/np.sqrt(abs(ev)));a=K[:,1];b=K[:,0]
    return[dict(kind='hyperbola',p=p,a=s*math.sqrt(rho)*(a+b)/2,
                 b=s*math.sqrt(rho)*(a-b)/2,anchors=list(pair),labels=list(labels),branch=s)for s in(-1,1)]

def finite_gates(N,v):
    return v[0]>1e-7 and abs(v@ETA@v-2)<3e-5 and all(C.compatible(g,4e-5)for g in N@ETA@v)

def proposals(N,rays,caplimit):
    edges,ends,sides=ledger(rays,N);points=np.array([e['point']for e in ends]);pool={};counts=Counter()
    pairs=sorted({tuple(e['sides'])for e in ends})
    def admit(f,t,recipe):
        if not math.isfinite(t)or t<=0:return
        v=C.evaluate(f,t)
        if not finite_gates(N,v):return
        values=points@ETA@v;cover=tuple(np.flatnonzero(values>2e-7).tolist())
        if not cover:return
        # A projected affine circuit is incompatible with a compact target.
        # Projection off the acute PD S6 adds nonnegative S6 coefficients,
        # so an exact positive kernel would also be a persistent ideal point
        # of the ambient chamber, incident to its retained pi5 pair.
        # This numerical near-null filter is conservative discovery only.
        g=N@ETA@v
        for i,j in itertools.combinations([i for i in sides if abs(g[i])<=2+1e-5],2):
            if abs(g[i])<1e-7 and abs(g[j])<1e-7:continue
            M=np.array([[N[i]@ETA@N[i],N[i]@ETA@N[j],g[i]],
                        [N[i]@ETA@N[j],N[j]@ETA@N[j],g[j]],[g[i],g[j],2.]])
            ee=np.linalg.eigvalsh(M)
            if abs(ee[0])<1e-6 and ee[1]>1e-6:
                counts['possible_affine_rejected']+=1;return
        key=tuple(np.round(v,7));score=(len(cover),float(sum(np.maximum(values,0))))
        rec=dict(vector=v,recipe=recipe,ideal_end_coverage=list(cover),score=score)
        if key not in pool:pool[key]=rec
    for pair in pairs:
        for labels in itertools.product(LABELS,repeat=2):
            for f in family_frames(N,pair,labels):
                counts['finite_endpoint_branches']+=1
                fixed=[];polys=[C.poly(f,np.array([-1.,0,0,0]))];bad=False
                for j,n in enumerate(N):
                    if j in C.S or j in pair:continue
                    cov=ETA@n;vary=max(abs(f['a']@cov),abs(f['b']@cov));val=f['p']@cov
                    if vary<1e-8:
                        if not C.compatible(val,4e-5):bad=True;break
                        fixed.append(j)
                    else:polys.append(C.poly(f,cov,2.))
                if not bad:
                    for lo,hi in C.intervals(polys,True):
                        cuts=sorted({t for point in points for t in C.roots(C.poly(f,ETA@point))if t>0 and(lo is None or t>lo+1e-8)and(hi is None or t<hi-1e-8)})
                        cuts=[lo]+cuts+[hi]
                        for a,b in zip(cuts,cuts[1:]):
                            ts=C.parameter_samples([a,b],True)
                            for t in ts[::max(1,len(ts)//3)]:
                                counts['continuous_samples']+=1
                                admit(f,t,dict(kind='two_contact_endpoint',anchors=pair,labels=labels,branch=f['branch'],parameter=t))
                # A third finite contact closes new loops. It need not lie
                # in the open all-other-dotted parameter domain above.
                for third in sides:
                    if third in pair or(max(pair)<55 and third<55):continue
                    for m in LABELS:
                        counts['third_contact_equations']+=1
                        q=C.poly(f,ETA@N[third],-product(m))
                        if max(abs(q))<1e-10:continue
                        for t in C.roots(q):admit(f,t,dict(kind='three_contact_continuation',anchors=list(pair)+[third],labels=list(labels)+[m],branch=f['branch'],parameter=t))
    rows=list(pool.values());counts['compatible_endpoint_candidates']=len(rows)
    # Keep coverage strata, so a single high numerical score cannot discard
    # every proposal for a less exposed old edge.
    groups={}
    for r in rows:groups.setdefault(tuple(r['ideal_end_coverage']),[]).append(r)
    for rr in groups.values():rr.sort(key=lambda r:(-r['score'][0],-r['score'][1],tuple(r['vector'])))
    selected=[];depth=0
    while len(selected)<caplimit:
        added=False
        for key in sorted(groups,key=lambda z:(-len(z),z)):
            if depth<len(groups[key]):selected.append(groups[key][depth]);added=True
            if len(selected)>=caplimit:break
        if not added:break
        depth+=1
    return selected,dict(counts),edges,ends

def measure(rays,N):
    edges,ends,sides=ledger(rays,N);counts=dict(Counter(r['kind']for r in rays))
    radius=max((float(r['v'][1:]@r['v'][1:])for r in rays),default=0.)
    return(len(ends),counts.get('unresolved_null',0),radius),dict(counts=counts,finite_edges=len(edges),open_edge_ends=len(ends),maximum_squared_klein_radius=radius)

def run(levels,width,caplimit):
    path=ROOT/f'results/{C.TAG}.json';d=json.loads(path.read_text());N=np.array(d['old_projected_normals'])
    rays=[dict(mask=sum(1<<i for i in r['active']),v=np.array(r['point']),kind=r['kind'])for r in d['actual_face_rays']]
    score,info=measure(rays,N);initial_edges,initial_ends,_=ledger(rays,N)
    initial_points=np.array([e['point']for e in initial_ends])
    beam=[dict(N=N,rays=rays,caps=[],measure=info,score=score)];history=[];best=[];compact=None;global_seen=set()
    for depth in range(levels):
        children=[]
        for si,state in enumerate(beam):
            props,counts,edges,ends=proposals(state['N'],state['rays'],caplimit)
            row=dict(depth=depth,state=si,starting_measure=state['measure'],counts=counts,selected_for_clipping=len(props));history.append(row)
            print('EDGE_STATE',depth,si,row,flush=True)
            for p in props:
                try:newrays=B.clip(state['rays'],p['vector'],len(state['N']))
                except ArithmeticError:counts['unresolved_clip']=counts.get('unresolved_clip',0)+1;continue
                newN=np.vstack([state['N'],p['vector']]);score,info=measure(newrays,newN)
                signature=tuple(sorted(tuple(np.round(c['vector'],7))for c in state['caps']+[p]))
                if signature in global_seen:continue
                global_seen.add(signature)
                caps=state['caps']+[p]
                covered=np.flatnonzero(np.max(initial_points@ETA@np.array([c['vector']for c in caps]).T,axis=1)>2e-7).tolist()
                info['initial_ends_cut']=covered
                score=(len(initial_ends)-len(covered),)+score
                child=dict(N=newN,rays=newrays,caps=caps,measure=info,score=score);children.append(child)
                if info['counts'].get('hyperideal_rays',0)==0 and info['counts'].get('unresolved_null',0)==0:
                    compact=child;break
            if compact:break
        children.sort(key=lambda s:(s['score'],tuple(tuple(c['vector'])for c in s['caps'])))
        best=children[:max(width,10)];beam=[];per=Counter()
        for s in children:
            key=tuple(s['measure']['initial_ends_cut'])
            if per[key]>=2:continue
            beam.append(s);per[key]+=1
            if len(beam)>=width:break
        print('EDGE_LEVEL',depth,'children',len(children),'best',[s['measure']for s in beam],flush=True)
        if compact or not beam:break
    def record(s):return dict(caps=s['caps'],measure=s['measure'],rays=[dict(active=[j for j in range(r['mask'].bit_length())if r['mask']>>j&1],point=r['v'],kind=r['kind'])for r in s['rays']])
    paths=[f'scripts/{TAG}.py',f'scripts/{C.TAG}.py',f'scripts/{B.TAG}.py',f'results/{C.TAG}.json',f'results/{C.TAG}_affine_gate.json']
    return C.clean(dict(status='BOUNDED_NUMERICAL_H3_EDGE_ENDPOINT_CONTINUATION',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        bounds=dict(levels=levels,beam_width=width,clip_limit_per_state=caplimit,labels=LABELS),initial_measure=info if not history else history[0]['starting_measure'],
        initial_edges=initial_edges,initial_open_ends=initial_ends,history=history,best_states=[record(s)for s in best],compact_proposal=None if compact is None else record(compact),
        scope='Numerical edge-driven bounded continuation. Every two-contact anchor frame and full endpoint triple must be positive definite; all current ambient contact gates remain checked. Third finite contacts may involve new walls, so this domain goes beyond the frozen old4329 bank. Parameter sampling, finite labels, beam truncation, numerical affine rejection and floating clipping are not an exact or exhaustive obstruction. No finite-volume or exact new wall claim.'))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--levels',type=int,default=6);p.add_argument('--width',type=int,default=6);p.add_argument('--caplimit',type=int,default=40);p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run(a.levels,a.width,a.caplimit);out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text());print('PASS saved edge continuation')
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],'compact',bool(d['compact_proposal']),flush=True)
