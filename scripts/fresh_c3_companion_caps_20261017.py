#!/usr/bin/env python3
"""Deletion-free companions for C(3), with complete intrinsic cone scoring.

Finite diagram domains are classified before the bounded numerical search.
High labels use explicit representative values; no sampled tail is excluded.
"""
import argparse,json
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
import fresh_dotted_c10_joint_scout_20261017 as J

H,B,M,N,S,P,ROOT=J.H,J.B,J.M,J.N,J.S,J.P,J.ROOT
TAG='fresh_c3_companion_caps_20261017'
C_ID=67
FRAMES=[([7,47,C_ID],[(0,1),(0,2)]),
        ([7,30,C_ID],[(0,2),(1,2)]),
        ([30,52,C_ID],[(0,1),(0,2)])]
HIGH=[6,8,12,24]

def setup():
    cur,labels,ns,h,profiles,primary=H.setup()
    f=np.array([ns[i]for i in [7,30,10]]);inv=np.linalg.inv(f@M@f.T)
    v=(inv@np.array([0.,0.,-1.]))@f;alpha=float(-v@M@v)
    rr,_=N.norm_roots(f,[0.,0.,-2.],2.)
    k=next((c-2*v)/np.sqrt(1+2*alpha)for _,c in rr if(c-2*v)@M@ns[55]<0)
    c=3*v+np.sqrt(1+4.5*alpha)*k
    ns[C_ID]=c
    return cur,labels,ns,h,profiles,primary

@lru_cache(None)
def finite_corner(i,j,oldbond,pi,mi,mj):
    # The H4 summand is orthogonal. The complete remaining five-node
    # diagram is A2 + two old roots + the proposed new root.
    g=[[2]*5 for _ in range(5)]
    for k in range(5):g[k][k]=1
    g[0][1]=g[1][0]=3
    for col,node in [(2,i),(3,j)]:
        if node==7:g[1][col]=g[col][1]=3
    for k,m in enumerate(N.G.profile_labels(pi)):g[k][4]=g[4][k]=m
    for a,b,m in [(2,3,oldbond),(2,4,mi),(3,4,mj)]:g[a][b]=g[b][a]=m
    return N.A.components(tuple(v for row in g for v in row))

def domains(ns):
    out=[]
    for frame,duties in FRAMES:
        bonds={(a,b):B.relation(float(ns[frame[a]]@M@ns[frame[b]]))for a,b in duties}
        assert all(isinstance(m,int)for m in bonds.values())
        rows=[]
        for pi,ms in product(range(7),product(range(2,7),repeat=3)):
            types=[finite_corner(frame[a],frame[b],bonds[a,b],pi,ms[a],ms[b])for a,b in duties]
            if any(x is None for x in types):continue
            rows.append(dict(profile=pi,labels=list(ms),high_positions=[i for i,m in enumerate(ms)if m==6],corner_types=types))
        out.append(dict(frame=frame,required_corner_pairs=[[frame[a],frame[b]]for a,b in duties],recipes=rows))
    return out

def run():
    cur,labels,ns,h,profiles,primary=setup();ids=labels+[C_ID];actual=cur['actual_sides']+[C_ID]
    domain=domains(ns);counts=Counter();blockers=Counter();unresolved=[];accepted=[];seen=set();best=[]
    def native(v,pi):
        return {i:float(v@M@ns[i]+(profiles[pi]['offsets'][labels.index(i)]if i!=C_ID else 0.))for i in ids}
    assert all(B.relation(float(ns[C_ID]@M@ns[i]))not in[None,'parallel_or_unresolved']for i in labels)
    assert float(ns[C_ID]@M@h)<0
    def score(extra):
        nd=dict(ns);nd.update({68+i:v for i,v in enumerate(extra)});sel=actual+list(range(68,68+len(extra)))
        bd=S.C.detailed_cone(nd,sel)
        ic=H.incidence_score([nd[i]for i in sel],sel,[float(nd[i]@M@nd[i])for i in sel])
        ok=not bd['warnings']and ic['status']=='numerically_separated_acute_incidence'
        ok=ok and bd['open_ends']==ic['open_ends']and bd['counts'].get('finite',0)==ic['finite']
        ok=ok and {tuple(sorted(r['active']))for r in bd['rays']if r['kind']=='finite'}=={tuple(sorted(r))for r in ic['finite_supports']}
        if not ok:return None,dict(boundary=bd,PD_incidence=ic)
        return dict(boundary=bd,PD_incidence=ic),None
    baseline,error=score([]);assert error is None
    assert baseline['boundary']['counts']=={'finite':27,'exterior':7}and baseline['boundary']['open_ends']==11
    for di,d in enumerate(domain):
        frame=d['frame'];f=[ns[i]for i in frame]
        for rec in d['recipes']:
            for hm in product(HIGH,repeat=len(rec['high_positions'])):
                ms=list(rec['labels'])
                for j,m in zip(rec['high_positions'],hm):ms[j]=m
                pi=rec['profile'];q=profiles[pi]['square']
                rhs=[-B.theta(m)-(profiles[pi]['offsets'][labels.index(i)]if i!=C_ID else 0.)for i,m in zip(frame,ms)]
                roots,info=N.norm_roots(f,rhs,q);counts['norm_equations']+=1
                if info['status'].startswith('unresolved'):
                    unresolved.append(dict(frame=frame,profile=pi,labels=ms,norm=info));continue
                for sg,v in roots:
                    counts['norm_branches']+=1
                    err=max(abs(float(v@M@v)-q),max(abs(float(v@M@w)-b)for w,b in zip(f,rhs)))
                    if err>1e-5:unresolved.append(dict(frame=frame,profile=pi,labels=ms,branch=sg,residual=err));continue
                    vals=native(v,pi)
                    bad=[(i,x)for i,x in vals.items()if B.relation(x)in[None,'parallel_or_unresolved']]
                    if bad:
                        blockers[di,tuple(i for i,_ in bad)]+=1
                        if len(bad)<=2:best.append(dict(frame=frame,profile=pi,labels=ms,branch=sg,blockers=bad))
                        continue
                    cp=float(v@M@h)+profiles[pi]['center_offset']
                    if cp>=-1e-8:counts['common_center_rejected']+=1;continue
                    key=(pi,)+tuple(np.round(v,8))
                    if key in seen:counts['duplicates']+=1;continue
                    seen.add(key);sc,error=score([v])
                    if error is not None:
                        unresolved.append(dict(frame=frame,profile=pi,labels=ms,branch=sg,stage='boundary',**error));continue
                    record=dict(index=len(accepted),frame=frame,profile=pi,labels=ms,branch=sg,normal=v.tolist(),
                      norm_residual=err,center_product=cp,native_products=sorted(vals.items()),
                      finite_contacts=[[i,B.relation(x)]for i,x in vals.items()if isinstance(B.relation(x),int)],**sc)
                    accepted.append(record);counts['coherent_companions']+=1
        print('C3_COMPANION_FRAME',frame,'accepted',len(accepted),'unresolved',len(unresolved),flush=True)
    pairs=[]
    for a,b in combinations(accepted,2):
        va,vb=np.array(a['normal']),np.array(b['normal'])
        ab=float(H.A.native(va,a['profile'],vb,b['profile'],profiles));rel=B.relation(ab)
        counts['companion_pairs_tested']+=1
        if rel in[None,'parallel_or_unresolved']:continue
        sc,error=score([va,vb])
        if error is not None:unresolved.append(dict(stage='pair_boundary',indices=[a['index'],b['index']],**error));continue
        pairs.append(dict(indices=[a['index'],b['index']],native_product=ab,relation=rel,**sc))
    paths=[Path(__file__),Path(J.__file__),Path(H.__file__),Path(N.__file__),Path(N.A.__file__),Path(S.C.__file__),
      H.PARENT,H.BOUNDARY,ROOT/'results/fresh_dotted_c10_domain_20261017.json',ROOT/'scripts/fresh_dotted_c10_domain_20261017.py']
    return dict(status='BOUNDED_DELETION_FREE_C3_COMPANION_DISCOVERY',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
      retained_control_walls=63,deleted_walls=[],old_finite_vertices_preserved_by_exact_C3=199,
      companion_preservation_scope='Conditional on exact coherent acute realization; all199 oldfinite vertices and oldcompact faces must persist.',
      C3_normal=ns[C_ID].tolist(),baseline=baseline,finite_diagram_domains=domain,high_label_samples=HIGH,
      counts=dict(counts),blockers=[dict(frame_index=k[0],walls=list(k[1]),count=n)for k,n in sorted(blockers.items())],
      near_native_recipes=best,companions=accepted,coherent_pairs=pairs,unresolved=unresolved,
      scope='All63 original control walls and C(3) are retained, with the old strict common center. Actual S6 facets use all7 inherited A2 profiles. For each mandatory finite-corner incidence the symbolic6 entry denotes an isolated I2(m), m>=6; its four tested values are only a bounded discovery sample. Every accepted wall/pair has all native contacts checked and a complete numerical target cone cross-checked by PD pair/triple incidence. No numerical geometry, all-label exclusion, ambient cone or finite-volume claim is promoted. Oldfinite/compact-face preservation follows conditionally for any exact coherent acute realization from the strict-center projection argument.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],'pairs',len(out['coherent_pairs']),'unresolved',len(out['unresolved']),flush=True)
