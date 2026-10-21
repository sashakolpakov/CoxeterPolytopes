#!/usr/bin/env python3
"""Changed four-side B-region closure at moving C(u) companion events.

One-cap and two-cap duties come from the complete numerical boundary.
Only high labels6/8/12/24 are sampled; no exact/no-tail claim is made.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
import fresh_cu_companion_contact_events_20261017 as U
import fresh_c3_companion_continuation_20261017 as E
C,H,B,M,N,S,P,ROOT=U.C,U.H,U.B,U.M,U.N,U.S,U.P,U.ROOT
TAG='fresh_cu_four_cycle_closure_20261017'
HIGH=[6,8,12,24]
FRAMES=[([30,52,47],[(0,1),(1,2)]),([30,52,68],[(0,1),(0,2)]),
        ([47,52,68],[(0,1),(0,2)]),([30,47,68],[(0,2),(1,2)])]

def run():
    src=ROOT/'results'/f'{U.TAG}.json';parent=json.loads(src.read_text())
    cur,labels,base,h,profiles,primary=C.setup()
    seeds=[(j,r)for j,r in enumerate(parent['coherent_events'])if r['boundary']['open_ends']==10]
    assert len(seeds)==17
    counts=Counter();blockers=Counter();records=[];unresolved=[];hist=Counter();improved=[]
    cycle=[30,52,47,68]
    off=lambda a,b:float(H.A.native(np.zeros(4),a,np.zeros(4),b,profiles))
    for seed_index,seed in seeds:
        fam=parent['moving_families'][seed['family']];nd=dict(base)
        nd[67]=np.array(seed['C_normal']);nd[68]=np.array(seed['normal'])
        pp={7:4,68:fam['profile']};ps=lambda i:pp.get(i,0)
        expected={tuple(sorted((i,j)))for i,j in zip(cycle,cycle[1:]+cycle[:1])}
        assert expected<={tuple(e['sides'])for e in seed['boundary']['edges']if e['kind']=='one_open'}
        def native(v,pi):
            return {i:float(v@M@nd[i])+(profiles[pi]['offsets'][labels.index(i)]if i in labels else off(pi,ps(i)))for i in nd}
        def score(extra):
            nn=dict(nd);nn.update({69+i:r['vector']for i,r in enumerate(extra)});sel=cur['actual_sides']+[67,68]+list(range(69,69+len(extra)))
            bd=S.C.detailed_cone(nn,sel)
            ic=H.incidence_score([nn[i]for i in sel],sel,[float(nn[i]@M@nn[i])for i in sel])
            good=not bd['warnings']and ic['status']=='numerically_separated_acute_incidence'
            good=good and bd['open_ends']==ic['open_ends']and bd['counts'].get('finite',0)==ic['finite']
            good=good and {tuple(sorted(r['active']))for r in bd['rays']if r['kind']=='finite'}=={tuple(sorted(r))for r in ic['finite_supports']}
            return dict(boundary=bd,PD_incidence=ic),good
        accepted=[];seen=set();domains=[]
        for ids,duties in FRAMES:
            bonds=[2]*9
            for a,b in duties:
                rel=B.relation(float(nd[ids[a]]@M@nd[ids[b]])+off(ps(ids[a]),ps(ids[b])))
                assert isinstance(rel,int);bonds[3*a+b]=rel
            domain=E.recipes(tuple(ps(i)for i in ids),tuple(bonds),tuple(duties))
            domains.append(dict(frame=ids,duties=duties,recipes=len(domain),high_recipes=sum(bool(z[2])for z in domain)))
            for pi,low,high in domain:
                for hh in product(HIGH,repeat=len(high)):
                    ms=list(low)
                    for i,m in zip(high,hh):ms[i]=m
                    rhs=[-B.theta(m)-off(pi,ps(i))for i,m in zip(ids,ms)]
                    roots,info=N.norm_roots([nd[i]for i in ids],rhs,profiles[pi]['square']);counts['norm_equations']+=1
                    if info['status'].startswith('unresolved'):
                        unresolved.append(dict(seed=seed_index,frame=ids,profile=pi,labels=ms,norm=info));continue
                    for sg,v in roots:
                        counts['norm_branches']+=1
                        err=max(abs(float(v@M@v)-profiles[pi]['square']),max(abs(float(v@M@nd[i])-r)for i,r in zip(ids,rhs)))
                        if err>2e-5:
                            unresolved.append(dict(seed=seed_index,frame=ids,profile=pi,labels=ms,branch=sg,stage='residual',error=err));continue
                        vals=native(v,pi);bad=[i for i,x in vals.items()if B.relation(x)in[None,'parallel_or_unresolved']]
                        if bad:blockers[tuple(bad)]+=1;continue
                        cp=float(v@M@h)+profiles[pi]['center_offset']
                        if cp>=-1e-8:counts['center_rejected']+=1;continue
                        key=(pi,)+tuple(np.round(v,8))
                        if key in seen:counts['duplicates']+=1;continue
                        seen.add(key)
                        finite=[]
                        for a,b in zip(cycle,cycle[1:]+cycle[:1]):
                            ab=B.relation(float(nd[a]@M@nd[b])+off(ps(a),ps(b)))
                            va,vb=B.relation(vals[a]),B.relation(vals[b])
                            finite.append(isinstance(va,int)and isinstance(vb,int)and P.corner3(ps(a),ps(b),pi,ab,va,vb)is not None)
                        r=dict(index=len(accepted),frame=ids,profile=pi,labels=ms,branch=sg,normal=v.tolist(),center_product=cp,norm_residual=err,
                          native_products=sorted(vals.items()),four_required_finite_corners=finite)
                        sc,good=score([dict(vector=v)])
                        if not good:unresolved.append(dict(seed=seed_index,stage='cone',**sc));continue
                        r.update(sc);accepted.append(r);counts['coherent_single_caps']+=1
                        bd=sc['boundary'];hist[(bd['counts']['finite'],bd['counts']['exterior'],bd['open_ends'])]+=1
                        if all(finite):
                            active=[z for z in bd['rays']if 69 in z['active']]
                            assert len(active)==4 and all(z['kind']=='finite'for z in active)
                            counts['whole_quadrilateral_closures']+=1
                        if bd['open_ends']<10:improved.append(dict(seed=seed_index,u=seed['u'],cap=r))
        pairs=[]
        for a,b in combinations(accepted,2):
            counts['cap_pairs_tested']+=1
            va,vb=np.array(a['normal']),np.array(b['normal'])
            ab=float(H.A.native(va,a['profile'],vb,b['profile'],profiles));rel=B.relation(ab)
            if rel in[None,'parallel_or_unresolved']:continue
            sc,good=score([dict(vector=va),dict(vector=vb)])
            if not good:unresolved.append(dict(seed=seed_index,stage='pair_cone',indices=[a['index'],b['index']],**sc));continue
            r=dict(indices=[a['index'],b['index']],native_product=ab,relation=rel,**sc);pairs.append(r);counts['coherent_cap_pairs']+=1
            bd=sc['boundary'];hist[(bd['counts']['finite'],bd['counts']['exterior'],bd['open_ends'])]+=1
            if bd['open_ends']<10:improved.append(dict(seed=seed_index,u=seed['u'],pair=r))
        records.append(dict(seed=seed_index,u=seed['u'],companion_profile=fam['profile'],companion_labels=fam['labels'],event=seed['origin'],domains=domains,caps=accepted,pairs=pairs))
        print('CU_FOUR_CYCLE',seed_index,seed['u'],'caps',len(accepted),'pairs',len(pairs),'improved',len(improved),'unresolved',len(unresolved),flush=True)
    paths=[Path(__file__),Path(U.__file__),src,Path(E.__file__),Path(P.__file__),H.PARENT,H.BOUNDARY]
    return dict(status='BOUNDED_CHANGED_CU_FOUR_CYCLE_CLOSURE',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
      retained_control_walls=63,deleted_walls=[],seed_count=len(seeds),boundary_cycle=cycle,high_label_representatives=HIGH,
      counts=dict(counts),records=records,improvements=improved,unresolved=unresolved,
      full_cone_histogram=[dict(finite=f,exterior=e,open_ends=i,count=n)for(f,e,i),n in sorted(hist.items())],
      blockers=[dict(walls=list(k),count=n)for k,n in sorted(blockers.items())],
      scope='Numerical whole-region continuation at17 changed-u native-coherent contact events. All original63walls,C(u),A are retained. Complete finite A2-profile classifications derive each three-contact frame; symbolic high labels are only sampled6/8/12/24. Every accepted cap/pair passes all native products, strict center, full target cone and independent PD incidence. No exact geometry, all-label exclusion, ambient cone or finite-volume claim. Deletion debit zero; old199finite vertices/compact faces persist conditionally in an exact acute realization.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],'improvements',len(out['improvements']),'unresolved',len(out['unresolved']),flush=True)
