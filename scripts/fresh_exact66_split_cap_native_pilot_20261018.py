#!/usr/bin/env python3
"""Two independently solved caps on the fixed exact66 seed; numerical pilot."""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations,product
import numpy as np
import fresh_exact65_escape33_cap_20261018 as P
ROOT=P.ROOT;D=P.D;TAG='fresh_exact66_split_cap_native_pilot_20261018'
K,L=78,79
TRIANGLES=[[30,52,K],[52,47,K],[47,K,L],[47,76,L],[76,30,L],[30,K,L]]

def clean(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    if isinstance(x,dict):return {str(k):clean(v)for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [clean(v)for v in x]
    return x

def run(preflight=False):
    prior=json.loads((ROOT/f'results/{P.TAG}.json').read_text())
    for path,h in prior['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==h,path
    assert prior['physical']['summary']['total_walls']==66
    assert prior['physical']['summary']['complete_finite_vertices']==205
    s=D.T.Scout();D.P.H.H.old_frame_interface(s)
    nd={int(i):np.array(v)for i,v in prior['numerical']['snapshot']['old_normals'].items()}
    chosen=[r for r in prior['numerical']['accepted']if r['profile']==0 and r['labels']==[2,2,2]]
    assert len(chosen)==1;nd[77]=np.array(chosen[0]['normal'])
    tt={75:0,76:0,77:0};assert len(nd)==60 and all(i not in nd for i in(K,L))
    banks=[];survivors={}
    for new,frame,bonds in[(K,[30,52,47],[4,2]),(L,[47,76,30],[2,2])]:
        assert all(s.ps(i,tt)==0 for i in frame)
        rows=[];accepted=[];blocking=Counter()
        for pi,ms in product(range(7),product(range(2,6),repeat=3)):
            duties=[D.T.P.corner3(0,0,pi,b,ms[j],ms[j+1])for j,b in enumerate(bonds)]
            row=dict(profile=pi,labels=ms,duties=duties,branches=[])
            if any(x is None for x in duties):row['status']='finite_duty_rejected';rows.append(row);continue
            rhs=[-D.T.B.theta(m)for m in ms]
            roots,info=D.T.N.norm_roots([nd[i]for i in frame],rhs,s.profiles[pi]['square'])
            row['norm']=info;row['status']=info['status']
            for sign,v in roots:
                prods,bad,cp=s.native(v,pi,nd,tt)
                native=not bad and cp<0
                branch=dict(branch=sign,normal=v.tolist(),products=sorted(prods.items()),blockers=[[i,prods[i]]for i in bad],center_product=cp,native=native)
                row['branches'].append(branch);blocking.update(bad)
                if native:accepted.append(dict(profile=pi,labels=ms,branch=sign,normal=v.tolist()))
            rows.append(row)
        survivors[new]=accepted
        banks.append(dict(new_wall=new,frame=frame,frame_bonds=bonds,recipes=rows,
          native_survivors=accepted,blocking_histogram=dict(blocking),
          summary=dict(complete_recipes=len(rows),finite_eligible=sum(r['status']!='finite_duty_rejected'for r in rows),
            norm_branches=sum(len(r['branches'])for r in rows),native_survivors=len(accepted))))
        print('BANK',new,banks[-1]['summary'],flush=True)
    joins=[];cones=[]
    for ki,krow in enumerate(survivors[K]):
        kv=np.array(krow['normal']);kt=krow['profile']
        for li,lrow in enumerate(survivors[L]):
            lv=np.array(lrow['normal']);lt=lrow['profile']
            vals,bad,cp=s.native(lv,lt,{K:kv},{K:kt})
            full={**nd,K:kv,L:lv};types={**tt,K:kt,L:lt}
            duties=[s.corner(t,full,types)for t in TRIANGLES]
            accepted=not bad and cp<0 and all(t is not None for t in duties)
            row=dict(K_index=ki,L_index=li,mutual_product=vals[K],mutual_bad=bool(bad),duties=duties,accepted=accepted);joins.append(row)
            if accepted and not preflight:
                cone=D.T.S.C.detailed_cone(full,list(full))
                faces={str(i):[r for r in cone['rays']if i in r['active']]for i in[7,30,33,47,52,60,75,76,77,K,L]}
                cones.append(dict(K_index=ki,L_index=li,cone=cone,faces=faces))
                print('CONE',ki,li,len(cone['rays']),cone['open_ends'],flush=True)
    summary=dict(K_native=len(survivors[K]),L_native=len(survivors[L]),Cartesian_pairs=len(joins),
      accepted_joins=sum(r['accepted']for r in joins),cones_checked=len(cones))
    print('JOIN',summary,flush=True)
    if preflight:return
    sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in[f'scripts/{TAG}.py',
      f'scripts/{P.TAG}.py',f'results/{P.TAG}.json',f'results/{P.TAG}.md',
      'results/fresh_exact66_quadrilateral_cap_obstruction_20261018.md']}
    return clean(dict(status='PASS_BOUNDED_EXACT66_SPLIT_CAP_NUMERICAL_PILOT',sources=sources,
      fixed_seed=dict(stem=P.TAG,exact_facets=66,exact_finite_vertices=205,
        wall77_control=chosen[0],projected_normals=nd,profile_types=tt,
        profiles=s.profiles,metric=D.T.M,center=s.center),
      topology=dict(triangles=TRIANGLES,boundary=[30,52,47,76],vertices=6,edges=11,faces=6),
      banks=banks,joins=joins,cones=cones,summary=summary,
      deleted_original_walls=[],old_finite_vertex_debit=0,
      scope=dict(numerical_new_normals=True,profiles=list(range(7)),labels=[2,3,4,5],
        two_independent_banks=True,all_declared_pairs_checked=True,
        high_label_families_untested=True,exact_new_normal_certified=False,new_finite_volume_example=False)))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--preflight',action='store_true');ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';assert a.preflight or a.verify or not path.exists()
    out=run(a.preflight)
    if out is not None:
        if a.verify:assert out==json.loads(path.read_text())
        else:path.write_text(json.dumps(out,indent=2)+'\n')
        print(out['status'],flush=True)
