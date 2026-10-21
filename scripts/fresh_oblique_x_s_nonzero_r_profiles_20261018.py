#!/usr/bin/env python3
"""All six remaining R profiles at the frozen X-S high12 parents.

Exact inherited finite-duty filters, then numerical norm/native/center
gates and full ledgers for every surviving branch. No label expansion.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import product
from pathlib import Path
import numpy as np
import fresh_oblique_x_s_diagonal_native_pilot_20261018 as P
ROOT=P.ROOT;T=P.T;C=P.C;Q=P.Q
TAG='fresh_oblique_x_s_nonzero_r_profiles_20261018'

def reconstruct(s,row):
    nd={i:v for i,v in s.old.items()if i not in C.DELETED};tt={};steps=[]
    for j,ch in enumerate(row['reconstruction']):
        i=ch['wall'];v,proof=C.solve_saved(s,ch['frame'],ch['profile'],ch['labels'],ch['branch'],
            s.old if j<2 else nd,{}if j<2 else tt,row['prefix_normals'][str(i)])
        nd[i]=v;tt[i]=ch['profile'];steps.append(dict(wall=i,**proof))
    assert len(row['pairs'])==1
    for key,wall,branchkey in(('N',73,'N_branch'),('S',72,'S_branch')):
        rec=row[key];b=next(b for b in rec['branches']if b['branch']==row['pairs'][0][branchkey]);assert b['admissible']
        v,proof=C.solve_saved(s,rec['frame'],rec['profile'],rec['labels'],b['branch'],nd,tt,b['normal'])
        nd[wall]=v;tt[wall]=rec['profile'];steps.append(dict(wall=wall,**proof))
    assert set(tt)=={67,68,69,72,73} and 70 not in nd and 71 not in nd
    for i in tt:
        _,bad,cp=s.native(nd[i],tt[i],{j:v for j,v in nd.items()if j!=i},tt);assert not bad and cp<-1e-8
    return nd,tt,steps

def run():
    old=json.loads((ROOT/f'results/{P.TAG}.json').read_text())
    for p,h in old['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    assert old['summary']['complete_disks']==8 and old['summary']['unresolved_norm_recipes']==0
    s=T.Scout();s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[];s.counts=Counter();s.blockers=Counter()
    assert P.H.H.old_frame_interface(s)==old['old_frame_interface']
    debit=T.H.deletion(s.primary,s.cur,C.DELETED);assert json.loads(json.dumps(debit))==old['exact_deletion_debit']
    top=P.topology();assert top==old['topology'];ids=[i for i in s.cur['actual_sides']if i not in C.DELETED]
    records=[];accepted=[];statuses=Counter();blockers=Counter()
    for prior in old['records']:
        nd,tt,steps=reconstruct(s,prior)
        row=dict(parent=prior['index'],origins=prior['origins'],reconstruction=steps,
            normals={str(i):nd[i].tolist()for i in tt},profiles=tt,recipes=[]);records.append(row)
        for pi,m in product(range(1,7),range(2,6)):
            duties=[T.P.corner3(0,0,pi,12,2,2),T.P.corner3(0,0,pi,2,2,m)]
            rec=dict(index=len(row['recipes']),profile=pi,frame=[69,73,72],labels=[2,2,m],
                required_corners=[[69,73,70],[69,72,70]],finite_duty_types=duties,norm=None,branches=[])
            row['recipes'].append(rec)
            if any(d is None for d in duties):rec['status']='finite_duty_rejected';statuses[rec['status']]+=1;continue
            rhs=[-T.B.theta(z)-s.off(pi,tt[i])for i,z in zip(rec['frame'],rec['labels'])]
            roots,info=T.N.norm_roots([nd[i]for i in rec['frame']],rhs,s.profiles[pi]['square']);rec['norm']=info
            rec['status']=info['status'];statuses[rec['status']]+=1
            for branch,v in roots:
                vals,bad,cp=s.native(v,pi,nd,tt)
                err=max(abs(float(v@T.M@v)-s.profiles[pi]['square']),max(abs(float(v@T.M@nd[i])-q)for i,q in zip(rec['frame'],rhs)))
                b=dict(branch=branch,normal=v.tolist(),native_products=sorted(vals.items()),blockers=[[i,vals[i]]for i in bad],
                    center_product=cp,norm_residual=err,native_admissible=bool(not bad and cp<-1e-8 and err<2e-5),
                    retained_corner_types=None,all_new_gates=None,accepted_index=None)
                rec['branches'].append(b);blockers[str(tuple(bad))]+=1
                if not b['native_admissible']:continue
                full={**nd,70:v};ft={**tt,70:pi};b['retained_corner_types']=[s.corner(f,full,ft)for f in rec['required_corners']]
                if any(z is None for z in b['retained_corner_types']):continue
                gg=[]
                for i in ft:
                    vv,bb,cc=s.native(full[i],ft[i],{j:w for j,w in full.items()if j!=i},ft)
                    gg.append(dict(wall=i,blockers=[[j,vv[j]]for j in bb],center_product=cc))
                b['all_new_gates']=gg
                if any(g['blockers']or g['center_product']>=-1e-8 for g in gg):continue
                corners=[s.corner(f,full,ft)for f in top['triangles']]
                score=Q.P.B.B.score(s,full,ft,ids+list(ft),'nonzero_R_profile',ids+list(ft))
                raw=T.S.C.detailed_cone(full,list(full))
                faces={f['wall']:f for f in score['whole_faces']}
                exact={str(w):faces[w]['all_rays_finite']and sorted(sorted(q['active'])for q in faces[w]['rays'])==sorted(sorted(f)for f in top['triangles']if w in f)for w in(67,68,69)}
                b['accepted_index']=len(accepted)
                accepted.append(dict(parent=prior['index'],recipe=rec['index'],profile=pi,RS_label=m,branch=branch,
                    corner_types=corners,all_disk_corners=all(z is not None for z in corners),exact_whole_interior_caps=exact,
                    whole_compact_faces={str(w):f['all_rays_finite']for w,f in faces.items()},full_projected_boundary=raw,**score))
        print('PARENT',prior['index'],'recipes',len(row['recipes']),'accepted',sum(a['parent']==prior['index']for a in accepted),flush=True)
    recipes=[r for p in records for r in p['recipes']];branches=[b for r in recipes for b in r['branches']]
    assert len(recipes)==48
    summary=dict(parents=2,origins=4,profiles=list(range(1,7)),RS_labels=list(range(2,6)),recipes=len(recipes),
        status_counts=dict(statuses),branches=len(branches),native_branches=sum(b['native_admissible']for b in branches),
        accepted=len(accepted),complete_disks=sum(a['all_disk_corners']and all(a['exact_whole_interior_caps'].values())for a in accepted),
        boundary_patterns=dict(Counter(str((a['score']is not None,a['finite_rays'],a['open_ends']))for a in accepted)),
        unresolved_norm_recipes=sum(r['status'].startswith('unresolved')for r in recipes),blocker_histogram=dict(blockers))
    paths={Path(__file__)}|{ROOT/p for p in old['sources']}
    paths.update(ROOT/f'{folder}/{P.TAG}.{ext}'for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    return dict(status='BOUNDED_X_S_NONZERO_R_PROFILE_COMPLETION',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in sorted(paths)},
        predecessor=P.TAG,deleted_original_walls=list(C.DELETED),exact_deletion_debit=debit,topology=top,
        old_frame_interface=old['old_frame_interface'],wall_count_breakdown=old['wall_count_breakdown'],
        numerical_snapshot=old['numerical_snapshot'],omitted_timelike_constraints=old['omitted_timelike_constraints'],
        records=records,accepted=accepted,boundary_ledgers=s.ledgers,unresolved_cones=s.unresolved,summary=summary,
        scope='All six nonzero R profiles and all four RS labels2..5 on precisely the two frozen high12 X-S parents, including every finite-duty rejection and both roots of every eligible norm problem. Original52 projected plus6white and all fixed-new native/center gates precede full cones for survivors. Exact inherited old19/12/1 debit; new normals and cones numerical. No all-high-label, all-parent or ambient finite-volume claim.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';assert a.verify or not path.exists();out=json.loads(json.dumps(run()))
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['summary'],flush=True)
