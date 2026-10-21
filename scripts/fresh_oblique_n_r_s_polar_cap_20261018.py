#!/usr/bin/env python3
"""Bounded orthogonal cap at the persistent N/R/S exterior region.

All twelve frozen native X-S configurations and seven white profiles.
Only genuine NRT/SRT corners are required; N/S are dotted. Native and
strict common-center gates precede any numerical cone computation.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import numpy as np
import fresh_oblique_x_s_nonzero_r_profiles_20261018 as E
ROOT=E.ROOT;P=E.P;T=E.T;C=E.C;Q=E.Q
TAG='fresh_oblique_n_r_s_polar_cap_20261018'
WALL=74

def topology():
    tri=P.topology()['triangles']+[[73,70,74],[72,70,74]]
    vs=sorted(set(sum(tri,[])));ee=Counter(tuple(sorted(e))for f in tri for e in combinations(f,2))
    boundary=[8,10,73,74,72,52]
    assert(len(vs),len(ee),len(tri))==(10,21,12)
    assert{e for e,n in ee.items()if n==1}=={tuple(sorted(e))for e in zip(boundary,boundary[1:]+boundary[:1])}
    assert set(ee.values())=={1,2}and(72,73)not in ee
    links={}
    for v in vs:
        le=[sorted(set(f)-{v})for f in tri if v in f];deg=Counter(i for e in le for i in e)
        seen={next(iter(deg))}
        while True:
            new=seen|{i for e in le if seen.intersection(e)for i in e}
            if new==seen:break
            seen=new
        assert seen==set(deg)
        assert sorted(deg.values())==([1,1]+[2]*(len(deg)-2)if v in boundary else[2]*len(deg))
        links[str(v)]=le
    return dict(triangles=tri,vertices=vs,boundary=boundary,links=links,V=10,E=21,F=12,euler=1)

def inputs(d0,d1):
    rows=[]
    for source,d in ((P.TAG,d0),(E.TAG,d1)):
        for ai,a in enumerate(d['accepted']):
            parent=next(r for r in d0['records']if r['index']==a['parent'])
            if source==P.TAG:
                pair=parent['pairs'][0]
                rec=next(r for r in pair['R']if r['labels'][2]==a['RS_label'])
                b=next(b for b in rec['branches']if b['branch']==a['R_branch'])
                assert b['admissible'];pi=0
            else:
                r=next(r for r in d['records']if r['parent']==a['parent'])
                rec=r['recipes'][a['recipe']]
                b=next(b for b in rec['branches']if b['branch']==a['branch'])
                assert b['native_admissible'];pi=a['profile']
            rows.append((source,ai,parent,rec,b,pi,a))
    assert len(rows)==12 and len({(p['index'],pi,r['labels'][2],b['branch'])for _,_,p,r,b,pi,_ in rows})==12
    return rows

def run():
    d0=json.loads((ROOT/f'results/{P.TAG}.json').read_text())
    d1=json.loads((ROOT/f'results/{E.TAG}.json').read_text())
    for d in(d0,d1):
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    s=T.Scout();s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[];s.counts=Counter();s.blockers=Counter()
    assert P.H.H.old_frame_interface(s)==d0['old_frame_interface']==d1['old_frame_interface']
    debit=T.H.deletion(s.primary,s.cur,C.DELETED)
    assert json.loads(json.dumps(debit))==d0['exact_deletion_debit']==d1['exact_deletion_debit']
    ids=[i for i in s.cur['actual_sides']if i not in C.DELETED]
    top=topology()
    records=[];accepted=[];counts=Counter();badsets=Counter()
    for source,ai,parent,rr,rb,rpi,oldcap in inputs(d0,d1):
        nd,tt,steps=E.reconstruct(s,parent)
        v,proof=C.solve_saved(s,rr['frame'],rpi,rr['labels'],rb['branch'],nd,tt,rb['normal'])
        nd[70]=v;tt[70]=rpi;steps.append(dict(wall=70,**proof))
        assert set(tt)=={67,68,69,70,72,73} and WALL not in nd and 71 not in nd
        oldg=[]
        for i in tt:
            vv,bb,cc=s.native(nd[i],tt[i],{j:w for j,w in nd.items()if j!=i},tt)
            assert not bb and cc<-1e-8
            oldg.append(dict(wall=i,blockers=[],center_product=cc))
        row=dict(index=len(records),source=source,source_accepted_index=ai,parent=parent['index'],origins=parent['origins'],
            R_profile=rpi,RS_label=rr['labels'][2],R_branch=rb['branch'],reconstruction=steps,
            normals={str(i):nd[i].tolist()for i in tt},profiles=tt,old_new_gates=oldg,
            baseline_finite_rays=oldcap['finite_rays'],baseline_open_ends=oldcap['open_ends'],recipes=[])
        records.append(row)
        for pi in range(7):
            duties=[T.P.corner3(0,rpi,pi,2,2,2),T.P.corner3(0,rpi,pi,rr['labels'][2],2,2)]
            rec=dict(index=pi,wall=WALL,profile=pi,frame=[73,70,72],labels=[2,2,2],
                required_corners=[[73,70,WALL],[72,70,WALL]],finite_duty_types=duties,norm=None,branches=[])
            row['recipes'].append(rec)
            if any(t is None for t in duties):
                rec['status']='finite_duty_rejected';counts[rec['status']]+=1;continue
            rhs=[-s.off(pi,tt[i])for i in rec['frame']]
            roots,info=T.N.norm_roots([nd[i]for i in rec['frame']],rhs,s.profiles[pi]['square'])
            rec['norm']=info;rec['status']=info['status'];counts[rec['status']]+=1
            for branch,v in roots:
                vals,bad,cp=s.native(v,pi,nd,tt)
                err=max(abs(float(v@T.M@v)-s.profiles[pi]['square']),max(abs(float(v@T.M@nd[i])-q)for i,q in zip(rec['frame'],rhs)))
                b=dict(branch=branch,normal=v.tolist(),native_products=sorted(vals.items()),blockers=[[i,vals[i]]for i in bad],
                    center_product=cp,norm_residual=err,native_admissible=bool(not bad and cp<-1e-8 and err<2e-5),
                    retained_corner_types=None,all_new_gates=None,accepted_index=None)
                rec['branches'].append(b);badsets[str(tuple(bad))]+=1
                if not b['native_admissible']:continue
                full={**nd,WALL:v};ft={**tt,WALL:pi}
                b['retained_corner_types']=[s.corner(f,full,ft)for f in rec['required_corners']]
                if any(t is None for t in b['retained_corner_types']):continue
                gates=[]
                for i in ft:
                    vv,bb,cc=s.native(full[i],ft[i],{j:w for j,w in full.items()if j!=i},ft)
                    gates.append(dict(wall=i,blockers=[[j,vv[j]]for j in bb],center_product=cc))
                b['all_new_gates']=gates
                if any(g['blockers']or g['center_product']>=-1e-8 for g in gates):continue
                score=Q.P.B.B.score(s,full,ft,ids+list(ft),'N_R_S_orthogonal_cap',ids+list(ft))
                raw=T.S.C.detailed_cone(full,list(full))
                before={tuple(sorted(r['active']))for r in oldcap['full_projected_boundary']['rays']if r['kind']=='finite'}
                after={tuple(sorted(r['active']))for r in raw['rays']if r['kind']=='finite'}
                oldext={tuple(sorted(r['active']))for r in oldcap['full_projected_boundary']['rays']if r['kind']=='exterior'}
                newext={tuple(sorted(r['active']))for r in raw['rays']if r['kind']=='exterior'}
                faces={f['wall']:f for f in score['whole_faces']}
                exact={str(w):faces[w]['all_rays_finite']and sorted(sorted(q['active'])for q in faces[w]['rays'])==sorted(sorted(f)for f in top['triangles']if w in f)for w in(67,68,69,70)}
                b['accepted_index']=len(accepted)
                accepted.append(dict(configuration=row['index'],profile=pi,branch=branch,
                    old_disk_corner_types=[s.corner(f,full,ft)for f in d0['topology']['triangles']],
                    complete_disk_corner_types=[s.corner(f,full,ft)for f in top['triangles']],exact_whole_interior_caps=exact,
                    full_projected_boundary=raw,finite_supports_removed=sorted(before-after),finite_supports_added=sorted(after-before),
                    exterior_supports_removed=sorted(oldext-newext),exterior_supports_added=sorted(newext-oldext),
                    target_exterior_remaining={str(list(f)):f in newext for f in((70,72,73),(52,72,73))},**score))
        print('CONFIGURATION',row['index'],'Rprofile',rpi,'RS',rr['labels'][2],'accepted',sum(a['configuration']==row['index']for a in accepted),flush=True)
    recipes=[r for p in records for r in p['recipes']];branches=[b for r in recipes for b in r['branches']]
    summary=dict(configurations=len(records),profiles=list(range(7)),recipes=len(recipes),status_counts=dict(counts),
        branches=len(branches),native_branches=sum(b['native_admissible']for b in branches),accepted=len(accepted),
        unresolved_norm_recipes=sum(r['status'].startswith('unresolved')for r in recipes),blocker_histogram=dict(badsets),
        cone_patterns=dict(Counter(str((a['score']is not None,a['finite_rays'],a['open_ends']))for a in accepted)))
    assert len(records)==12 and len(recipes)==84
    assert len(accepted)==28 and all(a['score']is not None and a['finite_rays']==25 and a['open_ends']==13 for a in accepted)
    assert all(not a['finite_supports_removed']and all(a['exact_whole_interior_caps'].values())and all(a['complete_disk_corner_types'])for a in accepted)
    paths={Path(__file__)}|{ROOT/p for d in(d0,d1)for p in d['sources']}
    paths.update(ROOT/f'{folder}/{stem}.{ext}'for stem in(P.TAG,E.TAG)for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    return dict(status='BOUNDED_N_R_S_ORTHOGONAL_POLAR_CAP',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in sorted(paths)},
        predecessors=[P.TAG,E.TAG],new_wall=WALL,deleted_original_walls=list(C.DELETED),exact_deletion_debit=debit,
        inherited_topology=d0['topology'],topology=top,old_frame_interface=d0['old_frame_interface'],numerical_snapshot=d0['numerical_snapshot'],
        wall_count_breakdown=d0['wall_count_breakdown'],omitted_timelike_constraints=d0['omitted_timelike_constraints'],
        records=records,accepted=accepted,boundary_ledgers=s.ledgers,unresolved_cones=s.unresolved,summary=summary,
        scope='All12 frozen native X-S configurations, seven T profiles and both roots in the physical N/R/S right-duty frame. Only genuine NRT/SRT finite duties, no NST corner. Full retained-original/mutual/strict-center gates precede cone and complete whole-face diagnostics. Old19/12/1 debit exact; new normals/cones numerical. No other T labels, arbitrary caps, exact new normals or finite-volume claim.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';assert a.verify or not path.exists();out=json.loads(json.dumps(run()))
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['summary'],flush=True)
