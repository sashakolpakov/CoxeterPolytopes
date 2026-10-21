#!/usr/bin/env python3
"""Changed X-S diagonal, one declared high-N/high-S pair, all low R-S duties.

Numerical native-first construction pilot. No complete high-label or
continuous family claim, exact new normal, or ambient-volume conclusion.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import numpy as np
import fresh_oblique_high_n_right_r_join_20261018 as H
ROOT=H.ROOT; T=H.T; C=H.C; F=H.F; Q=H.Q
TAG='fresh_oblique_x_s_diagonal_native_pilot_20261018'
N_LABELS=[3,2,12]
S_LABELS=[12,2,2]
LOW=(2,3,4,5)


def topology():
    old=F.topology();remove={tuple(sorted(x))for x in[[68,69,71],[69,70,71],[70,71,72],[68,71,72]]}
    tri=[list(x)for x in old['triangles']if tuple(sorted(x))not in remove]+[[68,69,72],[69,70,72]]
    vs=sorted({x for f in tri for x in f});edges=Counter(tuple(sorted(e))for f in tri for e in combinations(f,2))
    boundary=[8,10,73,70,72,52];be={tuple(sorted(e))for e in zip(boundary,boundary[1:]+boundary[:1])}
    assert(len(vs),len(edges),len(tri))==(9,18,10)and 71 not in vs
    assert {e for e,c in edges.items()if c==1}==be and set(edges.values())=={1,2}
    links={}
    for v in vs:
        le=[sorted(set(f)-{v})for f in tri if v in f];adj={i:[]for e in le for i in e}
        for a,b in le:adj[a].append(b);adj[b].append(a)
        seen=set();todo=[next(iter(adj))]
        while todo:
            x=todo.pop()
            if x in seen:continue
            seen.add(x);todo+=adj[x]
        assert seen==set(adj)
        degree=sorted(map(len,adj.values()))
        assert degree==([1,1]+[2]*(len(degree)-2)if v in boundary else[2]*len(degree))
        links[str(v)]=le
    return dict(triangles=tri,vertices=vs,boundary=boundary,links=links,V=9,E=18,F=10,euler=1)


def bank(s,nd,tt,wall,frame,labels):
    pi=0;rhs=[-T.B.theta(m)for m in labels]
    roots,info=T.N.norm_roots([nd[i]for i in frame],rhs,2.)
    row=dict(wall=wall,profile=pi,frame=frame,labels=labels,norm=info,branches=[])
    for branch,v in roots:
        vals,bad,cp=s.native(v,pi,nd,tt)
        err=max(abs(float(v@T.M@v)-2),max(abs(float(v@T.M@nd[i])-r)for i,r in zip(frame,rhs)))
        row['branches'].append(dict(branch=branch,normal=v.tolist(),native_products=sorted(vals.items()),
            blockers=[[i,vals[i]]for i in bad],center_product=cp,norm_residual=err,
            admissible=bool(not bad and cp<-1e-8 and err<2e-5)))
    return row


def run():
    prior=json.loads((ROOT/f'results/{H.TAG}.json').read_text())
    for p,h in prior['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    s=T.Scout();s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[];s.counts=Counter();s.blockers=Counter()
    interface=H.H.old_frame_interface(s);assert interface==prior['old_frame_interface']
    debit=T.H.deletion(s.primary,s.cur,C.DELETED);assert json.loads(json.dumps(debit))==prior['exact_deletion_debit']
    top=topology();oldids=[i for i in s.cur['actual_sides']if i not in C.DELETED]
    assert len(oldids)==13 and len(set(s.old)-set(C.DELETED))==52
    omitted=[i for i in s.old if i not in C.DELETED and i not in oldids]
    assert len(omitted)==39 and float(s.center@T.M@s.center)<0
    omitted_checks=[dict(wall=i,projected_square=float(s.old[i]@T.M@s.old[i]),
        projected_center_product=float(s.old[i]@T.M@s.center))for i in omitted]
    assert all(r['projected_square']<-1e-5 and r['projected_center_product']<-1e-5 for r in omitted_checks)
    control=Q.P.B.B.score(s,s.old,{},s.cur['actual_sides'],'original63_control',s.cur['actual_sides'])
    assert (control['finite_rays'],control['open_ends'])==(26,10)
    records=[];accepted=[]
    for old in prior['records']:
        nd={i:v for i,v in s.old.items()if i not in C.DELETED};tt={};steps=[]
        for j,ch in enumerate(old['reconstruction'][:3]):
            w=ch['wall'];v,proof=C.solve_saved(s,ch['frame'],ch['profile'],ch['labels'],ch['branch'],
                s.old if j<2 else nd,{}if j<2 else tt,old['prefix_normals'][str(w)])
            nd[w]=v;tt[w]=ch['profile'];steps.append(dict(wall=w,**proof))
        assert set(tt)=={67,68,69}and all(i not in nd for i in(70,71,72,73))
        assert s.corner([68,69,52],nd,tt)is None #52-X is not a required edge.
        row=dict(index=old['index'],origins=old['origins'],reconstruction=steps,
            prefix_normals={str(i):nd[i].tolist()for i in tt},prefix_profiles=tt,
            N=bank(s,nd,tt,73,[10,67,69],N_LABELS),S=bank(s,nd,tt,72,[52,68,69],S_LABELS),pairs=[])
        records.append(row)
        for n in row['N']['branches']:
            if not n['admissible']:continue
            for z in row['S']['branches']:
                if not z['admissible']:continue
                pair=dict(N_branch=n['branch'],S_branch=z['branch'],R=[])
                row['pairs'].append(pair)
                base={**nd,73:np.array(n['normal']),72:np.array(z['normal'])};bt={**tt,73:0,72:0}
                sn=float(base[72]@T.M@base[73]);pair['SN_product']=sn;pair['SN_relation']=T.G.relation(sn)
                gates=[]
                for i in bt:
                    vv,bb,cp=s.native(base[i],bt[i],{j:v for j,v in base.items()if j!=i},bt)
                    gates.append(dict(wall=i,blockers=[[j,vv[j]]for j in bb],center_product=cp))
                pair['prefix_all_new_gates']=gates
                if any(r['blockers']or r['center_product']>=-1e-8 for r in gates):continue
                for m in LOW:
                    assert T.P.corner3(0,0,0,12,2,2)is not None and T.P.corner3(0,0,0,2,2,m)is not None
                    rr=bank(s,base,bt,70,[69,73,72],[2,2,m]);pair['R'].append(rr)
                    for branch in rr['branches']:
                        if not branch['admissible']:continue
                        full={**base,70:np.array(branch['normal'])};ft={**bt,70:0}
                        allg=[]
                        for i in ft:
                            vv,bb,cp=s.native(full[i],ft[i],{j:v for j,v in full.items()if j!=i},ft)
                            allg.append(dict(wall=i,blockers=[[j,vv[j]]for j in bb],center_product=cp))
                        assert all(not r['blockers']and r['center_product']<-1e-8 for r in allg)
                        cs=[s.corner(f,full,ft)for f in top['triangles']]
                        score=Q.P.B.B.score(s,full,ft,oldids+list(ft),'X_S_diagonal_native_R',oldids+list(ft))
                        full_boundary=T.S.C.detailed_cone(full,list(full))
                        assert not full_boundary['warnings']
                        assert full_boundary['counts'].get('finite',0)==score['finite_rays']
                        assert full_boundary['open_ends']==score['open_ends']
                        bd=s.ledgers[score['score']]if score['score']is not None else score['cone_diagnostic'][0]['boundary']
                        assert {tuple(sorted(r['active']))for r in full_boundary['rays']if r['kind']=='finite'}=={tuple(sorted(r['active']))for r in bd['rays']if r['kind']=='finite'}
                        faces={f['wall']:f for f in score['whole_faces']}
                        exact={str(w):faces[w]['all_rays_finite']and sorted(sorted(x['active'])for x in faces[w]['rays'])==sorted(sorted(f)for f in top['triangles']if w in f)for w in(67,68,69)}
                        accepted.append(dict(parent=old['index'],N_branch=n['branch'],S_branch=z['branch'],RS_label=m,
                            R_branch=branch['branch'],all_new_gates=allg,corner_types=cs,all_disk_corners=all(x is not None for x in cs),
                            exact_whole_interior_caps=exact,whole_compact_faces={str(w):f['all_rays_finite']for w,f in faces.items()},
                            full_projected_boundary=full_boundary,**score))
                print('PARENT',old['index'],'native_NS',pair['SN_relation'],'Rbanks',len(pair['R']),flush=True)
    rs=[b for r in records for p in r['pairs']for b in p['R']]
    rbranches=[b for r in rs for b in r['branches']]
    summary=dict(parents=len(records),origins=sum(len(r['origins'])for r in records),N_recipes=2,S_recipes=2,
        N_branches=sum(len(r['N']['branches'])for r in records),S_branches=sum(len(r['S']['branches'])for r in records),
        joint_prefixes=sum(bool(p['R'])for r in records for p in r['pairs']),R_recipes=len(rs),R_branches=len(rbranches),
        native_R=sum(b['admissible']for b in rbranches),accepted_scored=len(accepted),
        complete_disks=sum(a['all_disk_corners']and all(a['exact_whole_interior_caps'].values())for a in accepted),
        unresolved_norm_recipes=sum(x['norm']['status'].startswith('unresolved')for r in records for x in[r['N'],r['S']])+sum(x['norm']['status'].startswith('unresolved')for x in rs))
    paths={Path(__file__),Path(H.__file__),Path(C.__file__)}|{ROOT/p for p in prior['sources']}
    paths.update(ROOT/f'{folder}/{H.TAG}.{ext}'for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    return dict(status='BOUNDED_X_S_DIAGONAL_HIGH12_NATIVE_PILOT',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in sorted(paths)},
        N_labels=N_LABELS,S_labels=S_LABELS,RS_labels=list(LOW),topology=top,deleted_original_walls=list(C.DELETED),
        exact_deletion_debit=debit,removed_new_wall=71,wall_count_breakdown=dict(projected_originals=52,fixed_original_whites=6,retained_old_actual_boundary_sides=13),
        omitted_timelike_constraints=omitted_checks,control=control,numerical_snapshot=prior['numerical_snapshot'],
        old_frame_interface=interface,records=records,accepted=accepted,boundary_ledgers=s.ledgers,summary=summary,
        scope='Changed X-S diagonal after deleting new M. Two frozen J/L/X parents, fixed N-X=S-52=12 and N-10=3; all four low R-S labels and both norm roots. Full original/new native and center gates before numerical full-cone/PD/whole-face ledgers. No continuous/high-label completeness, exact new normal, full H9 certificate or volume claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';assert args.verify or not path.exists()
    out=json.loads(json.dumps(run()))
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['summary'],flush=True)
