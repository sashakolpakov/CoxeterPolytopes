#!/usr/bin/env python3
"""Resolve high-RM necessary right contacts by a one-dimensional polar.

No RM label cutoff: for each declared N and every frozen low S with SM=2,
R has profile0 and is orthogonal to X,N,S. New normals remain numerical.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from pathlib import Path
import numpy as np
import fresh_oblique_high_n_right_r_join_20261018 as H
ROOT=H.ROOT;T=H.T;C=H.C;Q=H.Q;F=H.F
TAG='fresh_oblique_high_rm_right_polar_20261018'


def run():
    names=(H.TAG,H.D.TAG,H.I.TAG,H.F.TAG)
    docs={n:json.loads((ROOT/f'results/{n}.json').read_text())for n in names}
    for d in docs.values():
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    old=docs[H.TAG];assert old['summary']['native_N']==20 and old['summary']['finite_MRS']==0
    s=T.Scout();s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[];s.counts=Counter();s.blockers=Counter()
    interface=H.H.old_frame_interface(s);assert interface==old['old_frame_interface']
    debit=T.H.deletion(s.primary,s.cur,C.DELETED);assert json.loads(json.dumps(debit))==old['exact_deletion_debit']
    parents=H.parents(s,docs[H.I.TAG],docs[H.F.TAG]);top=F.topology();records=[];accepted=[];hist=Counter()
    for pre,row in zip(parents,old['records']):
        assert pre['index']==row['index']and pre['origins']==row['origins']
        base,bt=pre['normals'],pre['types']
        for nc in row['caps']:
            nv,ncheck=C.solve_saved(s,[10,67,69],0,nc['labels'],nc['branch'],base,bt,nc['normal'])
            nd={**base,73:nv};tt={**bt,73:0}
            selected=[z for z in nc['S_bank']['caps']if z['labels'][2]==2]
            rec=dict(index=len(records),parent=row['index'],parent_origins=pre['origins'],N_index=nc['index'],N_labels=nc['labels'],N_reconstruction=ncheck,
                prefix_normals={str(i):nd[i].tolist()for i in tt},prefix_profiles=tt,
                complete_S_bank_size=len(nc['S_bank']['caps']),selected_S_indices=[z['index']for z in selected],polars=[])
            records.append(rec)
            for sc in selected:
                sv,scheck=C.solve_saved(s,[52,68,71],sc['profile'],sc['labels'],sc['branch'],nd,tt,sc['normal'])
                full={**nd,72:sv};ft={**tt,72:sc['profile']}
                for i in ft:
                    _,bad,cp=s.native(full[i],ft[i],{j:v for j,v in full.items()if j!=i},ft);assert not bad and cp<-1e-8
                assert T.G.relation(s.dot(sv,sc['profile'],full[71],0))==2
                frame=[69,73,72];G=np.array([[float(full[i]@T.M@full[j])for j in frame]for i in frame]);eig=np.linalg.eigvalsh(G)
                inv=np.linalg.inv(G)if min(abs(eig))>1e-8 else None
                polar=dict(S_index=sc['index'],S_profile=sc['profile'],S_labels=sc['labels'],S_branch=sc['branch'],S_normal=sv.tolist(),S_reconstruction=scheck,
                    frame=frame,frame_Gram=G.tolist(),frame_eigenvalues=eig.tolist(),branches=[])
                rec['polars'].append(polar)
                if inv is None or(sum(eig>1e-8),sum(eig< -1e-8))!=(2,1):
                    polar['status']='unresolved_frame';s.unresolved.append(dict(record=rec['index'],S=sc['index'],frame=polar));continue
                # Physical R profile0 is orthogonal to every fixed white, so
                # its three zero products are also projected zero products.
                cov=np.array([float(full[i]@T.M@full[71])for i in frame]);q=float(full[71]@T.M@full[71]-cov@inv@cov)
                assert q>0
                polar.update(frame_inverse=inv.tolist(),M_covector=cov.tolist(),M_complement_square=q,forced_RM_square=2*q)
                roots,info=T.N.norm_roots([full[i]for i in frame],[0.,0.,0.],2.)
                polar['norm']=info;polar['status']=info['status']
                if info['status'].startswith('unresolved'):
                    s.unresolved.append(dict(record=rec['index'],S=sc['index'],norm=info));continue
                for branch,v in roots:
                    vals,bad,cp=s.native(v,0,full,ft);rm=vals[71];rel=T.G.relation(rm)
                    err=max(abs(float(v@T.M@v)-2),max(abs(float(v@T.M@full[i]))for i in frame))
                    assert abs(rm*rm-2*q)<2e-6
                    high=isinstance(rel,int)and rel>=6
                    b=dict(branch=branch,normal=v.tolist(),native_products=sorted(vals.items()),blockers=[[i,vals[i]]for i in bad],center_product=cp,
                        norm_residual=err,RM_product=rm,RM_relation=rel,finite_high_RM=high)
                    polar['branches'].append(b);hist[str(tuple(bad))]+=1
                    if err>2e-5:s.unresolved.append(dict(record=rec['index'],S=sc['index'],branch=branch,error=err));continue
                    if bad or cp>=-1e-8 or not high:continue
                    chamber={**full,70:v};ct={**ft,70:0}
                    assert all(T.G.relation(s.dot(v,0,chamber[i],s.ps(i,ct)))==2 for i in(69,73,72))
                    cs=[s.corner(f,chamber,ct)for f in top['triangles']]
                    score=Q.P.B.B.score(s,chamber,ct,[i for i in s.cur['actual_sides']if i not in C.DELETED]+list(ct),'right_polar_high_RM',[i for i in F.WALLS if i in chamber])
                    faces={f['wall']:f for f in score['whole_faces']}
                    exact={str(w):faces[w]['all_rays_finite']and sorted(sorted(z['active'])for z in faces[w]['rays'])==sorted(sorted(f)for f in top['triangles']if w in f)for w in(67,68,69,71)}
                    accepted.append(dict(record=rec['index'],S=sc['index'],**b,corner_types=cs,all_disk_corners=all(x is not None for x in cs),exact_whole_interior_caps=exact,**score))
            rec['complete_polar_digest']=sha256(json.dumps(rec['polars'],separators=(',',':')).encode()).hexdigest()
            print('N',rec['parent'],nc['labels'],'rightS',len(selected),flush=True)
    pp=[p for r in records for p in r['polars']];bb=[b for p in pp for b in p['branches']];future=[b for b in bb if b['center_product']< -1e-8]
    summary=dict(parent_origins=4,distinct_parents=2,N_geometries=len(records),complete_S_caps=old['summary']['S_caps'],right_S_polars=len(pp),
        branches=len(bb),future_branches=len(future),native_future_branches=sum(not b['blockers']for b in future),finite_high_RM_branches=sum(b['finite_high_RM']for b in bb),
        accepted=len(accepted),new_cones=len(accepted),unresolved=len(s.unresolved),minimum_abs_RM=min(abs(b['RM_product'])for b in bb),maximum_abs_RM=max(abs(b['RM_product'])for b in bb),
        future_blocker_histogram=dict(Counter(str(tuple(i for i,v in b['blockers']))for b in future)),complete_blocker_histogram=dict(hist))
    paths={Path(__file__),Path(H.__file__)};paths.update(ROOT/p for p in old['sources'])
    paths.update(ROOT/f'{folder}/{n}.{ext}'for n in names for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    return dict(status='BOUNDED_HIGH_RM_RIGHT_POLAR_DIAGNOSTIC',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in sorted(paths)},
        deleted_walls=list(C.DELETED),exact_deletion_debit=debit,old_frame_interface=interface,topology=top,numerical_snapshot=old['numerical_snapshot'],
        parent_reconstruction=[dict(index=p['index'],origins=p['origins'],instruction_key=p['instruction_key'],reconstruction=p['reconstruction'],removed=p['removed'])for p in parents],
        records=records,accepted=accepted,summary=summary,boundary_ledgers=s.ledgers,unresolved=s.unresolved,
        necessary_conditions=dict(N_profile=0,R_profile=0,RX=2,RN=2,SM=2,SR=2,RM_label_lower_bound=6,RM_squared_lower_bound=3),
        wall_count_breakdown=dict(projected_originals=52,fixed_original_whites=6,fixed_new_before_R=6,explicit_R_gates=58,total_fixed_facets=64),
        scope='Every frozen low S with SMright on20declared highN geometries; reconstruct physical anchors, solveboth profile0 R polar sheets withRX/RN/RSright. No RMlabel cutoff but allnew numerical frames; no allN/allS or exactambient claim. No branch reaches cone/wholeface stages.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();path=ROOT/f'results/{TAG}.json'
    assert a.verify or not path.exists();out=json.loads(json.dumps(run()))
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['summary'],flush=True)
