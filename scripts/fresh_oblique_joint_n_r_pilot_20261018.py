#!/usr/bin/env python3
"""Four-prefix joint N/R pilot after the failed fixed-R left connector.

N is selected only when its sole native blocker is the removed R. The new
R70 is prescribed from X/M/N with XMR/XNR, then MRS is checked separately.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import product
from pathlib import Path
import numpy as np
import fresh_oblique_left_boundary_connector_20261018 as F
ROOT=F.ROOT;R=F.R;A=F.A;O=F.O;C=F.C;T=F.T;Q=F.Q;H=F.H
TAG='fresh_oblique_joint_n_r_pilot_20261018';RID=70;LOW=(2,3,4,5)


def selections(old):
    candidates=[]
    for row in old['records']:
        for q in row['recipes']:
            for b in q['branches']:
                if [i for i,_ in b['blockers']]==[70] and b['center_product']<-1e-8:
                    assert q['profile']==0 and q['labels'] in ([3,2,3],[3,2,4],[3,2,5])
                    candidates.append(dict(prefix=row['index'],recipe=q['index'],branch=b['branch'],
                        profile=q['profile'],labels=q['labels'],normal=b['normal'],former_NR=b['blockers'][0][1]))
    assert len(candidates)==1100
    key=lambda p:(abs(p['former_NR']),p['prefix'],p['recipe'],p['branch'])
    chosen=[];used=set()
    for label in (3,4,5):
        q=min((p for p in candidates if p['labels'][2]==label and p['prefix'] not in used),key=key)
        chosen.append(q);used.add(q['prefix'])
    q=min((p for p in candidates if p['prefix'] not in used),key=key);chosen.append(q)
    assert len(chosen)==4 and len({p['prefix'] for p in chosen})==4
    return chosen


def recipes(s,nd,tt):
    frame=[69,71,73];duties=[(0,1),(0,2)];bonds=[2]*9
    for a,b in duties:
        m=T.G.relation(s.dot(nd[frame[a]],s.ps(frame[a],tt),nd[frame[b]],s.ps(frame[b],tt)))
        assert isinstance(m,int);bonds[3*a+b]=m
    sy=T.E.recipes(tuple(s.ps(i,tt) for i in frame),tuple(bonds),tuple(duties))
    rr=[(pi,list(ms)) for pi,ms,hi in sy if not hi and all(m in LOW for m in ms)]
    direct=[(pi,list(ms)) for pi in range(7) for ms in product(LOW,repeat=3)
        if all(T.P.corner3(s.ps(frame[a],tt),s.ps(frame[b],tt),pi,bonds[3*a+b],ms[a],ms[b]) is not None for a,b in duties)]
    assert sorted(rr)==sorted(direct)
    return frame,duties,bonds,sy,rr


def run():
    names=(F.TAG,R.TAG,A.TAG,O.TAG,C.TAG,Q.TAG,Q.P.TAG)
    docs={n:json.loads((ROOT/f'results/{n}.json').read_text()) for n in names}
    for d in docs.values():
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    old=docs[F.TAG];assert old['summary']['native_N']==0
    s=T.Scout();s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[];s.counts=Counter();s.blockers=Counter()
    assert all(i not in s.index for i in range(67,74))
    interface=H.old_frame_interface(s);assert interface==old['old_frame_interface']
    debit=T.H.deletion(s.primary,s.cur,C.DELETED);assert json.loads(json.dumps(debit))==old['exact_deletion_debit']
    top=F.topology();ids=[i for i in s.cur['actual_sides'] if i not in C.DELETED]+[67,68,69,71,72,73]
    score=lambda nd,tt,ii,stage:Q.P.B.B.score(s,nd,tt,ii,stage,[i for i in F.WALLS if i in nd])
    control=score(s.old,{},s.cur['actual_sides'],'original63_control');assert control['score'] is not None and control['open_ends']==10
    prefixes=F.reconstruct(s,docs);selected=selections(old);records=[];accepted=[]
    for sel in selected:
        pre=prefixes[sel['prefix']];nd=dict(pre['normals']);tt=dict(pre['types'])
        nv,ch=C.solve_saved(s,[10,67,69],sel['profile'],sel['labels'],sel['branch'],nd,tt,sel['normal'])
        nd[73]=nv;tt[73]=sel['profile'];former=nd.pop(70);formerpi=tt.pop(70)
        steps=[*pre['reconstruction'],dict(wall=73,**ch)]
        for i in tt:
            _,bad,center=s.native(nd[i],tt[i],{j:v for j,v in nd.items() if j!=i},tt);assert not bad and center<-1e-8
        frame,duties,bonds,sy,rr=recipes(s,nd,tt);before=score(nd,tt,ids,'N_inserted_R_removed')
        row=dict(index=len(records),selection=sel,reconstruction=steps,former_R=former.tolist(),former_R_profile=formerpi,
            prefix_normals={str(i):nd[i].tolist() for i in tt},prefix_profiles=tt,before=before,
            frame=frame,duties=duties,frame_profiles=[s.ps(i,tt) for i in frame],duty_bonds=bonds,
            symbolic_high_families=[list(z) for z in sy if z[2]],recipes=[],caps=[])
        records.append(row);seen=set();digest=sha256();hist=Counter()
        for pi,ms in rr:
            rhs=[-T.B.theta(m)-(s.profiles[pi]['offsets'][s.index[i]] if i in s.index else s.off(pi,tt[i])) for i,m in zip(frame,ms)]
            roots,info=T.N.norm_roots([nd[i] for i in frame],rhs,s.profiles[pi]['square'])
            rec=dict(index=len(row['recipes']),profile=pi,labels=ms,norm=info,branches=[]);row['recipes'].append(rec)
            if info['status'].startswith('unresolved'):
                s.unresolved.append(dict(prefix=row['index'],recipe=rec['index'],norm=info));continue
            for branch,v in roots:
                vals,bad,cp=s.native(v,pi,nd,tt)
                err=max(abs(float(v@T.M@v)-s.profiles[pi]['square']),max(abs(float(v@T.M@nd[i])-q) for i,q in zip(frame,rhs)))
                b=dict(branch=branch,normal=v.tolist(),native_products=sorted(vals.items()),blockers=[[i,vals[i]] for i in bad],center_product=cp,norm_residual=err)
                rec['branches'].append(b);hist[str(tuple(bad))]+=1
                digest.update(json.dumps([rec['index'],branch,pi,ms,b['native_products'],cp],separators=(',',':')).encode())
                if err>2e-5:s.unresolved.append(dict(prefix=row['index'],recipe=rec['index'],branch=branch,error=err));continue
                if bad or cp>=-1e-8:continue
                key=(pi,*np.round(v,8))
                if key in seen:b['duplicate']=True;continue
                seen.add(key);full={**nd,RID:v};ft={**tt,RID:pi}
                for i in ft:
                    _,fail,center=s.native(full[i],ft[i],{j:w for j,w in full.items() if j!=i},ft);assert not fail and center<-1e-8
                assert all(s.corner([frame[a],frame[b],RID],full,ft) is not None for a,b in duties)
                sc=score(full,ft,ids+[RID],'joint_N_R' );faces={f['wall']:f for f in sc['whole_faces']}
                cs=[s.corner(f,full,ft) for f in top['triangles']]
                capsets={str(w):faces[w]['all_rays_finite'] and sorted(sorted(q['active']) for q in faces[w]['rays'])==sorted(sorted(f) for f in top['triangles'] if w in f) for w in (67,68,69,71)}
                cap=dict(index=len(row['caps']),recipe=rec['index'],profile=pi,labels=ms,branch=branch,normal=v.tolist(),
                    native_products=sorted(vals.items()),center_product=cp,norm_residual=err,corner_types=cs,
                    all_disk_corners=all(z is not None for z in cs),deferred_M_R_S=s.corner([71,RID,72],full,ft),
                    exact_whole_interior_caps=capsets,whole_compact_faces={str(w):faces[w]['all_rays_finite'] for w in faces},**sc)
                b['cap']=cap['index'];row['caps'].append(cap);accepted.append(dict(prefix=row['index'],**cap))
        row['branch_histogram']=dict(hist);row['complete_branch_digest']=digest.hexdigest()
        print('PREFIX',row['index'],'recipes',len(rr),'nativeR',len(row['caps']),flush=True)
    summary=dict(forward_sole_R_branches=1100,tested_prefixes=len(records),
        recipes=sum(len(r['recipes']) for r in records),branches=sum(len(q['branches']) for r in records for q in r['recipes']),
        native_R=len(accepted),clean_R=sum(c['score'] is not None for c in accepted),
        completed_disks=sum(c['score'] is not None and c['all_disk_corners'] and all(c['exact_whole_interior_caps'].values()) for c in accepted),
        fewer_than_ten_ends=sum(c['score'] is not None and c['open_ends']<10 for c in accepted),unresolved=len(s.unresolved))
    paths={Path(__file__),Path(F.__file__)};paths.update(ROOT/p for p in old['sources'])
    paths.update(ROOT/f'{f}/{n}.{e}' for n in names for f,e in (('scripts','py'),('results','json'),('results','md')))
    return dict(status='BOUNDED_OBLIQUE_JOINT_N_R_PILOT',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        deleted_walls=list(C.DELETED),retained_old_walls=58,new_id=70,
        wall_count_breakdown=old['wall_count_breakdown'],exact_deletion_debit=debit,
        old_frame_interface=interface,topology=top,control=control,selections=selected,
        selection_rule='One smallest absolute formerN/R product per N/X label3,4,5, then next smallest with a distinct frozen prefix; ties by prefix/recipe/branch.',
        records=records,accepted=accepted,summary=summary,boundary_ledgers=s.ledgers,unresolved=s.unresolved,
        numerical_snapshot=old['numerical_snapshot'],
        scope='Four declared forward prefixes with a sole removed-R obstruction. Prescribe XMR/XNR in frame X/M/N, defer MRS. Full retained native/center, target cone/PD and whole-face gates; numerical normals only. High families and other1100-prefix choices not exhausted; no ambientH9 volume claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';assert a.verify or not path.exists()
    out=json.loads(json.dumps(run()))
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['summary'],flush=True)
