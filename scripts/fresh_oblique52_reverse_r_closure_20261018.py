#!/usr/bin/env python3
"""Reverse R70 from X/M/S, imposing XMR/MRS before old10JR/JXR duties.

All prefix deduplication uses complete reconstruction instructions. The
four-prefix pilot is selected by the former S/R gap, with distinct frozen
R/M prefix indices. New normals and target cone checks remain numerical.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import product
from pathlib import Path
import numpy as np
import fresh_oblique52_boundary_connector_20261018 as A

ROOT=A.ROOT;O=A.O;C=A.C;T=A.T;Q=A.Q;H=A.H
TAG='fresh_oblique52_reverse_r_closure_20261018';RID=70;LOW=(2,3,4,5)


def reconstruct(s,docs):
    previous=A.reconstruct(s,docs[O.TAG],docs[C.TAG],docs[Q.TAG],docs[Q.P.TAG])
    out=[];seen={}
    for k,cap in enumerate(docs[A.TAG]['accepted']):
        p=previous[cap['prefix']];nd=dict(p['normals']);tt=dict(p['types'])
        sv,check=C.solve_saved(s,[52,68,71],cap['profile'],cap['labels'],cap['branch'],nd,tt,cap['normal'])
        nd[72]=sv;tt[72]=cap['profile']
        steps=[*p['reconstruction'],dict(wall=72,**check)]
        key=json.dumps([[z['wall'],z['frame'],z['profile'],z['labels'],z['branch']] for z in steps],separators=(',',':'))
        old_r=nd.pop(RID);old_pi=tt.pop(RID)
        for i in tt:
            _,bad,center=s.native(nd[i],tt[i],{j:v for j,v in nd.items() if j!=i},tt)
            assert not bad and center<-1e-8
        old_product=s.dot(old_r,old_pi,sv,tt[72]);assert abs(old_product-dict(cap['native_products'])[RID])<2e-5
        origin=dict(accepted=k,prefix=cap['prefix'],S=cap['index'])
        if key in seen:
            pre=out[seen[key]];assert tt==pre['types'] and all(np.max(abs(nd[i]-pre['normals'][i]))<1e-9 for i in tt)
            pre['origins'].append(origin);continue
        seen[key]=len(out)
        out.append(dict(index=len(out),origins=[origin],recipe_key=key,reconstruction=steps,normals=nd,types=tt,
            original_R=old_r.tolist(),original_R_profile=old_pi,former_SR_product=old_product))
    assert sum(len(p['origins']) for p in out)==449
    return out


def choose_pilot(prefixes):
    chosen=[];old=set()
    for p in sorted(prefixes,key=lambda p:(-p['former_SR_product'],p['index'])):
        key=p['origins'][0]['prefix']
        if key in old:continue
        old.add(key);chosen.append(p['index'])
        if len(chosen)==4:break
    assert len(chosen)==4
    return chosen


def recipes(s,nd,tt):
    frame=[69,71,72];duties=[(0,1),(1,2)];bonds=[2]*9
    for a,b in duties:
        m=T.G.relation(s.dot(nd[frame[a]],s.ps(frame[a],tt),nd[frame[b]],s.ps(frame[b],tt)))
        assert isinstance(m,int);bonds[3*a+b]=m
    sy=T.E.recipes(tuple(s.ps(i,tt) for i in frame),tuple(bonds),tuple(duties))
    rr=[(pi,list(ms)) for pi,ms,hi in sy if not hi and all(m in LOW for m in ms)]
    direct=[(pi,list(ms)) for pi in range(7) for ms in product(LOW,repeat=3)
        if all(T.P.corner3(s.ps(frame[a],tt),s.ps(frame[b],tt),pi,bonds[3*a+b],ms[a],ms[b]) is not None for a,b in duties)]
    assert sorted(rr)==sorted(direct)
    return frame,duties,bonds,sy,rr


def run(pilot=False):
    names=(A.TAG,O.TAG,C.TAG,Q.TAG,Q.P.TAG);docs={n:json.loads((ROOT/f'results/{n}.json').read_text()) for n in names}
    for d in docs.values():
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    old=docs[A.TAG];assert old['summary']['native_S']==old['summary']['clean_S']==449
    s=T.Scout();s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[];s.counts=Counter();s.blockers=Counter()
    assert all(i not in s.index for i in range(67,73))
    interface=H.old_frame_interface(s);assert interface==old['old_frame_interface']
    debit=T.H.deletion(s.primary,s.cur,C.DELETED);assert json.loads(json.dumps(debit))==old['exact_deletion_debit']
    top=A.topology();ids=[i for i in s.cur['actual_sides'] if i not in C.DELETED]+[67,68,69,71,72]
    score=lambda nd,tt,ii,stage:Q.P.B.B.score(s,nd,tt,ii,stage,[i for i in A.WALLS if i in nd])
    control=score(s.old,{},s.cur['actual_sides'],'original63_control');assert control['score'] is not None and control['open_ends']==10
    prefixes=reconstruct(s,docs);pilot_ids=choose_pilot(prefixes)
    chosen=[prefixes[i] for i in pilot_ids] if pilot else prefixes
    records=[];accepted=[]
    for pre in chosen:
        nd,tt=pre['normals'],pre['types'];frame,duties,bonds,sy,rr=recipes(s,nd,tt)
        before=score(nd,tt,ids,'remove_previous_R')
        row=dict(index=pre['index'],origins=pre['origins'],recipe_key=pre['recipe_key'],reconstruction=pre['reconstruction'],
            prefix_normals={str(i):nd[i].tolist() for i in tt},prefix_profiles=tt,original_R=pre['original_R'],
            original_R_profile=pre['original_R_profile'],former_SR_product=pre['former_SR_product'],before=before,
            frame=frame,duties=duties,frame_profiles=[s.ps(i,tt) for i in frame],duty_bonds=bonds,
            symbolic_high_families=[list(z) for z in sy if z[2]],recipes=[],caps=[])
        records.append(row);seen=set();digest=sha256();hist=Counter()
        for pi,ms in rr:
            rhs=[-T.B.theta(m)-(s.profiles[pi]['offsets'][s.index[i]] if i in s.index else s.off(pi,tt[i])) for i,m in zip(frame,ms)]
            roots,info=T.N.norm_roots([nd[i] for i in frame],rhs,s.profiles[pi]['square'])
            rec=dict(index=len(row['recipes']),profile=pi,labels=ms,norm=info,branches=[]);row['recipes'].append(rec)
            if info['status'].startswith('unresolved'):
                s.unresolved.append(dict(prefix=pre['index'],recipe=rec['index'],norm=info));continue
            for branch,v in roots:
                vals,bad,cp=s.native(v,pi,nd,tt)
                err=max(abs(float(v@T.M@v)-s.profiles[pi]['square']),max(abs(float(v@T.M@nd[i])-q) for i,q in zip(frame,rhs)))
                b=dict(branch=branch,normal=v.tolist(),native_products=sorted(vals.items()),blockers=[[i,vals[i]] for i in bad],center_product=cp,norm_residual=err)
                rec['branches'].append(b);hist[str(tuple(bad))]+=1
                digest.update(json.dumps([rec['index'],branch,pi,ms,b['native_products'],cp],separators=(',',':')).encode())
                if err>2e-5:s.unresolved.append(dict(prefix=pre['index'],recipe=rec['index'],branch=branch,error=err));continue
                if bad or cp>=-1e-8:continue
                key=(pi,*np.round(v,8))
                if key in seen:b['duplicate']=True;continue
                seen.add(key);full={**nd,RID:v};ft={**tt,RID:pi}
                for i in ft:
                    _,fail,center=s.native(full[i],ft[i],{j:w for j,w in full.items() if j!=i},ft);assert not fail and center<-1e-8
                assert all(s.corner([frame[a],frame[b],RID],full,ft) is not None for a,b in duties)
                sc=score(full,ft,ids+[RID],'reverse_R_X_M_S');faces={f['wall']:f for f in sc['whole_faces']}
                cs=[s.corner(f,full,ft) for f in top['triangles']]
                capsets={str(w):faces[w]['all_rays_finite'] and sorted(sorted(q['active']) for q in faces[w]['rays'])==sorted(sorted(f) for f in top['triangles'] if w in f) for w in (67,68,69,71)}
                cap=dict(index=len(row['caps']),recipe=rec['index'],profile=pi,labels=ms,branch=branch,normal=v.tolist(),
                    native_products=sorted(vals.items()),center_product=cp,norm_residual=err,corner_types=cs,
                    all_disk_corners=all(z is not None for z in cs),deferred_10_J_R=s.corner([10,67,RID],full,ft),
                    deferred_J_X_R=s.corner([67,69,RID],full,ft),exact_whole_interior_caps=capsets,
                    whole_compact_faces={str(w):faces[w]['all_rays_finite'] for w in faces},**sc)
                b['cap']=cap['index'];row['caps'].append(cap);accepted.append(dict(prefix=pre['index'],**cap))
        row['branch_histogram']=dict(hist);row['complete_branch_digest']=digest.hexdigest()
        print('PREFIX',pre['index'],'recipes',len(rr),'nativeR',len(row['caps']),flush=True)
    summary=dict(frozen_pointed_configurations=449,distinct_recipe_prefixes=len(prefixes),tested_prefixes=len(records),
        recipes=sum(len(r['recipes']) for r in records),branches=sum(len(q['branches']) for r in records for q in r['recipes']),
        native_R=len(accepted),clean_R=sum(c['score'] is not None for c in accepted),
        completed_disks=sum(c['score'] is not None and c['all_disk_corners'] and all(c['exact_whole_interior_caps'].values()) for c in accepted),
        fewer_than_ten_ends=sum(c['score'] is not None and c['open_ends']<10 for c in accepted),unresolved=len(s.unresolved))
    paths={Path(__file__),Path(A.__file__)};paths.update(ROOT/p for p in old['sources'])
    paths.update(ROOT/f'{f}/{n}.{e}' for n in names for f,e in (('scripts','py'),('results','json'),('results','md')))
    return dict(status='PILOT_OBLIQUE52_REVERSE_R' if pilot else 'BOUNDED_OBLIQUE52_REVERSE_R',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        deleted_walls=list(C.DELETED),retained_old_walls=len(s.old)-len(C.DELETED)+6,reconstructed_id=RID,
        exact_deletion_debit=debit,old_frame_interface=interface,topology=top,control=control,
        deduplication='Identical complete six-step reconstruction instructions; proximity alone never removes a prefix.',
        pilot_selection=dict(rule='Closest former S/R product to -2, breaking ties by prefix index, with four distinct frozen R/M prefix indices.',indices=pilot_ids),
        prefix_inventory=[dict(index=p['index'],origins=p['origins'],recipe_key=p['recipe_key'],former_SR_product=p['former_SR_product']) for p in prefixes],
        records=records,accepted=accepted,summary=summary,boundary_ledgers=s.ledgers,unresolved=s.unresolved,
        numerical_snapshot=old['numerical_snapshot'],
        scope='Bounded complete seven-profile low-label2..5 reverse-R banks over selected or all frozen449 S caps. X/M/R and M/S/R are prescribed;10/J/R and J/X/R deferred. All58 retained old walls, mutual native/center gates, complete numerical cone/PD and whole-face checks. Exact19/12/1 old-certificate debit unchanged. No high-label/continuous/all-prefix theorem, exact new normals or ambientH9 volume claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--pilot',action='store_true');a=ap.parse_args();assert not(a.verify and a.pilot)
    p=ROOT/f'results/{TAG}{"_pilot" if a.pilot else ""}.json';assert a.verify or not p.exists()
    out=json.loads(json.dumps(run(a.pilot)))
    if a.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['summary'],flush=True)
