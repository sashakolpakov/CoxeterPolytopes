#!/usr/bin/env python3
"""Bounded changed-incidence pilot: release the old R--30 duty via J.

Reconstruct Q without the obsolete obtuse R, build J from (7,30,Q), then
reconstruct R from (L,P,Q) and the complementary (Q,J,7) frame. No new
normal or ambient volume is certified.
"""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

import numpy as np
import fresh_b_region_lr_diagonal_20261017 as F

N,D,ROOT,M=F.N,F.D,F.ROOT,F.M
TAG='fresh_b_region_split30_native_qr_20261018'


def topology():
    old=F.topology();removed=[{30,7,70},{30,69,70}]
    fs=[t for t in old['triangles']if set(t)not in removed]
    assert len(fs)==8
    fs += [[30,7,71],[7,70,71],[30,69,71],[69,70,71]]
    bc=old['boundary'];es=Counter(tuple(sorted(e))for t in fs for e in combinations(t,2))
    vs=set(sum(fs,[]));assert(len(vs),len(es),len(fs))==(9,20,12)
    assert (30,70)not in es and (68,69)not in es and es[67,70]==2
    assert set(es.values())=={1,2}
    assert {e for e,n in es.items()if n==1}=={tuple(sorted(e))for e in zip(bc,bc[1:]+bc[:1])}
    links=[]
    for v in sorted(vs):
        ee=[tuple(i for i in t if i!=v)for t in fs if v in t];dd=Counter(i for e in ee for i in e)
        seen={next(iter(dd))}
        while True:
            nxt=seen|{i for e in ee if set(e)&seen for i in e}
            if nxt==seen:break
            seen=nxt
        assert seen==set(dd)
        assert sorted(dd.values())==([1,1]+[2]*(len(dd)-2)if v in bc else[2]*len(dd))
        links.append(dict(vertex=v,edges=ee))
    return dict(boundary=bc,triangles=fs,vertices=9,edges=20,faces=12,
                removed_edge=[30,70],new_vertex=71,J_link=[30,7,70,69,30],links=links)


def run():
    prior=json.loads((ROOT/f'results/{F.TAG}.json').read_text())
    for p,h in prior['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    s=N.Scout();disk=topology()
    verts=s.primary['retained_old_finite_supports']+[r['support']for r in s.primary['new_finite_vertices']]
    debit=json.loads(json.dumps(D.debit(verts,D.compact_incidence(verts),[47])))
    assert debit==prior['exact_deletion_debit']
    base={i:v for i,v in s.old.items()if i!=47};ids=[i for i in s.cur['actual_sides']if i!=47]
    deleted=s.score(base,{},ids,'delete47');assert s.ledgers[deleted]['open_ends']==10
    counts=Counter();rejected=Counter();unresolved=[];prefixes={};recoveries=[]

    def cofactor_root(frame,pi,ms,sg,nd,ty,native_values=False):
        rhs=[-(float(m)if native_values else N.B.theta(m))-
             (s.profiles[pi]['offsets'][s.index[i]]if i in s.index else s.off(pi,ty[i]))
             for i,m in zip(frame,ms)]
        roots,status=F.F.roots([nd[i]for i in frame],rhs,s.profiles[pi]['square'])
        return next((v for sign,v in roots if sign==sg),None),status

    def prefix(event):
        k=event['index']
        if k in prefixes:return prefixes[k]
        lp=event['L_profile'];pk=event['P']
        lv,status=cofactor_root([12,52,7],lp,[0,0,event['u']],event['L_branch'],base,{},True)
        assert lv is not None and not status.startswith('unresolved'),(k,status)
        nd={**base,67:lv};ty={67:lp}
        pv,status=cofactor_root([7,12,67],pk[0],pk[1],pk[2],nd,ty)
        assert pv is not None and not status.startswith('unresolved'),(k,status)
        nd[68]=pv;ty[68]=pk[0]
        assert max(np.max(abs(nd[i]-np.array(event['normals'][str(i)])))for i in(67,68))<2e-5
        for i in(67,68):
            vals,bad,cp=s.native(nd[i],ty[i],{j:v for j,v in nd.items()if j!=i},ty)
            assert not bad and cp < -1e-8,(k,i,bad,cp)
        prefixes[k]=nd,ty
        return nd,ty

    # Reconstruct every formerly R-only rejected Q; no saved normal is
    # accepted as geometry. Old R is used only to reproduce the ranking gap.
    for source_index,rec in enumerate(prior['near_native_rejections']):
        if [i for i,v in rec['blockers']]!=[70]:continue
        e=prior['events'][rec['event']];nd,ty=prefix(e);pi=rec['profile'];frame=rec['frame']
        rhs=[-N.B.theta(m)-(s.profiles[pi]['offsets'][s.index[i]]if i in s.index else s.off(pi,ty[i]))
             for i,m in zip(frame,rec['labels'])]
        roots,info=N.N.norm_roots([nd[i]for i in frame],rhs,s.profiles[pi]['square'])
        assert not info['status'].startswith('unresolved'),info
        qv=next((v for sign,v in roots if sign==rec['branch']),None);assert qv is not None
        err=max(abs(float(qv@M@qv)-s.profiles[pi]['square']),
                max(abs(float(qv@M@nd[i])-r)for i,r in zip(frame,rhs)))
        assert err<2e-5
        vals,bad,cp=s.native(qv,pi,nd,ty)
        assert not bad and cp < -1e-8,(source_index,bad,cp)
        oldqr=s.dot(qv,pi,np.array(e['normals']['70']),e['R'][0])
        assert oldqr>0 and abs(oldqr-rec['blockers'][0][1])<2e-5
        recoveries.append(dict(index=len(recoveries),source_rejection=source_index,event=e['index'],
            frame=frame,profile=pi,labels=rec['labels'],branch=rec['branch'],normal=qv.tolist(),
            obsolete_QR_product=oldqr,center_product=cp,norm_residual=err))
    assert len(recoveries)==440
    # Bounded representative pilot: one least-obtuse old QR choice per
    # event and Q profile. Neither old R nor its labels constrain new R.
    chosen={}
    for q in recoveries:
        key=q['event'],q['profile']
        if key not in chosen or(q['obsolete_QR_product'],q['index'])<(chosen[key]['obsolete_QR_product'],chosen[key]['index']):chosen[key]=q
    selected=[q['index']for key,q in sorted(chosen.items())]
    assert len(selected)<=66
    qrows=[];jrows=[];rrows=[];improvements=[]
    for qi in selected:
        q=recoveries[qi];e=prior['events'][q['event']];nn,tt=prefix(e)
        nd={**nn,69:np.array(q['normal'])};ty={**tt,69:q['profile']};actual=ids+[67,68,69]
        sid=s.score(nd,ty,actual,'L_P_Q_without_old_R')
        if sid is None:continue
        assert all(s.corner(t,nd,ty)is not None for t in [[7,12,68],[12,52,67],[52,30,69],[12,68,67],[52,67,69]])
        qr=dict(recovery=qi,event=q['event'],score=sid,J_bank=None);qrows.append(qr)
        if s.ledgers[sid]['open_ends']<10:improvements.append(dict(stage='L_P_Q',record=len(qrows)-1,score=sid))
        jbank=s.bank([7,30,69],[(0,1),(1,2)],nd,ty,actual,71,'J_7_30_Q')
        qr['J_bank']=jbank
        for j in jbank['caps']:
            ndj={**nd,71:np.array(j['normal'])};tyj={**ty,71:j['profile']}
            assert s.corner([7,30,71],ndj,tyj)is not None and s.corner([30,69,71],ndj,tyj)is not None
            jr=dict(q_record=len(qrows)-1,J=j,R_banks=[]);jrows.append(jr)
            if s.ledgers[j['score']]['open_ends']<10:improvements.append(dict(stage='L_P_Q_J',record=len(jrows)-1,score=j['score']))
            # The first frame makes LR/QR native. The complementary frame
            # instead enforces the two J-side corners that all first-frame
            # pilot survivors missed. Neither imposes the absent outer chord.
            for frame,duties,stage in [([67,68,69],[(0,1),(0,2)],'R_L_P_Q'),
                                       ([69,71,7],[(0,1),(1,2)],'R_Q_J_7')]:
              rbank=s.bank(frame,duties,ndj,tyj,actual+[71],70,stage)
              jr['R_banks'].append(rbank)
              for r in rbank['caps']:
                full={**ndj,70:np.array(r['normal'])};ft={**tyj,70:r['profile']};sidr=r['score']
                corners=[s.corner(t,full,ft)for t in disk['triangles']];faces=[]
                for node in(67,68,69,70,71):
                    rays=[x for x in s.ledgers[sidr]['rays']if node in x['active']]
                    expected={tuple(sorted(t))for t in disk['triangles']if node in t}
                    faces.append(dict(node=node,finite_exactly_expected=all(x['kind']=='finite'for x in rays)and
                                      {tuple(sorted(x['active']))for x in rays}==expected))
                done=all(c is not None for c in corners)and all(f['finite_exactly_expected']for f in faces)
                rrows.append(dict(j_record=len(jrows)-1,bank_stage=stage,R=r,finite_corners=corners,whole_cap_faces=faces,completed_disk=done))
                if s.ledgers[sidr]['open_ends']<10:improvements.append(dict(stage='five_caps',record=len(rrows)-1,score=sidr,completed_disk=done))
        print('SPLIT30_Q',qi,'event',q['event'],'profile',q['profile'],'ends',s.ledgers[sid]['open_ends'],
              'J',len(jbank['caps']),'R_total',len(rrows),'improvements',len(improvements),flush=True)
    sources=dict(prior['sources'])
    for folder,ext in [('scripts','py'),('results','json'),('results','md')]:
        p=ROOT/f'{folder}/{F.TAG}.{ext}';sources[str(p.relative_to(ROOT))]=sha256(p.read_bytes()).hexdigest()
    sources[f'scripts/{TAG}.py']=sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(status='BOUNDED_SPLIT30_NATIVE_QR_CONSTRUCTION_PILOT',sources=sources,disk=disk,
        deleted_walls=[47],retained_old_wall_count=62,exact_deletion_debit=debit,
        control_score=s.control,deletion_score=deleted,reconstructed_Q_candidates=recoveries,
        selection_rule='For each savedLR event and Q profile, select the least positive obsoleteQR product; oldR is never retained.',
        selected_recoveries=selected,high_label_samples=N.HIGH,Q_records=qrows,J_records=jrows,R_records=rrows,
        counts=dict(reconstructed_Q=len(recoveries),selected_Q=len(selected),scored_Q=len(qrows),
                    native_J=len(jrows),native_R=len(rrows)),bank_counts=dict(s.counts),
        blockers=[dict(stage=a,walls=list(b),count=n)for(a,b),n in sorted(s.blockers.items())],
        near_native_failures=s.near,boundary_ledgers=s.ledgers,
        histograms=[dict(stage=a,finite=b,exterior=c,open_ends=d,count=n)for(a,b,c,d),n in sorted(s.hist.items())],
        improvements=improvements,unresolved=s.unresolved,
        scope='Changed stellar subdivision of30-R. Reconstruct Q without obsolete obtuseR; J takes7/30 and30/Q duties, then newR usesL/P/Q and the complementaryQ/J/7 frame. Exact old47 incidence debit and disk only; bounded representatives of22 numerical event prefixes, seven A2profiles, low2..5/high6,8,12,24 norm recipes. All62 retained oldwall/native newwall products, strict commoncenter, fullnumericalS6cone and independentPD finite-support/open-end agreement at acceptedstages. No exact newnormal, all-label/continuous exclusion, ambient finitevertex census or finitevolume assertion. Original63/199/10 exact control unchanged.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],'improvements',len(out['improvements']),'unresolved',len(out['unresolved']),flush=True)
