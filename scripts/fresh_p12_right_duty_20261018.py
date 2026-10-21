#!/usr/bin/env python3
"""Purposeful P12-right moving prefix pilot; numerical evidence only."""
import argparse,json
from collections import Counter
from hashlib import sha256
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
import fresh_qj_connector_20261018 as Q
C,N,S,ROOT,M=Q.C,Q.N,Q.S,Q.ROOT,Q.M
TAG='fresh_p12_right_duty_20261018'

def run():
    f=C.Family();s=f.scout
    prior=json.loads((ROOT/f'results/{Q.TAG}.json').read_text())
    for p,h in prior['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    grid=sorted(set(map(float,[2+x for x in np.geomspace(1e-5,2,121)]+[f.u0]+[v['u']for v in prior['seeds']])))
    s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[]
    control=s.score(s.old,{},s.cur['actual_sides'],'control63');delete=s.score(f.base,{},f.ids,'delete47')
    rec=dict(frame=[7,12,67],profile=0,labels=[3,2,2]);records=[];accepted=[];unresolved=[];hist=Counter()
    assert N.P.corner3(4,0,0,2,3,2)is not None and N.P.corner3(0,0,0,2,2,2)is not None
    for u in grid:
        nd=dict(f.base);ty={};lr,ls=f.roots(f.recipes[67],u,nd,ty)
        if not lr:
            if ls.startswith('unresolved'):unresolved.append(dict(u=u,stage='L',status=ls))
            hist['L_nonreal_or_singular']+=1;continue
        lv=next(v for sg,v in lr if sg==f.recipes[67]['sheet']);nd[67]=lv;ty[67]=0
        lvals,lbad,lcp=s.native(lv,0,f.base,{})
        roots,status=f.roots(rec,u,nd,ty)
        if not roots:
            if status.startswith('unresolved'):unresolved.append(dict(u=u,stage='P',status=status))
            hist['P_nonreal_or_singular']+=1;continue
        for sheet,pv in roots:
            vals,bad,cp=s.native(pv,0,nd,ty)
            row=dict(index=len(records),u=u,P_sheet=sheet,L_normal=lv.tolist(),P_normal=pv.tolist(),
                L_blockers=[[i,lvals[i]]for i in lbad],L_center=lcp,P_blockers=[[i,vals[i]]for i in bad],P_center=cp,
                P_native_products=sorted(vals.items()),P_norm_error=abs(float(pv@M@pv)-2),
                frame_error=max(abs(s.dot(pv,0,nd[i],ty.get(i,0))+N.B.theta(m))for i,m in zip(rec['frame'],rec['labels'])))
            assert row['P_norm_error']<2e-5 and row['frame_error']<2e-5
            if not lbad and lcp< -1e-8 and not bad and cp< -1e-8:
                full={**nd,68:pv};ft={**ty,68:0};score=s.score(full,ft,f.ids+[67,68],'changedP_prefix')
                if score is not None:
                    row['score']=score;accepted.append(row['index'])
            hist['P_sheet_candidates']+=1;hist['native_P']+=not bad and cp< -1e-8
            records.append(row)
    paths=[Path(__file__)]+[ROOT/f'{folder}/{Q.TAG}.{ext}'for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    return dict(status='NUMERICAL_P12_RIGHT_DUTY_PREFIX_PILOT',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        P_old_recipe=f.recipes[68],P_new_recipe=rec,L_recipe=f.recipes[67],u_grid=grid,
        deleted_walls=[47],retained_old_wall_count=62,exact_deletion_debit=prior['exact_deletion_debit'],
        control_score=control,deletion_score=delete,records=records,accepted=accepted,counts=dict(hist),
        boundary_ledgers=s.ledgers,unresolved=unresolved+s.unresolved,
        scope='Purposeful single native duty change P12:3 to2 in frame7/12/L, with fixedP7=3 andPL=2, profile0. Both P norm sheets sampled at saved event values and a bounded logarithmic ugrid. All retained oldwall contacts, L/P native products and strictcenter are required before any closure continuation. Completecone/PD agreement for accepted prefixes. Numerical only; no interval/all-u exclusion or newfinitevolume chamber.')

def closure():
    f=C.Family();f.recipes[68]['labels']=[3,2,2]
    s=f.scout;s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[]
    disk=Q.topology();grid=sorted(set([float(2+x)for x in np.geomspace(1e-5,2,121)]+[f.u0]))
    cache={};branches=[];brackets=[];rejected=[];events=[];improvements=[];unresolved=[];counts=Counter()
    def lp(u):
        if u in cache:return cache[u]
        nd=dict(f.base);ty={}
        for i in(67,68):
            rec=f.recipes[i];rr,status=f.roots(rec,u,nd,ty)
            if not rr:return None
            nd[i]=next(v for sg,v in rr if sg==rec['sheet']);ty[i]=0
        cache[u]=(nd,ty);return cache[u]
    def state(u,key):
        pair=lp(u)
        if pair is None:return None
        nd,ty=pair;nd=dict(nd);ty=dict(ty)
        k,qs,pi,m7,mP,rs,q30,q52=key
        for node,rec,sign in((69,dict(frame=[30,52,67],profile=0,labels=[q30,q52,k]),qs),
                              (70,dict(frame=[7,67,68],profile=pi,labels=[m7,2,mP]),rs)):
            rr,status=f.roots(rec,u,nd,ty)
            if not rr:return None
            nd[node]=next(v for sg,v in rr if sg==sign);ty[node]=rec['profile']
        return s.dot(nd[69],ty[69],nd[70],ty[70]),nd,ty
    labels=[3,4,12,24]
    for q30,q52 in((3,2),(2,2),(2,3)):
      for k in labels:
        recipes=[]
        for pi in range(7):
          for m7 in range(2,6):
           for mP in range(2,6):
            if all(v is not None for v in(N.P.corner3(0,0,pi,k,2,2),N.P.corner3(0,0,pi,2,2,mP),N.P.corner3(4,0,pi,3,m7,mP))):
                if(q30,q52)==(3,2)or(pi,m7,mP)==(0,3,2):recipes.append((pi,m7,mP))
        for pi,m7,mP in recipes:
         for qs in(-1,1):
          for rs in(-1,1):
            key=[k,qs,pi,m7,mP,rs,q30,q52];prev=None;values=[];local=[]
            for u in grid:
                ans=state(u,key)
                if ans is None:prev=None;counts['nonreal_or_solver_rejected_samples']+=1;continue
                value=ans[0];values.append([u,value]);counts['real_samples']+=1
                if prev is not None and prev[1]*value<0:local.append([prev[0],u,prev[1],value])
                prev=u,value
            branches.append(dict(key=key,real_samples=len(values),brackets=local,
                closest=min(values,key=lambda r:abs(r[1]))if values else None,
                value_range=[min((v for u,v in values),default=None),max((v for u,v in values),default=None)]))
            seen=[]
            for interval in local:
                bi=len(brackets);brackets.append(dict(key=key,interval=interval))
                try:
                    u=float(brentq(lambda x:state(x,key)[0],*interval[:2],xtol=2e-14,maxiter=100));gap,nd,ty=state(u,key)
                except(ValueError,TypeError,StopIteration,RuntimeError)as error:
                    unresolved.append(dict(stage='QR_root',bracket=bi,error=str(error)));continue
                if any(abs(u-v)<1e-8 for v in seen):counts['duplicate_event_roots']+=1;continue
                seen.append(u);native=[];badall=[]
                for i in(67,68,69,70):
                    vals,bad,cp=s.native(nd[i],ty[i],{a:v for a,v in nd.items()if a!=i},ty)
                    row=dict(node=i,profile=ty[i],normal=nd[i].tolist(),native_products=sorted(vals.items()),
                             blockers=[[j,vals[j]]for j in bad],center_product=cp,
                             norm_error=abs(float(nd[i]@M@nd[i])-s.profiles[ty[i]]['square']))
                    assert row['norm_error']<2e-5;native.append(row)
                    if bad or cp>=-1e-8:badall.append(row)
                if N.G.relation(gap)!=2:
                    unresolved.append(dict(stage='QR_native_spacing',bracket=bi,u=u,product=gap));continue
                if badall:
                    rejected.append(dict(bracket=bi,u=u,QR=gap,bad=badall));continue
                sid=s.score(nd,ty,f.ids+[67,68,69,70],'P12right_LPQR')
                if sid is None:continue
                event=dict(index=len(events),key=key,u=u,bracket=bi,QR=gap,native=native,score=sid)
                bank=s.bank([7,30,70],[(0,1),(0,2)],nd,ty,f.ids+[67,68,69,70],71,'P12right_J')
                event['J_bank']=bank;events.append(event)
                if s.ledgers[sid]['open_ends']<10:improvements.append(dict(stage='LPQR',event=event['index'],score=sid))
                event['K_checks']=[]
                for cap in bank['caps']:
                    if s.ledgers[cap['score']]['open_ends']<10:improvements.append(dict(stage='LPQRJ',event=event['index'],J=cap['index'],score=cap['score']))
                    full={**nd,71:np.array(cap['normal'])};ft={**ty,71:cap['profile']}
                    assert ty[70]==cap['profile']==0 and q30==3
                    rj=N.G.relation(s.dot(full[70],0,full[71],0));assert rj in(3,4)
                    domain=Q.domain(q30,rj);assert all(r['high_position']is None for r in domain)
                    mat=np.array([full[i]@M for i in(30,69,70,71)]);condition=float(np.linalg.cond(mat))
                    kc=dict(J=cap['index'],frame=[30,69,70,71],condition=condition,domain_key=[q30,rj],records=[],accepted=[])
                    event['K_checks'].append(kc)
                    if not np.isfinite(condition)or condition>1e10:
                        unresolved.append(dict(stage='changedK_frame',event=event['index'],J=cap['index'],condition=condition));continue
                    for di,rec in enumerate(domain):
                        kpi=rec['profile'];ms=rec['labels'];offsets=np.array([s.profiles[kpi]['offsets'][s.index[i]]if i in s.index else s.off(kpi,ft[i])for i in(30,69,70,71)])
                        rhs=-np.array([N.B.theta(m)for m in ms])-offsets;kv=np.linalg.solve(mat,rhs)
                        residual=float(kv@M@kv)-s.profiles[kpi]['square']
                        row=dict(domain_record=di,profile=kpi,labels=ms,norm_residual=residual);kc['records'].append(row)
                        if abs(residual)>=2e-7:continue
                        vals,bad,cp=s.native(kv,kpi,full,ft);row.update(native_products=sorted(vals.items()),blockers=[[i,vals[i]]for i in bad],center_product=cp)
                        if bad or cp>=-1e-8:continue
                        withk={**full,72:kv};kt={**ft,72:kpi};kid=s.score(withk,kt,f.ids+list(range(67,73)),'changedP_K')
                        if kid is None:continue
                        corners=[s.corner(t,withk,kt)for t in disk['triangles']];faces=[]
                        for node in range(67,73):
                            rays=[r for r in s.ledgers[kid]['rays']if node in r['active']];expected={tuple(sorted(t))for t in disk['triangles']if node in t}
                            faces.append(dict(node=node,finite_exactly_expected=all(r['kind']=='finite'for r in rays)and{tuple(sorted(r['active']))for r in rays}==expected))
                        kc['accepted'].append(dict(domain_record=di,normal=kv.tolist(),score=kid,finite_corners=corners,whole_faces=faces,
                            completed_disk=all(x is not None for x in corners)and all(r['finite_exactly_expected']for r in faces)))
                        if s.ledgers[kid]['open_ends']<10:improvements.append(dict(stage='LPQRJK',event=event['index'],J=cap['index'],recipe=di,score=kid))
                print('P12_RIGHT_EVENT',event['index'],key,u,'ends',s.ledgers[sid]['open_ends'],'J',len(bank['caps']),flush=True)
        print('P12_RIGHT_LABEL',[q30,q52,k],'branches',len(branches),'brackets',len(brackets),'events',len(events),flush=True)
    return dict(domain=dict(P_recipe=f.recipes[68],L_recipe=f.recipes[67],u_grid=grid,LQ_labels=labels,
            primary_Q_endpoints=[3,2],adaptive_Q_endpoints=[[2,2],[2,3]],adaptive_R_recipe=[0,3,2],
            R_L=2,R_Q=2,R7_RP_labels=[2,3,4,5],both_Q_R_sheets=True,J_frame=[7,30,70],J_duties=[[0,1],[0,2]],J_high_samples=N.HIGH),
        branches=branches,brackets=brackets,rejected_events=rejected,events=events,boundary_ledgers=s.ledgers,
        counts=dict(counts,branches=len(branches),brackets=len(brackets),native_events=len(events),native_J=sum(len(e['J_bank']['caps'])for e in events),K_recipes=sum(len(k['records'])for e in events for k in e['K_checks']),native_K=sum(len(k['accepted'])for e in events for k in e['K_checks'])),
        improvements=improvements,unresolved=unresolved+s.unresolved,
        scope='ChangedP12-right QR event continuation; representative LQ3/4/12/24 and saved finite recipe subdomain, all sheets. Only sign-changing brackets on finiteugrid solved; tangencies, narrowwindows and rejectednear-singular normframes remain unsearched. All retained native contacts/center and fullcone/PD support agreement required before acceptance. No all-u or all-label exclusion.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();p=ROOT/f'results/{TAG}.json'
    out=run();out['closure']=closure();out=json.loads(json.dumps(out))
    if a.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],len(out['records']),out['counts'],'accepted',len(out['accepted']),'unresolved',len(out['unresolved']),'closure',out['closure']['counts'],flush=True)
