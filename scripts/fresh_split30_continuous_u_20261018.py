#!/usr/bin/env python3
"""Continue the split30 disk after removing the obsolete LR event equation.

Numerical discovery only. Every continued root uses a fixed cofactor sheet;
all native contacts are rechecked rather than inherited across parameters.
"""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
import fresh_b_region_split30_native_qr_20261018 as S

F,N,D,ROOT,M=S.F,S.N,S.D,S.ROOT,S.M
TAG='fresh_split30_continuous_u_20261018'


class Family:
    def __init__(self):
        self.scout=N.Scout();s=self.scout
        self.prior=json.loads((ROOT/f'results/{S.TAG}.json').read_text())
        self.old=json.loads((ROOT/f'results/{F.TAG}.json').read_text())
        for name,h in self.prior['sources'].items():
            assert sha256((ROOT/name).read_bytes()).hexdigest()==h,name
        self.base={i:v for i,v in s.old.items()if i!=47}
        self.ids=[i for i in s.cur['actual_sides']if i!=47]
        r=next(r for r in self.prior['R_records']if r['bank_stage']=='R_Q_J_7')
        j=self.prior['J_records'][r['j_record']]
        qrow=self.prior['Q_records'][j['q_record']]
        q=self.prior['reconstructed_Q_candidates'][qrow['recovery']]
        e=self.old['events'][q['event']]
        self.u0=e['u'];self.event=e
        self.recipes={67:dict(frame=[12,52,7],profile=e['L_profile'],labels=[2,2,'u'],sheet=e['L_branch']),
                      68:dict(frame=[7,12,67],profile=e['P'][0],labels=e['P'][1],sheet=e['P'][2])}
        self.reference={67:e['normals']['67'],68:e['normals']['68'],69:q['normal'],71:j['J']['normal'],70:r['R']['normal']}
        for node,cap in[(69,q),(71,j['J']),(70,r['R'])]:
            self.recipes[node]=dict(frame=cap['frame'],profile=cap['profile'],labels=cap['labels'],sheet=None)
        self.order=[67,68,69,71,70]
        nd=dict(self.base);ty={}
        for node in self.order:
            rec=self.recipes[node]
            roots,status=self.roots(rec,self.u0,nd,ty)
            assert not status.startswith('unresolved') and roots,(node,status)
            distances=[(float(np.max(abs(v-np.array(self.reference[node])))),sign,v)for sign,v in roots]
            error,sign,v=min(distances,key=lambda row:row[0])
            assert error<2e-5,(node,error)
            if rec['sheet']is not None:assert sign==rec['sheet']
            rec['sheet']=sign;rec['seed_reconstruction_error']=error
            nd[node]=v;ty[node]=rec['profile']
        self.seed=self.state(self.u0)
        assert self.seed['all_native']

    def roots(self,rec,u,nd,ty):
        s=self.scout;pi=rec['profile']
        magnitudes=[u if m=='u'else N.B.theta(m)for m in rec['labels']]
        rhs=[-a-(s.profiles[pi]['offsets'][s.index[i]]if i in s.index else s.off(pi,ty[i]))
             for i,a in zip(rec['frame'],magnitudes)]
        return F.F.roots([nd[i]for i in rec['frame']],rhs,s.profiles[pi]['square'])

    def state(self,u,include_old_R=True):
        s=self.scout;nd=dict(self.base);ty={};rows=[]
        for node in (self.order if include_old_R else self.order[:-1]):
            rec=self.recipes[node];roots,status=self.roots(rec,u,nd,ty)
            if status.startswith('unresolved')or not roots:
                return dict(u=u,status=status,failed_node=node,all_native=False,rows=rows)
            v=next(v for sign,v in roots if sign==rec['sheet'])
            vals,bad,cp=s.native(v,rec['profile'],nd,ty)
            rows.append(dict(node=node,normal=v.tolist(),native_products=sorted(vals.items()),
                blockers=[[i,vals[i]]for i in bad],center_product=cp,
                finite_contacts=[[i,N.G.relation(x)]for i,x in vals.items()if isinstance(N.G.relation(x),int)]))
            nd[node]=v;ty[node]=rec['profile']
        out=dict(u=u,status='real',all_native=all(not r['blockers']and r['center_product']< -1e-8 for r in rows),
                    prefix_native=all(not r['blockers']and r['center_product']< -1e-8 for r in rows if r['node']!=70),rows=rows)
        if include_old_R:out.update(RL=s.dot(nd[70],ty[70],nd[67],ty[67]),RP=s.dot(nd[70],ty[70],nd[68],ty[68]))
        return out


def pilot():
    f=Family()
    print('SEED',f.u0,'RECIPES',f.recipes,flush=True)
    for du in[-1e-4,0,1e-4]:
        d=f.state(f.u0+du)
        print('PERTURB',du,'native',d['all_native'],'pairs',[d.get('RL'),d.get('RP')],
              'finite',[(r['node'],r['finite_contacts'])for r in d['rows']],
              'bad',[(r['node'],r['blockers'])for r in d['rows']if r['blockers']],flush=True)
    for u in[2.0001,2.001,2.005,2.01,2.025,2.05,2.075,2.1,2.125,2.15,2.175,2.2,2.25,2.3,2.4,2.5,2.75,3,4,6,10,16]:
        d=f.state(u)
        print('CONTINUE',u,d['status'],'native',d['all_native'],'prefix',d.get('prefix_native'),
              'pairs',[d.get('RL'),d.get('RP')],'bad',[(r['node'],r['blockers'])for r in d['rows']if r['blockers']],flush=True)


def completion_domain():
    # A label>=6 at LQ isolates R in the required finite LQR corner.
    # R7/RP>=6 conflict with P7=3; RJ>=6 conflicts with JQ=5.
    rows=[]
    for pi in range(7):
      for m7 in range(2,6):
       for mP in range(2,6):
        for mJ in range(2,6):
         if all(x is not None for x in [N.P.corner3(0,0,pi,24,2,2),
             N.P.corner3(0,0,pi,2,2,mP),N.P.corner3(4,0,pi,3,m7,mP),
             N.P.corner3(4,0,pi,2,m7,mJ),N.P.corner3(0,0,pi,5,2,mJ)]):
            rows.append([pi,m7,mP,mJ])
    return rows


def native_first_pilot():
    f=Family();s=f.scout;domain=completion_domain()
    print('COMPLETE_NECESSARY_HIGH_LQ_DOMAIN',domain,flush=True)
    grid=sorted(set([float(2+x)for x in np.geomspace(1e-5,2,101)]+[f.u0]))
    cache={u:f.state(u,False)for u in grid};summary=[];brackets=[];coherent=[]
    for pi,m7,mP in sorted({tuple(r[:3])for r in domain}):
      rec=dict(frame=[7,67,68],profile=pi,labels=[m7,2,mP])
      for sheet in[-1,1]:
        vals=[];previous=None;bad=Counter()
        for u in grid:
            state=cache[u]
            if state['status']!='real':previous=None;continue
            nd={**f.base,**{r['node']:np.array(r['normal'])for r in state['rows']if r['node']!=70}}
            ty={i:f.recipes[i]['profile']for i in nd if i>=67}
            rr,status=f.roots(rec,u,nd,ty)
            if not rr:previous=None;continue
            rv=next(v for sg,v in rr if sg==sheet)
            qr=s.dot(rv,pi,nd[69],ty[69]);jr=s.dot(rv,pi,nd[71],ty[71])
            vv,bb,cp=s.native(rv,pi,nd,ty)
            bad[tuple(bb)]+=1
            if previous is not None and previous[1]*qr<0:
                brackets.append(dict(profile=pi,R7=m7,RP=mP,sheet=sheet,lower=previous[0],upper=u))
            previous=u,qr
            vals.append([u,qr,jr,state['prefix_native'],not bb and cp < -1e-8,bb])
            if state['prefix_native']and not bb and cp < -1e-8:
                coherent.append(dict(u=u,profile=pi,R7=m7,RP=mP,sheet=sheet,RQ=qr,RJ=jr))
        summary.append(dict(profile=pi,R7=m7,RP=mP,sheet=sheet,real=len(vals),
            RQ_range=[min((r[1]for r in vals),default=None),max((r[1]for r in vals),default=None)],
            closest=min(vals,key=lambda r:abs(r[1]))if vals else None,
            blockers=[dict(walls=k,count=v)for k,v in bad.most_common(3)]))
        print('NATIVE_FIRST',summary[-1],flush=True)
    print('ROOT_BRACKETS',brackets,flush=True)
    print('COHERENT',len(coherent),'best',sorted(coherent,key=lambda r:abs(r['RQ']))[:12],flush=True)


def run():
    f=Family();s=f.scout;disk=S.topology();counts=Counter();unresolved=[]
    # The inherited constructor also scores an unrelated58/66 deletion.
    # Start this package's ledger anew with only its actual controls.
    s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[]
    s.control=s.score(s.old,{},s.cur['actual_sides'],'control63')
    deletion_score=s.score(f.base,{},f.ids,'delete47')
    assert s.ledgers[deletion_score]['counts']==dict(finite=24,exterior=6)
    assert s.ledgers[deletion_score]['open_ends']==10
    grid=sorted(set([float(2+x)for x in np.geomspace(1e-5,2,121)]+[f.u0]))
    labels=[2,3,4,5,6,8,12,24]
    cache={};events=[];branches=[];brackets=[];rejected=[];improvements=[]
    def lp(u):
        if u in cache:return cache[u]
        nd=dict(f.base);ty={};rows=[]
        for i in[67,68]:
            rec=f.recipes[i];roots,status=f.roots(rec,u,nd,ty)
            if status.startswith('unresolved')or not roots:return None
            v=next(v for sg,v in roots if sg==rec['sheet'])
            vals,bad,cp=s.native(v,rec['profile'],nd,ty)
            rows.append(dict(node=i,normal=v.tolist(),native_products=sorted(vals.items()),
                             blockers=bad,center_product=cp))
            nd[i]=v;ty[i]=rec['profile']
        cache[u]=(nd,ty,rows)
        return cache[u]
    def state(u,k,qsign,pi,m7,mP,rsign,q30,q52):
        p=lp(u)
        if p is None:return None
        nd,ty,rows=p;nd=dict(nd);ty=dict(ty)
        qrec=dict(frame=[30,52,67],profile=0,labels=[q30,q52,k])
        qr,status=f.roots(qrec,u,nd,ty)
        if status.startswith('unresolved')or not qr:return None
        qv=next(v for sg,v in qr if sg==qsign);nd[69]=qv;ty[69]=0
        rrec=dict(frame=[7,67,68],profile=pi,labels=[m7,2,mP])
        rr,status=f.roots(rrec,u,nd,ty)
        if status.startswith('unresolved')or not rr:return None
        rv=next(v for sg,v in rr if sg==rsign);nd[70]=rv;ty[70]=pi
        return s.dot(rv,pi,qv,0),nd,ty,rows
    configs=[(3,2,k)for k in labels]+[(q30,q52,k)for q30,q52 in[(2,2),(2,3)]for k in labels]
    for q30,q52,k in configs:
        assert N.P.corner3(0,0,0,4,q30,q52)is not None
        recipes=[]
        for pi in range(7):
          for m7 in range(2,6):
           for mP in range(2,6):
            if all(c is not None for c in [N.P.corner3(0,0,pi,k,2,2),
                  N.P.corner3(0,0,pi,2,2,mP),N.P.corner3(4,0,pi,3,m7,mP)]):
                if (q30,q52)==(3,2)or(pi,m7,mP)==(0,3,2):recipes.append((pi,m7,mP))
        for pi,m7,mP in recipes:
          for qs in[-1,1]:
           for rs in[-1,1]:
            key=[k,qs,pi,m7,mP,rs,q30,q52];previous=None;values=[];local=[]
            for u in grid:
                ans=state(u,*key)
                if ans is None:previous=None;counts['sample_nonreal_or_singular']+=1;continue
                value=ans[0];values.append([u,value]);counts['sample_real']+=1
                if previous is not None and previous[1]*value<0:
                    local.append([previous[0],u,previous[1],value])
                previous=u,value
            branches.append(dict(key=key,real_samples=len(values),
                value_range=[min((v for u,v in values),default=None),max((v for u,v in values),default=None)],
                closest=min(values,key=lambda r:abs(r[1]))if values else None,brackets=local))
            roots_seen=[]
            for bracket in local:
                bi=len(brackets);brackets.append(dict(key=key,interval=bracket))
                try:
                    u=float(brentq(lambda x:state(x,*key)[0],*bracket[:2],xtol=2e-14,maxiter=100))
                    ans=state(u,*key)
                except (ValueError,TypeError,StopIteration,RuntimeError)as error:
                    unresolved.append(dict(stage='root',bracket=bi,error=str(error)));continue
                if any(abs(u-v)<1e-8 for v in roots_seen):counts['duplicate_roots']+=1;continue
                roots_seen.append(u);gap,nd,ty,rows=ans
                native=[];badall=[]
                for i in[67,68,69,70]:
                    vals,bad,cp=s.native(nd[i],ty[i],{j:v for j,v in nd.items()if j!=i},ty)
                    native.append(dict(node=i,normal=nd[i].tolist(),profile=ty[i],
                        native_products=sorted(vals.items()),blockers=bad,center_product=cp,
                        norm_residual=abs(float(nd[i]@M@nd[i])-s.profiles[ty[i]]['square'])))
                    assert native[-1]['norm_residual']<2e-5
                    if bad or cp>=-1e-8:badall.append(dict(node=i,blockers=[[j,vals[j]]for j in bad],center_product=cp))
                if N.G.relation(gap)!=2:
                    unresolved.append(dict(stage='RQ_event_spacing',bracket=bi,u=u,gap=gap));continue
                if badall:
                    rejected.append(dict(bracket=bi,u=u,gap=gap,blockers=badall));continue
                sid=s.score(nd,ty,f.ids+[67,68,69,70],'continuous_LPQR')
                if sid is None:continue
                event=dict(index=len(events),bracket=bi,key=key,u=u,gap=gap,native=native,score=sid,J_bank=None,completed=[])
                events.append(event)
                if s.ledgers[sid]['open_ends']<10:improvements.append(dict(stage='LPQR',event=event['index'],score=sid))
                bank=s.bank([7,30,70],[(0,1),(0,2)],nd,ty,f.ids+[67,68,69,70],71,'J_after_continuous_LPQR')
                event['J_bank']=bank
                for j in bank['caps']:
                    full={**nd,71:np.array(j['normal'])};ft={**ty,71:j['profile']};jid=j['score']
                    corners=[s.corner(t,full,ft)for t in disk['triangles']];faces=[]
                    for i in[67,68,69,70,71]:
                        rays=[r for r in s.ledgers[jid]['rays']if i in r['active']]
                        expected={tuple(sorted(t))for t in disk['triangles']if i in t}
                        faces.append(dict(node=i,finite_exactly_expected=all(r['kind']=='finite'for r in rays)and
                            {tuple(sorted(r['active']))for r in rays}==expected))
                    done=all(c is not None for c in corners)and all(x['finite_exactly_expected']for x in faces)
                    event['completed'].append(dict(cap=j['index'],finite_corners=corners,whole_cap_faces=faces,completed_disk=done))
                    if s.ledgers[jid]['open_ends']<10:improvements.append(dict(stage='LPQRJ',event=event['index'],cap=j['index'],score=jid,completed_disk=done))
                print('CONTINUOUS_EVENT',event['index'],key,u,'ends',s.ledgers[sid]['open_ends'],
                      'J',len(bank['caps']),'completed',sum(x['completed_disk']for x in event['completed']),flush=True)
        print('CONTINUOUS_LABEL',[q30,q52,k],'branches',len(branches),'brackets',len(brackets),'events',len(events),flush=True)
    # Reproduce the actual numerical finite-contact freedom around the seed.
    freedom=[f.state(f.u0+du,False)for du in[-1e-4,0,1e-4]]
    assert all(r['all_native']for r in freedom)
    sources=dict(f.prior['sources'])
    for folder,ext in [('scripts','py'),('results','json'),('results','md')]:
        p=ROOT/f'{folder}/{S.TAG}.{ext}';sources[str(p.relative_to(ROOT))]=sha256(p.read_bytes()).hexdigest()
    sources[f'scripts/{TAG}.py']=sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(status='BOUNDED_CONTINUOUS_U_NATIVE_FIRST_SCOUT',sources=sources,
        disk=disk,deleted_walls=[47],retained_old_wall_count=62,exact_deletion_debit=f.prior['exact_deletion_debit'],
        control_score=s.control,deletion_score=deletion_score,original_seed_u=f.u0,seed_recipes=f.recipes,seed_freedom=freedom,
        necessary_high_LQ_with_old_J=completion_domain(),
        domain=dict(u_grid=grid,LQ_labels=labels,Q_profile=0,Q_frame=[30,52,67],Q_other_labels=[3,2],
          adaptive_Q_endpoint_assignments=[[2,2],[2,3]],adaptive_R_recipe=[0,3,2],
          R_frame=[7,67,68],R_L=2,R_Q=2,R7_RP_labels=[2,3,4,5],R_profiles=list(range(7)),
          R_finite_duties=['LQR','LPR','7PR'],both_Q_R_cofactor_sheets=True,
          J_frame=[7,30,70],J_duties=[[0,1],[0,2]],J_high_samples=N.HIGH),
        branches=branches,brackets=brackets,rejected_events=rejected,events=events,
        counts=dict(counts,branches=len(branches),brackets=len(brackets),native_events=len(events),
            native_J=sum(len(e['J_bank']['caps'])for e in events),completed_disks=sum(c['completed_disk']for e in events for c in e['completed'])),
        boundary_ledgers=s.ledgers,histograms=[dict(stage=a,finite=b,exterior=c,open_ends=d,count=n)for(a,b,c,d),n in sorted(s.hist.items())],
        J_bank_counts=dict(s.counts),J_near_failures=s.near,
        J_blockers=[dict(stage=a,walls=list(b),count=n)for(a,b),n in sorted(s.blockers.items())],
        improvements=improvements,unresolved=unresolved+s.unresolved,
        scope='Numerical continuation of one fixed L/P family, independent of obsolete oldR and oldLR-event u. Both Q/R cofactor sheets; R right to L/Q by imposed subdomain, necessary for LQ>=6. Strict inherited native label spacing and all62 oldwall/newwall/commoncenter gates. Complete numerical targetS6 cone and independentPD incidence for every accepted event/J. Sign-change roots on saved finite grid only; tangencies, narrow sign windows, singular limits, u outsidegrid and unsampled LQ/J labels remain open. No exact newnormal, all-parameter exclusion, full ambient vertex census or finite-volume conclusion.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--pilot',action='store_true');ap.add_argument('--native-first-pilot',action='store_true');ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    if args.pilot:pilot()
    elif args.native_first_pilot:native_first_pilot()
    else:
        out=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
        if args.verify:assert out==json.loads(p.read_text())
        else:p.write_text(json.dumps(out,separators=(',',':'))+'\n')
        print(out['status'],out['counts'],'improvements',len(out['improvements']),'unresolved',len(out['unresolved']),flush=True)
