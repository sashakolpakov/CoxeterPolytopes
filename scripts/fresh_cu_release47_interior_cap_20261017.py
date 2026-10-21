#!/usr/bin/env python3
"""Fourth-cap stellar subdivision of the failed L/R edge.

Bounded numerical representatives, chosen across every nonempty saved
(A-family,P-label,L-profile) group. K is built once per event from A/L/P,
then tested against every saved R at that event. No all-label claim.
"""
import argparse,json
from collections import Counter,defaultdict
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
import fresh_cu_release47_quadratic_events_20261017 as D
U,E,H,B,M,N,S,P,ROOT=D.U,D.E,D.H,D.B,D.M,D.N,D.S,D.P,D.ROOT
TAG='fresh_cu_release47_interior_cap_20261017'
HIGH=[6,8,12,24]

def disk():
    old=D.Q.disk_check();removed=[{68,69,71},{70,69,71}]
    fs=[f for f in old['triangles']if set(f)not in removed]
    assert len(fs)==7
    fs += [[68,69,72],[68,72,71],[70,69,72],[70,72,71]]
    es=Counter(tuple(sorted(e))for f in fs for e in combinations(f,2));vs=set(sum(fs,[]));bc=old['boundary']
    assert len(vs)==9 and len(es)==19 and len(fs)==11 and len(vs)-len(es)+len(fs)==1
    assert set(es.values())=={1,2}and(69,71)not in es
    assert{e for e,n in es.items()if n==1}=={tuple(sorted(e))for e in zip(bc,bc[1:]+bc[:1])}
    links=[]
    for v in sorted(vs):
        ee=[tuple(i for i in f if i!=v)for f in fs if v in f];dd=Counter(i for e in ee for i in e)
        seen={next(iter(dd))}
        while True:
            new=seen|{b for a,b in ee if a in seen}|{a for a,b in ee if b in seen}
            if new==seen:break
            seen=new
        assert seen==set(dd)
        assert sorted(dd.values())==([1,1]+[2]*(len(dd)-2)if v in bc else[2]*len(dd))
        links.append(dict(vertex=v,edges=ee))
    return dict(boundary=bc,triangles=fs,vertices=9,edges=19,faces=11,links=links,removed_edge=[69,71],new_interior_vertex=72)

def select(parent):
    events={e['index']:e for e in parent['events']};rr=defaultdict(list)
    for r in parent['third_caps']:rr[r['event']].append(r)
    groups=defaultdict(list)
    for eid,rs in rr.items():
        e=events[eid];pm=parent['P_seeds'][e['P']]['labels'][0];lp=parent['L_seeds'][e['L']]['profile']
        groups[e['family'],pm,lp].append(eid)
    out=[]
    for key,ids in sorted(groups.items()):
        eid=max(ids,key=lambda i:(max(dict(r['native_products'])[69]for r in rr[i]),-i))
        best=max(dict(r['native_products'])[69]for r in rr[eid])
        out.append(dict(group=list(key),available_events=len(ids),selected_event=eid,closest_RL=best,R_indices=[r['index']for r in rr[eid]]))
    assert len(out)==83 and sum(len(r['R_indices'])for r in out)==516
    return out,events,rr

def run():
    src=ROOT/'results'/f'{D.TAG}.json';parent=json.loads(src.read_text())
    selected,events,rbank=select(parent)
    cur,oldlabels,base,h,profiles,*_=U.setup();base.pop(47);base.pop(67)
    oldindex={i:j for j,i in enumerate(oldlabels)};actual=[i for i in cur['actual_sides']if i!=47]
    off=lambda a,b:float(H.A.native(np.zeros(4),a,np.zeros(4),b,profiles))
    def native(v,pi,nd,types):
        vals={i:float(v@M@w)+(profiles[pi]['offsets'][oldindex[i]]if i in oldindex else off(pi,types.get(i,0)))for i,w in nd.items()}
        bad=[i for i,x in vals.items()if B.relation(x)in[None,'parallel_or_unresolved']]
        return vals,bad,float(v@M@h)+profiles[pi]['center_offset']
    def score(nd,types):
        sel=actual+[i for i in sorted(nd)if i not in oldindex]
        bd=S.C.detailed_cone(nd,sel)
        sq=[profiles[types.get(i,0)]['square']if i not in oldindex else float(nd[i]@M@nd[i])for i in sel]
        ic=H.incidence_score([nd[i]for i in sel],sel,sq)
        good=not bd['warnings']and ic['status']=='numerically_separated_acute_incidence'
        good=good and bd['open_ends']==ic['open_ends']and bd['counts'].get('finite',0)==ic['finite']
        good=good and{tuple(sorted(r['active']))for r in bd['rays']if r['kind']=='finite'}=={tuple(sorted(r))for r in ic['finite_supports']}
        return dict(boundary=bd,PD_incidence=ic),good
    counts=Counter();blockers=Counter();unresolved=[];records=[];improvements=[];hist=Counter();near=[]
    for selection in selected:
        eid=selection['selected_event'];ev=events[eid];lr=parent['L_seeds'][ev['L']]
        types={int(i):p for i,p in ev['profiles'].items()};nd=dict(base);nd.update({int(i):np.array(v)for i,v in ev['normals'].items()})
        ap,lp,pp=[types[i]for i in[68,69,70]];al,pl=ev['AL'],lr['labels'][2]
        # Only the actual AL and LP boundary edges enter the mandatory
        # finite-corner classification; the A/P nonedge is not declared right.
        bonds=[2]*9;bonds[1]=al;bonds[5]=pl
        domain=E.recipes((ap,lp,pp),tuple(bonds),((0,1),(1,2)))
        assert abs(float(H.A.native(nd[68],ap,nd[69],lp,profiles))+B.theta(al))<2e-5
        assert abs(float(H.A.native(nd[70],pp,nd[69],lp,profiles))+B.theta(pl))<2e-5
        for i in[67,68,69,70]:
            _,bad,cp=native(nd[i],types.get(i,0),{j:w for j,w in nd.items()if j!=i},types)
            assert not bad and cp<-1e-8
        caps=[];joins=[];seen=set();closest=[]
        for kp,low,hi in domain:
            for hh in product(HIGH,repeat=len(hi)):
                ms=list(low)
                for j,m in zip(hi,hh):ms[j]=m
                ids=[68,69,70];rhs=[-B.theta(m)-off(kp,types[i])for i,m in zip(ids,ms)]
                vv,info=N.norm_roots([nd[i]for i in ids],rhs,profiles[kp]['square']);counts['K_equations']+=1
                if info['status'].startswith('unresolved'):
                    unresolved.append(dict(stage='K_norm',event=eid,profile=kp,labels=ms,norm=info));continue
                if not vv:counts['K_no_real_equations']+=1
                for sg,k in vv:
                    counts['K_branches']+=1
                    err=max(abs(float(k@M@k)-profiles[kp]['square']),max(abs(float(k@M@nd[i])-r)for i,r in zip(ids,rhs)))
                    if err>2e-5:unresolved.append(dict(stage='K_residual',event=eid,profile=kp,labels=ms,error=err));continue
                    vals,bad,cp=native(k,kp,nd,types)
                    if bad or cp>=-1e-8:
                        blockers['K',tuple(bad)]+=1
                        if len(bad)<=2 and len(near)<100:near.append(dict(event=eid,profile=kp,labels=ms,branch=sg,blockers=[[i,vals[i]]for i in bad],center_product=cp))
                        continue
                    key=(kp,)+tuple(np.round(k,8))
                    if key in seen:counts['K_duplicates']+=1;continue
                    seen.add(key);nn=dict(nd);nn[72]=k;ty=dict(types);ty[72]=kp
                    sc,ok=score(nn,ty)
                    if not ok:unresolved.append(dict(stage='K_cone',event=eid,**sc));continue
                    rec=dict(index=len(caps),profile=kp,labels=ms,branch=sg,normal=k.tolist(),native_products=sorted(vals.items()),center_product=cp,norm_residual=err,score=D.compact_score(sc))
                    caps.append(rec);counts['native_K_caps']+=1
                    bd=sc['boundary'];hist['K',bd['counts']['finite'],bd['counts']['exterior'],bd['open_ends']]+=1
                    if bd['open_ends']<10:improvements.append(dict(stage='K',event=eid,cap=rec,complete=sc))
                    for rr in rbank[eid]:
                        r=np.array(rr['normal']);rp=rr['profile'];kr=float(H.A.native(k,kp,r,rp,profiles));rel=B.relation(kr)
                        counts['KR_pairs']+=1
                        if rel in[None,'parallel_or_unresolved']:
                            blockers['KR',(rel,)]+=1
                            if len(closest)<10:closest.append(dict(K=rec['index'],R=rr['index'],native_product=kr,relation=rel))
                            continue
                        if isinstance(rel,int):counts['finite_KR_pairs']+=1
                        else:counts['dotted_KR_pairs']+=1
                        # Check all new-wall inequalities once again in the
                        # assembled chamber, rather than inferring them from
                        # separate records or from a projected-product test.
                        both=dict(nn);both[71]=r;bt=dict(ty);bt[71]=rp
                        _,kb,kc=native(k,kp,{i:w for i,w in both.items()if i!=72},bt)
                        _,rb,rc=native(r,rp,{i:w for i,w in both.items()if i!=71},bt)
                        assert not kb and not rb and kc<-1e-8 and rc<-1e-8
                        ar,pr=rr['labels'][1],rr['labels'][2]
                        complete=isinstance(rel,int)and P.corner3(ap,rp,kp,ar,ms[0],rel)is not None and P.corner3(pp,rp,kp,pr,ms[2],rel)is not None
                        sc2,good=score(both,bt)
                        if not good:unresolved.append(dict(stage='KR_cone',event=eid,K=rec['index'],R=rr['index'],**sc2));continue
                        krays=[v for v in sc2['boundary']['rays']if 72 in v['active']]
                        expected={tuple(sorted(z))for z in[[68,69,72],[68,71,72],[70,69,72],[70,71,72]]}
                        whole_k={tuple(sorted(v['active']))for v in krays}==expected and all(v['kind']=='finite'for v in krays)
                        if complete:assert whole_k
                        jr=dict(K=rec['index'],R=rr['index'],native_product=kr,relation=rel,remaining_required_corners_finite=complete,whole_compact_K_quad=whole_k,score=D.compact_score(sc2))
                        joins.append(jr);counts['coherent_KR_joins']+=1
                        if complete:counts['completed_four_cap_disks']+=1
                        bd=sc2['boundary'];hist['KR',bd['counts']['finite'],bd['counts']['exterior'],bd['open_ends']]+=1
                        if bd['open_ends']<10:improvements.append(dict(stage='KR',event=eid,K=rec,R=rr,join=jr,complete=sc2))
        records.append(dict(selection=selection,A_profile=ap,L_profile=lp,P_profile=pp,AL=al,LP=pl,
          symbolic_domain=[list(z)for z in domain],caps=caps,joins=joins,first_rejected_KR=closest))
        print('CU_INTERIOR_K',eid,'group',selection['group'],'recipes',len(domain),'K',len(caps),'joins',len(joins),'improved',len(improvements),'unresolved',len(unresolved),flush=True)
    paths=[Path(__file__),Path(D.__file__),src,Path(U.__file__),Path(E.__file__),Path(N.__file__),Path(P.__file__),H.PARENT,H.BOUNDARY]
    return dict(status='BOUNDED_RELEASE47_INTERIOR_FOURTH_CAP',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
      topology=disk(),deleted_walls=[47],retained_control_walls=62,exact_deletion_ledger=parent['exact_deletion_ledger'],
      selection_rule='One event per nonempty(Afamily,P7label,Lprofile)group having savedR; select maximal nativeR.L, then test ALL savedR on that event.',
      selected_events=len(selected),selected_A_families=sorted({r['group'][0]for r in selected}),selected_R_caps=sum(len(r['R_indices'])for r in selected),
      high_label_samples=HIGH,counts=dict(counts),records=records,improvements=improvements,unresolved=unresolved,near_native_failures=near,
      histograms=[dict(stage=s,finite=f,exterior=e,open_ends=i,count=n)for(s,f,e,i),n in sorted(hist.items())],
      blockers=[dict(stage=k[0],walls=list(k[1]),count=n)for k,n in sorted(blockers.items(),key=lambda z:str(z[0]))],
      scope='Numerical fourth-cap topology with LR removed, not a rerun of the failed LR contact.83 representatives from all83 nonempty(Afamily,P7label,Lprofile)groups;18 A families with savedR,62 changedP4/P5 events and21oldP3controls,516Rcaps. Seven A2profiles and both norm branches; isolated highKlabels only6/8/12/24. Complete native/commoncenter checks for all62old/newwalls and completeS6cone/independentPDsupports for every coherentK andKR pair, including dotted joins. Other events/highlabels, exact root certification, ambient finite/exterior census and finite volume remain open. No old frozen geometry is changed.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],'improvements',len(out['improvements']),'unresolved',len(out['unresolved']),flush=True)
