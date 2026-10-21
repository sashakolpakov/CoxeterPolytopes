#!/usr/bin/env python3
"""Bounded double-flip discovery: XY->8Z and CZ->XW.

The second flip removes Z's finite-C duty. No exact chamber or all-label
claim follows from this numerical discovery program.
"""
import argparse,json
from collections import Counter,defaultdict
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
import fresh_old8_triangle_flip_20261017 as D
O,P,F,T,N,S,H,B,M,ROOT=D.O,D.P,D.F,D.T,D.N,D.S,D.H,D.B,D.M,D.ROOT
TAG='fresh_old8_triangle_double_flip_20261017'

def topology():
    old=D.topology();tri=old['triangles'];remove=[frozenset(t)for t in[['X','C','Z'],['C','Z','W']]]
    assert all(sum(frozenset(t)==s for t in tri)==1 for s in remove)
    assert not any({'X','W'}<=set(t)for t in tri)
    tri=[t for t in tri if frozenset(t)not in remove]+[['X','C','W'],['X','Z','W']]
    boundary=old['boundary'];new=list('ARCXYBZW');vertices=boundary+new
    edges=Counter(tuple(sorted(e))for t in tri for e in combinations(t,2))
    be={tuple(sorted(e))for e in zip(boundary,boundary[1:]+boundary[:1])}
    assert len(tri)==len({frozenset(t)for t in tri})==19 and len(edges)==31
    assert {e for e,n in edges.items()if n==1}==be
    assert all(n==(1 if e in be else 2)for e,n in edges.items())
    for v in vertices:
        link=defaultdict(set)
        for t in tri:
            if v in t:
                a,b=[x for x in t if x!=v];link[a].add(b);link[b].add(a)
        todo=[next(iter(link))];seen=set()
        while todo:
            x=todo.pop()
            if x not in seen:seen.add(x);todo.extend(link[x]-seen)
        assert len(seen)==len(link)
        assert sorted(map(len,link.values()))==([1,1]+[2]*(len(link)-2)if v in boundary else[2]*len(link))
    assert Counter(sum(v in new for v in t)for t in tri)=={1:5,2:8,3:6}
    assert not any({'C','Z'}<=set(t)or{'X','Y'}<=set(t)for t in tri)
    return dict(old,triangles=tri,second_flip='CZ->XW',conditional_end_change=-5)

def run():
    top=topology();cur,labels,ns,h,pro,primary=H.setup();actual=[i for i in cur['actual_sides']if i not in S.DELETE]
    parentpath=ROOT/'results'/f'{D.TAG}.json';bpath=ROOT/'results/fresh_two_connector_coupled_disk_20261017.json'
    parent=json.loads(parentpath.read_text());bank=json.loads(bpath.read_text());acs=parent['AC_bank']
    counts=Counter();blockers=Counter();obstructions=Counter();W_stops=Counter();WZ_samples=[];unresolved=[];prefixes=[];complete=[];rs={};acute_samples=[]
    def nativeprod(v,pi,w,pj):return float(H.A.native(v,pi,w,pj,pro))
    def offset(pi,pj):return nativeprod(np.zeros(4),pi,np.zeros(4),pj)
    def dotted(v,pi,w,pj):return float(v@M@w)<-np.sqrt(pro[pi]['square']*pro[pj]['square'])-1e-7
    def nonincident(v,pi,allow):return all(float(v@M@ns[i])<-np.sqrt(pro[pi]['square']*float(ns[i]@M@ns[i]))-1e-7 for i in actual if i not in allow)
    def finite(v,pi,w,pj):return F.label(nativeprod(v,pi,w,pj))
    def pairchecks(vs,ps):
        out=[dict(pair=[i,j],product=nativeprod(vs[i],ps[i],vs[j],ps[j]))for i,j in combinations(range(len(vs)),2)]
        return out if all(B.relation(z['product'])not in[None,'parallel_or_unresolved']for z in out)else None
    def solve(stage,frame,rhs,pi,assigned=(),context=None):
        counts[stage+'_equations']+=1;roots,info=N.norm_roots(frame,rhs,pro[pi]['square'])
        if info['status'].startswith('unresolved'):
            unresolved.append(dict(stage=stage,context=context,info=info));return []
        out=[]
        for sg,v in roots:
            counts[stage+'_branches']+=1
            error=max(abs(float(v@M@v)-pro[pi]['square']),max(abs(float(v@M@w)-q)for w,q in zip(frame,rhs)))
            if error>1e-5:
                unresolved.append(dict(stage=stage,context=context,branch=sg,error=error));continue
            vals=T.native_values(v,pi,labels,ns,pro);bad=T.gates(vals,assigned)
            if bad:
                blockers[stage,bad[0][0]]+=1
                if stage=='Z':
                    acute=all(t<=1e-7 for t in vals.values());obstructions['acute_angle_gap'if acute else'obtuse']+=1
                    if acute and len(acute_samples)<20:acute_samples.append(dict(context=context,branch=sg,bad=bad))
                continue
            if v@M@h+pro[pi]['center_offset']>=-1e-8:counts[stage+'_center_reject']+=1;continue
            counts[stage+'_native']+=1;out.append((sg,v))
        return out
    def boundary(vs):
        nd=dict(ns);nd.update({67+i:v for i,v in enumerate(vs)});sel=actual+list(range(67,67+len(vs)))
        bd=S.C.detailed_cone(nd,sel);ic=H.incidence_score([nd[i]for i in sel],sel,[float(nd[i]@M@nd[i])for i in sel])
        if bd['warnings']or ic['status']!='numerically_separated_acute_incidence'or bd['open_ends']!=ic['open_ends']or bd['counts'].get('finite',0)!=ic['finite']:
            unresolved.append(dict(stage='boundary',bd=bd,ic=ic));return None
        actual_finite={tuple(sorted(r['active']))for r in bd['rays']if r['kind']=='finite'}
        if actual_finite!={tuple(sorted(s))for s in ic['finite_supports']}:
            unresolved.append(dict(stage='boundary_finite_supports',bd=bd,ic=ic));return None
        return bd,ic
    for mr in bank['R44_labels']:
        rr=solve('R',[ns[8],ns[44],ns[63]],[0.,-B.theta(mr),0.],0,[8,44,63],mr);assert len(rr)==1;rs[mr]=rr[0][1]
    for ia,ac in enumerate(acs):
        a,c=map(np.array,ac['normals'])
        for ix,xr in enumerate(parent['X_bank'].get(str(ia),[])):
            x=np.array(xr['normal']);xp=xr['profile'];m8x,maxx,mcx=xr['labels']
            for mr,r in rs.items():
                if not all(dotted(v,pi,r,0)for v,pi in[(a,0),(c,0),(x,xp)]):continue
                counts['coherent_ARCX']+=1
                for zp,m8z,mrz,mxz in product(range(7),D.Z8R,D.Z8R,D.LOW):
                    if P.corner3(0,xp,zp,m8x,m8z,mxz)is None:continue
                    if P.corner3(0,0,zp,2,m8z,2)is None or P.corner3(0,0,zp,2,mrz,2)is None:continue
                    rem=pro[zp]['square']-(B.theta(m8z)**2+B.theta(mrz)**2)/2
                    if abs(rem)<1e-9:counts['zero_polar_limit_skipped']+=1;continue
                    if rem>0:counts['positive_polar_frame_skipped']+=1;continue
                    rhs=[-B.theta(m8z),-B.theta(mrz),-B.theta(mxz)-offset(zp,xp)]
                    for zsg,z in solve('Z',[ns[8],r,x],rhs,zp,[8],[ia,ix,mr,zp,m8z,mrz,mxz]):
                        if not nonincident(z,zp,[8])or not dotted(z,zp,a,0)or not dotted(z,zp,c,0):continue
                        counts['Z_intrinsic']+=1
                        for ysg,y in solve('Y',[ns[8],r,z],[0.,0.,0.],0,[8],[ia,ix,mr,zp,m8z,mrz,mxz,zsg]):
                            if not nonincident(y,0,[8])or not all(dotted(y,0,v,pi)for v,pi in[(a,0),(c,0),(x,xp)]):continue
                            vs=[a,r,c,x,y,z];ps=[0,0,0,xp,0,zp];pairs=pairchecks(vs,ps)
                            if pairs is None:continue
                            out=boundary(vs)
                            if out is None:continue
                            bd,ic=out;ys={tuple(sorted(s))for s in ic['finite_supports']if 71 in s}
                            expected={tuple(sorted(s))for s in[[8,68,71],[8,71,72],[68,71,72]]}
                            if ys!=expected:
                                unresolved.append(dict(stage='Y_actual_triangle',actual=sorted(ys),expected=sorted(expected)));continue
                            ip=len(prefixes);prefixes.append(dict(AC=ia,X=ix,R44=mr,Z_profile=zp,Z_labels=[m8z,mrz,mxz],branches=[zsg,ysg],normals=[v.tolist()for v in vs],profiles=ps,native_pairs=pairs,boundary=bd,PD_incidence=ic,Y_actual_finite_corners=sorted(ys)))
                            counts['coherent_six_prefixes']+=1
                            for wp,mw30,mwc,mwx in product(range(7),D.LOW,D.LOW,D.LOW):
                                if P.corner3(0,0,wp,2,mw30,mwc)is None or P.corner3(0,xp,wp,mcx,mwc,mwx)is None:continue
                                wrhs=[-B.theta(mw30),-B.theta(mwc),-B.theta(mwx)-offset(wp,xp)]
                                for wsg,w in solve('W',[ns[30],c,x],wrhs,wp,[30],[ip,wp,mw30,mwc,mwx]):
                                    if not nonincident(w,wp,[30]):W_stops['old_intrinsic_nonincidence']+=1;continue
                                    other=[name for name,v,pi in[('A',a,0),('R',r,0),('Y',y,0)]if not dotted(w,wp,v,pi)]
                                    if other:W_stops['new_nonincidence_'+','.join(other)]+=1;continue
                                    counts['W_intrinsic']+=1
                                    mwz=finite(w,wp,z,zp)
                                    if mwz is None:
                                        value=nativeprod(w,wp,z,zp);W_stops['WZ_'+str(B.relation(value))]+=1
                                        WZ_samples.append(dict(prefix=ip,W_profile=wp,W_labels=[mw30,mwc,mwx],branch=wsg,product=value))
                                        continue
                                    if P.corner3(xp,zp,wp,mxz,mwx,mwz)is None:W_stops['nonfinite_XZW_corner']+=1;continue
                                    if pairchecks(vs+[w],ps+[wp])is None:continue
                                    counts['coherent_seven_prefixes']+=1
                                    for ib,br in enumerate(bank['B_bank']):
                                        if br['R44']!=mr:continue
                                        b=np.array(br['B_projection']);bp=br['B_profile'];m30b,m63b,mrb=br['B_labels']
                                        if not all(dotted(b,bp,v,pi)for v,pi in[(a,0),(c,0),(x,xp),(y,0)]):continue
                                        mbz=finite(b,bp,z,zp);mbw=finite(b,bp,w,wp)
                                        if mbz is None or mbw is None:continue
                                        if P.corner3(0,bp,zp,mrb,mrz,mbz)is None or P.corner3(0,bp,wp,m30b,mw30,mbw)is None or P.corner3(bp,zp,wp,mbz,mbw,mwz)is None:continue
                                        allvs=[a,r,c,x,y,b,z,w];allps=[0,0,0,xp,0,bp,zp,wp];checks=pairchecks(allvs,allps)
                                        if checks is None:continue
                                        out=boundary(allvs)
                                        if out is None:continue
                                        bd,ic=out;mp={s:int(s)for s in top['boundary']}|dict(zip(list('ARCXYBZW'),range(67,75)))
                                        supports={tuple(sorted(s))for s in ic['finite_supports']}
                                        assert all(tuple(sorted(mp[s]for s in tri))in supports for tri in top['triangles'])
                                        complete.append(dict(prefix=ip,B=ib,W_profile=wp,W_labels=[mw30,mwc,mwx,mwz,mbw],W_branch=wsg,W=w.tolist(),native_pairs=checks,boundary=bd,PD_incidence=ic))
        print('DOUBLE_FLIP',ia,dict(counts),'prefixes',len(prefixes),'complete',len(complete),flush=True)
    paths=[Path(__file__),Path(D.__file__),parentpath,bpath,H.PARENT,H.BOUNDARY]
    return dict(status='BOUNDED_XY_8Z_CZ_XW_DOUBLE_FLIP_DISCOVERY',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},topology=top,counts=dict(counts),blockers=[[list(k),v]for k,v in sorted(blockers.items())],Z_obstruction_types=dict(obstructions),acute_Z_samples=acute_samples,W_stops=dict(W_stops),WZ_samples=WZ_samples,prefixes=prefixes,complete=complete,unresolved=unresolved,
      domain=dict(AC_events=27,X_events=139,R44=bank['R44_labels'],Z_profiles=list(range(7)),Z8_ZR=D.Z8R,ZX=D.LOW,Y='rightS6+8+R+Z, both polar signs',W_profiles=list(range(7)),W_labels=D.LOW,B_bank=len(bank['B_bank'])),deleted_walls=S.DELETE,exact_old_deletion_cost=parent['exact_old_deletion_cost'],
      scope='Exact abstract double-flip disk; bounded numerical geometry only. The finite C-Z duty is removed. Both norm branches, all retained native contacts, common fixed center, intrinsic nonincidences, full new pairs, actual Y finite corner rays and complete numerical cone+independent PD census are required for accepted six-prefixes. All19 finite corners are required for a complete patch. Zero polar limits are excluded only from this recipe. No all-label exclusion, exact geometry promotion, or finite-volume result is claimed.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=json.loads(json.dumps(run()));p=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],'prefixes',len(out['prefixes']),'complete',len(out['complete']),'unresolved',len(out['unresolved']),flush=True)
