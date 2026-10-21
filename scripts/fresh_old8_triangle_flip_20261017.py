#!/usr/bin/env python3
"""Flip XY to old8-Z; construct Y's triangular facet jointly with Z.

Bounded numerical discovery with full native and intrinsic gates. Frozen
control and predecessor files are unchanged.
"""
import argparse,json
from collections import Counter,defaultdict
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
import fresh_c10_endpoint_disk_20261017 as E
O,P,F,T,N,S,H,B,M,ROOT=E.O,E.P,E.F,E.T,E.N,E.S,E.H,E.B,E.M,E.ROOT
TAG='fresh_old8_triangle_flip_20261017'
Z8R=[2,3,4,5,6,8,12]
LOW=list(range(2,7))

def topology():
    old=O.topology(True);tri=[list(t)for t in old['triangles']]
    rem=[frozenset(t)for t in[['8','X','Y'],['X','Y','Z']]]
    assert all(sum(frozenset(t)==s for t in tri)==1 for s in rem)
    assert not any({'8','Z'}<=set(t)for t in tri)
    tri=[t for t in tri if frozenset(t)not in rem]+[['8','X','Z'],['8','Z','Y']]
    boundary=old['boundary'];new=list('ARCXYBZW');vs=boundary+new
    edge=Counter(tuple(sorted(e))for t in tri for e in combinations(t,2));be={tuple(sorted(e))for e in zip(boundary,boundary[1:]+boundary[:1])}
    assert len(tri)==len({frozenset(t)for t in tri})==19
    assert {e for e,k in edge.items()if k==1}==be and all(k==(1 if e in be else 2)for e,k in edge.items())
    for v in vs:
        lk=defaultdict(set)
        for t in tri:
            if v in t:
                a,b=[w for w in t if w!=v];lk[a].add(b);lk[b].add(a)
        todo=[next(iter(lk))];seen=set()
        while todo:
            a=todo.pop()
            if a not in seen:seen.add(a);todo.extend(lk[a]-seen)
        assert len(seen)==len(lk)
        assert sorted(map(len,lk.values()))==([1,1]+[2]*(len(lk)-2)if v in boundary else[2]*len(lk))
    xyz=Counter(sum(v in new for v in t)for t in tri)
    assert xyz=={1:5,2:8,3:6}
    assert len(vs)==13 and len(edge)==31 and 13-31+19==1
    yn={w for t in tri if'Y'in t for w in t}-{'Y'};assert yn=={'8','R','Z'}
    return dict(boundary=boundary,triangles=tri,vertices=13,edges=31,triangles_count=19,x=5,y=8,z=6,
        old_new_edges=13,new_new_edges=13,conditional_end_change=2*(13+13)-3*19,Y_neighbors=sorted(yn),
        scope='Exact abstract disk; actual embedded geometry and retention of the outside30--44 face require separate gates.')

def run():
    top=topology();cur,labels,ns,h,pro,primary=H.setup();actual=[i for i in cur['actual_sides']if i not in S.DELETE]
    ep=ROOT/'results'/f'{E.TAG}.json';bpath=ROOT/'results/fresh_two_connector_coupled_disk_20261017.json'
    events=json.loads(ep.read_text());bank=json.loads(bpath.read_text());acs=events['AC_bank']
    counts=Counter();badcounts=Counter();obstruction_types=Counter();acute_samples=[];unresolved=[];xbank=defaultdict(list);yzbank=[];prefixes=[];complete=[];rs={}
    def nativeprod(v,pi,w,pj):return float(H.A.native(v,pi,w,pj,pro))
    def offset(pi,pj):return nativeprod(np.zeros(4),pi,np.zeros(4),pj)
    def native(v,pi,assigned=()):
        bad=T.gates(T.native_values(v,pi,labels,ns,pro),assigned)
        if bad:return [('old'+str(i),z)for i,z in bad]
        if v@M@h+pro[pi]['center_offset']>=-1e-8:return [('center',0.)]
        return []
    def dotted(v,pi,w,pj):return float(v@M@w)<-np.sqrt(pro[pi]['square']*pro[pj]['square'])-1e-7
    def nonincident(v,pi,allow):return all(float(v@M@ns[i])<-np.sqrt(pro[pi]['square']*float(ns[i]@M@ns[i]))-1e-7 for i in actual if i not in allow)
    def solve(name,frame,rhs,pi,assigned=(),context=None):
        counts[name+'_equations']+=1;rr,info=N.norm_roots(frame,rhs,pro[pi]['square'])
        if info['status'].startswith('unresolved'):
            unresolved.append(dict(stage=name,context=context,info=info));return []
        out=[]
        for sg,v in rr:
            counts[name+'_branches']+=1
            err=max(abs(float(v@M@v)-pro[pi]['square']),max(abs(float(v@M@w)-t)for w,t in zip(frame,rhs)))
            if err>1e-5:
                unresolved.append(dict(stage=name,context=context,branch=sg,error=err));continue
            bad=native(v,pi,assigned)
            if bad:
                badcounts[name,bad[0][0]]+=1
                if name=='Z':
                    vals=T.native_values(v,pi,labels,ns,pro)
                    obtuse=[(i,z)for i,z in vals.items()if z>1e-7]
                    obstruction_types['obtuse'if obtuse else'acute_angle_gap']+=1
                    if not obtuse and len(acute_samples)<20:
                        acute_samples.append(dict(context=context,branch=sg,bad=bad,normal=v.tolist(),center=float(v@M@h+pro[pi]['center_offset'])))
                continue
            counts[name+'_native']+=1;out.append((sg,v))
        return out
    def finite(v,pi,w,pj):return F.label(nativeprod(v,pi,w,pj))
    def allpairs(vs,ps):
        out=[dict(pair=[i,j],product=nativeprod(vs[i],ps[i],vs[j],ps[j]))for i,j in combinations(range(len(vs)),2)]
        return out if all(B.relation(r['product'])not in[None,'parallel_or_unresolved']for r in out)else None
    for mr in bank['R44_labels']:
        roots=solve('R',[ns[8],ns[44],ns[63]],[0.,-B.theta(mr),0.],0,[8,44,63],mr);assert len(roots)==1;rs[mr]=roots[0][1]
    # X has no Y duty after the flip. Rebuild its whole bounded A/C domain.
    for ia,ac in enumerate(acs):
        a,c=map(np.array,ac['normals']);mac=ac['AC']
        for xp,m8,ma,mc in product(range(7),LOW,LOW,LOW):
            if T.corner(xp,m8,ma,0,2)is None or T.corner(xp,ma,mc,0,mac)is None:continue
            for sg,x in solve('X',[ns[8],a,c],[-B.theta(m8),-B.theta(ma),-B.theta(mc)],xp,[8],[ia,xp,m8,ma,mc]):
                if not nonincident(x,xp,[8]):continue
                xbank[ia].append(dict(profile=xp,labels=[m8,ma,mc],branch=sg,normal=x.tolist()))
    cgroups=defaultdict(list)
    for ia,ac in enumerate(acs):cgroups[ac['C10']].append(ia)
    for c10,indices in sorted(cgroups.items()):
        c=np.array(acs[indices[0]]['normals'][1]);assert all(np.max(abs(c-np.array(acs[i]['normals'][1])))<1e-8 for i in indices)
        for mr,r in rs.items():
            for zp,m8z,mrz,mcz in product(range(7),Z8R,Z8R,LOW):
                if P.corner3(0,0,zp,2,m8z,2)is None or P.corner3(0,0,zp,2,mrz,2)is None or P.corner3(0,0,zp,2,mcz,2)is None:continue
                # The projected8/R span is positive with Gram2I. A negative
                # Schur complement makes8/R/Z Lorentzian, so its polar Y
                # has positive square. Zero boundary is left unresolved.
                rem=pro[zp]['square']-(B.theta(m8z)**2+B.theta(mrz)**2)/2
                if abs(rem)<1e-9:counts['zero_polar_limit_skipped']+=1;continue
                if rem>0:counts['positive_8RZ_frame_skipped']+=1;continue
                for zsg,z in solve('Z',[ns[8],r,c],[-B.theta(m8z),-B.theta(mrz),-B.theta(mcz)],zp,[8],[c10,mr,zp,m8z,mrz,mcz]):
                    if not nonincident(z,zp,[8]):continue
                    for ysg,y in solve('Y',[ns[8],r,z],[0.,0.,0.],0,[8],[c10,mr,zp,m8z,mrz,mcz,zsg]):
                        if not nonincident(y,0,[8])or not dotted(y,0,c,0):continue
                        if allpairs([r,c,z,y],[0,0,zp,0])is None:continue
                        counts['native_compact_Y_triangle_pairs']+=1;iz=len(yzbank)
                        yzbank.append(dict(C10=c10,R44=mr,Z_profile=zp,Z_labels=[m8z,mrz,mcz],branches=[zsg,ysg],Y=y.tolist(),Z=z.tolist(),projected_8RZ_schur=rem))
                        for ia in indices:
                            a=np.array(acs[ia]['normals'][0])
                            if not all(dotted(a,0,v,pi)for v,pi in[(r,0),(z,zp),(y,0)]):continue
                            for xr in xbank[ia]:
                                x=np.array(xr['normal']);xp=xr['profile'];m8x,maxx,mcx=xr['labels']
                                if not dotted(x,xp,r,0)or not dotted(x,xp,y,0):continue
                                mxz=finite(x,xp,z,zp)
                                if mxz is None:counts['XZ_no_finite_angle']+=1;continue
                                if P.corner3(0,xp,zp,m8x,m8z,mxz)is None or P.corner3(0,xp,zp,mcx,mcz,mxz)is None:continue
                                vs=[a,r,c,x,y,z];ps=[0,0,0,xp,0,zp];pairs=allpairs(vs,ps)
                                if pairs is None:continue
                                counts['coherent_six_prefixes']+=1;ip=len(prefixes)
                                nd=dict(ns);nd.update({67+i:v for i,v in enumerate(vs)});sel=actual+list(range(67,73))
                                bd=S.C.detailed_cone(nd,sel);ic=H.incidence_score([nd[i]for i in sel],sel,[float(nd[i]@M@nd[i])for i in sel])
                                if bd['warnings']or ic['status']!='numerically_separated_acute_incidence'or bd['open_ends']!=ic['open_ends']:
                                    unresolved.append(dict(stage='six_prefix_cone',AC=ia,YZ=iz,bd=bd,ic=ic));continue
                                prefixes.append(dict(AC=ia,YZ=iz,X=xr,XZ=mxz,normals=[v.tolist()for v in vs],profiles=ps,native_pairs=pairs,boundary=bd,PD_incidence=ic))
                                # Finish B and W only after this changed six-prefix passes.
                                for ib,br in enumerate(bank['B_bank']):
                                    if br['R44']!=mr:continue
                                    b=np.array(br['B_projection']);bp=br['B_profile'];m30b,m63b,mrb=br['B_labels']
                                    if not all(dotted(b,bp,v,pi)for v,pi in[(a,0),(c,0),(x,xp),(y,0)]):continue
                                    mbz=finite(b,bp,z,zp)
                                    if mbz is None or P.corner3(0,bp,zp,mrb,mrz,mbz)is None:continue
                                    if allpairs(vs+[b],ps+[bp])is None:continue
                                    counts['coherent_seven_prefixes']+=1
                                    for wp,mw30,mwc,mwb in product(range(7),LOW,LOW,LOW):
                                        if P.corner3(0,0,wp,2,mw30,mwc)is None or P.corner3(0,bp,wp,m30b,mw30,mwb)is None:continue
                                        rhs=[-B.theta(mw30),-B.theta(mwc),-B.theta(mwb)-offset(wp,bp)]
                                        for wsg,w in solve('W',[ns[30],c,b],rhs,wp,[30],[ip,ib,wp,mw30,mwc,mwb]):
                                            if not nonincident(w,wp,[30])or not all(dotted(w,wp,v,pi)for v,pi in[(a,0),(r,0),(x,xp),(y,0)]):continue
                                            mwz=finite(w,wp,z,zp)
                                            if mwz is None or P.corner3(bp,zp,wp,mbz,mwb,mwz)is None or P.corner3(0,zp,wp,mcz,mwc,mwz)is None:continue
                                            allvs=[a,r,c,x,y,b,z,w];allps=[0,0,0,xp,0,bp,zp,wp];checks=allpairs(allvs,allps)
                                            if checks is None:continue
                                            nd=dict(ns);nd.update({67+i:v for i,v in enumerate(allvs)});sel=actual+list(range(67,75))
                                            bd=S.C.detailed_cone(nd,sel);ic=H.incidence_score([nd[i]for i in sel],sel,[float(nd[i]@M@nd[i])for i in sel])
                                            if bd['warnings']or ic['status']!='numerically_separated_acute_incidence'or bd['open_ends']!=ic['open_ends']:
                                                unresolved.append(dict(stage='full_cone',prefix=ip,B=ib,bd=bd,ic=ic));continue
                                            mp={str(i):int(i) for i in top['boundary']}|dict(zip(list('ARCXYBZW'),range(67,75)))
                                            supports={tuple(sorted(v))for v in ic['finite_supports']}
                                            assert all(tuple(sorted(mp[i]for i in tri))in supports for tri in top['triangles'])
                                            complete.append(dict(prefix=ip,B=ib,W_profile=wp,W_labels=[mw30,mwc,mwb,mwz],W_branch=wsg,W=w.tolist(),native_pairs=checks,boundary=bd,PD_incidence=ic))
        print('OLD8_TRIANGLE',c10,dict(counts),'pairs',len(yzbank),'prefixes',len(prefixes),'complete',len(complete),flush=True)
    paths=[Path(__file__),Path(E.__file__),Path(O.__file__),Path(P.__file__),Path(T.__file__),Path(N.__file__),Path(H.__file__),ep,bpath,H.PARENT,H.BOUNDARY]
    return dict(status='BOUNDED_XY_TO_OLD8Z_TRIANGLE_FACET_DISCOVERY',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},topology=top,counts=dict(counts),blockers=[[list(k),v]for k,v in sorted(badcounts.items())],Z_obstruction_types=dict(obstruction_types),acute_Z_samples=acute_samples,AC_bank=acs,X_bank={str(k):v for k,v in xbank.items()},YZ_bank=yzbank,prefixes=prefixes,complete=complete,unresolved=unresolved,
      domain=dict(C10_AC='27 frozen native endpoint events',R44=bank['R44_labels'],Z8_ZR=Z8R,ZC=LOW,Y='zero profile; right8,R,Z',X_W_profiles=list(range(7)),X_W_labels=LOW,Z_profiles=list(range(7))),
      deleted_walls=S.DELETE,exact_old_deletion_cost=events['exact_old_deletion_cost'],
      scope='New abstractXY-to8Z flip, bounded numerical Y/Z triangle construction. Both normbranches, native contacts, fixedcenter, prescribed nonincidences and full newpair gates checked. Degenerate8/R/Z limit is outside the positive-polar triangle recipe, not rejected for alltopologies. No all-label coverage or exact promotion; every complete prefix gets numerical cone+independentPD checks before any improvement claim. No finite-volume example unless further exact ambient certification succeeds.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=json.loads(json.dumps(run()));p=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],'complete',len(out['complete']),'unresolved',len(out['unresolved']),flush=True)
