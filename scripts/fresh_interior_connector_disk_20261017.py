#!/usr/bin/env python3
"""Bounded changed-incidence seven-cap disk, with native-gate fail safety.

Subdivide X--B in the previously flipped disk.  The new interior Z removes
that contact, allowing X to be constructed at A/C instead of B/Y.
"""
import argparse,json
from collections import Counter,defaultdict
from hashlib import sha256
from functools import lru_cache
from itertools import combinations,product
from pathlib import Path
import numpy as np
import fresh_two_connector_coupled_disk_20261017 as F
T,N,S,H,B,M,ROOT=F.T,F.N,F.S,F.H,F.B,F.M,F.ROOT
TAG='fresh_interior_connector_disk_20261017'
ALPHABET=range(2,6)

def disk():
    boundary=['7','8','44','63','30'];new=list('ARCXYBZ')
    triangles=[['7','8','A'],['8','44','R'],['44','63','R'],['63','30','B'],['30','7','C'],
     ['8','A','X'],['8','X','Y'],['8','Y','R'],['63','R','B'],['30','B','C'],['7','C','A'],
     ['A','X','C'],['X','Z','C'],['Z','B','C'],['X','Y','Z'],['Y','Z','B'],['Y','R','B']]
    e=Counter(tuple(sorted(p))for t in triangles for p in combinations(t,2));be={tuple(sorted(p))for p in zip(boundary,boundary[1:]+boundary[:1])}
    assert {p for p,c in e.items()if c==1}==be and all(c==(1 if p in be else 2)for p,c in e.items())
    for v in boundary+new:
        lk=defaultdict(set)
        for t in triangles:
            if v in t:
                a,b=[i for i in t if i!=v];lk[a].add(b);lk[b].add(a)
        todo=[next(iter(lk))];seen=set()
        while todo:
            a=todo.pop()
            if a not in seen:seen.add(a);todo.extend(lk[a]-seen)
        assert len(seen)==len(lk)
        assert sorted(map(len,lk.values()))==([1,1]+[2]*(len(lk)-2)if v in boundary else[2]*len(lk))
    xyz=Counter(sum(v in new for v in t)for t in triangles)
    assert (len(boundary)+len(new),len(e),len(triangles))==(12,28,17)
    assert xyz=={1:5,2:6,3:6}
    return dict(boundary=boundary,triangles=triangles,vertices=12,edges=28,triangles_count=17,
      x=5,y=6,z=6,conditional_end_change=-5,
      scope='Exact abstract simple disk only. Its embedding on the selected cavity side and all geometric incidences require separate checks; the retained outside30--44 chord cancels only under that embedding.')

@lru_cache(None)
def corner3(p0,p1,p2,m01,m02,m12):
    g=[[2]*5 for _ in range(5)]
    for i in range(5):g[i][i]=1
    g[0][1]=g[1][0]=3
    for j,pi in enumerate([p0,p1,p2],2):
        for i,m in enumerate(F.G.profile_labels(pi)):g[i][j]=g[j][i]=m
    for i,j,m in [(2,3,m01),(2,4,m02),(3,4,m12)]:g[i][j]=g[j][i]=m
    return N.A.components(tuple(v for row in g for v in row))

def run():
    topology=disk();current,labels,ns,h,profiles,primary=H.setup();actual=[i for i in current['actual_sides']if i not in S.DELETE]
    prior=ROOT/'results/fresh_two_connector_prefix_witnesses_20261017.json';pr=json.loads(prior.read_text())
    a=np.array(pr['records'][0]['normals'][0]);assert not T.gates(T.native_values(a,0,labels,ns,profiles),[7,8])
    assert abs(a@M@ns[53]+1)<1e-7
    old=ROOT/'results/fresh_two_connector_coupled_disk_20261017.json';bank=json.loads(old.read_text())
    counts=Counter();fail=Counter();unresolved=[];seeds=[];complete=[];bc_range=[float('inf'),float('-inf')]
    def native(v,pi,assigned=()):
        if v@M@h+profiles[pi]['center_offset']>=-1e-8:return [('center',0.)]
        return T.gates(T.native_values(v,pi,labels,ns,profiles),assigned)
    def nonincident(v,pi,allow):
        q=profiles[pi]['square']
        return all(v@M@ns[i]<-np.sqrt(q*(ns[i]@M@ns[i]))-1e-7 for i in actual if i not in allow)
    def finite(v,w,pi,pj):return F.label(H.A.native(v,pi,w,pj,profiles))
    def dotted(v,w,pi,pj):return v@M@w<-np.sqrt(profiles[pi]['square']*profiles[pj]['square'])-1e-7
    rs={}
    for mr in bank['R44_labels']:
        rr,_=N.norm_roots([ns[8],ns[44],ns[63]],[0.,-B.theta(mr),0.],2.)
        rs[mr]=[v for sg,v in rr if not native(v,0,[8,44,63])][0]
    for mac in [4,5,6]:
        cr,ci=N.norm_roots([ns[7],ns[30],a],[0.,0.,-B.theta(mac)],2.);counts['C_equations']+=1
        if ci['status'].startswith('unresolved'):unresolved.append(['C',mac,ci]);continue
        for csg,c in cr:
            if native(c,0,[7,30])or not nonincident(c,0,[7,30]):continue
            counts['native_C']+=1
            for xp,m8,maxx,mcx in product(range(7),ALPHABET,ALPHABET,ALPHABET):
                if T.corner(xp,m8,maxx,0,2)is None or T.corner(xp,maxx,mcx,0,mac)is None:continue
                xr,xi=N.norm_roots([ns[8],a,c],[-B.theta(m8),-B.theta(maxx),-B.theta(mcx)],profiles[xp]['square']);counts['X_equations']+=1
                if xi['status'].startswith('unresolved'):unresolved.append(['X',mac,xp,m8,maxx,mcx,xi]);continue
                for xsg,x in xr:
                    bad=native(x,xp,[8])
                    if bad:fail[('X',str(bad[0][0]))]+=1;continue
                    if not nonincident(x,xp,[8]):counts['X_extra_intrinsic_contact']+=1;continue
                    counts['native_ACX']+=1
                    for mr,r in rs.items():
                        if not dotted(x,r,xp,0)or not dotted(c,r,0,0)or not dotted(a,r,0,0):continue
                        for yp,my8,myr,myx in product(range(7),ALPHABET,ALPHABET,ALPHABET):
                            if T.corner(yp,my8,myr,0,2)is None or T.corner(yp,my8,myx,xp,m8)is None:continue
                            contacts=[-B.theta(my8),-B.theta(myr),-B.theta(myx)]
                            # The last native correction is obtained directly from the saved seven-profile Gram.
                            contacts[2]-=H.A.native(np.zeros(4),yp,np.zeros(4),xp,profiles)
                            yr,yi=N.norm_roots([ns[8],r,x],contacts,profiles[yp]['square']);counts['Y_equations']+=1
                            if yi['status'].startswith('unresolved'):unresolved.append(['Y',mac,xp,m8,maxx,mcx,mr,yp,my8,myr,myx,yi]);continue
                            for ysg,y in yr:
                                bad=native(y,yp,[8])
                                if bad:fail[('Y',str(bad[0][0]))]+=1;continue
                                if not nonincident(y,yp,[8])or not dotted(y,a,yp,0)or not dotted(y,c,yp,0):continue
                                counts['native_ACXYR']+=1
                                seed=dict(AC=mac,X_profile=xp,X_labels=[m8,maxx,mcx],R44=mr,Y_profile=yp,Y_labels=[my8,myr,myx],normals=[v.tolist()for v in[a,c,x,y,r]])
                                ndseed=dict(ns);ndseed.update({67+i:v for i,v in enumerate([a,r,c,x,y])})
                                seed['boundary']=S.C.detailed_cone(ndseed,actual+list(range(67,72)))
                                seed['PD_incidence']=H.incidence_score([ndseed[i]for i in actual+list(range(67,72))],actual+list(range(67,72)),[float(ndseed[i]@M@ndseed[i])for i in actual+list(range(67,72))])
                                seed['native_pair_checks']=[dict(pair=[i,j],product=float(H.A.native(v,pi,w,pj,profiles)),relation=B.relation(float(H.A.native(v,pi,w,pj,profiles))))for i,(v,pi)in enumerate(zip([a,c,x,y,r],[0,0,xp,yp,0]))for j,(w,pj)in enumerate(zip([a,c,x,y,r],[0,0,xp,yp,0]))if i<j]
                                assert all(q['relation']is not None and q['relation']!='parallel_or_unresolved'for q in seed['native_pair_checks'])
                                assert seed['boundary']['warnings']==[] and seed['boundary']['open_ends']==seed['PD_incidence']['open_ends']
                                seeds.append(seed)
                                for ib,br in enumerate(bank['B_bank']):
                                    if br['R44']!=mr:continue
                                    bp=br['B_profile'];b=np.array(br['B_projection']);m30,m63,mrb=br['B_labels'];mcb=finite(c,b,0,bp);myb=finite(y,b,yp,bp)
                                    cb=float(c@M@b);counts['matching_B_tests']+=1;bc_range=[min(bc_range[0],cb),max(bc_range[1],cb)]
                                    if B.relation(cb)=='dotted':counts['matching_B_C_dotted']+=1
                                    if mcb is None or myb is None:continue
                                    if not dotted(a,b,0,bp)or not dotted(x,b,xp,bp):continue
                                    if T.corner(bp,m30,mcb,0,2)is None or T.corner(yp,myr,myb,bp,mrb)is None:continue
                                    counts['native_six_cap_prefix']+=1
                                    for zp,mzx,mzy,mzc in product(range(7),ALPHABET,ALPHABET,ALPHABET):
                                        if corner3(xp,yp,zp,myx,mzx,mzy)is None or corner3(xp,0,zp,mcx,mzx,mzc)is None:continue
                                        cc=[-B.theta(mzx)-H.A.native(np.zeros(4),zp,np.zeros(4),xp,profiles),-B.theta(mzy)-H.A.native(np.zeros(4),zp,np.zeros(4),yp,profiles),-B.theta(mzc)]
                                        zr,zi=N.norm_roots([x,y,c],cc,profiles[zp]['square']);counts['Z_equations']+=1
                                        if zi['status'].startswith('unresolved'):unresolved.append(['Z',len(seeds)-1,ib,zp,mzx,mzy,mzc,zi]);continue
                                        for zsg,v in zr:
                                            if native(v,zp)or not nonincident(v,zp,[]):continue
                                            mzb=finite(v,b,zp,bp)
                                            if mzb is None or not dotted(v,a,zp,0)or not dotted(v,r,zp,0):continue
                                            if corner3(0,bp,zp,mcb,mzc,mzb)is None or corner3(yp,bp,zp,myb,mzy,mzb)is None:continue
                                            counts['coherent_complete_disks']+=1
                                            nd=dict(ns);new=[a,r,c,x,y,b,v];nd.update({67+i:n for i,n in enumerate(new)})
                                            bd=S.C.detailed_cone(nd,actual+list(range(67,74)))
                                            inc=H.incidence_score([nd[i]for i in actual+list(range(67,74))],actual+list(range(67,74)),[float(nd[i]@M@nd[i])for i in actual+list(range(67,74))])
                                            complete.append(dict(seed=len(seeds)-1,B_index=ib,Z_profile=zp,Z_labels=[mzx,mzy,mzc,mzb],Z=v.tolist(),boundary=bd,PD_incidence=inc))
        print('INTERIOR_CONNECTOR',mac,dict(counts),flush=True)
    paths=[Path(__file__),Path(F.__file__),Path(T.__file__),Path(N.__file__),Path(H.__file__),prior,old,H.PARENT,H.BOUNDARY]
    return dict(status='BOUNDED_NUMERICAL_SEVEN_CAP_INTERIOR_CONNECTOR_DISK',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},topology=topology,counts=dict(counts),matching_B_C_product_range=bc_range,blockers=[[list(k),v]for k,v in sorted(fail.items())],unresolved=unresolved,seeds=seeds,complete=complete,deleted_walls=S.DELETE,exact_old_deletion_cost=pr['exact_old_deletion_cost'],domain=dict(A='Frozen native A with A53=pi/3, right S6+7+8',AC=[4,5,6],R44=bank['R44_labels'],B='Frozen422 native branches',X_Y_Z_profiles=list(range(7)),assigned_labels=[2,3,4,5]),scope='Bounded numerical equations only. Unresolved branches are recorded, never accepted. All native contacts/common-center and stipulated intrinsic nonincidences are checked. Nonzero-profile finite-corner classifier is conservative and full finite support gates are inherited. High labels and continuously varied A remain open; no exact geometry or finite-volume promotion.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=json.loads(json.dumps(run()));path=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],len(out['complete']),flush=True)
