#!/usr/bin/env python3
"""New finite old10 anchor for C, freeing A.53 to be dotted.

Bounded numerical discovery; no all-label or exact chamber assertion.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
import fresh_old30_connector_flip_20261017 as O
P,F,T,N,S,H,B,M,ROOT=O.P,O.F,O.T,O.N,O.S,O.H,O.B,O.M,O.ROOT
TAG='fresh_c10_endpoint_disk_20261017'
C10=[3,4,5,6,7,8,9,10,12,15,20,30,50]
AC=[2,3,4,5,6,7,8,10,15]
ALPHABET=range(2,7)

def run():
    current,labels,ns,h,pro,primary=H.setup();actual=[i for i in current['actual_sides']if i not in S.DELETE]
    pp=ROOT/'results'/f'{P.TAG}.json';bp=ROOT/'results/fresh_two_connector_coupled_disk_20261017.json'
    parent=json.loads(pp.read_text());bank=json.loads(bp.read_text())
    counts=Counter();blocks=Counter();unresolved=[];ACbank=[];prefixes=[];zbank=[];complete=[];rs={}
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
        counts[name+'_equations']+=1;roots,info=N.norm_roots(frame,rhs,pro[pi]['square'])
        if info['status'].startswith('unresolved'):
            unresolved.append(dict(stage=name,context=context,info=info));return []
        out=[]
        for sg,v in roots:
            counts[name+'_branches']+=1
            err=max(abs(float(v@M@v)-pro[pi]['square']),max(abs(float(v@M@w)-q)for w,q in zip(frame,rhs)))
            if err>1e-5:
                unresolved.append(dict(stage=name,context=context,branch=sg,error=err));continue
            bad=native(v,pi,assigned)
            if bad:blocks[name,bad[0][0]]+=1;continue
            counts[name+'_native']+=1;out.append((sg,v))
        return out
    def pairchecks(vs,ps):
        out=[dict(pair=[i,j],product=nativeprod(vs[i],ps[i],vs[j],ps[j]))for i,j in combinations(range(len(vs)),2)]
        return out if all(B.relation(r['product'])not in[None,'parallel_or_unresolved']for r in out)else None
    for mr in bank['R44_labels']:
        roots=solve('R',[ns[8],ns[44],ns[63]],[0.,-B.theta(mr),0.],0,[8,44,63],mr)
        assert len(roots)==1;rs[mr]=roots[0][1]
    for c10 in C10:
        for csg,c in solve('C',[ns[7],ns[30],ns[10]],[0.,0.,-B.theta(c10)],0,[7,30,10],c10):
            if not nonincident(c,0,[7,30]):continue
            for mac in AC:
                for asg,a in solve('A',[ns[7],ns[8],c],[0.,0.,-B.theta(mac)],0,[7,8],[c10,mac]):
                    if not nonincident(a,0,[7,8]):continue
                    ia=len(ACbank);ACbank.append(dict(C10=c10,AC=mac,branches=[csg,asg],normals=[a.tolist(),c.tolist()],A53=float(a@M@ns[53])))
                    for xp,m8,mxa,mxc in product(range(7),ALPHABET,ALPHABET,ALPHABET):
                        if T.corner(xp,m8,mxa,0,2)is None or T.corner(xp,mxa,mxc,0,mac)is None:continue
                        for xsg,x in solve('X',[ns[8],a,c],[-B.theta(m8),-B.theta(mxa),-B.theta(mxc)],xp,[8],[ia,xp,m8,mxa,mxc]):
                            if not nonincident(x,xp,[8]):continue
                            counts['native_ACX']+=1
                            for mr,r in rs.items():
                                if not all(dotted(v,pi,r,0)for v,pi in[(a,0),(c,0),(x,xp)]):continue
                                for yp,my8,myr,myx in product(range(7),ALPHABET,ALPHABET,ALPHABET):
                                    if T.corner(yp,my8,myr,0,2)is None or T.corner(yp,my8,myx,xp,m8)is None:continue
                                    rhs=[-B.theta(my8),-B.theta(myr),-B.theta(myx)-offset(yp,xp)]
                                    for ysg,y in solve('Y',[ns[8],r,x],rhs,yp,[8],[ia,xp,m8,mxa,mxc,mr,yp,my8,myr,myx]):
                                        if not nonincident(y,yp,[8])or not dotted(y,yp,a,0)or not dotted(y,yp,c,0):continue
                                        pairs=pairchecks([a,r,c,x,y],[0,0,0,xp,yp])
                                        if pairs is None:continue
                                        counts['coherent_five_prefixes']+=1;ip=len(prefixes)
                                        prefix=dict(AC_index=ia,R44=mr,X_profile=xp,X_labels=[m8,mxa,mxc],Y_profile=yp,Y_labels=[my8,myr,myx],branches=[xsg,ysg],normals=[v.tolist()for v in[a,r,c,x,y]],native_pairs=pairs)
                                        prefixes.append(prefix)
                                        for zp,mzx,mzy,mzr in product(range(7),ALPHABET,ALPHABET,ALPHABET):
                                            if P.corner3(xp,yp,zp,myx,mzx,mzy)is None or P.corner3(yp,0,zp,myr,mzy,mzr)is None:continue
                                            rhs=[-B.theta(mzx)-offset(zp,xp),-B.theta(mzy)-offset(zp,yp),-B.theta(mzr)]
                                            for zsg,z in solve('Z',[x,y,r],rhs,zp,(),[ip,zp,mzx,mzy,mzr]):
                                                if not nonincident(z,zp,[])or not dotted(z,zp,a,0):continue
                                                mzc=F.label(nativeprod(z,zp,c,0))
                                                if mzc is None:counts['Z_C_nonfinite']+=1;continue
                                                if P.corner3(xp,0,zp,mxc,mzx,mzc)is None:counts['Z_C_nonfinite_corner']+=1;continue
                                                npairs=pairchecks([a,r,c,x,y,z],[0,0,0,xp,yp,zp])
                                                if npairs is None:continue
                                                counts['coherent_six_prefixes']+=1;nd=dict(ns);nd.update({67+i:v for i,v in enumerate([a,r,c,x,y,z])});sel=actual+list(range(67,73))
                                                bd=S.C.detailed_cone(nd,sel);ic=H.incidence_score([nd[i]for i in sel],sel,[float(nd[i]@M@nd[i])for i in sel])
                                                if bd['warnings']or ic['status']!='numerically_separated_acute_incidence'or bd['open_ends']!=ic['open_ends']:
                                                    unresolved.append(dict(stage='six_prefix_boundary',prefix=ip,bd=bd,ic=ic));continue
                                                zbank.append(dict(prefix=ip,Z_profile=zp,Z_labels=[mzx,mzy,mzr,mzc],branch=zsg,Z=z.tolist(),native_pairs=npairs,boundary=bd,PD_incidence=ic))
        print('C10_ENDPOINT',c10,dict(counts),'Z',len(zbank),flush=True)
    histogram=Counter()
    for i,prefix in enumerate(prefixes):
        vs=list(map(np.array,prefix['normals']));nd=dict(ns);nd.update({67+j:v for j,v in enumerate(vs)});sel=actual+list(range(67,72))
        bd=S.C.detailed_cone(nd,sel);ic=H.incidence_score([nd[j]for j in sel],sel,[float(nd[j]@M@nd[j])for j in sel])
        if bd['warnings']or ic['status']!='numerically_separated_acute_incidence'or bd['open_ends']!=ic['open_ends']:
            unresolved.append(dict(stage='five_prefix_boundary',prefix=i,bd=bd,ic=ic));continue
        prefix['boundary']=bd;prefix['PD_incidence']=ic
        histogram[str((bd['counts'],bd['open_ends']))]+=1
        if(i+1)%25==0:print('C10_PREFIX_BOUNDARY',i+1,'of',len(prefixes),dict(histogram),flush=True)
    paths=[Path(__file__),Path(O.__file__),Path(P.__file__),Path(T.__file__),Path(N.__file__),Path(H.__file__),pp,bp,H.PARENT,H.BOUNDARY]
    return dict(status='BOUNDED_NATIVE_C10_ENDPOINT_RECONSTRUCTION',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},counts=dict(counts),blockers=[[list(k),v]for k,v in sorted(blocks.items())],AC_bank=ACbank,prefixes=prefixes,Z_bank=zbank,unresolved=unresolved,
      prefix_boundary_histogram=dict(histogram),domain=dict(C10=C10,AC=AC,R44=bank['R44_labels'],X_Y_Z_profiles=list(range(7)),assigned_X_Y_Z_labels=list(ALPHABET)),
      deleted_walls=S.DELETE,exact_old_deletion_cost=parent['exact_old_deletion_cost'],topology=O.topology(True),
      scope='New old10 endpoint contact moves C and releases A53 while retaining every oldwall outside the frozen deletion set. Finite event lists are bounded, including label6 samples; no all-label or continuous coverage. Every saved prefix passes native/common-center/nonincidence gates; any Z prefix receives full numerical cone and independentPD incidence. B/W and complete patch/ambient certification remain required. No new chamber or finite-volume promotion.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=json.loads(json.dumps(run()));p=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],'Z',len(out['Z_bank']),'unresolved',len(out['unresolved']),flush=True)
