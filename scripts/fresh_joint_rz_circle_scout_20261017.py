#!/usr/bin/env python3
"""Joint low-RZ discovery on the mathematically bounded R circle.

The R44 bound is exact; geometry acceptance here remains numerical.
The high-Z8 tail is explicitly outside this bounded scout.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
import fresh_old8_triangle_double_flip_20261017 as D
O,P,F,T,N,S,H,B,M,ROOT=D.O,D.P,D.F,D.T,D.N,D.S,D.H,D.B,D.M,D.ROOT
TAG='fresh_joint_rz_circle_scout_20261017'

def run():
    pp=ROOT/'results'/f'{D.TAG}.json';rp=ROOT/'results/fresh_r_circle_admissibility_20261017.json'
    parent=json.loads(pp.read_text());circle=json.loads(rp.read_text());cur,labels,ns,h,pro,_=H.setup()
    assert circle['exact_circle']['necessary_labels']==list(range(2,16))
    a,_,c,x,_,_=map(np.array,parent['prefixes'][0]['normals']);actual=[i for i in cur['actual_sides']if i not in S.DELETE]
    counts=Counter();blockers=Counter();unresolved=[];prefixes=[];complete=[]
    def nativeprod(v,pi,w,pj):return float(H.A.native(v,pi,w,pj,pro))
    def offset(pi,pj):return nativeprod(np.zeros(4),pi,np.zeros(4),pj)
    def dotted(v,pi,w,pj):return float(v@M@w)<-np.sqrt(pro[pi]['square']*pro[pj]['square'])-1e-7
    def nonincident(v,pi,allow):return all(float(v@M@ns[i])<-np.sqrt(pro[pi]['square']*float(ns[i]@M@ns[i]))-1e-7 for i in actual if i not in allow)
    def finite(v,pi,w,pj):return F.label(nativeprod(v,pi,w,pj))
    def pairs(vs,ps):
        out=[dict(pair=[i,j],product=nativeprod(vs[i],ps[i],vs[j],ps[j]))for i,j in combinations(range(len(vs)),2)]
        return out if all(B.relation(q['product'])not in[None,'parallel_or_unresolved']for q in out)else None
    def solve(stage,frame,rhs,pi,assigned=(),context=None):
        counts[stage+'_equations']+=1;roots,info=N.norm_roots(frame,rhs,pro[pi]['square'])
        if info['status'].startswith('unresolved'):
            unresolved.append(dict(stage=stage,context=context,info=info));return []
        out=[]
        for sg,v in roots:
            counts[stage+'_branches']+=1
            error=max(abs(float(v@M@v)-pro[pi]['square']),max(abs(float(v@M@w)-q)for w,q in zip(frame,rhs)))
            if error>1e-5:unresolved.append(dict(stage=stage,context=context,branch=sg,error=error));continue
            bad=T.gates(T.native_values(v,pi,labels,ns,pro),assigned)
            if bad:blockers[stage,bad[0][0]]+=1;continue
            if v@M@h+pro[pi]['center_offset']>=-1e-8:counts[stage+'_center_reject']+=1;continue
            counts[stage+'_native']+=1;out.append((sg,v))
        return out
    def boundary(vs):
        nd=dict(ns);nd.update({67+i:v for i,v in enumerate(vs)});sel=actual+list(range(67,67+len(vs)))
        bd=S.C.detailed_cone(nd,sel);ic=H.incidence_score([nd[i]for i in sel],sel,[float(nd[i]@M@nd[i])for i in sel])
        if bd['warnings']or ic['status']!='numerically_separated_acute_incidence'or bd['open_ends']!=ic['open_ends']or bd['counts'].get('finite',0)!=ic['finite']:
            unresolved.append(dict(stage='boundary',bd=bd,ic=ic));return None
        if {tuple(sorted(v['active']))for v in bd['rays']if v['kind']=='finite'}!={tuple(sorted(v))for v in ic['finite_supports']}:
            unresolved.append(dict(stage='finite_support_equality',bd=bd,ic=ic));return None
        return bd,ic
    for rr in circle['numerical_R_candidates']:
        mr=rr['label']
        if mr not in circle['exact_circle']['necessary_labels']or not rr['native_and_center_pass']:continue
        r=np.array(rr['normal'])
        if not nonincident(r,0,[8,44,63])or not all(dotted(r,0,v,0)for v in[a,c,x]):counts['R_prefix_nonincidence_reject']+=1;continue
        counts['coherent_ARCX']+=1
        for zp,m8z,mrz,mxz in product(range(7),range(2,6),range(2,6),range(2,6)):
            if P.corner3(0,0,zp,2,m8z,mxz)is None:continue
            if P.corner3(0,0,zp,2,m8z,2)is None or P.corner3(0,0,zp,2,mrz,2)is None:continue
            rem=pro[zp]['square']-(B.theta(m8z)**2+B.theta(mrz)**2)/2
            if abs(rem)<1e-9:counts['zero_polar_limit_skipped']+=1;continue
            if rem>0:counts['positive_polar_frame_skipped']+=1;continue
            for zsg,z in solve('Z',[ns[8],r,x],[-B.theta(m8z),-B.theta(mrz),-B.theta(mxz)],zp,[8],[mr,zp,m8z,mrz,mxz]):
                if not nonincident(z,zp,[8])or not dotted(z,zp,a,0)or not dotted(z,zp,c,0):counts['Z_nonincidence_reject']+=1;continue
                for ysg,y in solve('Y',[ns[8],r,z],[0.,0.,0.],0,[8],[mr,zp,m8z,mrz,mxz,zsg]):
                    if not nonincident(y,0,[8])or not all(dotted(y,0,v,0)for v in[a,c,x]):continue
                    vs=[a,r,c,x,y,z];ps=[0,0,0,0,0,zp];pchecks=pairs(vs,ps)
                    if pchecks is None:continue
                    bdic=boundary(vs)
                    if bdic is None:continue
                    bd,ic=bdic;yrays=[v for v in bd['rays']if 71 in v['active']]
                    if len(yrays)!=3 or any(v['kind']!='finite'for v in yrays):
                        unresolved.append(dict(stage='Y_actual_triangle',rays=yrays));continue
                    ip=len(prefixes);prefixes.append(dict(R44=mr,Z_profile=zp,Z_labels=[m8z,mrz,mxz],branches=[zsg,ysg],normals=[v.tolist()for v in vs],profiles=ps,native_pairs=pchecks,boundary=bd,PD_incidence=ic))
                    counts['coherent_six_prefixes']+=1
                    # The last two walls are searched only if the genuinely
                    # changed low-RZ six-prefix survives all old native gates.
                    wbank=[];bbank=[]
                    for wp,mw30,mwc,mwx in product(range(7),range(2,7),range(2,7),range(2,6)):
                        if P.corner3(0,0,wp,2,mw30,mwc)is None or P.corner3(0,0,wp,2,mwc,mwx)is None:continue
                        for wsg,w in solve('W',[ns[30],c,x],[-B.theta(mw30),-B.theta(mwc),-B.theta(mwx)],wp,[30],[ip,wp,mw30,mwc,mwx]):
                            if nonincident(w,wp,[30])and all(dotted(w,wp,v,0)for v in[a,r,y]):wbank.append((wp,mw30,mwc,mwx,wsg,w))
                    for bp,mb30,mb63,mbr in product(range(7),range(2,7),range(2,7),range(2,6)):
                        if P.corner3(0,0,bp,2,mb30,mb63)is None or P.corner3(0,0,bp,2,mb63,mbr)is None:continue
                        for bsg,b in solve('B',[ns[30],ns[63],r],[-B.theta(mb30),-B.theta(mb63),-B.theta(mbr)],bp,[30,63],[ip,bp,mb30,mb63,mbr]):
                            if nonincident(b,bp,[30,63])and all(dotted(b,bp,v,0)for v in[a,c,y]):bbank.append((bp,mb30,mb63,mbr,bsg,b))
                    for wr,br in product(wbank,bbank):
                        wp,mw30,mwc,mwx,wsg,w=wr;bp,mb30,mb63,mbr,bsg,b=br
                        mbz=finite(b,bp,z,zp);mbw=finite(b,bp,w,wp)
                        if mbz is None or mbw is None or P.corner3(0,zp,bp,mrz,mbr,mbz)is None or P.corner3(0,bp,wp,mb30,mw30,mbw)is None:continue
                        mwz=finite(w,wp,z,zp);mbx=finite(b,bp,x,0);kind=None
                        if dotted(b,bp,x,0)and mwz is not None and P.corner3(0,zp,wp,mxz,mwx,mwz)is not None and P.corner3(bp,zp,wp,mbz,mbw,mwz)is not None:kind='original_WZ'
                        if dotted(w,wp,z,zp)and mbx is not None and P.corner3(0,zp,bp,mxz,mbx,mbz)is not None and P.corner3(0,bp,wp,mbx,mwx,mbw)is not None:kind='flipped_XB'
                        if kind is None:continue
                        allvs=[a,r,c,x,y,b,z,w];allps=[0,0,0,0,0,bp,zp,wp];checks=pairs(allvs,allps)
                        if checks is None:continue
                        bdic=boundary(allvs)
                        if bdic is None:continue
                        bd,ic=bdic;complete.append(dict(prefix=ip,topology=kind,B_profile=bp,B_labels=[mb30,mb63,mbr,mbz,mbw,mbx],B_branch=bsg,B=b.tolist(),W_profile=wp,W_labels=[mw30,mwc,mwx,mwz],W_branch=wsg,W=w.tolist(),native_pairs=checks,boundary=bd,PD_incidence=ic))
        print('JOINT_RZ',mr,dict(counts),'prefixes',len(prefixes),'complete',len(complete),flush=True)
    paths=[Path(__file__),Path(D.__file__),pp,rp,ROOT/'scripts/fresh_r_circle_admissibility_20261017.py',H.PARENT,H.BOUNDARY]
    return dict(status='BOUNDED_JOINT_LOW_RZ_CIRCLE_DISCOVERY',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},counts=dict(counts),blockers=[[list(k),v]for k,v in sorted(blockers.items())],prefixes=prefixes,complete=complete,unresolved=unresolved,
      domain=dict(fixed_C10=9,fixed_AC=3,fixed_X_labels=[2,3,2],R44='complete necessary2..15 from exact native circle gate',Z_profiles=list(range(7)),Z8_ZR_ZX=list(range(2,6)),Y='rightS6+8+R+Z, both polar signs',remaining_B_W='bounded assigned labels2..6 if a six-prefix survives'),
      high_Z8_tail='OPEN: withZX=2 and zero A2 profile, Z8 may be any finite label>=6. This scout does not exclude that continuous tail. A highZX is incompatible with a finite8/X/Z corner unlessZ8=2, which cannot produce the required Lorentz8/R/Z frame while its Y-corners remain finite.',
      exact_old_deletion_cost=parent['exact_old_deletion_cost'],scope='The R44 circle bound is exact; this joint low-RZ/low-Z8 discovery and all geometry acceptance are numerical. It does not exclude the stated high-Z8 tail or all continuously varied A/C/X prefixes. Any completed chamber requires separate exact native, facet, full-boundary and volume certification. Frozen63/199/10 infinite-volume control remains unchanged.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=json.loads(json.dumps(run()));p=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],'complete',len(out['complete']),'unresolved',len(out['unresolved']),flush=True)
