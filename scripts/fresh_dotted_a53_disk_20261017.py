#!/usr/bin/env python3
"""Bounded one-parameter rebuilding of the old30 flipped disk.

A stays right to S6+7+8, while -A.53 varies and may be dotted. Rebuild
C,X,Y and the mandatory Z instead of reusing the twenty frozen vectors.
Numerical discovery only; no continuous-coverage or exact-existence claim.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
import fresh_old30_connector_flip_20261017 as O
P,F,T,N,S,H,B,M,ROOT=O.P,O.F,O.T,O.N,O.S,O.H,O.B,O.M,O.ROOT
TAG='fresh_dotted_a53_disk_20261017'

def run():
    current,labels,ns,h,pro,primary=H.setup();actual=[i for i in current['actual_sides']if i not in S.DELETE]
    pp=ROOT/'results'/f'{P.TAG}.json';gp=ROOT/'results/fresh_old30_connector_native_gate_20261017.json'
    bp=ROOT/'results/fresh_two_connector_coupled_disk_20261017.json'
    parent=json.loads(pp.read_text());gate=json.loads(gp.read_text());bank=json.loads(bp.read_text())
    counts=Counter();blocks=Counter();unresolved=[];accepted=[];rows=[];R={}
    def prod(v,pi,w,pj):return float(H.A.native(v,pi,w,pj,pro))
    def offset(pi,pj):return prod(np.zeros(4),pi,np.zeros(4),pj)
    def native(v,pi,assigned=()):
        bad=T.gates(T.native_values(v,pi,labels,ns,pro),assigned)
        if bad:return [('old'+str(i),z)for i,z in bad]
        if v@M@h+pro[pi]['center_offset']>=-1e-8:return [('center',0.)]
        return []
    def dot(v,pi,w,pj):return float(v@M@w)<-np.sqrt(pro[pi]['square']*pro[pj]['square'])-1e-7
    def nonincident(v,pi,allow):
        return all(float(v@M@ns[i])<-np.sqrt(pro[pi]['square']*float(ns[i]@M@ns[i]))-1e-7 for i in actual if i not in allow)
    def solve(name,frame,contact,pi,assigned=(),context=None):
        counts[name+'_equations']+=1
        roots,info=N.norm_roots(frame,contact,pro[pi]['square'])
        if info['status'].startswith('unresolved'):
            unresolved.append(dict(stage=name,context=context,info=info));return []
        out=[]
        for sg,v in roots:
            counts[name+'_branches']+=1
            err=max(abs(float(v@M@v)-pro[pi]['square']),max(abs(float(v@M@w)-c)for w,c in zip(frame,contact)))
            if err>1e-5:
                unresolved.append(dict(stage=name,context=context,branch=sg,maximum_equation_error=err));continue
            bad=native(v,pi,assigned)
            if bad:blocks[name,bad[0][0]]+=1;continue
            counts[name+'_native']+=1;out.append((sg,v))
        return out
    for mr in sorted({r['R44']for r in parent['seeds']}):
        roots=solve('R',[ns[8],ns[44],ns[63]],[0.,-B.theta(mr),0.],0,[8,44,63],mr)
        assert len(roots)==1;R[mr]=roots[0][1]
    parameters=[dict(name='control_pi3',a=1.)]+[dict(name='pi'+str(m),a=B.theta(m))for m in[4,5,6,10]]
    parameters += [dict(name='dotted_'+str(a),a=a)for a in[2.001,2.01,2.05,2.1,2.25,2.5,3.,4.,6.,8.,12.,16.,24.,32.]]
    for par in parameters:
        a53=par['a'];start=Counter(counts);before=len(accepted)
        aa=solve('A',[ns[7],ns[8],ns[53]],[0.,0.,-a53],0,[7,8],par['name'])
        for asg,a in aa:
            if not nonincident(a,0,[7,8]):counts['A_extra_intrinsic_contact']+=1;continue
            for seed,rec in enumerate(parent['seeds']):
                xp,yp=rec['X_profile'],rec['Y_profile'];m8,maxx,mcx=rec['X_labels'];my8,myr,myx=rec['Y_labels'];r=R[rec['R44']]
                if not dot(a,0,r,0):counts['A_R_nondotted']+=1;continue
                ctx=[par['name'],seed]
                cc=solve('C',[ns[7],ns[30],a],[0.,0.,-B.theta(rec['AC'])],0,[7,30],ctx)
                for csg,c in cc:
                    if not nonincident(c,0,[7,30])or not dot(c,0,r,0):continue
                    xx=solve('X',[ns[8],a,c],[-B.theta(m8),-B.theta(maxx),-B.theta(mcx)],xp,[8],ctx)
                    for xsg,x in xx:
                        if not nonincident(x,xp,[8])or not dot(x,xp,r,0):continue
                        yy=solve('Y',[ns[8],r,x],[-B.theta(my8),-B.theta(myr),-B.theta(myx)-offset(yp,xp)],yp,[8],ctx)
                        for ysg,y in yy:
                            if not nonincident(y,yp,[8])or not dot(y,yp,a,0)or not dot(y,yp,c,0):continue
                            counts['coherent_five_prefixes']+=1
                            for case in gate['records'][seed]['complete_Z_cases']:
                                zp=case['profile'];mx,my,mr=case['labels']
                                zz=solve('Z',[x,y,r],[-B.theta(mx)-offset(zp,xp),-B.theta(my)-offset(zp,yp),-B.theta(mr)],zp,(),ctx+[zp,mx,my,mr])
                                for zsg,z in zz:
                                    if not nonincident(z,zp,[])or not dot(z,zp,a,0):counts['Z_extra_intrinsic_contact']+=1;continue
                                    mzc=F.label(prod(z,zp,c,0))
                                    if mzc is None:counts['Z_C_nonfinite']+=1;continue
                                    if P.corner3(xp,0,zp,mcx,mx,mzc)is None:counts['Z_C_nonfinite_corner']+=1;continue
                                    counts['coherent_six_prefixes']+=1
                                    new=[a,r,c,x,y,z];pis=[0,0,0,xp,yp,zp]
                                    pairs=[dict(pair=[i,j],product=prod(new[i],pis[i],new[j],pis[j]))for i,j in combinations(range(6),2)]
                                    if any(B.relation(p['product'])in[None,'parallel_or_unresolved']for p in pairs):counts['all_pair_failure']+=1;continue
                                    accepted.append(dict(parameter=par,seed=seed,branches=[asg,csg,xsg,ysg,zsg],Z_profile=zp,Z_labels=[mx,my,mr,mzc],normals=[v.tolist()for v in new],profiles=pis,native_pairs=pairs))
        row=dict(parameter=par,counts=dict(counts-start),accepted=len(accepted)-before);rows.append(row)
        print('DOTTED_A53',par['name'],row['counts'],'accepted',row['accepted'],flush=True)
    paths=[Path(__file__),Path(O.__file__),Path(P.__file__),Path(T.__file__),Path(N.__file__),Path(H.__file__),pp,gp,bp,H.PARENT,H.BOUNDARY]
    return dict(status='BOUNDED_A53_PARAMETER_REBUILDING_OF_FLIPPED_DISK',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
      parameters=parameters,counts=dict(counts),rows=rows,blockers=[[list(k),v]for k,v in sorted(blocks.items())],accepted=accepted,unresolved=unresolved,
      deleted_walls=S.DELETE,exact_old_deletion_cost=parent['exact_old_deletion_cost'],topology=O.topology(True),
      domain='A rightS6+7+8, norm2, nativeA53=-a. One control value, four other finite angles and14 explicit dotted parameters. C,X,Y rebuilt from the20 saved contact recipes; R choices retained. Z uses295low+2high-placeholder recipes, with6 sampled in the latter, not all-label coverage.',
      scope='Bounded numerical discovery, not continuous interval coverage. Every retained native gate, fixed common center, equation residual and required intrinsic nonincidence is checked before saving. No new exact geometry, complete patch or finite-volume promotion; oldcontrol63/199/10 remains unchanged.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=json.loads(json.dumps(run()));p=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],'accepted',len(out['accepted']),'unresolved',len(out['unresolved']),flush=True)
