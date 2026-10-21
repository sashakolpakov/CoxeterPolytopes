#!/usr/bin/env python3
"""All finite profiles exclude one wall terminating both12/52 and30/52."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_nonarithmetic_replace47_close12_20261017 as C
P,U,I,D=C.P,C.U,C.I,C.D
ROOT=C.ROOT;TAG='fresh_nonarithmetic_replace47_edges1252_3052_gate_20261017'


def run():
    z=C.geometry();ns=z['normals'];pair=z['pair']
    frame=[ns[i]for i in [12,30,52]]
    gram=[[I(2),C.literal(z,12,30),-D.theta(6)],
          [C.literal(z,12,30),I(2),-D.theta(4)],[-D.theta(6),-D.theta(4),I(2)]]
    inverse,pivots=P.C53.inverse(gram)
    e7=[U.rnd(a+(b+2*c)/3)for a,b,c in zip(ns[7],ns[5],ns[6])]
    kernel_products=[I(),-D.theta(5),C.literal(z,7,52)]
    kernel_coeff=U.imv(inverse,kernel_products)
    kernel=[U.rnd(a-b)for a,b in zip(e7,P.combine(kernel_coeff,frame))]
    kernel_square=U.rnd(I(F(4,3))-U.idot(kernel_products,kernel_coeff));assert kernel_square.l>0
    profiles=[dict(node=None,label=2,contact30=m)for m in [2,3]]
    profiles +=[dict(node=node,label=label,contact30=m)for node in [5,6]for label in [3,4,5]for m in [2,3]]
    records=[];accepted_profiles=0;branch_count=0
    for profile in profiles:
        node,label,m30=profile['node'],profile['label'],profile['contact30']
        # A3 joined through a3 to the B2 end is B5. A4/5 attachment
        # instead gives two large bonds and fails the second endpoint.
        contacts={i:I()for i in C.S+[12,30,52]}
        if node is not None:contacts[node]=-D.theta(label)
        contacts[30]=-D.theta(m30)
        endpoint_data=[];rejected=None
        for endpoint in [[12,52],[30,52]]:
            ids=C.S+endpoint
            g=[[C.literal(z,i,j)for j in ids]+[contacts[i]]for i in ids]
            g.append([contacts[j]for j in ids]+[I(2)])
            # Use the old positive8x8 endpoint block; its scalar Schur
            # complement is a complete definiteness decision here.
            oldg=[row[:-1]for row in g[:-1]];oldpd=P.pd(oldg);gi,ip=P.C53.inverse(oldg)
            rhs=[contacts[i]for i in ids];schur=U.rnd(2-U.idot(rhs,U.imv(gi,rhs)))
            assert schur.l>0 or schur.h<0
            endpoint_data.append(dict(support=ids,old_PD_pivots=oldpd,inverse_pivots=ip,new_schur=schur.enc()))
            if schur.h<0:rejected=endpoint
        if rejected is not None:
            records.append(dict(profile=profile,rejected_endpoint=rejected,endpoint_checks=endpoint_data));continue
        accepted_profiles+=1
        w=[I()for _ in range(10)];white_square=I()
        if node is not None:
            c=D.theta(label)
            coeff=[-2*c/3,-c/3]if node==5 else[-c/3,-2*c/3]
            w=P.combine(coeff,[ns[5],ns[6]]);white_square=U.rnd(2*c*c/3)
        residual_products=[I(),-D.theta(m30),I()]
        coeff=U.imv(inverse,residual_products);p=P.combine(coeff,frame)
        residual=U.rnd(2-white_square-U.idot(residual_products,coeff));assert residual.l>0
        scale=U.sqrt(residual/kernel_square);branches=[]
        for sign in [-1,1]:
            n=[U.rnd(a+b+sign*scale*c)for a,b,c in zip(w,p,kernel)]
            # The retained4 and7 suffice; no center or47 is used.
            blockers={i:U.rnd(pair(n,ns[i]))for i in [4,7]}
            bad=next((i for i,q in blockers.items()if q.l>0),None)
            assert bad is not None,(profile,sign,{i:q.enc()for i,q in blockers.items()})
            branch_count+=1
            branches.append(dict(sign=sign,normal=P.enc(n),blocker=bad,product=blockers[bad].enc(),two_blocker_products={str(i):q.enc()for i,q in blockers.items()}))
        records.append(dict(profile=profile,endpoint_checks=endpoint_data,white_square=white_square.enc(),remaining_norm=residual.enc(),positive_scale=scale.enc(),branches=branches))
    assert len(records)==14 and accepted_profiles==10 and branch_count==20
    return dict(status='PASS_ALL_LABEL_ACTUAL_S6_FACET_TWO_EDGE_OBSTRUCTION_1252_3052',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in[Path(__file__),Path(C.__file__),ROOT/'results'/(C.TAG+'.json')]},
        support=C.S,frame=[12,30,52],frame_gram=[P.enc(r)for r in gram],frame_inverse_pivots=pivots,
        projected7=P.enc(e7),kernel=P.enc(kernel),kernel_square=kernel_square.enc(),profiles=records,
        profile_count=14,PD_profile_count=accepted_profiles,real_branches=branch_count,
        completeness='At a finite endpoint on S6+12+52=H4+A2+I2(6), the new root is rightH4,12,52 by the complete finite Coxeter classification; the A2 contacts are zero or exactly one bond3/4/5 at5or6. At the endpoint on S6+30+52=H4+A2+B2, right52 allows contact30 only2or3. This gives14 profiles. Exact Schur complements reject4; each of10 remaining profiles has exactly2 real norm branches in the4dimensional S6 complement after its unique white projection is removed. Every branch is strictly obtuse to retained4or7.',
        scope='No actual S6-facet wall can give finite endpoints to both12/52 and30/52 while retaining4and7, even if47 is deleted. All six nonzero A2 profiles are covered. This excludes that single-wall duty in a minimal two-wall split replacement, not all coordinated replacements or finite-volume completion.')


def main():
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
    path=ROOT/'results'/(TAG+'.json')
    if a.verify:assert json.loads(path.read_text())==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['PD_profile_count'],out['real_branches'],flush=True)

if __name__=='__main__':main()
