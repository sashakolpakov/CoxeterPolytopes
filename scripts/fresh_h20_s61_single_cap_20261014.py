#!/usr/bin/env python3
"""Complete two-anchor one-cap S61 diagnostic; exact outward interval exclusion.

This tests only a cap perpendicular to S61 and meeting BOTH old side walls96/97
at finite S61 endpoints. It does not exclude a multi-wall completion or other
support contacts. No frozen predecessor is changed.
"""
import argparse,json
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from fresh_h20_coupled_roles_20261014_field import Field,parse,enc
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_s61_single_cap_20261014'
FILES=['results/fresh_h20_e7_edge_cusp_repair_20261014_exact.json','results/fresh_h20_e7_bridge_face_20261014_exact.json','results/fresh_h20_bridge_edge_closure_20261014.json','results/fresh_h20_coupled_roles_20261014_exact.json']
BITS=384
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def pt(x):return(F(x),F(x))
def add(a,b):return(a[0]+b[0],a[1]+b[1])
def neg(a):return(-a[1],-a[0])
def sub(a,b):return add(a,neg(b))
def mul(a,b):
    x=[u*v for u in a for v in b];return(min(x),max(x))
def scale(a,q):return mul(a,pt(q))
def div(a,b):
    assert b[0]>0 or b[1]<0
    return mul(a,(1/b[1],1/b[0]))
def sqrt(a):
    assert a[0]>0
    d=1<<BITS;l=isqrt(a[0].numerator*d*d//a[0].denominator);h=isqrt(a[1].numerator*d*d//a[1].denominator)+1
    out=(F(l,d),F(h,d));assert out[0]**2<=a[0]and out[1]**2>a[1];return out

def run():
    old,bridge,edge,current=[json.loads((ROOT/p).read_text())for p in FILES]
    for d in[old,bridge,edge,current]:
        for p,h in d['sources'].items():assert sha(p)==h,p
    A=Field(old,bridge,edge);k=2;z=A.z[k]
    rec=edge['recipe'];assert int(rec['boost_integer'])==1
    assert rec['support']==current['S61_limitation']['support']
    X=tuple(parse(x,k)for x in rec['X']);W=tuple(parse(x,k)for x in rec['W']);a0=parse(rec['A'],k);w=parse(rec['W_norm'],k)
    R=tuple(parse(x,k)for x in bridge['new_normal']);N0=A.vlift(tuple(parse(x,1)for x in old['normals'][0]),1,k)
    om=tuple(parse(x,k)for x in rec['dual_weight']);p=tuple(parse(x,k)for x in rec['particular'])
    assert p==A.lin([om,R],[A.c(-1,k),A.c(F(-1,2),k)],k)
    assert A.dot(X,X,k)==A.scale(a0,-1,k)and A.dot(W,W,k)==w
    assert A.dot(R,R,k)==A.c(2,k)and A.dot(R,N0,k)==A.c(-1,k)
    assert all(A.dot(v,u,k)==z for v,u in[(X,W),(X,R),(W,R)])
    assert A.dot(om,om,k)==A.c(F(2,3),k)
    assert all(A.dot(om,v,k)==z for v in[X,W,R])
    delta=parse(edge['field']['new_positive_sqrt_radical'],k)
    assert delta==A.mul(a0,A.add(w,A.c(F(-5,6),k),k),k)
    for raw,sign in zip(edge['new_normals'],[1,-1]):
        target=tuple((A.mul(a0,A.add(p[j],A.scale(W[j],sign,k),k),k),X[j])for j in range(21))
        assert tuple(parse(x,3)for x in raw)==target
    # Therefore the serialized ambient side normals really have the stated
    # projected frame formula after division by the positive A.
    assert A.sign(a0,k)>0 and A.sign(A.add(w,A.c(-2,k),k),k)>0
    # The raw side normals, divided by positive A, project to these vectors:
    # L±=sqrt(w-5/6)*x-r/sqrt(2)±sqrt(w)*f, squares4/3.
    # The ambient normals have square2, since the omitted omega square is2/3.
    wb=A.bounds(w,k,BITS);ab=A.bounds(a0,k,BITS)
    root2=sqrt(pt(2));root5=sqrt(pt(5));rtw=sqrt(wb);rtA=sqrt(ab)
    kk=sqrt(sub(wb,pt(F(5,6))));BB=sub(scale(wb,2),pt(F(8,3)));assert BB[0]>0
    cos={2:pt(0),3:pt(F(1,2)),4:scale(root2,F(1,2)),5:scale(add(pt(1),root5),F(1,4))}
    # At m>=6, the full20 Schur pivot is <=2-3*(3/4)=-1/4.
    pivots={m:sub(pt(2),scale(mul(c,c),3))for m,c in cos.items()}
    assert all(v[0]>0 for v in pivots.values())
    gx=div(A.bounds(A.dot(X,N0,k),k,BITS),rtA)
    ge=neg(div(pt(1),root2))
    gf=div(A.bounds(A.dot(W,N0,k),k,BITS),rtw)
    rows=[]
    for m,cm in cos.items():
        for n,cn in cos.items():
            ss=add(cm,cn);dd=sub(cn,cm)
            # a*x+b*e+c*f; q=-a²+b²+c²=2.
            # a*k+b/sqrt2=cm+cn, c*sqrt(w)=cn-cm.
            DD=add(scale(mul(ss,ss),2),mul(BB,sub(pt(2),div(mul(dd,dd),wb))))
            assert DD[0]>0
            for branch in[-1,1]:
                aa=div(add(scale(mul(ss,kk),2),scale(sqrt(DD),branch)),BB)
                bb=mul(root2,sub(ss,mul(aa,kk)));cc=div(dd,rtw)
                gp95=mul(root2,bb)
                gp0=add(add(mul(aa,gx),mul(bb,ge)),mul(cc,gf))
                blocker=95 if gp95[0]>0 else 0
                witness=gp95 if blocker==95 else gp0
                assert witness[0]>0,(m,n,branch,blocker)
                # Interval identities are sanity checks; exact identities follow
                # from the displayed quadratic formula and symbolic algebra.
                norm=sub(add(mul(bb,bb),mul(cc,cc)),mul(aa,aa));assert norm[0]<=2<=norm[1]
                contact_plus=add(sub(neg(mul(aa,kk)),div(bb,root2)),mul(cc,rtw))
                contact_minus=sub(sub(neg(mul(aa,kk)),div(bb,root2)),mul(cc,rtw))
                for result,target in[(contact_plus,scale(cm,-2)),(contact_minus,scale(cn,-2))]:assert result[0]<=target[1]and target[0]<=result[1]
                rows.append(dict(labels=[m,n],branch=branch,frame_coeff_bounds=[aa,bb,cc],discriminant_bounds=DD,old0_pair_bounds=gp0,old95_pair_bounds=gp95,positive_blocker=blocker,positive_margin_lower=witness[0]))
    assert len(rows)==32 and sum(r['positive_blocker']==0 for r in rows)==16 and sum(r['positive_blocker']==95 for r in rows)==16
    paths=FILES+[f'scripts/{TAG}.py','scripts/fresh_h20_coupled_roles_20261014_field.py','scripts/fresh_h20_e7_causal_residuals_20261014_bridge_review.py']
    return enc(dict(status='PASS_COMPLETE_ONE_CAP_S61_QUADRILATERAL_GATE',sources={p:sha(p)for p in paths},interval_bits=BITS,retained_side_walls=[95,96,97],support=current['S61_limitation']['support'],frame=dict(X=X,A=a0,W=W,W_norm=w,r95=R,positive_normalization='x=X/sqrt(A),e=r95/sqrt(2),f=W/sqrt(W_norm)',old0_unitframe_pair_bounds=[gx,ge,gf]),complete_finite_contact_labels=[2,3,4,5],Schur_pivot_bounds=pivots,excluded_m_ge6_pivot_upper=F(-1,4),quadratic_formula=dict(k_squared='w-5/6',B='2w-8/3',s='cos(pi/m)+cos(pi/n)',d='cos(pi/n)-cos(pi/m)',D='2s²+B*(2-d²/w)',a='(2s*k±sqrt(D))/B',b='sqrt(2)*(s-a*k)',c='d/sqrt(w)'),candidates=rows,candidate_count=32,label_pair_count=16,positive_blocker_counts={'0':16,'95':16},minimum_positive_margin=min(r['positive_margin_lower']for r in rows),scope='All normals perpendicular toS61 that meet both retained96/97 sides at finite endpoints, up to positive normalization. Every conic branch is obtuse to retained0 or95. Thus no single added Coxeter wall can closeS61 as a compact quadrilateral. This does not exclude multiple new walls, nonright S61 contacts, or changed retained walls. No numerical acceptance or global nonexistence claim.'))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['candidate_count'],d['positive_blocker_counts'],flush=True)
