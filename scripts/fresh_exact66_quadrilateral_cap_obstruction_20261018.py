#!/usr/bin/env python3
"""Exact necessary norm obstruction to a single four-finite-corner cap."""
import argparse,json
from hashlib import sha256
from itertools import product
import fresh_exact65_escape33_cap_20261018 as P
A,F,D=P.A,P.F,P.D
ROOT=P.ROOT;TAG='fresh_exact66_quadrilateral_cap_obstruction_20261018';FRAME=[30,52,47,76]
V,add,sub,mul,neg,dot=A.V,A.add,A.sub,A.mul,A.neg,A.dot
enc=A.enc

def read(s):return json.loads((ROOT/f'results/{s}.json').read_text())

def calculate(corner):
    seed=read(A.TAG);prior=read(F.TAG);current=read(P.TAG)
    for d in(seed,prior,current):
        for path,h in d['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==h,path
    roots={int(i):list(map(V,x))for i,x in seed['physical_input']['roots'].items()}
    metric=[list(map(V,r))for r in seed['physical_input']['metric']]
    assert all(seed['white_covectors'][str(i)]==enc([V(0)]*6)for i in FRAME[:-1])
    assert prior['physical']['profile']==0 and prior['physical']['new_wall']==76
    identities=[];g=[]
    for i in FRAME:
        row=[]
        for j in FRAME:
            if i==j:v=V(2)
            elif 76 in(i,j):v=V(prior['physical']['native_products'][str(j if i==76 else i)])
            else:
                v,pr=A.refine(dot(roots[i],metric,roots[j]),i,j)
                if pr:identities.append(pr)
            row.append(v)
        g.append(row)
    zero=V(0);theta4=A.E.theta(4)
    assert g[0][1]==neg(theta4)and g[0][3]==g[1][2]==g[2][3]==zero
    assert g[0][2][1]<-2 and g[1][3][1]<-2
    inv,det=A.P.inverse(g);assert det[1]<0
    # The original projected metric has signature(3,1). Nonsingularity
    # makes these four independent walls a basis of that entire space.
    accepted=[];rejected=[]
    for pi,ms in product(range(7),product(range(2,6),repeat=4)):
        duties=[corner(pi,b,ms[i],ms[(i+1)%4])for i,b in enumerate([4,2,2,2])]
        row=dict(profile=pi,labels=ms,duties=duties)
        if any(x is None for x in duties):rejected.append(row);continue
        rhs=[neg(A.E.theta(m))for m in ms]
        profile_square=V(2)if pi==0 else sub(V(2),mul(V(A.F(2,3)),mul(A.E.theta(3+(pi-1)%3),A.E.theta(3+(pi-1)%3))))
        square=dot(rhs,inv,rhs);gap=sub(square,profile_square)
        assert gap[1]<0 or gap[0]>0,(pi,ms,enc(gap))
        row.update(rhs=rhs,projected_square=square,required_profile_square=profile_square,norm_gap=gap,
            sign=-1 if gap[1]<0 else 1);accepted.append(row)
    assert len(accepted)==84 and len(rejected)==1708
    highs=[];left=A.sqrt(V(3))
    for axis,opposite in[(2,0),(3,1)]:
        for m in[2,3]:
            e=[V(-int(i==axis))for i in range(4)];b=[neg(A.E.theta(m))if i==opposite else zero for i in range(4)]
            aa=dot(e,inv,e);bb=mul(V(2),dot(e,inv,b));cc=sub(dot(b,inv,b),V(2))
            atleft=add(add(mul(aa,V(3)),mul(bb,left)),cc)
            derivative=add(mul(V(2),mul(aa,left)),bb)
            assert aa[0]>0 and bb[0]>=0 and atleft[0]>0 and derivative[0]>0
            highs.append(dict(high_wall=FRAME[axis],opposite_wall=FRAME[opposite],opposite_label=m,
                profile=0,coefficients=[aa,bb,cc],theta_interval=[left[0],A.F(2)],
                value_at_sqrt3=atleft,derivative_at_sqrt3=derivative,strictly_increasing=True))
    return dict(frame=FRAME,cycle_labels=[4,2,2,2],dotted_diagonals=[[30,47],[52,76]],
        original_native_identity_refinements=identities,frame_Gram=g,frame_inverse=inv,frame_determinant=det,
        low=dict(profiles=list(range(7)),labels=[2,3,4,5],complete_recipes=1792,finite_duty_rejected=rejected,
            eligible=accepted,eligible_count=84,norm_excluded_count=84),
        high=highs,all_declared_profile_and_finite_label_caps_excluded=True)

def run():
    out=calculate(lambda pi,b,m,n:D.T.P.corner3(0,0,pi,b,m,n))
    sources={f'{folder}/{s}.{ext}':sha256((ROOT/f'{folder}/{s}.{ext}').read_bytes()).hexdigest()
        for s,pairs in[(TAG,[('scripts','py')]),(P.TAG,[('scripts','py'),('results','json'),('results','md')]),
            (F.TAG,[('scripts','py'),('results','json'),('results','md')]),(A.TAG,[('scripts','py'),('results','json'),('results','md')])]
        for folder,ext in pairs}
    return dict(status='PASS_EXACT66_FOUR_FINITE_DUTY_CAP_OBSTRUCTION',sources=sources,necessary_norm_certificate=out,
        scope=dict(fixed_exact66_parent=True,seven_profiles_only=True,all_finite_labels=True,
            ideal_or_dotted_duties_not_covered=True,other_incidence_or_parent_not_excluded=True,
            chamber_or_cone_search=False,finite_volume_completion_excluded=False,new_geometric_example=False),
        deletion_debit=dict(deleted_old_walls=[],old_finite_vertices_lost=0))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    p=ROOT/f'results/{TAG}.json';assert a.verify or not p.exists();out=enc(run())
    if a.verify:assert json.loads(p.read_text())==out
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],flush=True)
