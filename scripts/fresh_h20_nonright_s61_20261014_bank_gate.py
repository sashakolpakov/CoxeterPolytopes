#!/usr/bin/env python3
"""Exact common escape for the finite 847-recipe nonright-S61 bank.

Each root is defined by named support contacts, two named anchor contacts,
the norm-two equation and an isolated quadratic branch.  We do not promote
the numerical bank's remaining Coxeter fits to exact identities.  The common
escape therefore applies even if every recipe is admitted as a halfspace.
"""
import argparse,json
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from fresh_h20_coupled_roles_20261014_field import Field,parse,pd,B
from fresh_h20_s61_multicap_20261014_exact import I,fdot,pairing,FILES,BITS

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_nonright_s61_20261014_bank_gate'
BASE='results/fresh_h20_s61_multicap_20261014_exact.json'
BANK='results/fresh_h20_nonright_s61_20261014_scout.json'
POINT=[
 '13.33759289659472','25.81206850238294','39.48371047163468',
 '34.81222895686587','30.14074804607308','25.46926773925626',
 '20.79778803641542','16.12630893755058','11.45483044266173',
 '6.78335255174887','2.11187526481199','-2.55960141814890',
 '-7.23107749713380','-9.64212144420163','-12.05316478729347',
 '-13.84931928051131','-15.03058492385513','-16.21184996322297',
 '18.34312408799651','-3.64389568995840','3.45897476925484']
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def saved(x):
    """Smaller outward dyadic bounds; all decisions use 448-bit intervals."""
    den=1<<192
    return[str(F((x.l*den).__floor__(),den)),str(F((x.h*den).__ceil__(),den))]
def vecout(v):return[saved(x)for x in v]
def cost(c,J):return sum((c[i]*J[i][j]*c[j]for i in range(18)for j in range(18)),I(0))
def two_solve(M,b):
    det=M[0][0]*M[1][1]-M[0][1]*M[1][0]
    assert det.l>0 or det.h<0
    return[(b[0]*M[1][1]-M[0][1]*b[1])/det,
           (M[0][0]*b[1]-b[0]*M[1][0])/det],det
def run():
    old,bridge,edge,polar,cur=[json.loads((ROOT/p).read_text())for p in FILES]
    base=json.loads((ROOT/BASE).read_text());bank=json.loads((ROOT/BANK).read_text())
    for data in[old,bridge,edge,polar,cur,base,bank]:
        for p,h in data['sources'].items():assert sha(p)==h,p
    A=Field(old,bridge,edge);z=A.z
    def bound(x,k=3):return I(*A.bounds(x,k,BITS))
    M=[tuple(parse(x,1)for x in v)for v in old['normals']]
    N=[A.vlift(v,1)for v in M]+[A.vlift(tuple(parse(x,2)for x in bridge['new_normal']),2)]
    N +=[tuple(parse(x,3)for x in v)for v in edge['new_normals']]
    N +=[tuple(A.c(x)for x in v)for v in polar['new_normals']]
    norms=[A.c(q)for q in old['norms']]+[A.c(2)]+[parse(edge['new_normal_norm'],3)]*2+[A.c(q)for q in polar['new_norms']]
    N.append(A.vlift(tuple(parse(x,1)for x in cur['cap_R']['normal']),1))
    norms.append(A.lift(parse(cur['cap_R']['norm'],1),1));assert len(N)==101
    nd=[(bound(q)/2).sqrt()for q in norms]
    S=base['support'];assert S==bank['support'] and len(S)==18
    G=[[A.rational(A.dot(M[i],M[j],1))for j in S]for i in S]
    J=[[F(x)for x in row]for row in base['support_inverse']]
    assert pd(G) and B.inv(G)==J and all(x>=0 for row in J for x in row)
    omega=[A.lin([M[i]for i in S],[A.c(c,1)for c in row],1)for row in J]
    frame=[tuple(parse(x,2)for x in base['frame'][key])for key in['X','r95','W']]
    fsq=[parse(base['frame']['A'],2),A.c(2,2),parse(base['frame']['W_norm'],2)]
    assert [[A.dot(v,w,2)for w in frame]for v in frame]==[
        [A.scale(fsq[0],-1,2),z[2],z[2]],
        [z[2],fsq[1],z[2]],[z[2],z[2],fsq[2]]]
    assert all(A.dot(v,A.vlift(M[i],1,2),2)==z[2]for v in frame for i in S)
    den=[bound(q,2).sqrt()for q in fsq];frame3=[A.vlift(v,2)for v in frame]
    E=tuple(A.c(x)for x in POINT)
    qE=A.rational(A.dot(E,E),3);assert qE>0
    ep=[bound(A.dot(v,E))/d for v,d in zip(frame3,den)]
    we=[bound(A.dot(A.vlift(w,1),E))for w in omega]
    Eold=tuple(A.c(x)for x in base['strict_escape']['rational_point'])
    op=[bound(A.dot(v,Eold))/d for v,d in zip(frame3,den)]
    wo=[bound(A.dot(A.vlift(w,1),Eold))for w in omega]
    H=tuple(parse(x,1)for x in cur['strict_center']['point'])
    qH=A.dot(H,H,1);assert A.sign(qH,1)<0
    hd=bound(A.scale(qH,-1,1),1).sqrt()
    hp=[bound(A.dot(v,A.vlift(H,1)))/(d*hd)for v,d in zip(frame3,den)]
    wh=[bound(A.dot(w,H,1),1)/hd for w in omega]
    eps=F(base['strict_center']['epsilon'])
    C=[I(1),I(F(-1,2)),I(0)]
    qC=I(F(-3,4))+2*eps*pairing(C,hp)-eps*eps;assert qC.h<0
    EC=pairing(C,ep)+eps*bound(A.dot(E,A.vlift(H,1)))/hd;assert EC.h<0
    oldpairs=[bound(A.dot(E,n))/d for n,d in zip(N,nd)]
    hb=[tuple(parse(x,1)for x in cur['hex_basis'][key])for key in['C','U','V']]
    hq=[parse(x,1)for x in cur['hex_basis']['squared_lengths']]
    for co in cur['hexagon']['frame_coefficients']:
        oldpairs.append(sum((bound(A.lift(parse(c,0),0))*bound(A.dot(E,A.vlift(v,1)))/bound(A.lift(q,1)).sqrt()for c,v,q in zip(co,hb,hq)),I(0)))
    assert len(oldpairs)==107
    coeff=[[I(*p)for p in row]for row in base['chain']['frame_coefficient_bounds']]
    oldpairs +=[pairing(c,ep)for c in coeff]
    assert len(oldpairs)==133 and all(g.h<0 for g in oldpairs)
    # Only these five actual anchors occur; no rounded coefficient enters a
    # contact equation. Every anchor is positively normalized to norm two.
    anchors=sorted({i for r in bank['proposals']for i in r['anchors']})
    assert anchors==[0,65,84,92,95]
    P={i:[bound(A.dot(v,N[i]))/(d*nd[i])for v,d in zip(frame3,den)]for i in anchors}
    WP={i:[bound(A.dot(A.vlift(w,1),N[i]))/nd[i]for w in omega]for i in anchors}
    P[95]=[I(0),I(2).sqrt(),I(0)]
    assert all(A.dot(A.vlift(w,1),N[95])==z[3]for w in omega)
    WP[95]=[I(0)]*18
    labels={2:I(0),3:I(1),4:I(2).sqrt(),5:(1+I(5).sqrt())/2,6:I(3).sqrt()}
    expected=[([j],[m])for j in range(18)for m in[3,4,5,6]]
    expected +=[(list(s),[3,3])for s in combinations(range(18),2)]
    profiles={r['id']:r for r in bank['profiles']+bank['skipped']}
    assert len(expected)==bank['profiles_total']==len(profiles)==225
    for i,(ss,ls)in enumerate(expected):
        assert profiles[i]['support']==ss and profiles[i]['labels']==ls
        assert profiles[i]['walls']==[S[j]for j in ss]
    budgets=[];budget_counts={'cost_below_two':0,'cost_equal_two':0,'cost_above_two':0}
    for i,(ss,ls)in enumerate(expected):
        if len(ss)==2:
            aa=sum(J[j][k]for j in ss for k in ss);bb=F(0)
        else:
            j=ss[0];m=ls[0]
            aa=J[j][j]*{3:F(1),4:F(2),5:F(3,2),6:F(3)}[m]
            bb=J[j][j]/2 if m==5 else F(0)
        delta=I(aa-2)+bb*I(5).sqrt()
        if bb==0 and aa==2:kind='cost_equal_two'
        elif delta.h<0:kind='cost_below_two'
        else:assert delta.l>0;kind='cost_above_two'
        budget_counts[kind]+=1
        budgets.append(dict(profile=i,cost_Qsqrt5=[str(aa),str(bb)],classification=kind))
    assert len(bank['proposals'])==847
    records=[]
    for index,row in enumerate(bank['proposals']):
        pr=profiles[row['profile']];c=[I(0)]*18
        for j,m in zip(pr['support'],pr['labels']):c[j]=labels[m]
        assert all(x.l-F(1,10**12)<F(str(v))<x.h+F(1,10**12)for x,v in zip(c,row['contacts']))
        mass=2-cost(c,J);assert mass.h<0
        ii,jj=row['anchors'];L=[P[ii],P[jj]]
        rhs=[sum((a*b for a,b in zip(c,WP[i])),I(0))-labels[m]for i,m in zip([ii,jj],row['anchor_labels'])]
        cross=[L[0][1]*L[1][2]-L[0][2]*L[1][1],L[0][2]*L[1][0]-L[0][0]*L[1][2],L[0][0]*L[1][1]-L[0][1]*L[1][0]]
        free=max(range(3),key=lambda j:max(abs(cross[j].l),abs(cross[j].h)))
        other=[j for j in range(3)if j!=free]
        LL=[[L[i][j]for j in other]for i in range(2)]
        pp,det=two_solve(LL,rhs);vv,_=two_solve(LL,[-L[i][free]for i in range(2)])
        p=[I(0)]*3;v=[I(0)]*3;v[free]=I(1)
        for j,a,b in zip(other,pp,vv):p[j]=a;v[j]=b
        qa=fdot(v,v);qb=2*fdot(p,v);qc=fdot(p,p)-mass
        assert qa.l>0 or qa.h<0
        disc=qb*qb-4*qa*qc;assert disc.l>0
        roots=[]
        # The decimal proposal is solely a branch selector. The exact normal
        # is the unique algebraic solution in this rational box.
        box=[(F(str(x))-F(1,10000),F(str(x))+F(1,10000))for x in row['residual_frame']]
        rejected=[]
        for branch in[-1,1]:
            t=(-qb+branch*disc.sqrt())/(2*qa)
            y=[a+t*b for a,b in zip(p,v)]
            inside=all(l<x.l<=x.h<h for x,(l,h)in zip(y,box))
            if inside:roots.append((branch,y))
            else:
                assert any(x.h<l or x.l>h for x,(l,h)in zip(y,box))
                rejected.append(branch)
        assert len(roots)==1 and len(rejected)==1,(index,'branch isolation')
        branch,y=roots[0];assert y[0].l>0
        escape=pairing(y,ep)-sum((a*b for a,b in zip(c,we)),I(0));assert escape.h<0,(index,escape.out())
        oldcut=pairing(y,op)-sum((a*b for a,b in zip(c,wo)),I(0));assert oldcut.l>0
        center=fdot(y,C)+eps*(pairing(y,hp)-sum((a*b for a,b in zip(c,wh)),I(0)))
        assert center.h<0
        records.append(dict(index=index,profile=row['profile'],anchors=[ii,jj],anchor_labels=row['anchor_labels'],
            free_coordinate=free,branch=branch,selector_box=[[str(a),str(b)]for a,b in box],
            linear_determinant=saved(det),quadratic_coefficients=vecout([qa,qb,qc]),discriminant=saved(disc),
            residual_norm=saved(mass),residual_frame=vecout(y),escape_pair=saved(escape),
            former_escape_pair=saved(oldcut),common_center_pair=saved(center)))
        if (index+1)%100==0:print('exact bank roots',index+1,flush=True)
    # Full-real causal reduction: J>=0 implies w.N<=0 for every nonsupport
    # old wall, in particular all 29 polygon sides.
    # Thus y.N<=0. The inherited 29-gon cone is future timelike except zero.
    # The following exact positive spatial relation fixes time orientation.
    assert base['face']['cyclic_walls']==[95,96]+list(range(107,133))+[97]
    assert len(base['face']['vertices'])==29
    mid=coeff[12:14]
    assert mid[0][0].l>0 and mid[0][1].h<0 and mid[0][2].l>0 and mid[1][2].h<0
    assert all(x.h<0 for x in[I(*v['norm_bound'])for v in base['face']['vertices']])
    paths=FILES+[BASE,BANK,'scripts/fresh_h20_coupled_roles_20261014_field.py',
        'scripts/fresh_h20_s61_multicap_20261014_exact.py',f'scripts/{TAG}.py']
    return dict(status='EXACT_COMMON_ESCAPE_FOR_ALL_847_ALGEBRAIC_RECIPES',sources={p:sha(p)for p in paths},
        interval_bits=BITS,scope='Finite archived recipe bank, not all normals, contact profiles, anchor choices or continuous families. No exact claim for the bank\'s full Coxeter pair gates.',
        support=S,support_inverse_nonnegative=True,anchors=anchors,recipe_count=847,
        exact_profile_budgets=budgets,profile_budget_counts=budget_counts,
        escape=dict(rational_point=POINT,norm=str(qE),old133_pairs=vecout(oldpairs),frame_pairings=vecout(ep),
            common_center_pair=saved(EC),maximum_bank_pair_upper=str(max(F(r['escape_pair'][1])for r in records))),
        common_center=dict(inherited_recipe=base['strict_center']['recipe'],epsilon=str(eps),norm_bound=saved(qC),
            maximum_bank_pair_upper=str(max(F(r['common_center_pair'][1])for r in records))),
        records=records,
        all_real_residual_lemma=dict(polygon_vertices=29,spatial_positive_relation='n119+n120+sqrt(2)*sqrt(a^2+1)*n95=(2a,0,0)',
            conclusion='Any nonzero coherent acute added root r=-w+y has future timelike y, hence w^2>r^2. The y=0 case contradicts strict interior because Jc is nonnegative.',
            dependency='Frozen compact 29-gon and complete old pair acuteness; no bounded-label assumption.'))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert data==json.loads(path.read_text())
    else:path.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['recipe_count'],'norm',data['escape']['norm'],flush=True)
