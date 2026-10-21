#!/usr/bin/env python3
"""Exact five-wall completion of the persistent affineG2 cusp in H20.

Algebraic vectors are defined by projection, positive-square-root Gram-Schmidt
and a rational LDL factor. Interval arithmetic verifies strict inequalities;
exact Gram identities follow from the recorded construction equations.
"""
import argparse,json
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
from fresh_h20_coupled_roles_20261014_field import Field,parse,B
from fresh_h20_s61_multicap_20261014_exact import I,FILES,BITS
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_off_cusp_cap_20261016'
BASE='results/fresh_h20_s61_multicap_20261014_exact.json'
ARC='results/fresh_h20_dotted_arc_cap_20261016.json'
CAP='results/fresh_h20_complementary_cap_20261016.json'
REVIEW='results/fresh_h20_complementary_cap_20261016_review.json'
FRONT='results/fresh_h20_complementary_frontier_20261016.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def saved(x):
    den=1<<192
    return[str(F((x.l*den).__floor__(),den)),str(F((x.h*den).__ceil__(),den))]
def outs(v):return [saved(x)for x in v]
def has(x,t=0):return x.l<=t<=x.h

def run():
    old,bridge,edge,polar,cur=[json.loads((ROOT/p).read_text())for p in FILES]
    base,arc,cap,review,front=[json.loads((ROOT/p).read_text())for p in [BASE,ARC,CAP,REVIEW,FRONT]]
    for doc in [old,bridge,edge,polar,cur,base,arc,cap,review,front]:
        for p,h in doc['sources'].items():assert sha(p)==h,p
    assert review['sources'][CAP]==sha(CAP)
    A=Field(old,bridge,edge);Q=A.Q
    rows=[[(j,q)for j,q in enumerate(row)if q]for row in Q]
    def bd(x,k=3):return I(*A.bounds(x,k,BITS))
    def dot(v,w):return sum((v[i]*w[j]*q for i,row in enumerate(rows)for j,q in row),I(0))
    def lin(vs,cs):return [sum((v[j]*c for v,c in zip(vs,cs)),I(0))for j in range(21)]
    def plus(v,w):return [x+y for x,y in zip(v,w)]
    def scale(v,c):return [x*c for x in v]
    M=[tuple(parse(x,1)for x in v)for v in old['normals']]
    exact=[A.vlift(v,1)for v in M]+[A.vlift(tuple(parse(x,2)for x in bridge['new_normal']),2)]
    exact +=[tuple(parse(x,3)for x in v)for v in edge['new_normals']]
    exact +=[tuple(A.c(x)for x in v)for v in polar['new_normals']]
    exact.append(A.vlift(tuple(parse(x,1)for x in cur['cap_R']['normal']),1))
    norms=[A.c(q)for q in old['norms']]+[A.c(2)]+[parse(edge['new_normal_norm'],3)]*2+[A.c(q)for q in polar['new_norms']]+[A.lift(parse(cur['cap_R']['norm'],1),1)]
    N=[scale([bd(x)for x in v],1/(bd(q)/2).sqrt())for v,q in zip(exact,norms)]
    hb=[tuple(parse(x,1)for x in cur['hex_basis'][key])for key in ['C','U','V']]
    hq=[parse(x,1)for x in cur['hex_basis']['squared_lengths']]
    hv=[scale([bd(x,1)for x in v],1/bd(q,1).sqrt())for v,q in zip(hb,hq)]
    for co in cur['hexagon']['frame_coefficients']:N.append(lin(hv,[bd(A.lift(parse(c,0),0))for c in co]))
    assert len(N)==107
    rawframe=[tuple(parse(x,2)for x in base['frame'][key])for key in ['X','r95','W']]
    fsq=[parse(base['frame']['A'],2),A.c(2,2),parse(base['frame']['W_norm'],2)]
    frame=[scale([bd(x,2)for x in v],1/bd(q,2).sqrt())for v,q in zip(rawframe,fsq)]
    for co in base['chain']['frame_coefficient_bounds']:N.append(lin(frame,[I(*x)for x in co]))
    S=base['support'];J=[[F(x)for x in row]for row in base['support_inverse']]
    omega=[lin([N[j]for j in S],[I(x)for x in row])for row in J]
    r133=plus(lin(frame,[I(*x)for x in arc['definition']['frame_bounds']]),scale(omega[S.index(31)],-I(*arc['definition']['dual_coefficient_bounds'])))
    r134=plus(lin(frame,[I(*x)for x in cap['definition']['frame_bounds']]),scale(omega[S.index(5)],-1))
    N +=[r133,r134];assert len(N)==135
    Hraw=tuple(parse(x,1)for x in cur['strict_center']['point']);Hq=A.dot(Hraw,Hraw,1)
    Hunit=scale([bd(x,1)for x in Hraw],1/bd(A.scale(Hq,-1,1),1).sqrt())
    H=plus(lin(frame,[I(1),I(F(-1,2)),I(0)]),scale(Hunit,F(base['strict_center']['epsilon'])))
    assert dot(H,H).h<0 and all(dot(H,v).h<0 for v in N)
    nu=lin([N[134],N[0],N[1]],[I(3).sqrt(),I(2),I(1)])
    nh=dot(nu,H);assert nh.h<0 and has(dot(nu,nu))
    completion=json.loads((ROOT/'results/fresh_h20_affine_a17_completion_20261016.json').read_text())
    reviewed=json.loads((ROOT/'results/fresh_h20_affine_a17_completion_20261016_review.json').read_text())
    assert reviewed['sources']['results/fresh_h20_affine_a17_completion_20261016.json']==sha('results/fresh_h20_affine_a17_completion_20261016.json')
    N +=[[I(*x)for x in row]for row in completion['definition']['new_normal_bounds']]
    assert len(N)==140
    coef=[F(8),F(-3,2),F(-77,10)];d=F(14,5);j17=S.index(17);j31=S.index(31)
    r=plus(lin(frame,[I(x)for x in coef]),scale(omega[j17],-d))
    q=-coef[0]**2+coef[1]**2+coef[2]**2+d*d*J[j17][j17]
    assert q==F(73,50)>0 and has(dot(r,r),q)
    finite=set(S)-{17};pairs=[]
    for i,n in enumerate(N):
        g=dot(r,n)
        if i in finite:assert has(g);relation='right'
        else:assert g.h<0 and g.h*g.h>2*q,(i,saved(g));relation='dotted'
        pairs.append(dict(wall=i,relation=relation,pair_bounds=saved(g)))
    assert len(finite)==17
    centerpairs=[dot(n,H)for n in N]+[dot(r,H)];hq=dot(H,H)
    assert hq.h<0 and all(g.h<0 for g in centerpairs)
    facetnorms=[hq-g*g/(2 if i<140 else q)for i,g in enumerate(centerpairs)]
    assert all(x.h<0 for x in facetnorms)
    ef=list(map(F,cap['strict_escape']['frame_coefficients']));eb=list(map(F,cap['strict_escape']['dual_coefficients']))
    E1=plus(lin(frame,[I(x)for x in ef]),lin(omega,[I(x)for x in eb]))
    cut=-coef[0]*ef[0]+coef[1]*ef[1]+coef[2]*ef[2]-d*sum(J[j17][k]*eb[k]for k in range(18))
    assert cut==F(961,1000)>0 and dot(r,E1).l>0 and has(dot(r,E1),cut)
    neweb=eb[:];neweb[j17]=F(-3,10);neweb[j31]=F(-2,5)
    E2=plus(lin(frame,[I(x)for x in ef]),lin(omega,[I(x)for x in neweb]))
    eq=-ef[0]**2+ef[1]**2+ef[2]**2+sum(neweb[i]*J[i][j]*neweb[j]for i in range(18)for j in range(18))
    assert eq>0 and has(dot(E2,E2),eq)
    escape_pairs=[dot(E2,n)for n in N]+[dot(E2,r)]
    assert all(g.h<0 for g in escape_pairs)
    new_escape_pair=-coef[0]*ef[0]+coef[1]*ef[1]+coef[2]*ef[2]-d*sum(J[j17][k]*neweb[k]for k in range(18))
    assert new_escape_pair==F(-19,1000) and has(escape_pairs[-1],new_escape_pair)
    assert dot(E2,H).h<0
    nullpair=dot(r,nu);assert nullpair.h<0
    oldcusps=[]
    for raw in [old['completed_cusp']['point']]+[row['point']for row in old['preserved_old_bounded_cusps']]:
        cv=[bd(A.lift(parse(x,0),0))for x in raw];g=dot(r,cv);assert g.h<0;oldcusps.append(saved(g))
    src=[BASE,ARC,CAP,REVIEW,FRONT,'results/fresh_h20_affine_a17_completion_20261016.json',
         'results/fresh_h20_affine_a17_completion_20261016_review.json',f'scripts/{TAG}.py']+FILES+[
         'scripts/fresh_h20_coupled_roles_20261014_field.py','scripts/fresh_h20_s61_multicap_20261014_exact.py']
    return dict(status='PASS_EXACT_141_WALL_OFF_CUSP_E1_CUTTER',sources={p:sha(p)for p in src},
       definition=dict(formula='r=8*x-(3/2)*e-(77/10)*f-(14/5)*omega17',frame_coefficients=list(map(str,coef)),
           support_wall=17,dual_coefficient=str(d),norm=str(q),physical_normal_bounds=outs(r)),
       retained_walls=140,new_wall=140,actual_facets=141,pairs=pairs,right_neighbors=sorted(finite),
       relations=dict(right=17,dotted=123),all26_chain_contacts_dotted=True,
       common_center=dict(norm_bounds=saved(hq),all141_pair_bounds=outs(centerpairs),facet_projection_norm_bounds=outs(facetnorms)),
       E1_cut=dict(pair=str(cut),pair_bounds=saved(dot(r,E1)),old_escape=CAP),
       new_proper_cusp_preserved=dict(null_ray='sqrt3*n134+2*n0+n1',new_pair_bounds=saved(nullpair),
           link='affineG2+affineA17 unchanged because the new wall is strictly negative at nu'),
       preserved_three_other_cusp_pairs=oldcusps,
       complete_new_finite_vertices=0,known_total_finite_vertices_at_least=170,
       finite_census_proof='Only17 oldfinite-neighbor walls; allother123oldwalls strictlydotted. No20-root PDfinitecornersupport caninclude the newwall. Oldfiniteverticespersist by commoncenter andacuteness.',
       strict_escape=dict(formula='E2=x-e/100-(49/50)*f-(3/10)*omega17-(2/5)*omega31-eps*sum(other16omegas)',
           epsilon='1/100000',frame_coefficients=list(map(str,ef)),dual_coefficients=list(map(str,neweb)),
           norm=str(eq),all141_pair_bounds=outs(escape_pairs),new_cap_pair=str(new_escape_pair),
           center_pair_bounds=saved(dot(E2,H))),
       scope='Exact coherent141-wall Coxeter extension retaining the completedG2+A17cusp; cutsE1,all141facetsactual,oldfinitevertices and4propercuspspersist,no newfinitevertices. E2 is strictspacelike, so volume remains infinite. No full ambientcone census and no new dimension claim.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert data==json.loads(path.read_text())
    else:path.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['strict_escape']['norm'],flush=True)
