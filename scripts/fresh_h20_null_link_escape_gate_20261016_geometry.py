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
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_null_link_escape_gate_20261016_geometry'
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
    paths=[[3,4],[6,7,8,9,10],[12,13,14,15],[17],[31]]
    white=sum(paths,[]);assert white==[3,4,6,7,8,9,10,12,13,14,15,17,31]
    Uids=white+[0,1];U=[N[i]for i in Uids]
    GU=[[A.rational(A.dot(M[i],M[j],1),1)for j in Uids]for i in Uids];JU=B.inv(GU)
    assert all(has(dot(nu,v))for v in U)
    def orthU(v):return plus(v,lin(U,[-sum((I(JU[i][j])*dot(v,U[j])for j in range(15)),I(0))for i in range(15)]))
    completed=json.loads((ROOT/'results/fresh_h20_affine_a17_completion_20261016.json').read_text())
    free=[[I(*v)for v in row]for row in completed['definition']['free_frame_bounds']]
    assert len(free)==4
    for i,v in enumerate(free):
        for j,w in enumerate(free):assert has(dot(v,w),int(i==j))
        assert all(has(dot(v,w))for w in U+[nu])
    ez=list(map(F,cap['strict_escape']['frame_coefficients']));eb=list(map(F,cap['strict_escape']['dual_coefficients']))
    E=plus(lin(frame,[I(x)for x in ez]),lin(omega,[I(x)for x in eb]))
    incident=set(front['exact_null_frontier']['incident_walls'])
    def coeff(v):return dict(free=outs([dot(f,v)for f in free]),null=saved(dot(nu,v)),white=outs([dot(N[j],v)for j in white]))
    ns=[coeff(v)for v in N]
    # For the sole zero-profile/all-S-right case, y is the line orthogonal
    # to the timelike frame of134 and the residual oldroot0.
    y134=[I(*x)for x in cap['definition']['frame_bounds']]
    p0=[dot(v,N[0])for v in frame]
    z=y134[0]/y134[2]
    e=-(p0[0]+p0[2]*z)/p0[1]
    axis=[I(1),e,z]
    axis95=sum((axis[i]*dot(frame[i],N[95])for i in range(3)),I(0))
    axis107=sum((axis[i]*dot(frame[i],N[107])for i in range(3)),I(0))
    assert axis95.h<0 and axis107.l>0
    files=[BASE,ARC,CAP,REVIEW,FRONT,'results/fresh_h20_affine_a17_completion_20261016.json',
           'results/fresh_h20_null_link_profiles_20261016.json',f'scripts/{TAG}.py']+FILES+[
        'scripts/fresh_h20_coupled_roles_20261014_field.py','scripts/fresh_h20_s61_multicap_20261014_exact.py']
    return dict(status='PASS_EXACT_NULL_LINK_OPTIMIZATION_COEFFICIENT_ENCLOSURES',
       sources={p:sha(p)for p in files},white=white,incident=sorted(incident),normals=ns,
       center=coeff(H),escape=coeff(E),positive_Gram=[[str(x)for x in row[:13]]for row in GU[:13]],
       zero_profile_all_S_right=dict(axis_frame=outs(axis),pair95=saved(axis95),pair107=saved(axis107),
        proof='RighttoS18 forces pureframe;right134givesz/t=y134_t/y134_f;right0fixese/t. If t0 thenwholevector0. Axis95negativeandaxis107positive force t0 byacuteness, hence no norm2wall.'),
       scope='Exact symbolic orthonormal quotient parameterization inherited from A17 construction; rigorous rational coefficient enclosures. No optimizer conclusion in this coefficient file.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert data==json.loads(path.read_text())
    else:path.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],flush=True)
