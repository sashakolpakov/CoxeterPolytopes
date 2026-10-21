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
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_affine_a17_completion_20261016'
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
    hbar=orthU(H);assert dot(hbar,hbar).h<0
    k=scale(hbar,-1/nh);eta=plus(k,scale(nu,dot(k,k)/2))
    assert has(dot(nu,eta),-1) and has(dot(eta,eta))
    # A canonical positive orthonormal4-frame in Uperp cap nuperp cap etaprep.
    free=[];pivots=[]
    for coordinate in range(21):
        v=[I(int(i==coordinate))for i in range(21)];v=orthU(v)
        v=plus(v,lin([nu,eta],[dot(v,eta),dot(v,nu)]))
        for w in free:v=plus(v,scale(w,-dot(v,w)))
        q=dot(v,v)
        if q.l<=F(1,10**30):continue
        free.append(scale(v,1/q.sqrt()));pivots.append(dict(coordinate=coordinate,norm_bounds=saved(q)))
        if len(free)==4:break
    assert len(free)==4
    for i,v in enumerate(free):
        assert all(has(dot(v,w))for w in U+[nu,eta])
        for j,w in enumerate(free):assert has(dot(v,w),int(i==j))
    GW=[row[:13]for row in GU[:13]];JW=B.inv(GW)
    contact=[[F(0)]*13 for _ in range(5)];ends=[]
    for i in range(5):
        ends.append([paths[i][-1],paths[(i+1)%5][0]])
        for w in ends[-1]:contact[i][white.index(w)]=-1
    pc=[[sum(contact[i][k]*JW[k][j]for k in range(13))for j in range(13)]for i in range(5)]
    schur=[[F(2*(i==j))-sum(pc[i][k]*contact[j][k]for k in range(13))for j in range(5)]for i in range(5)]
    assert all(sum(row)==0 for row in schur)
    weights=[F(1,len(p)+1)for p in paths]
    for i in range(5):
        for j in range(5):
            val=weights[i]+weights[(i+1)%5]if i==j else -weights[i]if j==(i-1)%5 else -weights[(i+1)%5]if j==(i+1)%5 else F(0)
            assert schur[i][j]==val
    L=[[F(int(i==j))for j in range(4)]for i in range(4)];DD=[]
    for i in range(4):
        d=schur[i][i]-sum(L[i][k]*L[i][k]*DD[k]for k in range(i));assert d>0;DD.append(d)
        for j in range(i+1,4):L[j][i]=(schur[j][i]-sum(L[j][k]*L[i][k]*DD[k]for k in range(i)))/d
    zz=[lin(free,[I(L[i][j])*I(DD[j]).sqrt()for j in range(4)])for i in range(4)]
    zz.append(lin(zz,[I(-1)]*4))
    unshifted=[plus(lin([N[j]for j in white],[I(x)for x in pc[i]]),zz[i])for i in range(5)]
    assert all(has(x)for x in lin(unshifted+[N[j]for j in white],[I(1)]*18))
    for i,v in enumerate(unshifted):
        assert has(dot(v,nu))
        for j,w in enumerate(unshifted):assert has(dot(v,w),2*int(i==j))
        for k,w in enumerate(white):assert has(dot(v,N[w]),contact[i][k])
        assert all(has(dot(v,N[w]))for w in [134,0,1])
    incident=set(front['exact_null_frontier']['incident_walls']);assert len(incident)==16
    nus=[dot(nu,v)for v in N]
    assert all(has(g)if i in incident else g.h<0 for i,g in enumerate(nus))
    shifted=[];heights=[];constraints=[]
    for i,v in enumerate(unshifted):
        # Rational integerstrictmargin: lambda greater than every dotted threshold.
        thresholds=[((dot(v,N[j])+2)/(-nus[j]),j)for j in range(135)if j not in incident]
        thresholds.append((dot(v,H)/(-nh),'center'))
        largest,j=max(thresholds,key=lambda q:q[0].h)
        height=max(1,largest.h.__ceil__()+1)
        shifted.append(plus(v,scale(nu,height)));heights.append(height)
        constraints.append(dict(connector=135+i,endpoints=ends[i],height=height,
            maximum_dotted_threshold=saved(largest),controlling_wall=j,
            all_dotted_thresholds=[dict(wall=j,bounds=saved(q))for q,j in thresholds]))
    pairs=[]
    for i,v in enumerate(shifted):
        for j,w in enumerate(N):
            g=dot(v,w)
            if j in incident:
                expected=-1 if j in ends[i]else 0;assert has(g,expected)
                relation='pi/3'if expected==-1 else 'right'
            else:assert g.h<-2,(i,j,saved(g));relation='dotted'
            pairs.append(dict(new_wall=135+i,old_wall=j,relation=relation,pair_bounds=saved(g)))
        for j in range(i):assert has(dot(v,shifted[j]))
    centerpairs=[dot(v,H)for v in N+shifted];hq=dot(H,H)
    assert all(g.h<0 for g in centerpairs)
    facetnorms=[hq-g*g/2 for g in centerpairs];assert all(q.h<0 for q in facetnorms)
    # Full18cycle, with old13nodes andthefive newconnectors.
    cycle=[]
    for i,p in enumerate(paths):cycle+=p+[135+i]
    assert len(cycle)==18
    full=N+shifted
    for i,v in enumerate(cycle):
        for j,w in enumerate(cycle):
            expected=2 if i==j else -1 if (i-j)%18 in [1,17]else 0
            assert has(dot(full[v],full[w]),expected)
    nullsum=lin([full[i]for i in cycle],[I(1)]*18)
    assert all(has(x-y*sum(heights))for x,y in zip(nullsum,nu))
    ez=list(map(F,cap['strict_escape']['frame_coefficients']));eb=list(map(F,cap['strict_escape']['dual_coefficients']))
    E=plus(lin(frame,[I(x)for x in ez]),lin(omega,[I(x)for x in eb]))
    en=dot(E,nu);assert en.h<0
    ep=[dot(E,v)for v in shifted]
    assert all(g.h<0 for g in ep)
    cuts=[]
    for i,v in enumerate(unshifted):
        cutheight=dot(v,E)/(-en)
        cuts.append(dict(connector=135+i,E1_cut_height_upper=saved(cutheight),
            chosen_height=heights[i],chosen_pair=saved(ep[i])))
    cusps=[]
    for raw in [old['completed_cusp']['point']]+[r['point']for r in old['preserved_old_bounded_cusps']]:
        cv=[bd(A.lift(parse(x,0),0))for x in raw]
        vals=[dot(v,cv)for v in shifted];assert all(g.h<0 for g in vals);cusps.append(outs(vals))
    print('PASS140facets,685newpairgates,affineG2+A17cusp,E1survives',heights,flush=True)
    src=[BASE,ARC,CAP,REVIEW,FRONT,f'scripts/{TAG}.py']+FILES+[
        'scripts/fresh_h20_coupled_roles_20261014_field.py','scripts/fresh_h20_s61_multicap_20261014_exact.py']
    return dict(status='PASS_EXACT_140_WALL_AFFINE_G2_A17_CUSP_COMPLETION',sources={p:sha(p)for p in src},
        interval_bits=BITS,retained_walls=135,new_walls=list(range(135,140)),actual_facets=140,
        definition=dict(old_paths=paths,connector_endpoints=ends,white_basis=white,positive_basis=Uids,
            positive_Gram=[[str(x)for x in row]for row in GU],
            null_ray='nu=sqrt3*n134+2*n0+n1',
            transverse_recipe='hbar=H-proj_U(H);k=hbar/(-nu.H);eta=k+(k^2/2)*nu; project coordinatevectors offU thenadd(v.eta)*nu+(v.nu)*eta; orderedGramSchmidt',
            free_pivots=pivots,free_frame_bounds=[outs(v)for v in free],
            contact_matrix=[[str(x)for x in row]for row in contact],old_projection_coefficients=[[str(x)for x in row]for row in pc],
            Schur_matrix=[[str(x)for x in row]for row in schur],cycle_weights=list(map(str,weights)),
            Schur_LDL_L=[[str(x)for x in row]for row in L],Schur_LDL_D=list(map(str,DD)),
            free_remainders='z_i=sum_j L_ij*sqrt(D_j)*free_j fori<4;z_4=-sum_i<4 z_i',
            normals='n135+i=sum_j projection_ij*oldwhite_j+z_i+height_i*nu',heights=heights,
            new_normal_bounds=[outs(v)for v in shifted]),
        height_constraints=constraints,new_old_pairs=pairs,new_new_right_pairs=10,
        common_center=dict(all140_pair_bounds=outs(centerpairs),norm_bounds=saved(hq),facet_projection_norm_bounds=outs(facetnorms)),
        complete_new_finite_vertices=0,known_total_finite_vertices_at_least=170,
        finite_census_proof='Newrootsare dotted to everyoldwalloutside21-node affineG2+A17star. Any20-root finitecornersubsetofthatstar omitsatmostone node, leavingone fullaffinecomponent. Hence noPD20support involvinganewwall.',
        new_proper_cusp=dict(incident_walls=sorted([0,1,134]+cycle),components=['affineG2','affineA17'],
            Euclidean_rank=19,horosphere_dimension=19,A17_cycle=cycle,
            A17_null_sum_multiplier=sum(heights),affine_G2_kernel=['sqrt3','2','1'],
            future_center_pair=saved(nh),bounded_link='productofaffineG2alcove andaffineA17simplex',
            exact_incidence=True,old_nonincident_pairs=outs([g for i,g in enumerate(nus)if i not in incident])),
        preserved_three_cusp_new_pairs=cusps,remaining_escape=dict(dependency=CAP,norm=cap['strict_escape']['norm'],
            new_five_pairs=outs(ep),nu_pair=saved(en),height_cut_diagnostic=cuts),
        scope='Exactcoherent140-wall localproper-cusp completion. All140facetsactual;oldfinitevertices andthreepropercuspspreserved;completezero newfinitevertexcensus. E1survivesstrictly,soinfinitevolume. Onecanonicalpathorder/freeframe;nooptimal-heightorall-orientationexclusionclaim.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],flush=True)
