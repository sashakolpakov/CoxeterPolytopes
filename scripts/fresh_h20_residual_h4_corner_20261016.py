#!/usr/bin/env python3
"""Exact three-cap cluster completing an A1^4 residual H4 corner in H20.

All four residual normals and their unique finite corner have rational
coefficients in the orthonormal frame plus the two isolated support duals.
Outward intervals check the inherited physical old-wall inequalities.
"""
import argparse,json
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
from fresh_h20_coupled_roles_20261014_field import Field,parse,B
from fresh_h20_s61_multicap_20261014_exact import I,FILES,BITS
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_residual_h4_corner_20261016'
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
    PREV='results/fresh_h20_off_cusp_cap_20261016.json'
    PREV_REVIEW='results/fresh_h20_off_cusp_cap_20261016_review.json'
    previous=json.loads((ROOT/PREV).read_text());pr=json.loads((ROOT/PREV_REVIEW).read_text())
    assert previous['status']=='PASS_EXACT_141_WALL_OFF_CUSP_E1_CUTTER'
    assert pr['sources'][PREV]==sha(PREV)
    for doc in [previous,pr]:
        for p,h in doc['sources'].items():assert sha(p)==h,p
    j17=S.index(17);j31=S.index(31)
    pc=list(map(F,previous['definition']['frame_coefficients']))
    pd=F(previous['definition']['dual_coefficient']);pq=F(previous['definition']['norm'])
    assert pc==[F(8),F(-3,2),F(-77,10)] and pd==F(14,5) and pq==F(73,50)
    r140=plus(lin(frame,[I(x)for x in pc]),scale(omega[j17],-pd))
    N.append(r140);norms_old=[F(2)]*140+[pq];assert len(N)==141
    MIX='results/fresh_h20_mixed_isolated_cap_20261016.json'
    mixed=json.loads((ROOT/MIX).read_text())
    assert mixed['status']=='PASS_EXACT_142_WALL_MIXED_ISOLATED_E2_CUTTER'
    for p,h in mixed['sources'].items():assert sha(p)==h,p
    g5=[F(-1),F(1),F(1),F(1,2),F(1,2)]
    def fd5(v,w):return sum(g*a*b for g,a,b in zip(g5,v,w))
    def vec5(v):return plus(lin(frame,[I(x)for x in v[:3]]),lin([omega[j17],omega[j31]],[I(x)for x in v[3:]]))
    recipes=[
        [F(17),F(-3,2),F(-67,4),F(-2),F(-3)],
        [F(13),F(0),F(-64,5),F(-21,10),F(-3)],
        [F(2727,100),F(-773,300),F(-267,10),F(-5),F(-5)],
        [F(276851,39846),F(-8739,13282),F(-135452,19923),F(-1),F(-3,2)]]
    norms4=[fd5(v,v)for v in recipes]
    assert norms4==[F(5,16),F(309,200),F(19717,22500),F(6567095,1058469144)]
    gram4=[[fd5(v,w)for w in recipes]for v in recipes]
    assert all(gram4[i][j]==(norms4[i]if i==j else 0)for i in range(4)for j in range(4))
    assert mixed['definition']['frame_coefficients']==list(map(str,recipes[0][:3]))
    assert mixed['definition']['dual_coefficients']==['2','3']
    N.append(vec5(recipes[0]));norms_old.append(norms4[0]);assert len(N)==142
    new=[vec5(v)for v in recipes[1:]]
    S16=set(S)-{17,31};pairs=[];finite_by_new=[]
    for k,r in enumerate(new):
        q=norms4[k+1];assert has(dot(r,r),q)
        finite=S16|{141}|set(range(142,142+k))
        if k==0:finite|={95}
        finite_by_new.append(sorted(finite))
        for i,n in enumerate(N):
            g=dot(r,n)
            if i in finite:assert has(g);relation='right'
            else:assert g.h<0 and g.h*g.h>norms_old[i]*q,(142+k,i,saved(g));relation='dotted'
            pairs.append(dict(walls=[i,142+k],norms=[str(norms_old[i]),str(q)],relation=relation,pair_bounds=saved(g)))
        N.append(r);norms_old.append(q)
    assert len(N)==145 and len(pairs)==429
    assert sum(p['relation']=='right'for p in pairs)==55
    # Exact right relations follow dual orthogonality, the displayed Gram4,
    # and the zero e coefficient of u (old95 is sqrt2*e).
    centerpairs=[dot(n,H)for n in N];hq=dot(H,H)
    assert hq.h<0 and all(g.h<0 for g in centerpairs)
    facetnorms=[hq-g*g/nq for nq,g in zip(norms_old,centerpairs)]
    assert all(x.h<0 for x in facetnorms)
    corner5=[F(1),F(-5971425,64223509),F(-314433802,321117545),F(-46580508,321117545),F(-67243924,321117545)]
    cq=fd5(corner5,corner5)
    assert cq==F(-10669433158276,103116477706827025)<0
    assert all(fd5(corner5,r)==0 for r in recipes)
    C=vec5(corner5);incident=sorted(S16|{141,142,143,144})
    assert len(incident)==20 and has(dot(C,C),cq)
    cornerpairs=[dot(C,n)for n in N]
    for i,g in enumerate(cornerpairs):
        if i in incident:assert has(g)
        else:assert g.h<0,(i,saved(g))
    assert dot(C,H).h<0
    Gs=[[F(x)for x in row]for row in base['support_gram']]
    p16=sorted(S16);GG=[[Gs[S.index(i)][S.index(j)]for j in p16]for i in p16]
    G20=[row+[F(0)]*4 for row in GG]+[[F(0)]*16+row for row in gram4]
    A20=[row[:]for row in G20];pivots=[]
    for k in range(20):
        pivot=A20[k][k];assert pivot>0;pivots.append(str(pivot))
        for i in range(k+1,20):
            for j in range(k+1,20):A20[i][j]-=A20[i][k]*A20[k][j]/pivot
    # The complete finite-neighbor list of any new cap forces this one corner.
    neighbors={142:S16|{95,141,143,144},143:S16|{141,142,144},144:S16|{141,142,143}}
    from itertools import combinations
    supports=[]
    for i,nn in neighbors.items():
        for co in combinations(sorted(nn),19):
            if 95 in co and any(k in co for k in [141,143,144]):continue
            support=tuple(sorted(co+(i,)));assert support==tuple(incident);supports.append(list(support))
    assert len(supports)==3
    # A strict escape survives the complete cluster.
    ef=list(map(F,mixed['strict_escape']['frame_coefficients']))
    eb=list(map(F,mixed['strict_escape']['dual_coefficients']))
    E3=plus(lin(frame,[I(x)for x in ef]),lin(omega,[I(x)for x in eb]))
    e3chart=ef+[eb[j17],eb[j31]]
    cuts=[fd5(v,e3chart)for v in recipes[1:]]
    assert cuts[0]==F(21,250)>0 and cuts[1]==F(1403,30000)>0
    neweb=eb[:];neweb[j31]=F(-1,20)
    E4=plus(lin(frame,[I(x)for x in ef]),lin(omega,[I(x)for x in neweb]))
    eq=-ef[0]**2+ef[1]**2+ef[2]**2+sum(neweb[i]*J[i][j]*neweb[j]for i in range(18)for j in range(18))
    assert eq==F(33750061,5000000000)>0 and has(dot(E4,E4),eq)
    ep=[dot(E4,n)for n in N];assert all(g.h<0 for g in ep)
    e4chart=ef+[neweb[j17],neweb[j31]]
    newep=[fd5(v,e4chart)for v in recipes[1:]]
    assert all(x<0 for x in newep) and all(has(ep[142+i],x)for i,x in enumerate(newep))
    assert dot(E4,H).h<0
    # A whole rational rectangle of strict escapes remains, not only E4.
    abox=[F(29,100),F(31,100)];bbox=[F(3,100),F(3,50)];rect=[]
    for aa in abox:
        for bb in bbox:
            cb=neweb[:];cb[j17]=-aa;cb[j31]=-bb
            EE=plus(lin(frame,[I(x)for x in ef]),lin(omega,[I(x)for x in cb]))
            gp=[dot(EE,n)for n in N];ch=dot(EE,H)
            assert all(x.h<0 for x in gp) and ch.h<0
            rect.append(dict(a=str(aa),b=str(bb),all145_pair_bounds=outs(gp),center_pair_bounds=saved(ch)))
    minrect=eq+(abox[0]**2-F(3,10)**2+bbox[0]**2-F(1,20)**2)/2
    assert minrect==F(15000061,5000000000)>0
    nullpairs=[dot(r,nu)for r in new];assert all(g.h<0 for g in nullpairs)
    oldcusps=[]
    for raw in [old['completed_cusp']['point']]+[row['point']for row in old['preserved_old_bounded_cusps']]:
        cv=[bd(A.lift(parse(x,0),0))for x in raw];gs=[dot(r,cv)for r in new]
        assert all(g.h<0 for g in gs);oldcusps.append(outs(gs))
    src=[BASE,ARC,CAP,REVIEW,FRONT,PREV,PREV_REVIEW,MIX,
         'results/fresh_h20_affine_a17_completion_20261016.json',
         'results/fresh_h20_affine_a17_completion_20261016_review.json',f'scripts/{TAG}.py']+FILES+[
         'scripts/fresh_h20_coupled_roles_20261014_field.py','scripts/fresh_h20_s61_multicap_20261014_exact.py']
    return dict(status='PASS_EXACT_145_WALL_RESIDUAL_H4_FINITE_CORNER_CLUSTER',sources={p:sha(p)for p in src},
       definition=dict(chart='coefficients in x,e,f,omega17,omega31; chartGram=diag(-1,1,1,1/2,1/2)',
           four_walls=[141,142,143,144],four_recipes=[list(map(str,v))for v in recipes],
           four_norms=list(map(str,norms4)),four_gram=[list(map(str,row))for row in gram4],
           new_normal_bounds=[outs(r)for r in new]),
       retained_walls=142,new_walls=[142,143,144],actual_facets=145,pairs=pairs,
       relations=dict(right=55,dotted=374),all26_chain_contacts_dotted=True,
       new_finite_neighbor_lists={str(k):sorted(v)for k,v in neighbors.items()},
       common_center=dict(norm_bounds=saved(hq),all145_pair_bounds=outs(centerpairs),facet_projection_norm_bounds=outs(facetnorms)),
       new_finite_corner=dict(incident=incident,chart=list(map(str,corner5)),norm=str(cq),all145_pair_bounds=outs(cornerpairs),
           center_pair_bounds=saved(dot(C,H)),gram=[list(map(str,row))for row in G20],ldl_pivots=pivots),
       complete_new_finite_vertices=1,known_total_finite_vertices_at_least=171,
       finite_census_proof='Complete finite-neighbor supports of each newwall give exactlyS16+walls141,142,143,144; any supportincludingold95 with141/143/144 is notPD because thatpairis dotted. Alloldfiniteverticespersist. No complete oldambientcensus.',
       new_corner_finite_edges_to_old_vertices=0,
       new_corner_edge_statement='The newcorner containsfour newlyintroduced finitevertex walls141..144; everyoldfinitevertex hasallthese inequalitiesstrict, so nooldfinitevertexshares19incidentwalls. This doesnot classifyideal endpoints.',
       E3_cut_pairs=list(map(str,cuts)),
       residual_escape_rectangle=dict(a_bounds=list(map(str,abox)),b_bounds=list(map(str,bbox)),corners=rect,norm_minimum=str(minrect),proof='Everywallpairing andcenterpairing isaffine ina,b andstrictat4corners; norm isconstant+(a^2+b^2)/2, minimizedat lowerpositiveendpoints.'),
       completed_cusp_preserved=dict(null_ray='sqrt3*n134+2*n0+n1',new_pair_bounds=outs(nullpairs)),
       preserved_three_other_cusp_pairs=oldcusps,
       strict_escape=dict(formula='E4=x-e/100-(49/50)*f-(3/10)*omega17-(1/20)*omega31-eps*sum(other16omegas)',
           epsilon='1/100000',frame_coefficients=list(map(str,ef)),dual_coefficients=list(map(str,neweb)),
           norm=str(eq),all145_pair_bounds=outs(ep),new_cap_pairs=list(map(str,newep)),center_pair_bounds=saved(dot(E4,H))),
       scope='Exact coherent145-wall Coxeter extension withonecomplete newfinite H20cornerS16+A1^4; all145facetsactual, oldfinitevertices and4propercuspspersist. E4strictspacelike impliesinfinitevolume. No full oldcone census or finitevolume/newdimension claim.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert data==json.loads(path.read_text())
    else:path.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['strict_escape']['norm'],flush=True)
