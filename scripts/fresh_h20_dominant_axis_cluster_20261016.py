#!/usr/bin/env python3
"""Exact compatible two-axis cluster after the dominant142 cap.

156 actual facets and four new finite corners; strict E6 survives.
"""
import argparse,json
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
from fresh_h20_coupled_roles_20261014_field import Field,parse,B
from fresh_h20_s61_multicap_20261014_exact import I,FILES,BITS
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_dominant_axis_cluster_20261016'
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
    CLUSTER='results/fresh_h20_residual_h4_corner_20261016.json'
    cluster=json.loads((ROOT/CLUSTER).read_text())
    assert cluster['status']=='PASS_EXACT_145_WALL_RESIDUAL_H4_FINITE_CORNER_CLUSTER'
    for p,h in cluster['sources'].items():assert sha(p)==h,p
    kk=-cq;qs=norms4[0];qw=norms4[3]
    Avec=scale(N[141],-(I(2)/qs).sqrt())
    Bvec=plus(scale(C,(I(2)/kk).sqrt()),scale(N[144],-(I(2)/qw).sqrt()))
    aa=dot(Avec,N[140]);bb=dot(Bvec,N[140]);target=-(I(pq)/2).sqrt()
    assert aa.l>0 and bb.h<0
    zz=(target-aa)/bb;assert zz.l>1
    rr=plus(Avec,scale(Bvec,zz));assert has(dot(rr,rr),2)
    right=S16|{142,143};newpairs=[]
    for i,nn in enumerate(N):
        gg=dot(rr,nn)
        if i in right:assert has(gg);rel='right'
        elif i==140:assert has(gg-target);rel='pi/3'
        elif i==141:assert has(gg+(I(2*qs)).sqrt());rel='parallel'
        else:assert gg.h<0 and gg.h*gg.h>2*norms_old[i],(i,saved(gg));rel='dotted'
        newpairs.append(dict(wall=i,norm=str(norms_old[i]),relation=rel,pair_bounds=saved(gg)))
    assert len(right)==18 and sum(p['relation']=='dotted'for p in newpairs)==125
    oldcenterpairs=centerpairs[:];centerpairs.append(dot(rr,H));assert centerpairs[-1].h<0
    allnorms=norms_old+[F(2)];facetnorms=[hq-g*g/nq for nq,g in zip(allnorms,centerpairs)]
    assert all(x.h<0 for x in facetnorms)
    # Exact new ideal ray from the parallel pair; its complete incidence is checked.
    nu2=plus(scale(C,1/I(kk).sqrt()),scale(N[144],-1/I(qw).sqrt()))
    assert has(dot(nu2,nu2)) and dot(nu2,H).h<0
    idealincident=sorted(S16|{141,142,143,145})
    npairs=[dot(nu2,n)for n in N]+[dot(nu2,rr)]
    for i,gg in enumerate(npairs):
        if i in idealincident:assert has(gg)
        else:assert gg.h<0,(i,saved(gg))
    assert len(idealincident)==20
    ef=list(map(F,cluster['strict_escape']['frame_coefficients']));eb=list(map(F,cluster['strict_escape']['dual_coefficients']))
    E4=plus(lin(frame,[I(x)for x in ef]),lin(omega,[I(x)for x in eb]))
    cut=dot(rr,E4);assert cut.l>0
    eb5=eb[:];eb5[j31]=F(-1,25)
    E5=plus(lin(frame,[I(x)for x in ef]),lin(omega,[I(x)for x in eb5]))
    eq=-ef[0]**2+ef[1]**2+ef[2]**2+sum(eb5[i]*J[i][j]*eb5[j]for i in range(18)for j in range(18))
    assert eq==F(31500061,5000000000)>0 and has(dot(E5,E5),eq)
    ep=[dot(E5,n)for n in N]+[dot(E5,rr)];assert all(x.h<0 for x in ep)
    assert dot(E5,H).h<0
    nupair=dot(rr,nu);assert nupair.h<0
    oldcusps=[]
    for raw in [old['completed_cusp']['point']]+[row['point']for row in old['preserved_old_bounded_cusps']]:
        cv=[bd(A.lift(parse(x,0),0))for x in raw];gg=dot(rr,cv);assert gg.h<0;oldcusps.append(saved(gg))
    # Theonly19 finite oldneighbors include the dottedpair140/142.
    finite=sorted(right|{140});assert len(finite)==19
    assert dot(N[140],N[142]).h<0 and dot(N[140],N[142]).h**2>norms_old[140]*norms_old[142]
    PARTNER='results/fresh_h20_h4_parallel_partner_20261016.json'
    partner=json.loads((ROOT/PARTNER).read_text())
    assert partner['status']=='PASS_EXACT_146_WALL_FINITE_PARALLEL_E4_CUTTER'
    for p,h in partner['sources'].items():assert sha(p)==h,p
    N.append(rr);norms_old.append(F(2));assert len(N)==146 and len(norms_old)==146
    spherical=[list(range(1,11)),list(range(12,16)),[61,87],[142],[143]]
    assert set(sum(spherical,[]))==S16|{142,143}
    newnorms=[F(2),F(2),F(2),norms4[1],norms4[2]]
    unshifted=[lin([N[j]for j in comp],[I(-1)]*len(comp))for comp in spherical]
    oldproper=[nu]
    for raw in [old['completed_cusp']['point']]+[row['point']for row in old['preserved_old_bounded_cusps']]:
        oldproper.append([bd(A.lift(parse(x,0),0))for x in raw])
    other=[j for j in range(146)if j not in idealincident]
    assert len(other)==126
    nu2pairs=[dot(nu2,n)for n in N]
    for j in other:assert nu2pairs[j].h<0
    special=[('center',H),('escapeE5',E5)]+[(f'proper_cusp_{j}',v)for j,v in enumerate(oldproper)]
    for name,v in special:assert dot(nu2,v).h<0,name
    heights=[];constraints=[];added=[]
    for i,(pp,nq)in enumerate(zip(unshifted,newnorms)):
        assert has(dot(pp,pp),nq) and has(dot(pp,nu2))
        bounds=[]
        for j in other:
            bound=(dot(pp,N[j])+(I(nq*norms_old[j])).sqrt())/(-nu2pairs[j])
            bounds.append(dict(kind='dotted_old_wall',index=j,strict_height_lower_bounds=saved(bound)))
        for name,v in special:
            bound=dot(pp,v)/(-dot(nu2,v));bounds.append(dict(kind=name,strict_height_lower_bounds=saved(bound)))
        hh=max(1,1+max(F(row['strict_height_lower_bounds'][1]).__floor__()for row in bounds))
        heights.append(hh);constraints.append(bounds);added.append(plus(pp,scale(nu2,hh)))
    newpairs=[]
    for i,(rnew,nq)in enumerate(zip(added,newnorms)):
        assert has(dot(rnew,rnew),nq)
        comp=spherical[i]
        for j,nn in enumerate(N):
            g=dot(rnew,nn)
            if j in idealincident:
                if len(comp)==1 and j==comp[0]:
                    assert has(g,-nq);rel='parallel'
                elif len(comp)>1 and j in [comp[0],comp[-1]]:
                    assert nq==norms_old[j]==2 and has(g,-1);rel='pi/3'
                else:assert has(g);rel='right'
            else:assert g.h<0 and g.h*g.h>nq*norms_old[j],(146+i,j,saved(g));rel='dotted'
            newpairs.append(dict(walls=[j,146+i],norms=[str(norms_old[j]),str(nq)],relation=rel,pair_bounds=saved(g)))
        for j in range(i):
            g=dot(rnew,added[j]);assert has(g)
            newpairs.append(dict(walls=[146+j,146+i],norms=[str(newnorms[j]),str(nq)],relation='right',pair_bounds=saved(g)))
    assert len(newpairs)==740
    relcounts={name:sum(row['relation']==name for row in newpairs)for name in ['right','pi/3','parallel','dotted']}
    assert relcounts==dict(right=102,**{'pi/3':6},parallel=2,dotted=630)
    FULL=N+added;fullnorms=norms_old+newnorms
    hp=[dot(v,H)for v in FULL];assert len(hp)==151 and all(g.h<0 for g in hp)
    facetnorms=[hq-g*g/nq for nq,g in zip(fullnorms,hp)];assert all(x.h<0 for x in facetnorms)
    allcusp=sorted(idealincident+list(range(146,151)))
    components=[comp+[146+i]for i,comp in enumerate(spherical)]+[[141,145]]
    assert len(allcusp)==25 and set(sum(components,[]))==set(allcusp)
    # Normalize each root to square2 only for the exact diagram certificate.
    expected=[[F(0)]*25 for _ in range(25)];pos={wall:i for i,wall in enumerate(allcusp)}
    for comp in components:
        for i,a in enumerate(comp):
            for j,b in enumerate(comp):
                val=2 if i==j else -2 if len(comp)==2 else -1 if (i-j)%len(comp)in[1,len(comp)-1]else 0
                expected[pos[a]][pos[b]]=F(val)
    for i,a in enumerate(allcusp):
        assert has(dot(nu2,FULL[a]))
        for j,b in enumerate(allcusp):
            actual=2*dot(FULL[a],FULL[b])/(I(fullnorms[a]*fullnorms[b])).sqrt()
            assert has(actual,expected[i][j]),(a,b)
    DD=[row[:]for row in expected];pivs=[]
    for i in range(25):
        pp=DD[i][i];assert pp>=0
        if pp==0:assert all(DD[i][j]==0 for j in range(i,25))
        else:
            for j in range(i+1,25):
                for k in range(i+1,25):DD[j][k]-=DD[j][i]*DD[i][k]/pp
        pivs.append(str(pp))
    assert sum(F(v)>0 for v in pivs)==19 and sum(F(v)==0 for v in pivs)==6
    nullmultipliers=[]
    for i,comp in enumerate(components):
        coeff=[(I(2)/fullnorms[j]).sqrt()for j in comp]
        nullsum=lin([FULL[j]for j in comp],coeff)
        mult=I(heights[i])*(I(2)/newnorms[i]).sqrt()if i<5 else I(2).sqrt()*zz
        assert mult.l>0 and all(has(a-b)for a,b in zip(nullsum,scale(nu2,mult)))
        nullmultipliers.append(saved(mult))
    cuspvals=[dot(nu2,v)for v in FULL]
    for i,g in enumerate(cuspvals):
        if i in allcusp:assert has(g)
        else:assert g.h<0
    preserved=[[dot(v,cp)for v in added]for cp in oldproper]
    assert all(g.h<0 for row in preserved for g in row)
    escapes=[dot(E5,v)for v in FULL];assert all(g.h<0 for g in escapes)
    assert dot(E5,H).h<0
    SECOND='results/fresh_h20_second_cusp_completion_20261016.json'
    second=json.loads((ROOT/SECOND).read_text())
    assert second['status']=='PASS_EXACT_151_WALL_SECOND_PROPER_CUSP_COMPLETION'
    for p,h in second['sources'].items():assert sha(p)==h,p
    U=[FULL[j]for j in [95,142,145]]
    aa=dot(U[0],U[2]);assert aa.h<-2
    den=I(4)-aa*aa;assert den.h<0
    JU=[[I(2)/den,I(0),-aa/den],[I(0),I(1)/norms4[1],I(0)],[-aa/den,I(0),I(2)/den]]
    def projU(w):
        dots=[dot(w,v)for v in U]
        cs=[sum((j*x for j,x in zip(row,dots)),I(0))for row in JU]
        return lin(U,cs)
    vv=plus(FULL[150],scale(projU(FULL[150]),-1));vnorm=dot(vv,vv);assert vnorm.l>0
    raw=omega[j31];kv=plus(raw,scale(projU(raw),-1));kv=plus(kv,scale(vv,-dot(kv,vv)/vnorm))
    knorm=dot(kv,kv);assert knorm.l>0
    cosine=(I(2)+I(2).sqrt()).sqrt()/2
    target=-(I(2*norms4[2])).sqrt()*cosine
    pp=scale(vv,target/vnorm);pnorm=target*target/vnorm
    lam=((2-pnorm)/knorm).sqrt();assert lam.l>0
    trials=[];accepted=[]
    for sg in [-1,1]:
        cand=plus(pp,scale(kv,sg*lam));pairs=[];bad=[]
        for j,nn in enumerate(FULL):
            gg=dot(cand,nn)
            if j in S16|{95,142,145}:good=has(gg);rel='right'
            elif j==150:good=has(gg-target);rel='pi/8'
            else:good=gg.h<0 and gg.h*gg.h>2*fullnorms[j];rel='dotted'
            if not good:bad.append(j)
            pairs.append(dict(wall=j,old_norm=str(fullnorms[j]),relation=rel,pair_bounds=saved(gg)))
        hh=dot(cand,H);cut=dot(cand,E5)
        trials.append(dict(sign=sg,failed_gates=bad,center_pair_bounds=saved(hh),E5_pair_bounds=saved(cut)))
        if not bad and hh.h<0 and cut.l>0:accepted.append((sg,cand,pairs,hh,cut))
    assert len(accepted)==1
    sg,r151,ppairs,hnew,cut=accepted[0];assert has(dot(r151,r151),2)
    counts={rel:sum(x['relation']==rel for x in ppairs)for rel in ['right','pi/8','dotted']}
    assert counts=={'right':19,'pi/8':1,'dotted':131}
    FULL2=FULL+[r151];norms2=fullnorms+[F(2)]
    hp2=[dot(H,v)for v in FULL2];assert all(g.h<0 for g in hp2)
    facets2=[hq-g*g/nq for g,nq in zip(hp2,norms2)];assert all(g.h<0 for g in facets2)
    cusp2=[dot(r151,cp)for cp in oldproper+[nu2]];assert all(g.h<0 for g in cusp2)
    # Unique new finite corner: S16 plus142,145,150,new.
    cvec=H[:]
    for j in S16:cvec=plus(cvec,scale(omega[S.index(j)],-dot(H,FULL[j])))
    cvec=plus(cvec,scale(FULL[142],-dot(H,FULL[142])/norms4[1]))
    cvec=plus(cvec,scale(FULL[145],-dot(H,FULL[145])/2))
    dc=I(2*norms4[2])-target*target;assert dc.l>0
    h150=dot(H,FULL[150]);hr=dot(H,r151)
    cvec=plus(cvec,scale(FULL[150],-(2*h150-target*hr)/dc))
    cvec=plus(cvec,scale(r151,-(norms4[2]*hr-target*h150)/dc))
    cnorm=dot(cvec,cvec);assert cnorm.h<0 and dot(cvec,H).h<0
    cinc=sorted(S16|{142,145,150,151});assert len(cinc)==20
    cpairs=[dot(cvec,v)for v in FULL2]
    for j,gg in enumerate(cpairs):
        if j in cinc:assert has(gg)
        else:assert gg.h<0,(j,saved(gg))
    neighbors=sorted(S16|{95,142,145,150});assert len(neighbors)==20
    rejects=[]
    for omit in neighbors:
        remaining=set(neighbors)-{omit}
        if omit==95:continue
        if 95 in remaining and 145 in remaining:witness=[95,145]
        elif 95 in remaining and 150 in remaining:witness=[95,150]
        else:raise AssertionError(omit)
        x,y=witness;gg=dot(FULL[x],FULL[y]);assert gg.h<0 and gg.h*gg.h>fullnorms[x]*fullnorms[y]
        rejects.append(dict(omitted=omit,dotted_pair=witness,pair_bounds=saved(gg)))
    assert len(rejects)==19
    # A conservative rational escape, to be refined only within this fresh source.
    eb6=eb5[:];eb6[j17]=F(-7,25);eb6[j31]=F(-1,20)
    E6=plus(lin(frame,[I(x)for x in ef]),lin(omega,[I(x)for x in eb6]))
    eq6=-ef[0]**2+ef[1]**2+ef[2]**2+sum(eb6[i]*J[i][j]*eb6[j]for i in range(18)for j in range(18))
    assert eq6>0 and has(dot(E6,E6),eq6)
    esc6=[dot(E6,v)for v in FULL2];assert all(g.h<0 for g in esc6),[(j,saved(x))for j,x in enumerate(esc6)if x.h>=0]
    assert dot(E6,H).h<0
    PI8='results/fresh_h20_second_cusp_pi8_cap_20261016.json'
    pi8=json.loads((ROOT/PI8).read_text())
    assert pi8['status']=='PASS_EXACT_ALTERNATIVE_152_WALL_PI8_CUSP_EDGE_CAP'
    for p,h in pi8['sources'].items():assert sha(p)==h,p
    k8=-cnorm;assert k8.l>0
    r152=lin([FULL[142],FULL[145],cvec],[-I(1)/I(norms4[1]).sqrt(),-I(3).sqrt()/2,I(1)/(2*k8).sqrt()])
    assert has(dot(r152,r152),2)
    target142=-I(norms4[1]).sqrt();target145=-I(3).sqrt()
    npairs=[]
    for j,nn in enumerate(FULL2):
        gg=dot(r152,nn)
        if j in S16|{150,151}:assert has(gg);rel='right'
        elif j==142:assert has(gg-target142);rel='pi/4'
        elif j==145:assert has(gg-target145);rel='pi/6'
        else:assert gg.h<0 and gg.h*gg.h>2*norms2[j],(j,saved(gg));rel='dotted'
        npairs.append(dict(wall=j,old_norm=str(norms2[j]),relation=rel,pair_bounds=saved(gg)))
    rc={rel:sum(x['relation']==rel for x in npairs)for rel in ['right','pi/4','pi/6','dotted']}
    assert rc=={'right':18,'pi/4':1,'pi/6':1,'dotted':132}
    FULL3=FULL2+[r152];norms3=norms2+[F(2)]
    hp3=[dot(H,v)for v in FULL3];assert all(g.h<0 for g in hp3)
    f3=[hq-g*g/nq for g,nq in zip(hp3,norms3)];assert all(v.h<0 for v in f3)
    cusps3=[dot(r152,cp)for cp in oldproper+[nu2]];assert all(g.h<0 for g in cusps3)
    # H18 is the projection of H perpendicular to S16 and roots150/151.
    h142=dot(H,FULL[142]);h145=dot(H,FULL[145]);h152=dot(H,r152)
    H18=lin([cvec,FULL[142],FULL[145]],[I(1),h142/norms4[1],h145/2])
    assert dot(H18,H18).h<0
    face_support=sorted(S16|{150,151});assert len(face_support)==18
    for j in face_support:assert has(dot(H18,FULL3[j]))
    corners=[]
    for oldid,nq,targ in [(142,norms4[1],target142),(145,F(2),target145)]:
        ho=dot(H,FULL3[oldid]);det=2*nq-targ*targ;assert det.l>0
        cc=lin([H18,FULL3[oldid],r152],[I(1),-(2*ho-targ*h152)/det,-(nq*h152-targ*ho)/det])
        ci=sorted(set(face_support)|{oldid,152});assert len(ci)==20
        qq=dot(cc,cc);assert qq.h<0 and dot(cc,H).h<0
        allpairs=[dot(cc,nn)for nn in FULL3]
        for j,gg in enumerate(allpairs):
            if j in ci:assert has(gg)
            else:assert gg.h<0,(oldid,j,saved(gg))
        corners.append(dict(side_pair=[oldid,152],incident=ci,norm_bounds=saved(qq),point_bounds=outs(cc),
           all153_pair_bounds=outs(allpairs),center_pair_bounds=saved(dot(cc,H)),two_root_determinant_bounds=saved(det)))
    neigh=sorted(S16|{142,145,150,151});assert len(neigh)==20
    rejects=[]
    for omit in neigh:
        if omit in [142,145]:continue
        assert {142,145}<=set(neigh)-{omit}
        rejects.append(dict(omitted=omit,indefinite_subdiagram=[142,145,152],normalized_Gram_determinant='-2'))
    assert len(rejects)==18
    esc7=[dot(E6,nn)for nn in FULL3];assert all(gg.h<0 for gg in esc7)
    TRIANGLE='results/fresh_h20_pi8_triangle_face_20261016.json'
    triangle=json.loads((ROOT/TRIANGLE).read_text())
    assert triangle['status']=='PASS_EXACT_153_WALL_COMPACT_246_TRIANGLE_FACE'
    for p,h in triangle['sources'].items():assert sha(p)==h,p
    qu=norms4[1];gamma=(I(1)+I(5).sqrt())/4
    ccboost=(qu*gamma+5*(25-qu+qu*gamma*gamma).sqrt())/(25-qu)
    ssboost=I(qu).sqrt()*(ccboost+gamma)/5
    assert ccboost.l>1 and ssboost.l>0 and has(ccboost*ccboost-ssboost*ssboost,1)
    assert has((25-qu)*ccboost*ccboost-2*qu*gamma*ccboost-(25+qu*gamma*gamma))
    r153=lin([FULL3[142],C],[-(I(2)/qu).sqrt()*ccboost,(I(2)/kk).sqrt()*ssboost])
    assert has(dot(r153,r153),2)
    target149=-(I(2*qu)).sqrt()*gamma
    right=S16|{141,143,144};np4=[]
    for j,nn in enumerate(FULL3):
        gg=dot(r153,nn)
        if j in right:assert has(gg);rel='right'
        elif j==149:assert has(gg-target149);rel='pi/5'
        else:assert gg.h<0 and gg.h*gg.h>2*norms3[j],(j,saved(gg));rel='dotted'
        np4.append(dict(wall=j,old_norm=str(norms3[j]),relation=rel,pair_bounds=saved(gg)))
    rel4={rel:sum(p['relation']==rel for p in np4)for rel in ['right','pi/5','dotted']}
    assert rel4=={'right':19,'pi/5':1,'dotted':133}
    FULL4=FULL3+[r153];norms4all=norms3+[F(2)]
    hp4=[dot(H,v)for v in FULL4];assert all(gg.h<0 for gg in hp4)
    fac4=[hq-gg*gg/nq for gg,nq in zip(hp4,norms4all)];assert all(gg.h<0 for gg in fac4)
    cp4=[dot(r153,cp)for cp in oldproper+[nu2]];assert all(gg.h<0 for gg in cp4)
    # H18 is the old center projection perpendicular to S16+141+143.
    H18b=lin([C,FULL3[142],FULL3[144]],[-dot(H,C)/kk,dot(H,FULL3[142])/qu,dot(H,FULL3[144])/qw])
    support18=sorted(S16|{141,143});assert len(support18)==18
    assert dot(H18b,H18b).h<0
    for j in support18:assert has(dot(H18b,FULL4[j]))
    corners4=[];hr=dot(H,r153)
    for oldid,nq,targ in [(144,qw,I(0)),(149,qu,target149)]:
        ho=dot(H,FULL4[oldid]);det=2*nq-targ*targ;assert det.l>0
        point=lin([H18b,FULL4[oldid],r153],[I(1),-(2*ho-targ*hr)/det,-(nq*hr-targ*ho)/det])
        incident=sorted(set(support18)|{oldid,153});assert len(incident)==20
        pn=dot(point,point);ph=dot(point,H);assert pn.h<0 and ph.h<0
        pp4=[dot(point,nn)for nn in FULL4]
        for j,gg in enumerate(pp4):
            if j in incident:assert has(gg)
            else:assert gg.h<0,(oldid,j,saved(gg))
        corners4.append(dict(side_pair=[oldid,153],incident=incident,norm_bounds=saved(pn),point_bounds=outs(point),all154_pair_bounds=outs(pp4),center_pair_bounds=saved(ph),two_root_determinant_bounds=saved(det)))
    neighbors4=sorted(right|{149});assert len(neighbors4)==20
    dotw=dot(FULL3[144],FULL3[149]);assert dotw.h<0 and dotw.h*dotw.h>qw*qu
    rejects4=[]
    for omit in neighbors4:
        if omit in [144,149]:continue
        assert {144,149}<=set(neighbors4)-{omit}
        rejects4.append(dict(omitted=omit,dotted_pair=[144,149],pair_bounds=saved(dotw)))
    assert len(rejects4)==18
    esc8=[dot(E6,nn)for nn in FULL4];assert all(gg.h<0 for gg in esc8)
    aaxis=ssboost/ccboost;assert 0<aaxis.l<aaxis.h<1
    DOMINANT='results/fresh_h20_dominant142_cap_20261016.json'
    dominant=json.loads((ROOT/DOMINANT).read_text())
    assert dominant['status']=='PASS_EXACT_154_WALL_DOMINANT_AXIS_PI5_CAP'
    for p,h in dominant['sources'].items():assert sha(p)==h,p
    sc=I(1)/(2*ssboost);cc=(1+sc*sc).sqrt();assert sc.l>1
    assert has(sc*ssboost,F(1,2))
    axisroots=[lin([FULL4[axis],C],[-(I(2)/norms3[axis]).sqrt()*cc,(I(2)/kk).sqrt()*sc])for axis in [141,143]]
    assert has(dot(axisroots[0],axisroots[1])+2*sc*sc) and dot(axisroots[0],axisroots[1]).h<-2
    FULL6=FULL4[:];norms6=norms4all[:];newpairs6=[];newcorners6=[];newdefinitions=[];newrejects=[]
    for rootindex,(axis,root) in enumerate(zip([141,143],axisroots),154):
        assert has(dot(root,root),2)
        right6=S16|({141,142,143,144}-{axis});assert len(right6)==19
        for j,nn in enumerate(FULL6):
            gg=dot(root,nn)
            if j in right6:assert has(gg);rel='right'
            elif j==153:assert has(gg,-1);rel='pi/3'
            else:assert gg.h<0 and gg.h*gg.h>2*norms6[j],(rootindex,j,saved(gg));rel='dotted'
            newpairs6.append(dict(walls=[j,rootindex],old_norm=str(norms6[j]),relation=rel,pair_bounds=saved(gg)))
        FULL6.append(root);norms6.append(F(2))
        # Project the common center to the rank18 support, leaving C and the two relevant old axes.
        sup18=sorted(S16|({141,142,143,144}-{axis,142}));assert len(sup18)==18
        ha=dot(H,FULL6[axis]);hu=dot(H,FULL6[142]);hroot=dot(H,root)
        h18=lin([C,FULL6[axis],FULL6[142]],[-dot(H,C)/kk,ha/norms6[axis],hu/qu])
        assert dot(h18,h18).h<0
        for j in sup18:assert has(dot(h18,FULL6[j]))
        cn=[]
        for oldid,nq,targ in [(142,qu,I(0)),(153,F(2),I(-1))]:
            ho=dot(H,FULL6[oldid]);det=2*nq-targ*targ;assert det.l>0
            point=lin([h18,FULL6[oldid],root],[I(1),-(2*ho-targ*hroot)/det,-(nq*hroot-targ*ho)/det])
            inc=sorted(set(sup18)|{oldid,rootindex});assert len(inc)==20
            pn=dot(point,point);ph=dot(point,H);assert pn.h<0 and ph.h<0
            pp=[dot(point,n)for n in FULL4+axisroots]
            for j,gg in enumerate(pp):
                if j in inc:assert has(gg)
                else:assert gg.h<0,(rootindex,oldid,j,saved(gg))
            cn.append(dict(side_pair=[oldid,rootindex],incident=inc,norm_bounds=saved(pn),point_bounds=outs(point),all156_pair_bounds=outs(pp),center_pair_bounds=saved(ph),two_root_determinant_bounds=saved(det)))
        newcorners6.extend(cn)
        neighbors6=sorted(right6|{153});assert len(neighbors6)==20
        dotw=dot(FULL6[142],FULL6[153]);assert dotw.h<0 and dotw.h*dotw.h>2*qu
        rejects=[]
        for omit in neighbors6:
            if omit in [142,153]:continue
            assert {142,153}<=set(neighbors6)-{omit}
            rejects.append(dict(omitted=omit,dotted_pair=[142,153],pair_bounds=saved(dotw)))
        assert len(rejects)==18
        newrejects.append(dict(new_wall=rootindex,finite_neighbors=neighbors6,rejected_supports=rejects))
        newdefinitions.append(dict(new_wall=rootindex,axis=axis,recipe='sqrt2*(-c_star*n_axis/sqrt(q_axis)+s_star*C/sqrt(-C^2))',norm='2',physical_normal_bounds=outs(root)))
    rc6={rel:sum(v['relation']==rel for v in newpairs6)for rel in ['right','pi/3','dotted']}
    assert len(newpairs6)==309 and rc6=={'right':38,'pi/3':2,'dotted':269}
    hp6=[dot(H,v)for v in FULL6];assert all(gg.h<0 for gg in hp6)
    fac6=[hq-gg*gg/nq for gg,nq in zip(hp6,norms6)];assert all(gg.h<0 for gg in fac6)
    cusp6=[[dot(root,cp)for cp in oldproper+[nu2]]for root in axisroots];assert all(gg.h<0 for row in cusp6 for gg in row)
    esc9=[dot(E6,n)for n in FULL6];assert all(gg.h<0 for gg in esc9)
    cutoff=sc/cc;assert 0<cutoff.l<cutoff.h<1
    src=[BASE,ARC,CAP,REVIEW,FRONT,PREV,PREV_REVIEW,MIX,CLUSTER,PARTNER,SECOND,PI8,TRIANGLE,DOMINANT,
         'results/fresh_h20_affine_a17_completion_20261016.json',
         'results/fresh_h20_affine_a17_completion_20261016_review.json',f'scripts/{TAG}.py']+FILES+[
         'scripts/fresh_h20_coupled_roles_20261014_field.py','scripts/fresh_h20_s61_multicap_20261014_exact.py']
    return dict(status='PASS_EXACT_156_WALL_COORDINATED_DOMINANT_AXIS_CLUSTER',sources={p:sha(p)for p in src},
       definition=dict(s_star='1/(2*s153)',c_star='sqrt(1+s_star^2)',s_star_bounds=saved(sc),c_star_bounds=saved(cc),
          exact_norm_identity='c_star^2-s_star^2=1',exact_pi3_identity='r154.r153=r155.r153=-2*s_star*s153=-1',
          exact_mutual_identity='r154.r155=-2*s_star^2<-2',mutual_pair_bounds=saved(dot(axisroots[0],axisroots[1])),
          new_roots=newdefinitions,residual_chart_inequalities=['z141<=s_star/c_star','z143<=s_star/c_star'],axis_cutoff_bounds=saved(cutoff)),
       retained_walls=154,new_walls=[154,155],actual_facets=156,new_pairs=newpairs6,relations=rc6,
       common_center=dict(norm_bounds=saved(hq),all156_pair_bounds=outs(hp6),facet_projection_norm_bounds=outs(fac6)),
       preserved_five_proper_cusp_new_pairs=[outs(row)for row in cusp6],new_finite_corners=newcorners6,
       complete_new_finite_vertices=4,known_total_finite_vertices_at_least=180,finite_neighbor_censuses=newrejects,
       strict_escape=dict(dependency=PI8,formula=pi8['strict_escape']['formula'],norm=str(eq6),all156_pair_bounds=outs(esc9),center_pair_bounds=saved(dot(E6,H))),
       scope='Exactcompatibletwo-capcluster,156actualfacets,4complete newfinitevertices,5propercuspsandcompact246trianglepreserved. Cutsbothremaininglargeaxes. StrictE6survives,soinfinitevolume;nofullambientfinite/exteriorcensusorvolumeclosureclaim.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();data=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert data==json.loads(path.read_text())
    else:path.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['relations'],flush=True)
