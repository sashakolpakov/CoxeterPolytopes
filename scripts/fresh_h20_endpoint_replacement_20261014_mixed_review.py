#!/usr/bin/env python3
"""Separate FLINT rational-matrix and rational-interval mixed-orbit audit.

No producer module is imported.  The inherited 395-wall subsystem's mutual
pair certificate is retained; all orbit transfers used here are rechecked.
"""
import argparse,json,math
from collections import Counter
from hashlib import sha256
from itertools import product
from pathlib import Path
from flint import fmpq as Q, fmpq_mat as M

ROOT=Path(__file__).resolve().parents[1]
BASE='fresh_h20_endpoint_replacement_20261014_mixed'
TAG=BASE+'_review'
PREC=10**85

def sha(p): return sha256((ROOT/p).read_bytes()).hexdigest()
def read(p): return json.loads((ROOT/p).read_text())
def sign(a,b=Q(0)):
    # Sign of a+b sqrt(7), independent of the producer's pair arithmetic.
    if b==0: return int(a>0)-int(a<0)
    if a==0: return int(b>0)-int(b<0)
    if a>0 and b>0:return 1
    if a<0 and b<0:return -1
    return (int(a*a>7*b*b)-int(a*a<7*b*b))*(1 if a>0 else -1)
def p(x):return(Q(x),Q(x))
def ia(x,y):return(x[0]+y[0],x[1]+y[1])
def im(x,y):
    a=[s*t for s in x for t in y];return min(a),max(a)
def minus(x,y):return(x[0]-y[1],x[1]-y[0])
def sqroot(x):
    out=[]
    for i,a in enumerate(x):
        assert a>=0
        n=math.isqrt(int(a.numerator)*PREC**2//int(a.denominator))
        lo,hi=Q(n,PREC),Q(n+1,PREC)
        assert lo*lo<=a<hi*hi
        out.append((lo,hi)[i])
    return tuple(out)
S7=sqroot(p(7))
def kval(a,b):return ia(p(a),im(p(b),S7))
def mat(rows):return M([[Q(x)for x in row]for row in rows])
def dot(x,y,G):return (M([x])*G*M([y]).transpose())[0,0]
def enc(x):
    if isinstance(x,Q):return str(x)
    if isinstance(x,dict):return {str(k):enc(v)for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [enc(y)for y in x]
    return x

def run():
    source='scripts/'+BASE+'.py';cert='results/'+BASE+'.json'
    oldpath='results/fresh_h20_balanced_integer_exact_20261014.json'
    d=read(cert);old=read(oldpath)
    for a in(d,old):
        for f,h in a['sources'].items():assert sha(f)==h,f
    G=M(19,19)
    blocks=[(0,6),(6,6),(12,7)]
    for start,n in blocks:
        for i in range(n):
            G[start+i,start+i]=1
            if i+1<n:G[start+i,start+i+1]=G[start+i+1,start+i]=Q(-1,2)
    # A direct LDL calculation checks signature(20,1) after the hyperbolic plane.
    gg=M(G);pivots=[]
    for i in range(19):
        pivot=gg[i,i];assert pivot>0;pivots.append(pivot)
        for j in range(i+1,19):
            for k in range(i+1,19):gg[j,k]-=gg[j,i]*gg[i,k]/pivot
    roots=mat(old['walls']);seed=[[roots[i,j]for j in range(19)]for i in range(22)]
    assert all(roots[i,19]==1 and roots[i,20]==0 for i in range(22))
    factors=[];wallstart=0
    for start,n in blocks:
        B=M([[G[start+i,start+j]for j in range(n)]for i in range(n)])
        inv=B.inv();vv=[]
        expected=[[-Q(1)]*n]+[[Q(i==j)for j in range(n)]for i in range(n)]
        assert [s[start:start+n]for s in seed[wallstart:wallstart+n+1]]==expected
        for omitted in range(n+1):
            # Solve all retained affine-facet equations, not the producer's
            # closed-form A_n inverse formula.
            rhs=M([[Q(1)-(n+1 if j+1==omitted else 0)]for j in range(n)])
            v=inv*rhs;vv.append([Q(0)]*start+[v[j,0]for j in range(n)]+[Q(0)]*(19-start-n))
        factors.append(vv);wallstart+=n+1
    local=list(product(range(7),range(7),range(8)))
    Z=[[sum(factors[f][a[f]][i]for f in range(3))for i in range(19)]for a in local]
    assert all(dot(z,z,G)==196 for z in Z)
    # Full cyclic Gram invariance, including each factor's origin; hence
    # simultaneous shifts transfer all blended-center distance tests.
    cyclic_tests=0
    for vv in factors:
        F=M(vv);K=F*G*F.transpose();n=len(vv)
        for shift in range(n):
            for i in range(n):
                for j in range(n):
                    assert K[i,j]==K[(i+shift)%n,(j+shift)%n];cyclic_tests+=1
    C0=[];C1=[]
    for a,saved in zip(local,d['centers']):
        base=Z[local.index(a)]
        delta=[sum(factors[f][(a[f]+1)%7][i]-factors[f][a[f]][i]for f in (0,1))for i in range(19)]
        c0=[base[i]+2*delta[i]for i in range(19)]
        c1=[Q(-4,7)*delta[i]for i in range(19)]
        assert [[str(x),str(y)]for x,y in zip(c0,c1)]==saved
        C0.append(c0);C1.append(c1)
    A,B=M(C0),M(C1);ZZ=M(Z);SS=M(seed)
    AA=A*G*A.transpose()+7*B*G*B.transpose()
    AB=A*G*B.transpose()+B*G*A.transpose()
    N0,N1=Q(916),Q(-288);Q0,Q1=Q(308),Q(-112)
    assert all(AA[i,i]==N0 and AB[i,i]==N1 for i in range(392))
    assert sign(Q0,Q1)>0 and sign(14,-4)>0
    assert (14**2+7*(-4)**2,2*14*(-4))==(Q0,Q1)
    assert (Q0-28*14+84,Q1-28*(-4))==(0,0)
    assert d['radius']==['14','-4'] and d['radius_square']==[str(Q0),str(Q1)]
    assert d['center_square']==[str(N0),str(N1)]
    # Separate endpoint expression rho/R=sqrt(c8^2+17/7)-c8.
    c8=sqroot(im(p(Q(1,4)),ia(p(2),sqroot(p(2)))))
    c82=im(c8,c8);mu=minus(sqroot(ia(c82,p(Q(17,7)))),c8)
    ri=kval(14,-4);rho=im(ri,mu);rho2=im(rho,rho)
    assert rho[0]>0 and rho2[1]<21
    assert Q(d['endpoint_radius_interval'][0])<=rho[0]<=rho[1]<=Q(d['endpoint_radius_interval'][1])
    # Polynomial identity in s,c, with s^2=c^2+17/7:
    # (s-c)^2+2c(s-c)+1-24/7=s^2-c^2-17/7=0.
    import sympy as sp
    s,c=sp.symbols('s c')
    identity=sp.expand((s-c)**2+2*c*(s-c)+1-sp.Rational(24,7))
    assert sp.rem(identity,s*s-c*c-sp.Rational(17,7),s)==0
    assert im(ia(rho,p(8)),ia(rho,p(8)))[1]<196
    assert kval(N0,N1)[0]>im(ia(ri,p(8)),ia(ri,p(8)))[1]
    # Literal seed contacts for every mixed cap, independent of the cyclic proof.
    SG0=A*G*SS.transpose();SG1=B*G*SS.transpose();seedcounts=Counter()
    for i in range(392):
        cnt=Counter()
        for j in range(22):
            u,v=1-SG0[i,j],-SG1[i,j]
            if u==0 and v==0:kind='right'
            elif u==14 and v==-4:kind='parallel'
            else:assert sign(u-14,v+4)>0;kind='dotted'
            cnt[kind]+=1
        assert cnt=={'right':17,'parallel':2,'dotted':3};seedcounts.update(cnt)
    # Complete mixed-mixed pair matrix (not just reference distances).
    paircounts=Counter();neighbors=[]
    for i in range(392):
        for j in range(i):
            dist0,dist1=2*(N0-AA[i,j]),2*(N1-AB[i,j])
            if (dist0,dist1)==(3*Q0,3*Q1):kind='pi3'
            else:assert sign(dist0-4*Q0,dist1-4*Q1)>0;kind='dotted'
            paircounts[kind]+=1
    assert dict(paircounts)==d['mixed_mutual_pairs']
    for j in range(1,392):
        if (2*(N0-AA[0,j]),2*(N1-AB[0,j]))==(3*Q0,3*Q1):neighbors.append(j)
    assert neighbors==[8,48,56,336]
    AZ=A*G*ZZ.transpose();BZ=B*G*ZZ.transpose()
    limit=im(ia(ri,rho),ia(ri,rho));de=[]
    for j in range(392):
        a,b=N0+196-2*AZ[0,j],N1-2*BZ[0,j];de.append((a,b))
        if j==0:assert (a,b)==(Q(24,7)*Q0,Q(24,7)*Q1)
        else:assert kval(a,b)[0]>limit[1]
        assert sign(a-Q0,b-Q1)>0
    assert d['mixed_endpoint_pairs']=={'pi8':392,'dotted':392*391}
    # The five nonseed finite neighbors are mutually dotted.
    for i in neighbors:
        assert kval(N0+196-2*AZ[i,0],N1-2*BZ[i,0])[0]>limit[1]
        for j in neighbors:
            if i!=j:assert sign(2*(N0-AA[i,j])-4*Q0,2*(N1-AB[i,j])-4*Q1)>0
    assert d['finite_neighbor_star']['maximum_finite_clique_size']==19
    # Inherited middle/motif balls are strictly inside B8; both new orbits
    # are outside B8.  This gives every omitted cross-pair strictly dotted.
    oldballs=[]
    for i in range(415,810):
        w=[roots[i,j]for j in range(19)];nn=dot(w,w,G);q=nn-2*roots[i,19]
        assert nn<36 and 0<q<=4 and -200-roots[i,19]<0
        oldballs.append((w,q))
    assert len(oldballs)==395
    # H=(0,200,1), norm -400, is strictly inside all1202 walls.  Acute
    # pairings make each orthogonal facet projection timelike and strict.
    assert sign(N0-Q0+400,N1-Q1)>0
    assert 196-rho2[1]+400>0
    assert len(seed)+1+len(oldballs)+784==d['coherent_wall_count']==1202
    # Exact preserved endpoint: P_z=(z,(196+rho^2)/2,1).
    # Nineteen seed zeros are independent.  All mixed caps are strict by
    # de>Q; other endpoint caps are strict since their centers differ.
    endpoint_seed_count=0
    for a,z in zip(local,Z):
        gaps=[1-dot(z,w,G)for w in seed]
        retained=[i for i,g in enumerate(gaps)if g==0]
        assert len(retained)==19 and min(gaps)==0
        K=M([seed[i]for i in retained])*G*M([seed[i]for i in retained]).transpose()
        assert K.det()>0
        endpoint_seed_count+=1
    assert min(dot([u-v for u,v in zip(Z[0],z)],[u-v for u,v in zip(Z[0],z)],G)for z in Z[1:])==84
    assert 84>4*rho2[1]
    # Two representative tangent cusps, with full incidence tested against
    # both new orbits; old balls are excluded by the B8 separation.
    cusps=[]
    for seed_id in(1,8):
        w=seed[seed_id]
        y0=[C0[0][i]+14*w[i]for i in range(19)]
        y1=[C1[0][i]-4*w[i]for i in range(19)]
        yy0=dot(y0,y0,G)+7*dot(y1,y1,G);yy1=2*dot(y0,y1,G)
        incident=[]
        for j,n in enumerate(seed):
            a,b=1-dot(y0,n,G),-dot(y1,n,G)
            assert sign(a,b)>=0
            if a==0 and b==0:incident.append(j)
        assert len(incident)==17 and seed_id in incident
        other=[j for j in incident if j!=seed_id]
        assert all(dot(seed[seed_id],seed[j],G)==0 for j in other)
        O=M([seed[j]for j in other]);assert (O*G*O.transpose()).rank()==16
        for i in range(392):
            d0=[y0[j]-C0[i][j]for j in range(19)];d1=[y1[j]-C1[i][j]for j in range(19)]
            power0=dot(d0,d0,G)+7*dot(d1,d1,G)-Q0
            power1=2*dot(d0,d1,G)-Q1
            assert (power0==power1==0) if i==0 else sign(power0,power1)>0
            e0=[y0[j]-Z[i][j]for j in range(19)]
            assert kval(dot(e0,e0,G)+7*dot(y1,y1,G),2*dot(e0,y1,G))[0]>rho2[1]
        cusps.append(dict(seed=seed_id,other_seed_incidence=other,quotient_rank=17))
    # Full finite old rays are all distinct under the cyclic group; the
    # two parallel facets have different tangent points on each mixed ball.
    tangent_keys=set()
    for i,a in enumerate(local):
        for sid in((a[0]+1)%7,7+(a[1]+1)%7):
            w=seed[sid];key=tuple((C0[i][j]+14*w[j],C1[i][j]-4*w[j])for j in range(19))
            tangent_keys.add(key)
    assert len(tangent_keys)==784
    # Independent strict escape replay, all walls (and radial relaxation).
    saved=d['strict_escape'];x=list(map(Q,saved['spatial']));D=Q(saved['D'])
    assert x==[Q(99,200)*(u+v)for u,v in zip(Z[0],Z[1])]
    xx=dot(x,x,G);assert xx-2*D==1 and -200-D<0
    assert min(1-dot(x,s,G)for s in seed)==Q(1,100)
    mixpowers=[]
    for i in range(392):
        d0=[x[j]-C0[i][j]for j in range(19)]
        a=dot(d0,d0,G)+7*dot(C1[i],C1[i],G)-Q0
        b=-2*dot(d0,C1[i],G)-Q1
        assert sign(a-1,b)>0;mixpowers.append((a,b))
    mp=min(mixpowers,key=lambda ab:float(ab[0])+math.sqrt(7)*float(ab[1]))
    assert all(sign(a-mp[0],b-mp[1])>=0 for a,b in mixpowers)
    assert [str(u)for u in mp]==saved['mixed_power']
    ed=min(dot([u-v for u,v in zip(x,z)],[u-v for u,v in zip(x,z)],G)for z in Z)
    assert ed-rho2[1]>1
    op=min(dot([u-v for u,v in zip(x,c)],[u-v for u,v in zip(x,c)],G)-q for c,q in oldballs)
    rp=min(xx-4-(dot(x,z,G)-4)**2/192 for z in Z)
    assert op==Q(saved['old_middle_power']) and rp==Q(saved['all_real_radial_power'])
    assert min(op,rp,xx-4)>1
    oldx=[Q(99,400)*sum(Z[j][i]for j in(0,8,56,64))for i in range(19)]
    ox=[oldx[j]-C0[0][j]for j in range(19)]
    oldpower=(dot(ox,ox,G)+7*dot(C1[0],C1[0],G)-Q0,
              -2*dot(ox,C1[0],G)-Q1)
    assert oldpower==(Q(2270077,5000),Q(-176)) and sign(*oldpower)<0
    # A7 adjacent centers differ by square112, hence the extra R/8
    # shift has square7Q/4 and gives the positive pair Q/8.
    assert dot([u-v for u,v in zip(Z[0],Z[1])],
               [u-v for u,v in zip(Z[0],Z[1])],G)==112
    assert (2-Q(7,4))/2==Q(1,8)
    paths=[source,cert,oldpath,'scripts/'+TAG+'.py']
    return enc(dict(status='PASS_SEPARATE_RATIONAL_MATRIX_MIXED_GEOMETRY_REVIEW',
        sources={p:sha(p)for p in paths},signature=[20,1],spatial_positive_pivots=pivots,
        coherent_facets=1202,cyclic_gram_identities=cyclic_tests,
        literal_mixed_centers=392,literal_mixed_seed_pairs=dict(seedcounts),
        complete_mixed_mutual_pairs=dict(paircounts),reference_mixed_endpoint_pairs=392,
        transferred_mixed_endpoint_pairs=392**2,endpoint_radius_interval=rho,
        mixed_endpoint_pi8_identity='(rho/R)^2+2cos(pi/8)*(rho/R)=17/7',
        actual_facets_by_strict_center_projection=1202,preserved_finite_endpoints=endpoint_seed_count,
        maximum_finite_clique_containing_mixed_cap=19,new_finite_vertices_involving_mixed_cap=0,
        representative_forced_cusps=cusps,distinct_transported_forced_cusps=784,
        strict_escape_norm=1,strict_escape_center_pair=-200-D,
        mixed_escape_power=mp,old_middle_escape_power=op,all_real_radial_power=rp,
        previous_square_escape_ball_power=oldpower,
        scope='Separate FLINT rational-matrix arithmetic and independently rounded rational radical bounds; no producer imports. Full mixed-mixed and seed pair matrices are checked. Mixed-endpoint pairs use independently verified cyclic Gram transfer. The 395 inherited middle/motif mutual pair certificate is inherited, not re-enumerated. Finite-neighbor argument excludes new finite vertices involving mixed caps; this is not a complete vertex/cusp census. Exact strict escape proves infinite volume.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=run()
    path=ROOT/'results'/f'{TAG}.json'
    if a.verify:assert read(str(path.relative_to(ROOT)))==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['coherent_facets'],flush=True)
