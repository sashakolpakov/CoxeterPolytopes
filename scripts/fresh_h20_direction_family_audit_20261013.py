#!/usr/bin/env python3
"""Independent exact census and all-height obstruction for old vertex directions.

Only standard-library Fraction arithmetic and the earlier independent
Qsqrt6 quadratic-extension checker are used. No producer/numeric imports.
"""
import argparse
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts import fresh_h20_targeted_escape_20261013_second_independent as E

ROOT=E.ROOT
BASE=E.BASE
SHORT='data/enumerations/fresh_h20_targeted_escape_20261013_short_vertices.json'
LONG='data/enumerations/fresh_h20_targeted_escape_20261013_long_vertices.json'
BANK='data/enumerations/fresh_h20_targeted_escape_20261013_vertex_bank.json'
GEOM='data/enumerations/fresh_h20_targeted_escape_20261013_bank_geometry.json'
THIRD='data/enumerations/fresh_h20_targeted_escape_20261013_third.json'
OUT='results/fresh_h20_direction_family_audit_20261013.json'
ZERO=(Q(0),Q(0));ONE=(Q(1),Q(0))
def qa(x,y):return(x[0]+y[0],x[1]+y[1])
def qn(x):return(-x[0],-x[1])
def qm(x,y):return(x[0]*y[0]+10*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def qd(x,y):
    den=y[0]**2-10*y[1]**2;assert den
    z=qm(x,(y[0],-y[1]));return(z[0]/den,z[1]/den)
def qs(x):
    a,b=x;sa=(a>0)-(a<0);sb=(b>0)-(b<0)
    if not sb:return sa
    if not sa or sa==sb:return sb
    d=a*a-10*b*b;return sa*((d>0)-(d<0))
def qr(x):return(Q(x),Q(0))
def qsum(xs):
    out=ZERO
    for x in xs:out=qa(out,x)
    return out
def qdot(a,b):return qsum(qm(x,y)for x,y in zip(a,b,strict=True))
def dot(a,b):return sum((x*y for x,y in zip(a,b,strict=True)),Q(0))
def mv(M,x):return [dot(row,x)for row in M]
def enc(x):return list(map(str,x))
def solve(M,b):
    A=[list(row)+[v]for row,v in zip(M,b,strict=True)];n=len(A);assert all(len(row)==n+1 for row in A)
    for j in range(n):
        i=next(i for i in range(j,n)if A[i][j]!=ZERO);A[i],A[j]=A[j],A[i]
        p=A[j][j];A[j]=[qd(x,p)for x in A[j]]
        for i in range(n):
            if i==j:continue
            t=A[i][j]
            if t!=ZERO:A[i]=[qa(x,qn(qm(t,y)))for x,y in zip(A[i],A[j],strict=True)]
    return [row[-1]for row in A]

def census():
    b=E.read(BASE);short=E.read(SHORT);long=E.read(LONG)
    assert short['source_sha256']==long['source_sha256']==E.sha(BASE)
    assert all(Q(x[1])==0 for row in b['form']for x in row)
    M=[[Q(x[0])for x in row[:20]]for row in b['form'][:20]]
    assert len(b['form'])==21 and all(Q(b['form'][20][i][0])==Q(b['form'][i][20][0])==0 for i in range(20))
    assert Q(b['form'][20][20][0])==3 and b['strict_center'][-1]==['0','0']
    R=[[Q(x[0])for x in r[:20]]for r in b['old_normals']]
    assert len(R)==36 and all(Q(x[1])==0 for row in b['old_normals']for x in row)
    assert b['normals'][:24]==b['old_normals'][:24]and all(row[-1]==['0','0']for row in b['old_normals'])
    cov=[mv(M,r)for r in R];G=[[dot(a,c)for c in cov]for a in R]
    assert all(G[i][i]==(2 if i<24 else 6)for i in range(36))
    assert all(G[i][j]in(0,-1)for i in range(24)for j in range(i))
    assert all(G[i][j]<=-6 for i in range(24,36)for j in range(24,i))
    assert all(G[i][j]in(0,-3)for i in range(24,36)for j in range(24))
    adj=[{j for j in range(24)if i!=j and G[i][j]==-1}for i in range(24)]
    @lru_cache(None)
    def pd(ids):
        # A cycle's all-one vector has square<=0 in its principal Gram.
        # A forest is decided by exact leaf Schur elimination.
        S=set(ids);nbr={i:adj[i]&S for i in ids};seen=set()
        for root in ids:
            if root in seen:continue
            stack=[(root,-1)]
            while stack:
                i,parent=stack.pop()
                if i in seen:return False
                seen.add(i);stack.extend((j,i)for j in nbr[i]if j!=parent)
        diag={i:Q(2)for i in ids};left=set(ids)
        while left:
            i=min(j for j in left if len(nbr[j])<=1)
            if diag[i]<=0:return False
            for j in nbr[i]:diag[j]-=1/diag[i];nbr[j].remove(i)
            left.remove(i)
        return True
    accepted={ids for ids in combinations(range(24),19)if pd(ids)}
    expected={tuple(r['incident_short_walls'])for r in short['vertices']}
    assert len(expected)==len(short['vertices'])==2046 and accepted==expected
    counts=Counter();longaccepted=set();candidates=0
    for L in range(24,36):
        N=[j for j in range(24)if G[L][j]==-3];O=[j for j in range(24)if G[L][j]==0];assert len(N)==3
        # Two long neighbors have negative3x3 Gram determinant. One neighbor
        # joined to another selected short wall gives zero/negative determinant.
        for p in(0,-1):assert 24-36-6*p*p+18*p<0
        assert 6*(4-1)-9*2==0
        possibilities=[tuple(s)for s in combinations(O,18)]
        for a in N:possibilities.extend(tuple(sorted((a,)+s))for s in combinations([j for j in O if G[a][j]==0],17))
        assert len(possibilities)==len(set(possibilities))
        for S in possibilities:
            candidates+=1
            if pd(S):longaccepted.add(S+(L,));counts[L]+=1
    assert candidates==22116 and len(longaccepted)==6300
    assert longaccepted=={tuple(r['incident_walls'])for r in long['vertices']}
    assert len(long['vertices'])==6300 and set(counts.values())=={525}
    h=[Q(x[0])for x in b['strict_center'][:20]];hc=mv(M,h)
    for row in short['vertices']+long['vertices']:
        ids=row.get('incident_short_walls',row.get('incident_walls'));x=list(map(Q,row['projection']));coeff=list(map(Q,row['gram_solution']))
        assert all(x[k]==h[k]-sum((a*R[i][k]for a,i in zip(coeff,ids,strict=True)),Q(0))for k in range(20))
        ps=[dot(x,c)for c in cov]
        assert all(p==0 if i in ids else p<0 for i,p in enumerate(ps))
        q=dot(x,mv(M,x));assert q==Q(row['norm'])<0 and dot(x,hc)<0
    # For every long vertex, its incident norm6 wall was split into opposite
    # height coordinates 2epsilon and -epsilon. Acuteness forces t=0, leaving
    # the negative norm of its timelike old projection.
    for L in range(24,36):
        r,s=b['normals'][24+2*(L-24):26+2*(L-24)]
        assert all(Q(r[i][0])==Q(s[i][0])==R[L][i]and Q(r[i][1])==Q(s[i][1])==0 for i in range(20))
        z,zp=Q(r[-1][0]),Q(s[-1][0]);assert z==-2*zp and z*zp<0
    return b,M,R,G,short,long,{'short_subsets':42504,'short_directions':2046,'long_reduced_subsets':candidates,'long_directions':6300,
        'all_finite_directions':8346,'long_directions_no_spacelike_acute_height':6300,'pd_method':'Cycle null-vector rejection; exact rational leaf Schur elimination on forests. Long incidence reduces to an isolated I2(6) block plus the short forest.'}

def target_and_bounds(b,M,short):
    bank=E.read(BANK);entry=bank['caps'][916]
    assert len(bank['caps'])==1593 and(entry['vertex_index'],entry['anchor'],entry['label'])==(1400,26,4)
    fullM=[[Q(x[0])for x in row]for row in b['form']]
    R=[[Q(x[0])for x in r]for r in b['normals'][:48]];assert all(Q(x[1])==0 for r in b['normals'][:48]for x in r)
    cov=[mv(fullM,r)for r in R]
    rows=short['vertices'];x=[Q(v)/55 for v in rows[1400]['projection']]+[Q(0)];assert dot(x,mv(fullM,x))==-2
    t=(Q(-4),Q(1));anchor=list(map(qr,x));anchor[-1]=t
    active=[0,2,4,5,7,9,10,11,12,13,14,15,16,17,19,20,21,36,42]
    equations=[list(map(qr,cov[i]))for i in active]
    equations.append([qsum(qm(qr(a),z)for a,z in zip(row,anchor,strict=True))for row in fullM])
    equations.append([qr(int(i==20))for i in range(21)])
    w=solve(equations,[ZERO]*20+[qr(-1)]);assert all(qdot(row,w)==ZERO for row in equations[:20])
    proposed=E.read(GEOM)['whole_bank_geometry']['point'];scale=-Q(str(proposed[-1]));assert scale>0
    for exact,numeric in zip(w,proposed,strict=True):
        delta=qa(exact,qn(qr(Q(str(numeric))/scale)))
        assert qs(qa(delta,qr(Q(1,10**8))))>0 and qs(qa(qr(Q(1,10**8)),qn(delta)))>0
    wc=[qsum(qm(qr(a),z)for a,z in zip(row,w,strict=True))for row in fullM]
    norm=qdot(w,wc);assert qs(norm)>0
    bounds=[];possible=[]
    for i,row in enumerate(rows):
        v=list(map(Q,row['projection']))+[Q(0)];a=[dot(v,c)for c in cov];z=[c[-1]for c in cov]
        assert all(a[j]<=0 for j in range(48)if z[j]==0)
        lo=max(-aa/zz for aa,zz in zip(a,z,strict=True)if zz<0);hi=min(-aa/zz for aa,zz in zip(a,z,strict=True)if zz>0)
        assert lo<=hi
        intercept=qsum(qm(qr(xx),yy)for xx,yy in zip(v,wc,strict=True));maximum=qa(intercept,qr(-3*lo))
        if qs(maximum)>0:possible.append(i)
        bounds.append({'direction':i,'height_min':str(lo),'height_max':str(hi),'max_target_pairing_under_acuteness':enc(maximum)})
    assert possible==[1400,1511]
    return fullM,w,norm,bounds,{'active_old_walls':active,'extra_active_mirror':{'direction':1400,'height':enc(t),'bank_entry':916},'point':[enc(x)for x in w],'norm':enc(norm),'potential_cut_directions':possible}

def exceptions(b,M,w,short):
    cov=[[sum((Q(a[0])*Q(z[0])for a,z in zip(row,r,strict=True)),Q(0))for row in b['form']]for r in b['normals'][:48]]
    wc=[qsum(qm(qr(a),z)for a,z in zip(row,w,strict=True))for row in M]
    x=[Q(v)/55 for v in short['vertices'][1400]['projection']]+[Q(0)]
    xx=dot(x,mv(M,x));assert xx==-2 and dot(x,cov[26])==-6 and cov[26][-1]==-6 and dot(bare(b['normals'][26]),cov[26])==18
    intercept=qsum(qm(qr(a),z)for a,z in zip(x,wc,strict=True));assert intercept==qm(qr(3),(Q(-4),Q(1)))
    # A cutting t lies in[-1,sqrt10-4). The wall26 ratio
    # F=2(1+t)^2/(-2+3t²) is strictly increasing there: derivative has
    # the sign of -(t+1)(3t+2), positive for -1<t<-sqrt(2/3).
    # Thus F is in[0,1/2). Coxeter choices are0 and1/4 only.
    assert qs(qa((Q(-4),Q(1)),qr(1)))>0
    t0=(Q(-4),Q(1));q0=qa(qr(-2),qm(qr(3),qm(t0,t0)))
    assert qs(q0)>0 and qs(qa(t0,qr(Q(2,3))))<0
    assert qm(qr(4),qm(qa(ONE,t0),qa(ONE,t0)))==q0
    # Right height t=-1: wall49 has ratio(5-2sqrt6)/6 in(0,1/4).
    f=E.Field((Q(14),Q(0)));zero=f.rational(0)
    r49=tuple(map(f.base,b['normals'][49]));fm=tuple(tuple(f.base(z)for z in row)for row in b['form']);c49=f.mv(fm,r49)
    A49=f.dot(tuple(f.rational(v)for v in x),c49);Z49=c49[-1];N49=f.dot(r49,c49)
    assert A49==f.base((Q(-2),Q(-1,2)))and Z49==f.rational(-3)and N49==f.rational(3)
    def qratio_delta(t,c):
        p=E.add(A49,f.mul(Z49,t));q=E.add(f.rational(-2),f.mul(f.rational(3),f.mul(t,t)))
        return E.sub(f.mul(p,p),f.mul(f.rational(c),f.mul(N49,q)))
    tr=f.rational(-1);assert f.sign(qratio_delta(tr,0))>0 and f.sign(qratio_delta(tr,Q(1,4)))<0
    # pi3 height solves5t²+16t+10=0; the second root is<-1.
    t3=(Q(-8,5),Q(0),Q(1,5),Q(0));assert f.sign(E.add(t3,f.rational(1)))>0
    assert f.sign(E.add(f.rational(-2),f.mul(f.rational(3),f.mul(t3,t3))))>0
    assert E.add(E.add(f.mul(f.rational(5),f.mul(t3,t3)),f.mul(f.rational(16),t3)),f.rational(10))==zero
    assert f.sign(qratio_delta(t3,Q(5,6)))>0 and f.sign(qratio_delta(t3,Q(17,20)))<0
    # cos²(pi7) is the largest root of64X³-80X²+24X-1.
    # P(5/6)>0 and P'>0 on[5/6,1], so it is<5/6;
    # cos²(pi8)=(2+sqrt2)/4>17/20 since sqrt2>7/5.
    assert 64*Q(5,6)**3-80*Q(5,6)**2+24*Q(5,6)-1>0
    assert 192*Q(5,6)**2-160*Q(5,6)+24>0 and 384*Q(5,6)-160>0 and Q(7,5)**2<2
    # The second direction needs only a coarser rational target bound.
    y=[Q(v)/29 for v in short['vertices'][1511]['projection']]+[Q(0)]
    assert dot(y,mv(M,y))==-6 and dot(y,cov[26])==-9
    int2=qsum(qm(qr(a),z)for a,z in zip(y,wc,strict=True));cut2=qd(int2,qr(3))
    assert qs(qa(cut2,qr(Q(3,2))))>0 and qs(qa(cut2,qr(Q(29,20))))<0
    # For-3/2<t<-29/20, F26=2(t+3/2)^2/(-6+3t²) lies(0,1/4).
    # A safe bound is numerator<1/200, denominator>123/400.
    assert Q(1,200)/Q(123,400)<Q(1,4)
    assert dot(y,cov[3])==-1 and cov[3][-1]==0
    assert Q(1)/(2*(-6+3*Q(3,2)**2))==Q(2,3)
    return {'direction1400':{'primitive_projection':list(map(str,x)),'norm':-2,'acute_lower_height':'-1','cut_upper_height':enc((Q(-4),Q(1))),
        'wall26_squared_cosine':'2(1+t)^2/(-2+3t^2)','possible_cutting_coxeter_heights':['-1','(-8+sqrt14)/5'],
        'right_height_wall49_squared_cosine':'(5-2sqrt6)/6 in(0,1/4)',
        'pi3_height_wall49_squared_cosine_bounds':['5/6','17/20'],
        'pi3_height_exclusion':'cos²(pi7)<5/6<ratio<17/20<cos²(pi8)'},
        'direction1511':{'primitive_projection':list(map(str,y)),'norm':-6,'acute_lower_height':'-3/2','exact_cut_upper_height':enc(cut2),
        'rational_cut_upper_bound':'-29/20','wall26_squared_cosine':'2(t+3/2)^2/(-6+3t^2)',
        'only_possible_cutting_coxeter_height':'-3/2','short3_squared_cosine':'2/3',
        'exclusion':'4cos²(pi/m) is an algebraic integer; the rational value8/3 is not.'}}

def bare(row):
    assert all(Q(x[1])==0 for x in row)
    return list(Q(x[0])for x in row)

def strict_witness(b,w):
    f=E.Field((Q(10),Q(0)));zero=f.rational(0)
    M=tuple(tuple(f.base(x)for x in row)for row in b['form']);h=tuple(map(f.base,b['strict_center']));hc=f.mv(M,h)
    assert f.sign(f.dot(h,hc))<0
    target=tuple((a,Q(0),b,Q(0))for a,b in w);tc=f.mv(M,target)
    old=[tuple(map(f.base,r))for r in b['normals']]+[tuple(map(f.base,r['vector']))for r in E.read(E.POLAR)['cuts']]
    assert len(old)==82
    extras=E.read(E.CAPS)['caps']+[{'field':E.read(E.FIRST)['field'],'vector':E.read(E.FIRST)['normal']},{'field':E.read(E.SOURCE)['field'],'vector':E.read(E.SOURCE)['normal']}]
    def checks(v,strict):
        vc=f.mv(M,v);ps=[]
        for r in old:
            p=f.dot(r,vc);assert f.sign(p)<0 if strict else f.sign(p)<=0;ps.append(enc(p))
        for row in extras:
            d=f.base(row['field']['radicand']);rr=[(f.base(x[0]),f.base(x[1]))for x in row['vector']]
            p=(f.dot(tuple(x[0]for x in rr),vc),f.dot(tuple(x[1]for x in rr),vc))
            assert f.extsign(p,d)<0 if strict else f.extsign(p,d)<=0;ps.append([enc(x)for x in p])
        return ps
    closed=checks(target,False);assert f.sign(f.dot(target,tc))>0 and f.sign(f.dot(target,hc))<0
    checks(h,True);epsilon=Q(1,1048576);v=tuple(E.add(x,f.mul(f.rational(epsilon),y))for x,y in zip(target,h,strict=True));vc=f.mv(M,v)
    assert f.sign(f.dot(v,vc))>0 and f.sign(f.dot(v,hc))<0
    ps=checks(v,True)
    return {'unperturbed_old87_pairings':closed,'perturbation':str(epsilon),'strict_point':[enc(x)for x in v],
            'field_basis':['1','sqrt6','sqrt10','sqrt60'],'norm':enc(f.dot(v,vc)),'center_pairing':enc(f.dot(v,hc)),
            'old87_strict_pairings':ps,'universal_direction_family_strictness':'Every allowed cap has nonpositive pairing with the unperturbed target and strictly negative pairing with the old center. The positive rational perturbation makes all such inequalities strict simultaneously.'}

def third_comparison(b,strict):
    data=E.read(THIRD);f=E.Field((Q(10),Q(0)));zero=f.rational(0)
    assert data['field']['radicand']==['63','0'];d=f.base(data['field']['radicand'])
    M=tuple(tuple(f.base(x)for x in row)for row in b['form'])
    r=[(f.base(x[0]),f.base(x[1]))for x in data['normal']]
    rc=tuple(zip(f.mv(M,tuple(x[0]for x in r)),f.mv(M,tuple(x[1]for x in r)),strict=True))
    norm=(zero,zero)
    for a,c in zip(r,rc,strict=True):
        p=f.extmul(a,c,d);norm=tuple(E.add(x,y)for x,y in zip(norm,p,strict=True))
    assert norm==(f.rational(2),zero)
    w=tuple(tuple(map(Q,x))for x in strict['strict_point']);wc=f.mv(M,w)
    pairing=(f.dot(tuple(x[0]for x in r),wc),f.dot(tuple(x[1]for x in r),wc))
    assert f.extsign(pairing,d)>0
    old_pairs=[];zeros=[]
    for i,old in enumerate(b['old_normals']):
        c=f.mv(M,tuple(map(f.base,old)))
        p=(f.dot(tuple(x[0]for x in r),c),f.dot(tuple(x[1]for x in r),c));s=f.extsign(p,d)
        assert s<=0
        if s==0:zeros.append(i)
        old_pairs.append({'old':i,'pairing':[enc(x)for x in p],'sign':s})
    assert zeros==data['short_orthogonal_support']and len(zeros)==18 and all(i<24 for i in zeros)
    return {'source':THIRD,'source_sha256':E.sha(THIRD),'norm':[enc(x)for x in norm],
            'pairing_with_universal_strict_escape':[enc(x)for x in pairing],'pairing_sign':1,
            'field_basis':'Coefficients in(1,sqrt6,sqrt10,sqrt60), then extension by positive sqrt63.',
            'old_spatial_projection_zero_walls':zeros,'all36_original_wall_pairings':old_pairs,
            'outside_direction_family':'The old spatial component annihilates exactly18short walls and no long wall, whereas each old finite-vertex direction annihilates19old walls. The new cap therefore leaves the excluded family and cuts its exact universal witness.',
            'scope':'Independent norm,36old-spatial incidences and target cut. The new cap all87 Coxeter relation certificate is an existing separate dependency, not repeated here.'}

def run():
    b,M,R,G,short,long,counts=census();fullM,w,q,bounds,target=target_and_bounds(b,M,short)
    cases=exceptions(b,fullM,w,short);strict=strict_witness(b,w);comparison=third_comparison(b,strict)
    deps=[BASE,SHORT,LONG,E.POLAR,E.CAPS,E.FIRST,E.SOURCE,BANK,GEOM,THIRD]
    return {'status':'PASS_ALL_REAL_HEIGHTS_FOR_8346_OLD_FINITE_VERTEX_DIRECTIONS','source_sha256':{p:E.sha(p)for p in deps},
        'census':counts,'target':target,'all2046_acute_height_bounds':bounds,'two_exceptional_directions':cases,'universal_strict_escape':strict,'new88th_cap_comparison':comparison,
        'family':'Positive multiples of(x,t), where x is an oldH19 finite-vertex ray and t is arbitrary real. Require positive norm and Coxeter compatibility with the retained87 walls. All finite labels and all dotted or parallel relations are allowed.',
        'scope':'No finite-volume polytope completion retaining the87-wall chamber can use finitely many caps drawn only from this direction family: one strict future spacelike escape survives every allowable real-height cap. No uniform-neighborhood claim about an infinite cap intersection is needed. This does not exclude18-short directions, changed old spatial directions, or replacement of retained walls. The numerical1593-bank and continuous interval scouts are motivation only; their sampled coverage is not used in the theorem.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
    if a.verify:assert E.read(OUT)==out
    else:(ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['census'])
