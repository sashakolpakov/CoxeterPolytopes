#!/usr/bin/env python3
"""Independent rational enclosure proof, no project/numerical imports.

Uses exact midpoint inverses and a Neumann-series error bound, rather than
the producer's interval Gauss-Jordan inverse. No normal/cone reconstruction.
"""
from fractions import Fraction as Q
from hashlib import sha256
from itertools import product
import json
from math import isqrt
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_compact10_cross_bridge_high_interval_20261018'
D=10**28


def box(x):return tuple(map(Q,x)) if isinstance(x,(list,tuple)) else (Q(x),Q(x))
def out(l,h):
    assert l<=h
    return Q((l.numerator*D)//l.denominator,D),Q(-((-h.numerator*D)//h.denominator),D)
def plus(x,y):return out(x[0]+y[0],x[1]+y[1])
def minus(x,y):return out(x[0]-y[1],x[1]-y[0])
def times(x,y):
    z=[a*b for a in x for b in y];return out(min(z),max(z))
def quotient(x,y):
    assert y[0]*y[1]>0
    return times(x,(1/y[1],1/y[0]))
def sigma(xs):
    z=box(0)
    for x in xs:z=plus(z,x)
    return z
def root(x):
    assert x[0]>=0
    a=isqrt(x[0].numerator*D*D//x[0].denominator)
    b=isqrt(x[1].numerator*D*D//x[1].denominator)
    if Q(b,D)**2<x[1]:b+=1
    assert Q(a,D)**2<=x[0] and Q(b,D)**2>=x[1]
    return Q(a,D),Q(b,D)
def bilinear(x,g,y):return sigma(times(times(a,g[i][j]),b) for i,a in enumerate(x) for j,b in enumerate(y))


def midpoint_inverse_enclosure(g):
    n=len(g)
    # Intersect symmetric entries, since the actual Gram is symmetric.
    g=[row[:] for row in g]
    for i in range(n):
        for j in range(i):
            g[i][j]=g[j][i]=(max(g[i][j][0],g[j][i][0]),min(g[i][j][1],g[j][i][1]))
            assert g[i][j][0]<=g[i][j][1]
    m=[[(x[0]+x[1])/2 for x in row] for row in g]
    a=[row[:]+[Q(i==j) for j in range(n)] for i,row in enumerate(m)];pivots=[]
    for i in range(n):
        p=a[i][i];assert p;pivots.append(p)
        a[i]=[x/p for x in a[i]]
        for j in range(n):
            if j==i:continue
            t=a[j][i];a[j]=[x-t*y for x,y in zip(a[j],a[i])]
    c=[row[n:] for row in a]
    assert all(sum(m[i][k]*c[k][j] for k in range(n))==Q(i==j) for i in range(n) for j in range(n))
    cnorm=max(sum(abs(x) for x in row) for row in c)
    error=max(sum(max(abs(g[i][j][0]-m[i][j]),abs(g[i][j][1]-m[i][j])) for j in range(n)) for i in range(n))
    eta=cnorm*error;assert 0<=eta<1
    delta=cnorm*eta/(1-eta)
    # ||(M+E)^(-1)-M^(-1)||∞ <= ||M^(-1)||∞ eta/(1-eta).
    # This also prevents an eigenvalue crossing along M+tE,0<=t<=1.
    return [[out(x-delta,x+delta) for x in row] for row in c],pivots


def check():
    saved=json.loads((ROOT/f'results/{TAG}.json').read_text())
    assert saved['status']=='PASS_EXACT_INTERVAL_HIGH_BRIDGE_EXCLUSION'
    for p,h in saved['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    primary=json.loads((ROOT/'results/fresh_nonarithmetic_release49_62_replacement_20261017.json').read_text())
    prior=json.loads((ROOT/'results/fresh_compact10_cross_bridge_20261018.json').read_text())
    white=[0,1,2,3,5,6];ids=white+[8,10,33,7]
    center=list(map(box,primary['common_center']));metric=[list(map(box,row)) for row in primary['computational_gram']]
    normals={}
    for r in primary['facet_witnesses']:
        w=r['wall']
        if w in ids:
            cp=box(primary['center_products'][str(w)]);assert cp[1]<0
            normals[w]=[quotient(times(box(2),minus(x,y)),cp) for x,y in zip(center,map(box,r['point']))]
    assert set(normals)==set(ids)
    gram=[[bilinear(normals[i],metric,normals[j]) for j in ids] for i in ids]
    assert all(gram[i][i][0]<=2<=gram[i][i][1] for i in range(10))
    wi,wp=midpoint_inverse_enclosure([r[:6] for r in gram[:6]]);assert all(x>0 for x in wp)
    residual=[[minus(gram[i][j],bilinear(gram[i][:6],wi,[gram[k][j] for k in range(6)])) for j in range(6,10)] for i in range(6,10)]
    gi,gp=midpoint_inverse_enclosure([r[:3] for r in residual[:3]])
    assert gp[0]>0 and gp[1]>0 and gp[2]<0
    c=[residual[i][3] for i in range(3)]
    v=[sigma(times(gi[i][j],c[j]) for j in range(3)) for i in range(3)]
    s7=minus(residual[3][3],bilinear(c,gi,c));assert s7[0]>0
    original=saved['projection'];assert original['walls']==ids and original['residual_walls']==[8,10,33,7]
    for own,row in zip(residual,original['residual_gram']):
        for a,b in zip(own,map(box,row)):assert max(a[0],b[0])<=min(a[1],b[1])
    assert max(s7[0],box(original['seven_residual_square'])[0])<=min(s7[1],box(original['seven_residual_square'])[1])
    sqrt2,sqrt3,sqrt5=[root(box(n)) for n in (2,3,5)]
    phi=quotient(plus(box(1),sqrt5),box(2));theta={2:box(0),3:box(1),4:sqrt2,5:phi,6:(sqrt3[0],Q(2))}
    profiles=[[box(0)]*6]
    for pos in (4,5):
        for m in (3,4,5):
            rho=[box(0)]*6;rho[pos]=minus(box(0),theta[m]);profiles.append(rho)
    high=[(0,[a,2,b],[i for i,m in enumerate((a,2,b)) if m==6]) for a,b in product(range(2,7),repeat=2) if max(a,b)==6]
    assert json.loads(json.dumps(high))==saved['symbolic_high_families']
    assert len(saved['continuous_boxes'])==9 and len(saved['low_point_domains'])==42
    low=[(r['profile'],r['labels']) for r in prior['bank']['recipes']]
    assert [(r['profile'],r['labels']) for r in saved['low_point_domains']]==low
    margins=[];counts={'positive':0,'native_gap':0}
    for typ,domains in [('high',high),('low',[(p,m,[]) for p,m in low])]:
        for i,(pi,labels,slots) in enumerate(domains):
            rr=saved['continuous_boxes'][i] if typ=='high' else saved['low_point_domains'][i]
            assert rr['index']==i and rr['profile']==pi
            assert rr['symbolic_labels' if typ=='high' else 'labels']==labels
            rho=profiles[pi]
            offsets=[bilinear(rho,wi,[gram[k][6+j] for k in range(6)]) for j in range(4)]
            square=minus(box(2),bilinear(rho,wi,rho))
            rhs=[minus(minus(box(0),theta[m]),offsets[j]) for j,m in enumerate(labels)]
            rem=minus(square,bilinear(rhs,gi,rhs));assert rem[0]>0
            linear=plus(sigma(times(rhs[j],v[j]) for j in range(3)),offsets[3])
            radical=root(times(rem,s7));values=[minus(linear,radical),plus(linear,radical)]
            for j,value in enumerate(values):
                if value[0]>0:
                    margins.append(value[0]);counts['positive']+=1;kind='positive';gap=None
                else:
                    ok=[]
                    for m in (2,3,4,5):
                        lower=minus(box(0),theta[m+1] if m<5 else sqrt3)
                        upper=minus(box(0),theta[m])
                        left,right=value[0]-lower[1],upper[0]-value[1]
                        if min(left,right)>0:ok.append((m,left,right))
                    assert len(ok)==1,(typ,i,j,value)
                    m,left,right=ok[0];margins.extend([left,right]);counts['native_gap']+=1;kind='native_gap';gap=[m,m+1]
                if typ=='high':
                    assert kind==('native_gap' if j==0 else 'positive')
                    if j==0:assert gap==rr['excluded_between_native_labels']
                    old=box(rr['minus_native7' if j==0 else 'plus_native7'])
                else:
                    br=rr['branches'][j];assert br['sign']==(-1 if j==0 else 1) and br['exclusion']['kind']==kind
                    if gap:assert br['exclusion']['labels']==gap
                    old=box(br['native7'])
                assert max(value[0],old[0])<=min(value[1],old[1])
    assert counts==dict(positive=54,native_gap=48) and min(margins)>0
    assert saved['covered_high_families']==9 and saved['covered_high_norm_signs']==18
    assert saved['covered_low_recipes']==42 and saved['covered_low_norm_signs']==84
    assert saved['total_symbolic_domains']==51 and saved['total_norm_sign_domains']==102 and saved['two_parameter_boxes']==1
    assert saved['new_normals']==saved['new_cone_checks']==0 and saved['prior_control_open_ends']==10
    assert saved['exact_deletion_debit']==prior['exact_deletion_debit'] and saved['deleted_walls']==[30,44]
    supports=primary['retained_old_finite_supports']+[r['support'] for r in primary['new_finite_vertices']]
    assert len(supports)==199
    lost=[v for v in supports if set(v)&{30,44}]
    assert lost==saved['exact_deletion_debit']['lost_finite_supports'] and len(lost)==16 and sum(8 in v for v in lost)==2
    assert saved['exact_deletion_debit']['compact_faces_by_dimension']=={'2':12,'3':1}
    print('PASS_INDEPENDENT_BRIDGE_INTERVAL_PROOF',dict(high_boxes=9,low_points=42,sign_domains=102,
        positive=counts['positive'],native_gaps=counts['native_gap'],minimum_separation=str(min(margins)),
        deletion=[16,12,1],new_normals=0,new_cones=0),flush=True)


if __name__=='__main__':check()
