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
TAG='fresh_compact10_three_chain_gate_20261018'
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
    assert saved['status']=='PASS_THREE_CHAIN_FIRST_CAP_EXACT_GATE'
    for p,h in saved['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    primary=json.loads((ROOT/'results/fresh_nonarithmetic_release49_62_replacement_20261017.json').read_text())
    prior=json.loads((ROOT/'results/fresh_moving33_compact10_bridge_20261018.json').read_text())
    white=[0,1,2,3,5,6];sides=[4,7,8,10,44,63,65];ids=white+sides
    h=list(map(box,primary['common_center']));metric=[list(map(box,r)) for r in primary['computational_gram']]
    normals={}
    for r in primary['facet_witnesses']:
        w=r['wall']
        if w in ids:
            cp=box(primary['center_products'][str(w)]);assert cp[1]<0
            normals[w]=[quotient(times(box(2),minus(x,y)),cp) for x,y in zip(h,map(box,r['point']))]
    assert set(normals)==set(ids)
    gram=[[bilinear(normals[a],metric,normals[b]) for b in ids] for a in ids]
    wi,wp=midpoint_inverse_enclosure([r[:6] for r in gram[:6]]);assert all(x>0 for x in wp)
    cov={a:gram[ids.index(a)][:6] for a in sides}
    res={a:{b:minus(gram[ids.index(a)][ids.index(b)],bilinear(cov[a],wi,cov[b])) for b in sides} for a in sides}
    assert saved['projection']['walls']==ids
    for a in sides:
        for b in sides:
            own=res[a][b];old=box(saved['projection']['residual'][str(a)][str(b)])
            assert max(own[0],old[0])<=min(own[1],old[1])
    for a,expected in {8:[0]*6,10:[0,0,0,0,-1,0],44:[0]*6,63:[0]*6,65:[0]*6}.items():
        assert all(x[0]<=v<=x[1] for x,v in zip(cov[a],expected))
    for a,b in ((8,10),(8,65),(63,65),(8,44),(10,44),(63,44)):
        x=gram[ids.index(a)][ids.index(b)];assert x[0]<=0<=x[1]

    def make_frame(fs,gate):
        inverse,pivots=midpoint_inverse_enclosure([[res[a][b] for b in fs] for a in fs])
        assert pivots[0]>0 and pivots[1]>0 and pivots[2]<0
        c=[res[a][gate] for a in fs]
        v=[sigma(times(inverse[i][j],c[j]) for j in range(3)) for i in range(3)]
        sq=minus(res[gate][gate],bilinear(c,inverse,c));assert sq[0]>0
        return fs,gate,inverse,v,sq

    def evaluate(fr,wc,ms,continuous=False):
        fs,gate,inv,v,sq=fr
        offs=[bilinear(wc,wi,cov[a]) for a in fs+[gate]]
        rhs=[minus(minus(box(0),x),o) for x,o in zip(ms,offs)]
        rem=minus(minus(box(2),bilinear(wc,wi,wc)),bilinear(rhs,inv,rhs));assert rem[1]>=0
        if continuous:rem=(max(Q(),rem[0]),rem[1])
        else:assert rem[0]>0
        lin=plus(sigma(times(x,y) for x,y in zip(rhs,v)),offs[-1]);rad=root(times(rem,sq))
        return [minus(lin,rad),plus(lin,rad)]

    fr4=make_frame([8,10,65],4);fr63=make_frame([8,10,65],63)
    cross=minus(res[4][63],bilinear([res[a][4] for a in (8,10,65)],fr4[2],[res[a][63] for a in (8,10,65)]))
    assert cross[1]<0 and box(saved['gate4_63_residual_cross'])[1]<0
    theta={2:box(0),3:box(1),4:root(box(2)),5:quotient(plus(box(1),root(box(5))),box(2))}
    profiles=[[box(0)]*6]
    for pos in (4,5):
        for m in (3,4,5):
            w=[box(0)]*6;w[pos]=minus(box(0),theta[m]);profiles.append(w)
    rows=saved['low_recipes'];assert len(rows)==46
    domain=[(r['profile'],tuple(r['labels'])) for r in rows];assert len(set(domain))==46
    assert all(0<=p<7 and len(ms)==3 and set(ms)<={2,3,4,5} for p,ms in domain)
    count={'positive':0,'native_gap':0}
    for i,r in enumerate(rows):
        assert r['index']==i and len(r['exclusions'])==2
        p,ms=domain[i];vs={4:evaluate(fr4,profiles[p],[theta[m] for m in ms]),63:evaluate(fr63,profiles[p],[theta[m] for m in ms])}
        for sign,proof in enumerate(r['exclusions']):
            gate=proof['gate'];si=sign if gate==63 else 1-sign
            assert gate in vs and proof['canonical63_sign']==sign and proof['gate_sign']==si
            value=vs[gate][si];old=box(r['evaluations'][str(gate)]['branches'][si])
            assert max(old[0],value[0])<=min(old[1],value[1])
            if proof['kind']=='positive':assert value[0]>0
            else:assert proof['kind']=='native_gap' and proof['labels']==[2,3] and -1<value[0]<=value[1]<0
            count[proof['kind']]+=1
    assert count==saved['low_branch_decisions']
    high=[[0,[2,m,6],[2]] for m in range(2,6)]+[[0,[6,2,2],[0]]]
    assert high==saved['symbolic_high_families']
    values=evaluate(fr63,profiles[0],[(root(box(3))[0],Q(2)),box(0),box(0)],True)
    assert all(v[0]>0 for v in values)
    for value,old in zip(values,map(box,saved['continuous_high8']['branches'])):assert max(value[0],old[0])<=min(value[1],old[1])
    fr65=make_frame([8,10,63],65)
    t12=quotient(plus(root(box(2)),root(box(6))),box(2))
    t11=Q(2)-Q(3,11)**2+Q(3,11)**4/12;t13=Q(2)-Q(22,91)**2
    assert t11==Q(saved['theta11_rational_upper']) and t13==Q(saved['theta13_rational_lower'])
    assert len(saved['high65'])==4
    for m,r in zip(range(2,6),saved['high65']):
        assert r['label10']==m
        minus_value,plus_value=evaluate(fr65,profiles[0],[box(0),theta[m],box(0)])
        assert plus_value[0]>0
        if m==2:assert r['minus_branch']['kind']=='old44_only'
        else:
            lo,hi=(t12[1],t13) if m==3 else (t11,t12[0])
            assert lo<-minus_value[1]<=-minus_value[0]<hi
            assert r['minus_branch']['labels']==([12,13] if m==3 else [11,12])
    supports=primary['retained_old_finite_supports']+[r['support'] for r in primary['new_finite_vertices']]
    assert len(supports)==199
    lost=[v for v in supports if set(v)&{30,33,44}]
    assert len(lost)==18 and sum(8 in v for v in lost)==2
    assert lost==saved['exact_deletion_debit']['lost_finite_supports']
    assert saved['exact_deletion_debit']==prior['exact_deletion_debit'] and saved['exact_deletion_debit']['compact_faces_by_dimension']=={'2':12,'3':1}
    assert saved['complete_low_recipes']==46 and saved['complete_low_signs']==92 and saved['high_families']==5
    assert saved['changed_native_caps']==saved['new_normals']==saved['new_cones']==0 and saved['old44_control']
    print('PASS_INDEPENDENT_THREE_CHAIN_GATE',dict(low_signs=92,decisions=count,continuous_high8=1,native_first_high65=4,deletion=[18,12,1]),flush=True)


if __name__=='__main__':check()


