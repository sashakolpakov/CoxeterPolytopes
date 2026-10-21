#!/usr/bin/env python3
"""Standalone rational-interval check of the fixed exact66 quadrilateral cap.

No project imports. Directed arithmetic follows the frozen saved-checker design;
adjugate reconstruction and a separate Gaussian inverse certify the necessary
norm signs. The finite catalogue is rebuilt explicitly, including a symbolic
high-bond category. The original physical scaffold/profile theorem is inherited.
"""
from fractions import Fraction as F
from hashlib import sha256
from itertools import permutations, product
from functools import lru_cache
from math import isqrt
from pathlib import Path
import json
import math
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_exact66_quadrilateral_cap_obstruction_20261018'
SEED='fresh_remote60_exact_original_extension_20261018'
PRIOR='fresh_exact64_escape7_cap_20261018'
CURRENT='fresh_exact65_escape33_cap_20261018'
SCALE=10**70

def iv(x):
    return tuple(map(F, x)) if isinstance(x, (list, tuple)) else (F(x), F(x))

def roundout(v):
    a, b = v
    assert a <= b
    return F(a.numerator*SCALE//a.denominator, SCALE), F(-((-b.numerator*SCALE)//b.denominator), SCALE)

def add(a, b): return roundout((a[0]+b[0], a[1]+b[1]))

def neg(a): return -a[1], -a[0]

def sub(a, b): return add(a, neg(b))

def mul(a, b):
    z = [x*y for x in a for y in b]
    return roundout((min(z), max(z)))

def div(a, b):
    assert b[0]*b[1] > 0
    return mul(a, (1/b[1], 1/b[0]))

def total(values):
    out = iv(0)
    for v in values: out = add(out, v)
    return out

def sqrt(a):
    assert a[0] >= 0
    lo = isqrt(a[0].numerator*SCALE*SCALE//a[0].denominator)
    hi = isqrt(a[1].numerator*SCALE*SCALE//a[1].denominator)
    if F(hi, SCALE)**2 < a[1]: hi += 1
    return F(lo, SCALE), F(hi, SCALE)

def dot(x, g, y): return total(mul(mul(a, g[i][j]), b) for i, a in enumerate(x) for j, b in enumerate(y))

def enc(x):
    if isinstance(x, F): return str(x)
    if isinstance(x, dict): return {str(k): enc(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [enc(v) for v in x]
    return x

def determinant(g):
    n, terms = len(g), []
    for p in permutations(range(n)):
        value = iv((-1)**sum(p[i] > p[j] for i in range(n) for j in range(i+1, n)))
        for i in range(n): value = mul(value, g[i][p[i]])
        terms.append(value)
    return total(terms)

def inverse(g):
    n, d = len(g), determinant(g)
    assert d[0]*d[1] > 0
    return [[div(mul(iv((-1)**(i+j)), determinant([[g[a][b] for b in range(n) if b != i]
        for a in range(n) if a != j])), d) for j in range(n)] for i in range(n)], d

def gauss(g):
    n = len(g)
    a = [r[:]+[iv(int(i == j)) for j in range(n)] for i, r in enumerate(g)]
    pivots = []
    for i in range(n):
        p = a[i][i]
        assert p[0]*p[1] > 0
        pivots.append(p)
        a[i] = [div(x, p) for x in a[i]]
        for j in range(n):
            if i != j:
                q = a[j][i]
                a[j] = [sub(x, mul(q, y)) for x, y in zip(a[j], a[i])]
    return [r[n:] for r in a], pivots

def theta(m):
    if m in (2, 3): return iv(m-2)
    if m == 4: return sqrt(iv(2))
    if m == 5: return div(add(iv(1), sqrt(iv(5))), iv(2))
    if m == 6: return sqrt(iv(3))
    assert isinstance(m,int) and m>=7
    def atan(x):
        z = sum((-1)**k*x**(2*k+1)/F(2*k+1) for k in range(100))
        return z, z+x**201/F(201)
    pi = sub(mul(iv(16), atan(F(1, 5))), mul(iv(4), atan(F(1, 239))))
    x = div(pi, iv(m))
    xx, term, value = mul(x, x), iv(1), iv(1)
    for k in range(1, 60):
        term = div(mul(neg(xx), term), iv((2*k-1)*(2*k)))
        value = add(value, term)
    e = xx[1]*max(abs(t) for t in term)/F(119*120)
    return mul(iv(2), (value[0]-e, value[1]+e))

def read(stem):return json.loads((ROOT/f'results/{stem}.json').read_text())

@lru_cache(None)
def graph_key(flat):
    n=math.isqrt(len(flat));g=[flat[i*n:(i+1)*n]for i in range(n)]
    return min(tuple(g[p[i]][p[j]]for i in range(n)for j in range(i+1,n))for p in permutations(range(n)))

@lru_cache(None)
def finite_catalogue():
    out={}
    def add_graph(name,n,edges):
        g=[[2]*n for _ in range(n)]
        for i in range(n):g[i][i]=1
        for i,j,m in edges:g[i][j]=g[j][i]=m
        out[n,graph_key(tuple(x for r in g for x in r))]=name
    for n in range(1,6):add_graph('A'+str(n),n,[(i,i+1,3)for i in range(n-1)])
    for n in range(2,6):add_graph('B'+str(n),n,[(i,i+1,4 if i==0 else 3)for i in range(n-1)])
    for n in(4,5):add_graph('D'+str(n),n,[(i,i+1,3)for i in range(n-3)]+[(n-3,n-2,3),(n-3,n-1,3)])
    add_graph('F4',4,[(0,1,3),(1,2,4),(2,3,3)])
    for n in(3,4):add_graph('H'+str(n),n,[(i,i+1,5 if i==0 else 3)for i in range(n-1)])
    add_graph('I2(5)',2,[(0,1,5)])
    # Category6 denotes any finite label >=6: it must be an isolated I2.
    add_graph('I2(6)',2,[(0,1,6)])
    return out

def finite_duty(pi,b,m,n):
    g=[[2]*5 for _ in range(5)]
    for i in range(5):g[i][i]=1
    g[0][1]=g[1][0]=3
    if pi:g[0 if pi<=3 else 1][4]=g[4][0 if pi<=3 else 1]=3+(pi-1)%3
    for i,j,l in[(2,3,b),(2,4,m),(3,4,n)]:g[i][j]=g[j][i]=l
    unseen=set(range(5));names=[]
    while unseen:
        todo=[min(unseen)];comp=set()
        while todo:
            i=todo.pop()
            if i in comp:continue
            comp.add(i);todo.extend(j for j in unseen if j!=i and g[i][j]!=2)
        unseen-=comp;vs=sorted(comp);key=len(vs),graph_key(tuple(g[i][j]for i in vs for j in vs))
        name=finite_catalogue().get(key)
        if name is None:return None
        names.append(name)
    return sorted(names)

def normalized(names):
    if names is None:return None
    return sorted({'A2':'I2(3)','B2':'I2(4)'}.get(x,x) for x in names)

def native_intersection(raw, proof):
    assert proof['raw']==enc(raw)
    label=proof['identity'];exact=neg(theta(label))
    if label==2:
        assert -1<raw[0]<=0<=raw[1]<1
    else:
        assert raw[0]>-theta(label+1)[0] and raw[1]<-theta(label-1)[1]
    z=max(raw[0],exact[0]),min(raw[1],exact[1]);assert z[0]<=z[1]
    assert proof['identity_interval']==enc(exact) and proof['intersection']==enc(z)
    return z

def main():
    saved=read(TAG);d=saved['necessary_norm_certificate']
    assert saved['status']=='PASS_EXACT66_FOUR_FINITE_DUTY_CAP_OBSTRUCTION'
    assert len(saved['sources'])==10
    for path,h in saved['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==h,path
    seed,prior,current=read(SEED),read(PRIOR),read(CURRENT)
    for dep in (seed,prior,current):
        for path,h in dep['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==h,path
    assert current['physical']['summary']['total_walls']==66
    assert current['physical']['summary']['complete_finite_vertices']==205
    assert current['physical']['summary']['infinite_volume_proved'] is True
    frame=[30,52,47,76];assert d['frame']==frame
    assert d['cycle_labels']==[4,2,2,2] and d['dotted_diagonals']==[[30,47],[52,76]]
    assert all(seed['white_covectors'][str(i)]==enc([iv(0)]*6) for i in frame[:-1])
    assert prior['physical']['profile']==0 and prior['physical']['new_wall']==76
    roots={int(i):list(map(iv,v)) for i,v in seed['physical_input']['roots'].items()}
    metric=[list(map(iv,r)) for r in seed['physical_input']['metric']]
    proofs={tuple(p['pair']):p for p in d['original_native_identity_refinements']}
    assert len(proofs)==len(d['original_native_identity_refinements'])==4
    assert set(proofs)=={(30,52),(52,30),(52,47),(47,52)}
    gram=[]
    for i in frame:
        row=[]
        for j in frame:
            if i==j:z=iv(2)
            elif 76 in (i,j):z=iv(prior['physical']['native_products'][str(j if i==76 else i)])
            else:
                z=dot(roots[i],metric,roots[j])
                if (i,j) in proofs:z=native_intersection(z,proofs[i,j])
                else:assert z[1]<-2
            row.append(z)
        gram.append(row)
    assert enc(gram)==d['frame_Gram']
    assert gram[0][1]==neg(theta(4))
    assert gram[0][3]==gram[1][2]==gram[2][3]==iv(0)
    assert gram[0][2][1]<-2 and gram[1][3][1]<-2
    inv,det=inverse(gram);assert det[1]<0
    assert enc(inv)==d['frame_inverse'] and enc(det)==d['frame_determinant']
    other,pivots=gauss(gram)
    assert sum(p[0]>0 for p in pivots)==3 and sum(p[1]<0 for p in pivots)==1
    low=d['low'];assert low['profiles']==list(range(7)) and low['labels']==[2,3,4,5]
    assert low['complete_recipes']==1792 and low['eligible_count']==low['norm_excluded_count']==84
    assert len(low['eligible'])==84 and len(low['finite_duty_rejected'])==1708
    rows=low['eligible']+low['finite_duty_rejected']
    bykey={(r['profile'],tuple(r['labels'])):r for r in rows};assert len(bykey)==1792
    accepted=0
    for pi,ms in product(range(7),product(range(2,6),repeat=4)):
        r=bykey[pi,ms]
        duties=[normalized(finite_duty(pi,b,ms[i],ms[(i+1)%4])) for i,b in enumerate([4,2,2,2])]
        assert duties==r['duties']
        if None in duties:
            assert set(r)=={'profile','labels','duties'}
            continue
        accepted+=1
        rhs=[neg(theta(m)) for m in ms]
        q=iv(2) if pi==0 else sub(iv(2),mul(iv(F(2,3)),mul(theta(3+(pi-1)%3),theta(3+(pi-1)%3))))
        sq=dot(rhs,inv,rhs);gap=sub(sq,q);othergap=sub(dot(rhs,other,rhs),q)
        assert r['rhs']==enc(rhs) and r['projected_square']==enc(sq)
        assert r['required_profile_square']==enc(q) and r['norm_gap']==enc(gap)
        sign=-1 if gap[1]<0 else 1
        assert gap[1]<0 or gap[0]>0
        assert r['sign']==sign and (othergap[1]<0 if sign<0 else othergap[0]>0)
    assert accepted==84
    # Complete symbolic high-label category: classification depends only on
    # isolation of every I2(m), m>=6. No high-label numerical sampling is used.
    symbolic=[]
    for pi,ms in product(range(7),product(range(2,7),repeat=4)):
        if 6 not in ms:continue
        if all(finite_duty(pi,b,ms[i],ms[(i+1)%4]) is not None for i,b in enumerate([4,2,2,2])):
            symbolic.append((pi,ms))
    expect=[(0,(2,2,6,2)),(0,(3,2,6,2)),(0,(2,2,2,6)),(0,(2,3,2,6))]
    assert set(symbolic)==set(expect) and len(symbolic)==4
    assert len(d['high'])==4
    expected_order=[(2,0,2),(2,0,3),(3,1,2),(3,1,3)]
    root3=sqrt(iv(3))
    for r,(axis,opp,m) in zip(d['high'],expected_order):
        assert r['profile']==0 and r['high_wall']==frame[axis] and r['opposite_wall']==frame[opp] and r['opposite_label']==m
        e=[iv(-int(i==axis)) for i in range(4)]
        b=[neg(theta(m)) if i==opp else iv(0) for i in range(4)]
        aa=dot(e,inv,e);bb=mul(iv(2),dot(e,inv,b));cc=sub(dot(b,inv,b),iv(2))
        val=add(add(mul(aa,iv(3)),mul(bb,root3)),cc)
        deriv=add(mul(iv(2),mul(aa,root3)),bb)
        assert r['coefficients']==enc([aa,bb,cc]) and r['theta_interval']==enc([root3[0],F(2)])
        assert r['value_at_sqrt3']==enc(val) and r['derivative_at_sqrt3']==enc(deriv)
        assert aa[0]>0 and bb[0]>=0 and val[0]>0 and deriv[0]>0 and r['strictly_increasing'] is True
        a2=dot(e,other,e);b2=mul(iv(2),dot(e,other,b));c2=sub(dot(b,other,b),iv(2))
        assert a2[0]>0 and b2[0]>=0 and add(add(mul(a2,iv(3)),mul(b2,root3)),c2)[0]>0
    assert d['all_declared_profile_and_finite_label_caps_excluded'] is True
    assert saved['deletion_debit']=={'deleted_old_walls':[],'old_finite_vertices_lost':0}
    assert saved['scope']==dict(fixed_exact66_parent=True,seven_profiles_only=True,all_finite_labels=True,
        ideal_or_dotted_duties_not_covered=True,other_incidence_or_parent_not_excluded=True,
        chamber_or_cone_search=False,finite_volume_completion_excluded=False,new_geometric_example=False)
    print(json.dumps(dict(status='PASS',low_recipes=1792,finite_duty_rejects=1708,norm_exclusions=84,
        symbolic_high_families=4,high_parameters='complete closed [sqrt(3),2]',new_walls=0,
        scope='fixed-parent four-finite-duty norm obstruction; historical scaffold/profile theorem inherited')))

if __name__=='__main__':main()
