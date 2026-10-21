#!/usr/bin/env python3
"""Independent rational construction check, plus algebra in the written gap proof.

No producer imports. Published small-facet completeness and the geometric
opposite-cusp argument remain mathematical inputs, not a census replay.
"""
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PATH=ROOT/'results/fresh_ideal_ra4_double24_volume_20261018.json'

def vec(v):return tuple(Q(x) for x in v)
def inner(x,y):return sum(x[i]*y[i] for i in range(1,len(x)))-x[0]*y[0]
def eliminate(a):
    a=[list(map(Q,row)) for row in a]
    if not a:return a,[]
    piv=[];i=0
    for j in reversed(range(len(a[0]))):
        p=next((k for k in range(i,len(a)) if a[k][j]),None)
        if p is None:continue
        a[i],a[p]=a[p],a[i];z=a[i][j];a[i]=[x/z for x in a[i]]
        for k in range(len(a)):
            if k!=i:
                z=a[k][j];a[k]=[x-z*y for x,y in zip(a[k],a[i])]
        piv.append(j);i+=1
        if i==len(a):break
    return a,piv
def rk(a):return len(eliminate(a)[1])
def solution(a,b):
    # The augmented constant is first, so reverse elimination pivots spatial
    # columns before it. Four spatial pivots mean the affine solution is unique.
    rows,piv=eliminate([[y,*r] for r,y in zip(a,b)])
    if len(piv)!=4 or 0 in piv:return None
    ans=[Q(0)]*4
    for r,j in zip(rows,piv):ans[j-1]=r[0]
    return tuple(ans)

def check_poly(p):
    ns=list(map(vec,p['walls']));vs=list(map(vec,p['vertices']))
    assert len(ns)==len(set(ns)) and len(vs)==len(set(vs))
    assert all(inner(n,n)==1 for n in ns)
    assert all(v[0]==1 and inner(v,v)==0 for v in vs)
    assert all(inner(n,v)<=0 for n in ns for v in vs)
    assert all(inner(n,m)<=0 for n,m in combinations(ns,2))
    fs=[frozenset(i for i,v in enumerate(vs) if inner(n,v)==0) for n in ns]
    faces={frozenset(range(len(vs)))}
    changed=True
    while changed:
        more={a&b for a in faces for b in fs};more.discard(frozenset())
        changed=not more<=faces;faces|=more
    bydim={d:[] for d in range(5)}
    for f in faces:
        v=vs[min(f)]
        d=rk([[vs[i][j]-v[j] for j in range(1,5)] for i in f])
        bydim[d].append(f)
    assert len(bydim[4])==1
    for d in range(4):
        expected={frozenset(x) for x in p['faces'][str(d)]}
        assert len(expected)==len(p['faces'][str(d)]) and expected==set(bydim[d])
    fv=[len(bydim[d]) for d in range(4)]
    assert fv==p['f_vector'] and set(fs)==set(bydim[3])
    hist=Counter(map(len,bydim[2]))
    assert dict(hist)=={int(k):v for k,v in p['ridge_histogram'].items()}
    assert fv[1]==4*fv[0] and fv[2]==3*fv[0]+fv[3]
    assert sum(k*v for k,v in hist.items())==12*fv[0]
    assert p['ridge_defect']==sum((k-3)*v for k,v in hist.items())==3*(fv[0]-fv[3])
    assert Q(p['volume_coefficient_pi_squared'])==Q(fv[0]-fv[3]+4,3)
    assert Q(p['reflection_orbifold_euler_characteristic'])==Q(fv[0]-fv[3]+4,4)
    for e in bydim[1]:assert len(e)==2 and sum(e<=f for f in fs)==3
    for f in bydim[2]:
        ids=[i for i,h in enumerate(fs) if f<=h]
        assert len(ids)==2 and inner(ns[ids[0]],ns[ids[1]])==0
    assert len(p['cusp_incidence'])==len(vs)
    for i,row in enumerate(p['cusp_incidence']):
        ids=[j for j,f in enumerate(fs) if i in f]
        assert row['vertex']==i and row['facets']==ids and len(ids)==6
        pairs=[(j,k) for j,k in combinations(ids,2) if inner(ns[j],ns[k])==-1]
        assert [list(p) for p in pairs]==row['parallel_pairs']
        assert sorted(x for p in pairs for x in p)==ids
        assert all(inner(ns[j],ns[k]) in (0,-1) for j,k in combinations(ids,2))
        assert sum(i in e for e in bydim[1])==8 and sum(i in f for f in bydim[2])==12
    assert len(p['facet_cusp_shapes'])==len(fs)
    types=Counter()
    for fi,f in enumerate(fs):
        row=p['facet_cusp_shapes'][fi]
        assert row['facet']==fi and row['vertices']==len(f) and len(row['cusps'])==len(f)
        nsquare=0
        for vi,saved in zip(sorted(f),row['cusps']):
            ids=[j for j,h in enumerate(fs) if j!=fi and vi in h and f&h in bydim[2]]
            assert len(ids)==4
            pairs=[(j,k) for j,k in combinations(ids,2) if inner(ns[j],ns[k])==-1]
            assert len(pairs)==2 and sorted(x for p in pairs for x in p)==sorted(ids)
            sums=[tuple(ns[j][l]+ns[k][l] for l in range(5)) for j,k in pairs]
            assert all(w==tuple(w[0]*x for x in vs[vi]) for w in sums)
            ratio=sums[0][0]/sums[1][0];assert ratio>0
            ratio=max(ratio,1/ratio)
            assert saved=={'vertex':vi,'aspect_ratio':str(ratio)}
            nsquare+=ratio!=1
        assert row['nonsquare_cusps']==nsquare
        if len(f)==6:assert nsquare==0
        if len(f)==9:assert nsquare==3
        edges=[e for e in bydim[1] if e<=f];ridges=[r for r in bydim[2] if r<=f]
        types[(len(f),len(edges),len(ridges),tuple(sorted(Counter(map(len,ridges)).items())))]+=1
    claimed=Counter()
    for t in p['facet_types']:
        key=(t['vertices'],t['edges'],t['ridges'],tuple(sorted((int(k),v) for k,v in t['ridge_histogram'].items())))
        claimed[key]+=t['count']
    assert types==claimed
    return ns,vs

def pmul(a,b):
    out=[Q(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    while len(out)>1 and not out[-1]:out.pop()
    return tuple(out)
def padd(a,b):
    out=[Q(0)]*max(len(a),len(b))
    for i,x in enumerate(a):out[i]+=x
    for i,x in enumerate(b):out[i]+=x
    while len(out)>1 and not out[-1]:out.pop()
    return tuple(out)
class RF:
    def __init__(self,n,d=(1,)):
        self.n=(Q(n),) if isinstance(n,int) else tuple(map(Q,n));self.d=tuple(map(Q,d))
    def __add__(self,o):
        if not isinstance(o,RF):o=RF(o)
        return RF(padd(pmul(self.n,o.d),pmul(o.n,self.d)),pmul(self.d,o.d))
    __radd__=__add__
    def __neg__(self):return RF(tuple(-x for x in self.n),self.d)
    def __sub__(self,o):return self+-RF(o) if not isinstance(o,RF) else self+-o
    def __rsub__(self,o):return -self+o
    def __mul__(self,o):
        if not isinstance(o,RF):o=RF(o)
        return RF(pmul(self.n,o.n),pmul(self.d,o.d))
    __rmul__=__mul__
    def __truediv__(self,o):
        if not isinstance(o,RF):o=RF(o)
        return RF(pmul(self.n,o.d),pmul(self.d,o.n))
    def __rtruediv__(self,o):return RF(o)/self
    def __eq__(self,o):
        if not isinstance(o,RF):o=RF(o)
        return pmul(self.n,o.d)==pmul(o.n,self.d)

def proof_algebra():
    c=RF((0,1));h2=(1-c)/(1+c);r2=2*c/(1+c);b2=2/(c*(1+c));mu2=2*(1-c)/c
    identities=[h2+r2==1,mu2*r2==4*h2,b2*c*c==r2,
      -1+b2+h2==mu2,-1+b2*c-h2==0,
      -1+b2*(2*c*c-1)+h2==-mu2,
      (1+h2)/(2*r2)==1/(2*c)]
    assert all(identities)
    # Exhaustive symbolic cases in the written proof, not a facet census.
    assert 2*(2+9)==22 and 3*(1+3+3)==21 and 4*(1+4)==20 and 4*5==20
    assert min(w for w in range(0,30,3) if 2*w>=20)==12
    return len(identities)

def main():
    data=json.loads(PATH.read_text())
    assert data['source_sha256']==sha256((ROOT/'scripts/fresh_ideal_ra4_double24_volume_20261018.py').read_bytes()).hexdigest()
    ns,vs=check_poly(data['original']);ds,dv=check_poly(data['double'])
    expected_walls=set()
    for ij in combinations(range(4),2):
        for signs in product((-1,1),repeat=2):
            x=[Q(0)]*4
            for i,s in zip(ij,signs):x[i]=Q(s)
            expected_walls.add((Q(1),*x))
    assert set(ns)==expected_walls
    found=set();bases=0
    for ids in combinations(range(24),4):
        x=solution([ns[i][1:] for i in ids],[1]*4)
        if x is not None and all(sum(a*b for a,b in zip(x,n[1:]))<=1 for n in ns):
            found.add((Q(1),*x));bases+=1
    assert found==set(vs) and bases==data['feasible_active_bases']
    n=vec(data['gluing_wall']);assert n in ns and inner(n,n)==1
    def refl(x):return tuple(x[i]-2*inner(x,n)*n[i] for i in range(5))
    rv={tuple(y/refl(x)[0] for y in refl(x)) for x in vs}
    assert set(dv)==set(vs)|rv
    assert set(ds)==(set(ns)-{n})|{refl(x) for x in ns if x!=n}
    assert set(map(vec,data['shared_cusps']))==set(vs)&rv and len(data['shared_cusps'])==6
    assert data['original']['f_vector']==[24,96,96,24]
    assert data['double']['f_vector']==[42,168,164,38]
    assert data['double']['ridge_histogram']=={'3':152,'4':12}
    assert data['intermediate_ridge_defects']==[3,6,9]
    assert data['possible_intermediate_volume_coefficients_pi_squared']==['5/3','2','7/3']
    assert data['second_minimum_proved'] is False and data['equality_classification_proved'] is False
    count=proof_algebra()
    print(json.dumps({'status':'PASS','construction_vertices':[24,42],'active_four_subsets':10626,
      'antiprism_rational_identities':count,'double_volume_pi_squared':'8/3',
      'scope':'exact construction and proof algebra; published facet classification and geometric counting are mathematical inputs'}))

if __name__=='__main__':main()
