#!/usr/bin/env python3
"""Standalone exact saved cone/polar checks; no producer imports.

Reconstructs every rank-nine active-set kernel, oriented ray, contact gate,
five-scale passive-lift test, conjugate Schur bound and one-pole successor.
Original union geometry is frozen. This checks extreme projected directions,
not arbitrary acute completion, and never sets unknown passive contacts to zero.
"""
import json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_raw23_union_polar_completion_preflight_20261018'
UNION='fresh_raw23_same_family_union_support_pilot_20261018'
def read(t):return json.loads((ROOT/f'results/{t}.json').read_text())
def bind(d):
    for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
def fields(xs):return[Q(*x)for x in xs]
def tuplekey(xs):return tuple((x.a,x.b)for x in xs)
def dot(x,g,y):return sum((x[i]*g[i][j]*y[j]for i in range(len(x))for j in range(len(y))),Q())
class Q:
    def __init__(self,a=0,b=0):self.a,self.b=F(a),F(b)
    def __add__(self,o):
        if not isinstance(o,Q):o=Q(o)
        return Q(self.a+o.a,self.b+o.b)
    __radd__=__add__
    def __neg__(self):return Q(-self.a,-self.b)
    def __sub__(self,o):return self+-o
    def __mul__(self,o):
        if not isinstance(o,Q):o=Q(o)
        return Q(self.a*o.a+5*self.b*o.b,self.a*o.b+self.b*o.a)
    __rmul__=__mul__
    def __truediv__(self,o):
        if not isinstance(o,Q):o=Q(o)
        n=o.a*o.a-5*o.b*o.b;assert n
        return self*Q(o.a/n,-o.b/n)
    def sign(self):
        if not self.b:return (self.a>0)-(self.a<0)
        if not self.a:return (self.b>0)-(self.b<0)
        if self.a*self.b>0:return (self.a>0)-(self.a<0)
        n=self.a*self.a-5*self.b*self.b
        return ((self.a>0)-(self.a<0))*((n>0)-(n<0))
    def __eq__(self,o):
        if not isinstance(o,Q):o=Q(o)
        return self.a==o.a and self.b==o.b

def mat(x):return [[Q(*y)for y in row]for row in x]
def principal(g,ids):return [[g[i][j]for j in ids]for i in ids]
def positive(g):
    a=[r[:]for r in g]
    for k in range(len(a)):
        if a[k][k].sign()<=0:return False
        for i in range(k+1,len(a)):
            for j in range(k+1,len(a)):a[i][j]=a[i][j]-a[i][k]*a[k][j]/a[k][k]
    return True

def inverse(g):
    n=len(g);a=[row[:]+[Q(int(i==j))for j in range(n)]for i,row in enumerate(g)]
    for i in range(n):
        k=next(k for k in range(i,n)if a[k][i]!=0);a[i],a[k]=a[k],a[i]
        z=a[i][i];a[i]=[x/z for x in a[i]]
        for k in range(n):
            if k!=i:
                z=a[k][i];a[k]=[x-z*y for x,y in zip(a[k],a[i])]
    return [row[n:]for row in a]


def kernel(rows):
    a=[r[:]for r in rows];p=[];i=0
    # Work from the last column, independent of the producer's free column.
    for j in reversed(range(10)):
        k=next((k for k in range(i,len(a))if a[k][j]!=0),None)
        if k is None:continue
        a[i],a[k]=a[k],a[i];v=a[i][j];a[i]=[x/v for x in a[i]]
        for k in range(len(a)):
            if k!=i:
                v=a[k][j];a[k]=[x-v*y for x,y in zip(a[k],a[i])]
        p.append(j);i+=1
        if i==len(a):break
    if i!=9:return i,None
    free=next(j for j in range(10)if j not in p);v=[Q(int(j==free))for j in range(10)]
    for i,j in enumerate(p):v[j]=-a[i][free]
    assert all(sum((x*y for x,y in zip(r,v)),Q())==0 for r in rows)
    return 9,v


def check_cone(c,rows,g):
    n=len(rows);assert c['wall_count']==n and c['dimension']==10 and mat(c['normal_rows'])==rows
    rays=c['rays'];assert [r['index']for r in rays]==list(range(len(rays)))
    keys={tuplekey(fields(r['products'])):r['index']for r in rays};assert len(keys)==len(rays)
    for r in rays:
        x=fields(r['coordinates']);prod=[sum((a*b for a,b in zip(row,x)),Q())for row in rows]
        assert len(x)==10 and prod==fields(r['products'])and all(z.sign()<=0 for z in prod)
        assert next(z for z in prod if z!=0)==-Q(1)
        norm=dot(x,g,x);assert norm==Q(*r['norm'])and norm.sign()==r['norm_sign']
        assert r['zero_facets']==[i for i,z in enumerate(prod)if z==0]
    domains=list(combinations(range(n),9));assert len(c['complete_active_sets'])==len(domains)
    observed=set()
    for active,record in zip(domains,c['complete_active_sets']):
        assert record['active']==list(active)
        rank,x=kernel([rows[i]for i in active]);assert record['rank']==rank
        if x is None:assert record['ray']is None;continue
        prod=[sum((a*b for a,b in zip(row,x)),Q())for row in rows]
        if all(z.sign()>=0 for z in prod):x=[-z for z in x];prod=[-z for z in prod]
        elif not all(z.sign()<=0 for z in prod):assert record['ray']is None;continue
        factor=-next(z for z in prod if z!=0);prod=[z/factor for z in prod]
        k=tuplekey(prod);assert k in keys and record['ray']==keys[k];observed.add(keys[k])
    assert observed==set(range(len(rays)))
    finite=[fields(r['coordinates'])for r in rays if r['norm_sign']<0]
    assert all(dot(a,g,b).sign()<0 for a in finite for b in finite)
    center=[sum((r[i]for r in finite),Q())for i in range(10)]
    assert center==fields(c['strict_timelike_center'])and dot(center,g,center).sign()<0
    products=[sum((a*b for a,b in zip(row,center)),Q())for row in rows]
    assert products==fields(c['center_products'])and all(z.sign()<0 for z in products)
    assert c['ray_counts']=={str(s):sum(r['norm_sign']==s for r in rays)for s in(-1,0,1)}
    return len(domains)


def gate(q):
    assert q.sign()>=0
    finite=(Q(),Q(1),Q(2),Q('3/2','1/2'),Q(3),Q('5/2','1/2'))
    if(q-4).sign()<0 and q not in finite:return'non_Coxeter_finite_angle'
    if(2*q.a).denominator!=1 or(q.a*q.a-5*q.b*q.b).denominator!=1:return'nonintegral_two_cycle'
    z=Q(q.a,-q.b)
    if z.sign()<0 or(z-4).sign()>0:return'conjugate_pair_minor'
    return None


def check_poles(p,c,g,sc):
    finite=[r for r in c['rays']if r['norm_sign']<0];outside=[r for r in c['rays']if r['norm_sign']>0]
    assert len(p['candidates'])==p['all_spacelike_extreme_rays_tested']==len(outside)
    assert not p['arbitrary_spacelike_combinations_tested']and not p['higher_cycle_arithmetic_sufficiency_claimed']
    assert p['arbitrary_passive_lift_relaxation_checked']and not p['unknown_passive_contacts_assumed_zero']
    passed=[]
    for r,test in zip(outside,p['candidates']):
        assert test['ray']==r['index'];x=fields(r['coordinates']);norm=Q(*r['norm'])
        fd=[dot(x,g,fields(z['coordinates']))for z in finite]
        assert fd==fields(test['finite_vertex_products'])
        preserve=all(v.sign()<0 for v in fd);assert preserve==test['preserves_all_finite_vertices']
        sq=[2*v*v/(norm*s)for v,s in zip(fields(r['products']),sc)]
        assert len(sq)==len(r['products'])==len(sc)and sq==fields(test['physical_squared_contacts'])
        gg=[gate(s)for s in sq];assert gg==test['gate_failures']
        accepted=preserve and not any(gg);assert accepted==test['passes_necessary_polar_gates']
        if accepted:passed.append(r['index'])
        cuts=[z['index']for z in outside if dot(x,g,fields(z['coordinates'])).sign()>0]
        assert cuts==test['cut_spacelike_rays']
        anchor=next(i for i,z in enumerate(sq)if z.sign()>0 and(z-4).sign()<0)
        assert test['passive_lift_anchor_wall']==anchor and Q(*test['passive_lift_anchor_square'])==sq[anchor]
        native_squares=(Q(1),Q(2),Q('3/2','1/2'),Q(3),Q('5/2','1/2'))
        assert len(test['passive_lift_tests'])==5;surviving=[]
        for native,rec in zip(native_squares,test['passive_lift_tests']):
            t=native/sq[anchor];within=(t-1).sign()<=0
            tc=Q(t.a,-t.b);conjugate_bound=tc.sign()>0 and(tc-1).sign()<=0
            assert Q(*rec['anchor_native_square'])==native and Q(*rec['projection_norm_fraction'])==t
            assert Q(*rec['passive_squared_norm'])==2*(Q(1)-t)and rec['within_physical_norm_bound']==within
            assert Q(*rec['conjugate_projection_norm_fraction'])==tc and rec['within_conjugate_known_span_bound']==conjugate_bound
            if within:
                ss=[t*z for z in sq];gg=[gate(z)for z in ss]
                assert fields(rec['scaled_known_contact_squares'])==ss and rec['gate_failures']==gg
                accept=conjugate_bound and preserve and not any(gg)
            else:
                assert rec['scaled_known_contact_squares']is None and rec['gate_failures']is None;accept=False
            assert rec['passes_necessary_lift_gates']==accept
            if accept:surviving.append(t)
        assert surviving==fields(test['passing_projection_norm_fractions'])==([Q(1)]if accepted else[])
    assert passed==p['passing_rays'];return len(outside)


def main():
    d,old=read(TAG),read(UNION);bind(d);bind(old)
    assert d['status']=='PASS_RAW23_EXTREME_RAY_POLAR_METHOD_PREFLIGHT'
    assert not d['old_matrix_changed']and d['optimizer_invocations']==0 and d['exported_columns']==[]
    assert not d['arbitrary_acute_completion_excluded']and not d['sequential_extreme_ray_polar_completion_found']
    assert not d['geometric_realization_claimed']and d['native_pair_tests_are_necessary_only']
    assert len(d['cases'])==len(old['cases'])==3
    active=poles=0
    for c,o in zip(d['cases'],old['cases']):
        assert c['kinds']==o['kinds']and c['relative_index']==o['relative_index']
        g0=mat(o['scaled_gram']);sc=fields(o['squared_scales']);basis=c['original_basis']
        assert c['original_scaled_gram']==o['scaled_gram']and c['original_squared_scales']==o['squared_scales']
        assert len(basis)==len(set(basis))==10 and all(0<=i<12 for i in basis)
        g=principal(g0,basis);assert g==mat(c['metric']);inverse(g)
        assert positive([[Q(x.a,-x.b)for x in row]for row in g])
        rows=[[g0[i][j]for j in basis]for i in range(12)]
        active+=check_cone(c['initial_cone'],rows,g);poles+=check_poles(c['initial_polar_tests'],c['initial_cone'],g,sc)
        assert c['initial_cone']['ray_counts']=={'-1':18,'0':0,'1':32}
        assert {tuple(r['zero_facets'])for r in c['initial_cone']['rays']if r['norm_sign']<0}==set(map(tuple,o['all_positive_definite_nine_supports']))
        passed=c['initial_polar_tests']['passing_rays'];nextstage=c['after_unique_polar']
        if c['kinds']==[1,6]:assert not passed and nextstage is None;continue
        assert len(passed)==1 and nextstage['added_initial_ray']==passed[0]
        r=c['initial_cone']['rays'][passed[0]];x=fields(r['coordinates']);norm=Q(*r['norm'])
        assert x==fields(nextstage['added_normal_coordinates'])and norm/2==Q(*nextstage['added_squared_scale'])
        test=next(t for t in c['initial_polar_tests']['candidates']if t['ray']==passed[0])
        assert len(test['cut_spacelike_rays'])==32
        newrow=[sum((x[i]*g[i][j]for i in range(10)),Q())for j in range(10)]
        active+=check_cone(nextstage['cone'],rows+[newrow],g)
        poles+=check_poles(nextstage['polar_tests'],nextstage['cone'],g,sc+[norm/2])
        assert nextstage['cone']['ray_counts']=={'-1':52,'0':0,'1':30}and not nextstage['polar_tests']['passing_rays']
    assert active==d['initial_support_sets']+d['additional_support_sets']==2090
    assert poles==d['initial_spacelike_poles']+d['subsequent_spacelike_poles']==156
    assert d['finite_passive_lift_scale_candidates']==5*poles==780 and d['passing_nonzero_passive_lifts']==0
    print('PASS_STANDALONE_RAW23_POLAR_METHOD_SAVED',dict(active_sets=active,polar_candidates=poles,passive_lift_scales=5*poles,arbitrary_completion_excluded=False),flush=True)


if __name__=='__main__':main()
