#!/usr/bin/env python3
"""Saved exact3D face checks, with no producer imports or native125 search.

Frozen rays/native-field hypotheses remain inputs. Exact Q(sqrt5) arithmetic
uses the earlier saved-checker design; rank pivots run in reverse columns.
Stationary bounds are checked by their defining linear equations.
"""
import json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from collections import Counter
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_raw23_three_face_anchor_preflight_20261018'
OLD='fresh_raw23_union_polar_completion_preflight_20261018'
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
        if not self.b:return(self.a>0)-(self.a<0)
        if not self.a:return(self.b>0)-(self.b<0)
        if self.a*self.b>0:return(self.a>0)-(self.a<0)
        n=self.a*self.a-5*self.b*self.b
        return((self.a>0)-(self.a<0))*((n>0)-(n<0))
    def __eq__(self,o):
        if not isinstance(o,Q):o=Q(o)
        return self.a==o.a and self.b==o.b

def vec(xs):return[Q(*x)for x in xs]
def mat(xs):return[vec(r)for r in xs]
def read(t):return json.loads((ROOT/f'results/{t}.json').read_text())
def bind(d):
    for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p

def rank(rows):
    a=[r[:]for r in rows]
    if not a:return 0
    piv=0
    for col in reversed(range(len(a[0]))):
        at=next((i for i in range(piv,len(a))if a[i][col]!=0),None)
        if at is None:continue
        a[piv],a[at]=a[at],a[piv];q=a[piv][col];a[piv]=[x/q for x in a[piv]]
        for i in range(piv+1,len(a)):
            q=a[i][col];a[i]=[x-q*y for x,y in zip(a[i],a[piv])]
        piv+=1
        if piv==len(a):break
    return piv

def bil(x,g,y):return sum((x[i]*g[i][j]*y[j]for i in range(len(x))for j in range(len(y))),Q())
def pd(g):
    a=[r[:]for r in g]
    for k in range(len(a)):
        if a[k][k].sign()<=0:return False
        for i in range(k+1,len(a)):
            for j in range(k+1,len(a)):a[i][j]=a[i][j]-a[i][k]*a[k][j]/a[k][k]
    return True

def bound(record,h,p,s):
    expected=[list(ids)for n in(1,2,3)for ids in combinations(range(3),n)]
    assert[rr['support']for rr in record['supports']]==expected
    admissible=[]
    for rr,ids in zip(record['supports'],expected):
        z=vec(rr['coefficients'])
        assert len(z)==len(ids)
        assert all(sum((h[i][j]*z[k]for k,j in enumerate(ids)),Q())==p[i]for i in ids)
        positive=all(x.sign()>0 for x in z);assert positive==rr['strictly_positive']
        if positive:
            v=2*sum((p[i]*z[k]for k,i in enumerate(ids)),Q())/s
            assert Q(*rr['bound'])==v;admissible.append((ids,v))
        else:assert rr['bound']is None
    assert admissible
    maximum=admissible[0][1]
    for _,v in admissible:
        if(v-maximum).sign()>0:maximum=v
    maximizers=[ids for ids,v in admissible if v==maximum]
    interior=[0,1,2]in maximizers
    assert Q(*record['maximum_square'])==maximum and record['maximizing_supports']==maximizers
    assert record['attained_in_relative_interior']==interior
    assert record['strictly_below_four_in_relative_interior']==((maximum-4).sign()<0 or maximum==4 and not interior)
    return maximum

def check():
    data=read(TAG);old=read(OLD);bind(data);bind(old)
    assert data['status']=='PASS_RAW23_THREE_FACE_ANCHOR_PREFLIGHT'
    assert len(data['cases'])==len(old['cases'])==3
    total=Counter()
    for case,original in zip(data['cases'],old['cases']):
        assert case['kinds']==original['kinds']
        cone=original['initial_cone'];rays=cone['rays'];byid={r['index']:r for r in rays}
        assert len(rays)==len(byid)==50
        xyz={i:vec(r['coordinates'])for i,r in byid.items()};rows=mat(cone['normal_rows'])
        metric=mat(original['metric']);sc=vec(original['original_squared_scales'])
        space=[i for i,r in byid.items()if r['norm_sign']>0];assert len(space)==32
        expected=list(combinations(space,3));assert len(case['triple_decisions'])==len(expected)==4960
        grouped={};rankcache={}
        for rr,ids in zip(case['triple_decisions'],expected):
            assert rr['rays']==list(ids)
            common=tuple(sorted(set.intersection(*(set(byid[i]['zero_facets'])for i in ids))))
            if common not in rankcache:rankcache[common]=rank([rows[i]for i in common])
            r=rankcache[common];independent=rank([xyz[i]for i in ids])==3
            found=r==7 and independent
            assert rr['common_active_walls']==list(common)and rr['common_rank']==r
            assert rr['linearly_independent']==independent and rr['three_dimensional_face']==found
            if found:grouped.setdefault(common,[]).append(list(ids))
        assert[rr['active_walls']for rr in case['faces']]==[list(k)for k in sorted(grouped)]
        counts=Counter()
        for face in case['faces']:
            active=tuple(face['active_walls']);full=[i for i in byid if set(active)<=set(byid[i]['zero_facets'])]
            assert full==face['all_extreme_rays']
            assert sorted(set.intersection(*(set(byid[i]['zero_facets'])for i in full)))==list(active)
            assert rank([xyz[i]for i in full])==3
            assert face['spacelike_generating_triples']==grouped[active]
            assert face['basis_rays']==grouped[active][0]
            ids=face['basis_rays'];h=[[bil(xyz[i],metric,xyz[j])for j in ids]for i in ids]
            assert h==mat(face['plane_gram'])and pd(h)and face['plane_inertia']==[3,0,0]
            assert face['positive_definite_plane']and face['all_rays_spacelike']
            assert face['extreme_ray_signs']==[byid[i]['norm_sign']for i in full]
            assert all(byid[i]['norm_sign']>0 for i in full)
            assert face['simplicial']==(len(full)==3)and len(full)in(3,4)
            triples=[list(t)for t in combinations(full,3)if rank([xyz[i]for i in t])==3]
            assert[rr['rays']for rr in face['independent_triple_cover']]==triples
            grams={}
            for rr,tt in zip(face['independent_triple_cover'],triples):
                g=[[bil(xyz[i],metric,xyz[j])for j in tt]for i in tt]
                assert g==mat(rr['gram'])and pd(g);grams[tuple(tt)]=g
            contacts=[[-Q(*byid[i]['products'][w])for i in full]for w in range(len(rows))]
            assert all(x.sign()>=0 for p in contacts for x in p)
            nonzero=[w for w,p in enumerate(contacts)if any(x!=0 for x in p)]
            assert[rr['wall']for rr in face['full_cone_contact_bounds']]==nonzero
            excludes=[]
            for rr in face['full_cone_contact_bounds']:
                w=rr['wall'];assert vec(rr['contact_on_all_rays'])==contacts[w]
                assert[ss['rays']for ss in rr['subcones']]==triples
                values=[]
                for ss,tt in zip(rr['subcones'],triples):
                    pp=[-Q(*byid[i]['products'][w])for i in tt]
                    zero=all(x==0 for x in pp);assert ss['zero_contact']==zero
                    if zero:assert Q(*ss['maximum_square'])==0;values.append(Q())
                    else:values.append(bound(ss,grams[tuple(tt)],pp,sc[w]))
                maximum=values[0]
                for v in values:
                    if(v-maximum).sign()>0:maximum=v
                assert Q(*rr['maximum_square'])==maximum
                below=(maximum-1).sign()<0;assert rr['strictly_below_one']==below
                if below:excludes.append(w)
            assert excludes==face['native_excluding_walls']and bool(excludes)==face['relative_interior_native_excluded']
            if face['simplicial']:
                assert mat(face['contact_magnitudes'])==contacts
                assert[rr['wall']for rr in face['anchor_bounds']]==nonzero
                eligible=[]
                for rr in face['anchor_bounds']:
                    w=rr['wall'];assert vec(rr['contact'])==contacts[w]
                    bound(rr,h,contacts[w],sc[w])
                    if rr['strictly_below_four_in_relative_interior']:eligible.append(w)
                assert face['eligible_anchor_walls']==eligible
                assert face['eligible_anchor_rank']==rank([contacts[w]for w in eligible])==3
                anchors=next(list(tt)for tt in combinations(eligible,3)if rank([contacts[w]for w in tt])==3)
                assert face['anchors']==anchors and face['anchor_status']=='three_independent_finite_anchors'
            else:
                assert face['anchor_status']==('nonsimplicial_strict_native_exclusion'if excludes else'nonsimplicial_no_strict_native_exclusion')
            counts[face['anchor_status']]+=1
        counts.update(triple_domain=4960,found_triples=sum(len(x)for x in grouped.values()),distinct_faces=len(grouped),
            all_spacelike_faces=len(grouped),simplicial_faces=sum(f['simplicial']for f in case['faces']),PD_faces=len(grouped),
            native_excluded_faces=sum(f['relative_interior_native_excluded']for f in case['faces']),
            independent_subcone_count=sum(len(f['independent_triple_cover'])for f in case['faces']))
        assert dict(counts)==case['counts'];total.update(counts)
    assert dict(total)==data['summary']
    assert total['native_excluded_faces']==total['distinct_faces']==552
    assert total['simplicial_faces']==192 and total['nonsimplicial_strict_native_exclusion']==360
    assert total['independent_subcone_count']==1632
    assert not data['native_candidate_enumeration_performed']and not data['nonsimplicial_faces_treated_as_unresolved']
    assert data['all_PD_faces_use_complete_triple_cover']and data['strict_less_than_one_only_for_full_cone_exclusion']
    assert data['all_entirely_spacelike_three_faces_covered']and not data['all_three_faces_regardless_of_ray_signs_claimed']
    assert data['arbitrary_passive_lifts_relaxed']and not data['unknown_passive_contacts_assumed_zero']
    assert data['exported_columns']==[]and data['optimizer_invocations']==0
    assert all(data[k]is False for k in('geometric_realization_claimed','new_normal_claimed','pair_exclusion_claimed'))
    print('PASS_RAW23_THREE_FACE_SAVED',dict(total))
if __name__=='__main__':check()
