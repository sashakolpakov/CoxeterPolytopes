#!/usr/bin/env python3
"""Saved exact4D face checks, with no producer imports or native625 search.

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
TAG='fresh_raw23_four_face_contact_preflight_20261018'
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

def inertia(g):
    # Congruence elimination, including zero-diagonal indefinite2x2 blocks.
    a=[r[:]for r in g];positive=negative=0
    while a:
        n=len(a);k=next((i for i in range(n)if a[i][i]!=0),None)
        if k is not None:
            order=[k]+[i for i in range(n)if i!=k];a=[[a[i][j]for j in order]for i in order]
            z=a[0][0]
            if z.sign()>0:positive+=1
            else:negative+=1
            a=[[a[i][j]-a[i][0]*a[0][j]/z for j in range(1,n)]for i in range(1,n)]
            continue
        pair=next(((i,j)for i in range(n)for j in range(i+1,n)if a[i][j]!=0),None)
        if pair is None:return[positive,negative,n]
        i,j=pair;order=[i,j]+[k for k in range(n)if k not in pair]
        a=[[a[i][j]for j in order]for i in order];z=a[0][1]
        positive+=1;negative+=1
        a=[[a[i][j]-(a[i][0]*a[1][j]+a[i][1]*a[0][j])/z for j in range(2,n)]for i in range(2,n)]
    return[positive,negative,0]

def check():
    data=read(TAG);old=read(OLD);bind(data);bind(old)
    assert data['status']=='PASS_RAW23_FOUR_FACE_CONTACT_PREFLIGHT'
    assert len(data['cases'])==len(old['cases'])==3
    total=Counter()
    for case,o in zip(data['cases'],old['cases']):
        assert case['kinds']==o['kinds'];cone=o['initial_cone'];rays=cone['rays'];byid={r['index']:r for r in rays}
        assert len(byid)==len(rays)==50
        xyz={i:vec(r['coordinates'])for i,r in byid.items()};metric=mat(o['metric']);rows=mat(cone['normal_rows']);sc=vec(o['original_squared_scales'])
        assert len(rows)==len(sc)==12 and len(metric)==10
        contacts={i:[-q for q in vec(r['products'])]for i,r in byid.items()}
        assert all(q.sign()>=0 for p in contacts.values()for q in p)
        rcache={};gcache={}
        def rayrank(ids):
            ids=tuple(ids)
            if ids not in rcache:rcache[ids]=rank([xyz[i]for i in ids])
            return rcache[ids]
        def gram(ids):
            result=[]
            for i in ids:
                rr=[]
                for j in ids:
                    key=tuple(sorted((i,j)))
                    if key not in gcache:gcache[key]=bil(xyz[i],metric,xyz[j])
                    rr.append(gcache[key])
                result.append(rr)
            return result
        expected=list(combinations(range(12),6));assert len(case['active_decisions'])==len(expected)==924
        grouped={};fulls={}
        for rr,active in zip(case['active_decisions'],expected):
            assert rr['active']==list(active)
            arank=rank([rows[i]for i in active]);full=tuple(i for i in byid if set(active)<=set(byid[i]['zero_facets']))
            dimension=rayrank(full);found=arank==6 and dimension==4
            assert rr['active_rank']==arank and rr['all_extreme_rays']==list(full)and rr['ray_span_rank']==dimension and rr['four_face']==found
            if not found:continue
            common=tuple(sorted(set.intersection(*(set(byid[i]['zero_facets'])for i in full))))
            assert rank([rows[i]for i in common])==6 and rr['closed_active_walls']==list(common)
            grouped.setdefault(common,[]).append(list(active))
            if common in fulls:assert fulls[common]==full
            fulls[common]=full
        assert[rr['active_walls']for rr in case['faces']]==[list(k)for k in sorted(grouped)]
        supportrows={tuple(rr['rays']):rr for rr in case['stationary_supports']}
        assert len(supportrows)==len(case['stationary_supports'])
        assert list(supportrows)==sorted(supportrows,key=lambda k:(len(k),k))
        used=set();counts=Counter()
        for face in case['faces']:
            active=tuple(face['active_walls']);full=fulls[active]
            assert face['all_extreme_rays']==list(full)and face['active_bases']==grouped[active]and face['dimension']==4
            basis=next(tt for tt in combinations(full,4)if rayrank(tt)==4)
            assert face['basis_rays']==list(basis)
            h=gram(basis);sig=inertia(h);assert mat(face['plane_gram'])==h and face['plane_inertia']==sig
            signs=[byid[i]['norm_sign']for i in full];allspace=all(s>0 for s in signs);positive=sig==[4,0,0]
            assert face['ray_signs']==signs and face['all_rays_spacelike']==allspace and face['positive_definite_plane']==positive
            assert face['simplicial']==(len(full)==4)
            if not(allspace and positive):
                assert face['status']=='outside_PD_spacelike_bound'and face['support_ids']==face['contact_bounds']==face['excluding_walls']==[]
                assert not face['relative_interior_excluded'];counts[face['status']]+=1;continue
            subsets=[ids for size in range(1,5)for ids in combinations(full,size)if rayrank(ids)==size]
            assert face['support_ids']==[list(t)for t in subsets];used.update(subsets)
            fullcontacts=[[contacts[i][w]for i in full]for w in range(12)]
            nonzero=[w for w,p in enumerate(fullcontacts)if any(x!=0 for x in p)]
            assert[rr['wall']for rr in face['contact_bounds']]==nonzero
            excludes=[]
            for rr in face['contact_bounds']:
                w=rr['wall'];assert vec(rr['contact_on_all_rays'])==fullcontacts[w]
                candidates=[]
                for ids in subsets:
                    ss=supportrows[ids]['stationary'][w]
                    if ss['strictly_positive']:candidates.append((ids,Q(*ss['candidate_square'])))
                assert candidates
                maximum=candidates[0][1]
                for ids,v in candidates:
                    if(v-maximum).sign()>0:maximum=v
                winners=[list(ids)for ids,v in candidates if v==maximum]
                assert Q(*rr['maximum_square'])==maximum and rr['maximizing_supports']==winners
                below=(maximum-1).sign()<0;assert rr['strictly_below_one']==below
                if below:excludes.append(w)
            assert face['excluding_walls']==excludes and face['relative_interior_excluded']==bool(excludes)
            assert face['status']==('strict_contact_exclusion'if excludes else'PD_no_strict_contact_exclusion')
            counts[face['status']]+=1
        assert used==set(supportrows)
        for ids,rr in supportrows.items():
            assert 1<=len(ids)<=4 and rayrank(ids)==len(ids)
            h=gram(ids);assert pd(h)and mat(rr['gram'])==h
            assert[ss['wall']for ss in rr['stationary']]==list(range(12))
            for ss in rr['stationary']:
                w=ss['wall'];p=[contacts[i][w]for i in ids];z=vec(ss['coefficients'])
                assert vec(ss['contact'])==p and len(z)==len(ids)
                assert all(sum((h[i][j]*z[j]for j in range(len(ids))),Q())==p[i]for i in range(len(ids)))
                admissible=all(x.sign()>0 for x in z);assert ss['strictly_positive']==admissible
                if admissible:assert Q(*ss['candidate_square'])==2*sum((a*b for a,b in zip(p,z)),Q())/sc[w]
                else:assert ss['candidate_square']is None
        counts.update(active_six_sets=924,found_active_bases=sum(len(x)for x in grouped.values()),distinct_faces=len(grouped),
            all_spacelike_faces=sum(r['all_rays_spacelike']for r in case['faces']),PD_faces=sum(r['positive_definite_plane']for r in case['faces']),
            simplicial_faces=sum(r['simplicial']for r in case['faces']),excluded_faces=sum(r['relative_interior_excluded']for r in case['faces']),
            unique_stationary_supports=len(used),support_occurrences=sum(len(r['support_ids'])for r in case['faces']))
        assert dict(counts)==case['counts'];total.update(counts)
        print('CHECKED_CASE',case['kinds'],dict(counts),flush=True)
    assert dict(total)==data['summary']
    assert total['active_six_sets']==2772 and total['distinct_faces']==2100
    assert total['all_spacelike_faces']==total['PD_faces']==total['excluded_faces']==540
    assert total['outside_PD_spacelike_bound']==1560
    assert total['unique_stationary_supports']==19968 and total['support_occurrences']==37344
    assert data['all_original_four_faces_classified']and data['all_PD_entirely_spacelike_faces_bounded']
    assert data['arbitrary_real_candidate_coefficients']and data['arbitrary_passive_lifts_relaxed']
    assert all(data[k]is False for k in('candidate_ground_field_restriction','unknown_passive_contacts_assumed_zero',
        'native_tuple_enumeration_performed','geometric_realization_claimed','new_normal_claimed','pair_exclusion_claimed'))
    assert data['exported_columns']==[]and data['optimizer_invocations']==0
    print('PASS_RAW23_FOUR_FACE_SAVED',dict(total))
if __name__=='__main__':check()
