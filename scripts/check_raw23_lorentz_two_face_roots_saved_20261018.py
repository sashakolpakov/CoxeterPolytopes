#!/usr/bin/env python3
"""Saved mixed-face roots, native contacts and exact infinite-volume escapes.

No producer imports. Uses the frozen standalone polar helper's Q(sqrt(5))
arithmetic and reverse-column cone checker. A separate reversed-radicand
multiquadratic implementation rechecks every finite-anchor tuple. Historical
ray completeness, native field catalogue and actual-source incidence remain
bound inputs; changed fourteen-wall cones are fully checked here.
"""
import json
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations,product
from math import isqrt
import check_raw23_union_polar_saved_20261018 as G
Q=G.Q;ROOT=G.ROOT
TAG='fresh_raw23_lorentz_two_face_native_roots_20261018'
def read(s):return json.loads((ROOT/f'results/{s}.json').read_text())
def vec(xs):return[Q(*x)for x in xs]
def add(x,y):return[a+b for a,b in zip(x,y)]
def scale(x,q):return[a*q for a in x]
def products(rows,x):return[sum((a*b for a,b in zip(row,x)),Q())for row in rows]
def rank(rows):
    a=[r[:]for r in rows];n=0
    if not a:return 0
    for j in reversed(range(len(a[0]))):
        k=next((k for k in range(n,len(a))if a[k][j]!=0),None)
        if k is None:continue
        a[n],a[k]=a[k],a[n];v=a[n][j];a[n]=[x/v for x in a[n]]
        for k in range(n+1,len(a)):
            v=a[k][j];a[k]=[x-v*y for x,y in zip(a[k],a[n])]
        n+=1
        if n==len(a):break
    return n
def inertia(g):
    a=[r[:]for r in g];pos=neg=zero=0
    while a:
        n=len(a);k=next((i for i in range(n)if a[i][i]!=0),None)
        if k is not None:
            ids=[k]+[i for i in range(n)if i!=k];a=[[a[i][j]for j in ids]for i in ids]
            s=a[0][0];pos+=s.sign()>0;neg+=s.sign()<0
            a=[[a[i][j]-a[i][0]*a[0][j]/s for j in range(1,n)]for i in range(1,n)]
        else:
            pair=next(((i,j)for i,j in combinations(range(n),2)if a[i][j]!=0),None)
            if pair is None:zero+=n;break
            i,j=pair;ids=[i,j]+[k for k in range(n)if k not in pair];a=[[a[i][j]for j in ids]for i in ids];z=a[0][1]
            pos+=1;neg+=1
            a=[[a[i][j]-(a[i][0]*a[1][j]+a[i][1]*a[0][j])/z for j in range(2,n)]for i in range(2,n)]
    return[pos,neg,zero]
def ratroot(q):
    if q<0:return None
    a,b=isqrt(q.numerator),isqrt(q.denominator)
    return Fraction(a,b)if a*a==q.numerator and b*b==q.denominator else None
def root(q):
    if q.sign()<0:return None
    if q==0:return Q()
    if q.b==0:
        z=ratroot(q.a/5)
        if z is not None:return Q(0,z)
    d=ratroot(q.a*q.a-5*q.b*q.b)
    if d is None:return None
    for u in((q.a-d)/2,(q.a+d)/2):
        z=ratroot(u)
        if z is None or z==0:continue
        w=Q(z,q.b/(2*z))
        if w*w==q:return w if w.sign()>0 else-w
    return None
class Extension:
    def __init__(self,radicands):
        self.basis=[];self.powers=[Q(1)]
        for r in reversed(radicands):
            assert r.sign()>0
            if any(root(r/p)is not None for p in self.powers):continue
            self.basis.append(r);self.powers +=[p*r for p in self.powers]
        self.n=len(self.powers);assert self.n<=4
    def scalar(self,z):return[z if isinstance(z,Q)else Q(z)]+[Q()]*(self.n-1)
    def sqrt(self,q):
        for i,p in enumerate(self.powers):
            z=root(q/p)
            if z is not None:
                x=self.scalar(0);x[i]=z;return x
        raise AssertionError('missing square class')
    def mul(self,x,y):
        out=self.scalar(0)
        for i,a in enumerate(x):
            for j,b in enumerate(y):out[i^j]=out[i^j]+a*b*self.powers[i&j]
        return out
    def sign(self,x):
        def recurse(z,k):
            if k==0:return z[0].sign()
            half=1<<(k-1);a,b=z[:half],z[half:];sa,sb=recurse(a,k-1),recurse(b,k-1)
            if not sa:return sb
            if not sb or sa==sb:return sa
            aa=self.mul(a+[Q()]*(self.n-half),a+[Q()]*(self.n-half))[:half]
            bb=self.mul(b+[Q()]*(self.n-half),b+[Q()]*(self.n-half))[:half]
            return sa*recurse([u-self.basis[k-1]*v for u,v in zip(aa,bb)],k-1)
        return recurse(x,len(self.basis))
FINITE=[Q(1),Q(2),Q('3/2','1/2'),Q(3),Q('5/2','1/2')]

def check_faces(c,o):
    rays=o['initial_cone']['rays'];rows=G.mat(o['initial_cone']['normal_rows']);g=G.mat(o['metric'])
    assert len(rays)==50 and len(c['two_faces'])==225
    by={tuple(f['rays']):f for f in c['two_faces']};assert len(by)==225;observed=set()
    for i,j in combinations(range(50),2):
        common=sorted(set(rays[i]['zero_facets'])&set(rays[j]['zero_facets']))
        if len(common)<8 or rank([rows[k]for k in common])!=8:continue
        full=[k for k,r in enumerate(rays)if set(common)<=set(r['zero_facets'])];assert full==[i,j]
        f=by[i,j];observed.add((i,j));assert f['active']==common
        xyz=[vec(rays[k]['coordinates'])for k in(i,j)]
        gg=[[G.dot(x,g,y)for y in xyz]for x in xyz];assert gg==G.mat(f['gram'])
        pp=[[-x for x in vec(rays[k]['products'])]for k in(i,j)]
        assert pp==G.mat(f['contacts'])and f['nonzero']==[k for k in range(12)if any(r[k]!=0 for r in pp)]
        assert f['signs']==[rays[k]['norm_sign']for k in(i,j)]and f['inertia']==inertia(gg)
    assert observed==set(by)
    assert Counter(tuple(f['signs'])for f in by.values())=={(1,1):112,(-1,1):64,(-1,-1):49}
    center=vec(o['initial_cone']['strict_timelike_center'])
    assert all(G.dot(vec(r['coordinates']),g,center).sign()<0 for r in rays if r['norm_sign']<0)
    return by

def check_anchors(c,o,faces):
    saved=c['finite_anchor_domain'];tests={(tuple(x['rays']),tuple(x['anchors'])):x for x in saved['tests']}
    assert len(tests)==len(saved['tests'])==320
    survivors={(tuple(x['rays']),tuple(x['anchors']),tuple(x['finite_square_indices'])):x for x in saved['survivors']}
    assert len(survivors)==len(saved['survivors'])
    sc=vec(o['original_squared_scales']);counts=Counter();fields={};observed=set();found=set()
    for ids,f in faces.items():
        if f['signs']!=[-1,1]:continue
        p=G.mat(f['contacts']);h=G.mat(f['gram']);assert len(f['nonzero'])==4;paircount=0
        for ia,ib in combinations(f['nonzero'],2):
            det=p[0][ia]*p[1][ib]-p[0][ib]*p[1][ia]
            if det==0:continue
            paircount+=1;key=(ids,(ia,ib));observed.add(key);rec=tests[key];assert len(rec['decisions'])==25
            for n,(fi,fj)in enumerate(product(range(5),repeat=2)):
                qa,qb=sc[ia]*FINITE[fi],sc[ib]*FINITE[fj];fk=(qa.a,qa.b,qb.a,qb.b)
                if fk not in fields:fields[fk]=Extension([qa,qb])
                k=fields[fk];ca,cb=k.sqrt(qa),k.sqrt(qb)
                a=scale(add(scale(ca,p[1][ib]),scale(cb,-p[1][ia])),Q(1)/det)
                b=scale(add(scale(cb,p[0][ia]),scale(ca,-p[0][ib])),Q(1)/det)
                counts['tuples']+=1
                if k.sign(a)<=0 or k.sign(b)<=0:status='outside_relative_interior'
                else:
                    norm=add(add(scale(k.mul(a,a),h[0][0]),scale(k.mul(b,b),h[1][1])),scale(k.mul(a,b),2*h[0][1]))
                    if k.sign(norm)<=0:status='causal_projection'
                    elif k.sign(add(norm,k.scalar(-2)))>0:status='norm_above2'
                    else:
                        contact=[add(scale(a,p[0][j]),scale(b,p[1][j]))for j in range(12)]
                        square=[scale(k.mul(v,v),Q(1)/s)for v,s in zip(contact,sc)]
                        if any(any(v!=0 for v in z[1:])for z in square):status='nonbase_two_cycle'
                        elif any(G.gate(z[0])for z in square):status='native_failure'
                        else:
                            status='survivor';assert not k.basis and norm==[Q(2)]
                            sk=(ids,(ia,ib),(fi,fj));s=survivors[sk];found.add(sk)
                            assert vec(s['coefficients'])==[a[0],b[0]]and Q(*s['norm'])==2
                            assert vec(s['physical_squares'])==[z[0]for z in square]
                assert status==rec['decisions'][n];counts[status]+=1
        assert paircount==5
    assert observed==set(tests)and found==set(survivors)and dict(counts)==saved['counts']
    return counts

def cycle_check(gg,sc,records):
    by={tuple(r['path']):Q(*r['value'])for r in records};assert len(by)==len(records)==93;seen=set()
    n=len(gg);adj={i:[j for j in range(n)if i!=j and gg[i][j]!=0]for i in range(n)}
    for start in range(n):
        def dfs(path):
            i=path[-1]
            for j in adj[i]:
                if j==start and len(path)>=3 and path[1]<path[-1]:
                    z=Q(1)
                    for a,b in zip(path,path[1:]+path[:1]):z=z*gg[a][b]/sc[a]
                    key=tuple(path);assert z==by[key]
                    assert(2*z.a).denominator==1 and(z.a*z.a-5*z.b*z.b).denominator==1;seen.add(key)
                elif j>start and j not in path:dfs(path+[j])
        dfs([start])
    assert seen==set(by)

def check_roots_joint(c,o):
    g=G.mat(o['metric']);rows=G.mat(o['initial_cone']['normal_rows']);rays=o['initial_cone']['rays'];sc=vec(o['original_squared_scales'])
    roots=c['new_roots'];assert len(roots)==len(c['finite_anchor_domain']['survivors'])
    if not roots:assert c['joint']is None;return 0
    assert len(roots)==2;xs=[]
    for rec,s in zip(roots,c['finite_anchor_domain']['survivors']):
        assert rec['source']==s;a,b=vec(s['coefficients']);u,v=[vec(rays[i]['coordinates'])for i in s['rays']]
        x=add(scale(u,a),scale(v,b));xs.append(x);assert x==vec(rec['coordinates'])and G.dot(x,g,x)==2
        pp=products(rows,x);assert pp==vec(rec['raw_contacts'])and all(z.sign()<=0 for z in pp)
        assert[v*v/s for v,s in zip(pp,sc)]==vec(s['physical_squares'])
        fp=[G.dot(x,g,vec(r['coordinates']))for r in rays if r['norm_sign']<0]
        assert fp==vec(rec['finite_vertex_products'])and len(fp)==18 and all(z.sign()<0 for z in fp)
        cuts=[r['index']for r in rays if r['norm_sign']>0 and G.dot(x,g,vec(r['coordinates'])).sign()>0]
        assert cuts==rec['cut_original_spacelike_rays']and len(cuts)==32
        pole=rays[o['initial_polar_tests']['passing_rays'][0]];v=G.dot(x,g,vec(pole['coordinates']))
        assert v.sign()==rec['old_pole_product_sign']==1 and 2*v*v/Q(*pole['norm'])==Q(*rec['old_pole_physical_product'])==1
    jo=c['joint'];co=jo['cone'];oldg=G.mat(o['original_scaled_gram']);sc +=[Q(1),Q(1)]
    newrows=[[sum((g[i][j]*x[j]for j in range(10)),Q())for i in range(10)]for x in xs]
    gg=[r+[products([row],x)[0]for x in xs]for r,row in zip(oldg,rows)]
    for i,x in enumerate(xs):gg.append([gg[j][12+i]for j in range(12)]+[G.dot(x,g,y)for y in xs])
    assert gg==G.mat(jo['gram'])and sc==vec(jo['squared_scales'])and gg[12][13]==Q(-1,-1)
    assert all(gg[i][j]==gg[j][i]for i,j in product(range(14),repeat=2))and all(gg[i][i]==2*sc[i]for i in range(14))
    assert inertia(gg)==jo['physical_inertia']==[9,1,4]
    assert inertia([[Q(z.a,-z.b)for z in r]for r in gg])==jo['conjugate_inertia']==[10,0,4]
    assert len(jo['native_pair_checks'])==91
    for (i,j),r in zip(combinations(range(14),2),jo['native_pair_checks']):
        z=gg[i][j]*gg[i][j]/sc[i]/sc[j]
        assert r['walls']==[i,j]and Q(*r['physical_square'])==z and gg[i][j].sign()<=0 and G.gate(z)is None
    cycle_check(gg,sc,jo['simple_cycles']);rows+=newrows
    active=G.check_cone(co,rows,g);assert active==2002 and co['ray_counts']=={'-1':34,'0':0,'1':72}
    assert len(jo['facet_centers'])==14
    for i,r in enumerate(jo['facet_centers']):
        ids=[v['index']for v in co['rays']if i in v['zero_facets']];fin=[j for j in ids if co['rays'][j]['norm_sign']<0]
        assert r['wall']==i and r['ray_indices']==ids and r['finite_ray_indices']==fin
        assert rank([vec(co['rays'][j]['coordinates'])for j in ids])==r['span_rank']==9
        x=[sum((Q(*co['rays'][j]['coordinates'][k])for j in fin),Q())for k in range(10)]
        z=[sum((Q(*co['rays'][j]['coordinates'][k])for j in ids),Q())for k in range(10)]
        eps=Q(*r['epsilon']);assert eps.sign()>0;w=add(x,scale(z,eps));pp=products(rows,w)
        assert w==vec(r['relative_interior_point'])and pp==vec(r['products'])and G.dot(w,g,w)==Q(*r['norm'])
        assert Q(*r['norm']).sign()<0 and pp[i]==0 and all(v.sign()<0 for j,v in enumerate(pp)if j!=i)
        assert G.dot(w,g,vec(co['strict_timelike_center'])).sign()<0
    e=jo['escape'];r=co['rays'][e['ray_index']];assert r['norm_sign']>0
    assert e['ray_coordinates']==r['coordinates']and e['ray_norm']==r['norm']and e['ray_products']==r['products']
    center=vec(co['strict_timelike_center']);eps=Q(*e['epsilon']);assert eps.sign()>0
    y=add(vec(r['coordinates']),scale(center,eps));pp=products(rows,y);yn=G.dot(y,g,y)
    assert y==vec(e['strict_spacelike_point'])and pp==vec(e['strict_products'])and all(z.sign()<0 for z in pp)
    assert yn==Q(*e['strict_norm'])and yn.sign()>0
    delta=add(y,scale(center,-1));poly=[G.dot(center,g,center),2*G.dot(center,g,delta),G.dot(delta,g,delta)]
    assert poly==vec(e['segment_null_polynomial'])and poly[0]==Q(*e['center_norm'])
    assert poly[0].sign()<0 and sum(poly,Q()).sign()>0
    assert jo['infinite_hyperbolic_volume']and jo['all_fourteen_walls_are_facets']
    return active

def main():
    d=read(TAG);old=read(G.TAG);G.bind(d);G.bind(old)
    assert d['status']=='PASS_RAW23_MIXED_TWO_FACE_NATIVE_ROOTS_AND_INFINITE_JOINT_CONES'
    assert len(d['cases'])==len(old['cases'])==3;counts=Counter();active=0
    for c,o in zip(d['cases'],old['cases']):
        assert c['kinds']==o['kinds']and c['relative_index']==o['relative_index']
        ff=check_faces(c,o);counts.update(check_anchors(c,o,ff));active+=check_roots_joint(c,o)
        print('CHECKED',c['kinds'],flush=True)
    assert dict(counts)==d['counts']and counts['tuples']==d['finite_anchor_tuples']==24000 and counts['survivor']==d['new_norm2_roots']==4
    assert(d['two_face_count'],d['mixed_faces'],d['timelike_pair_faces'])==(675,192,147)
    assert(active,d['joint_cones'],d['joint_active_sets'],d['joint_spacelike_rays'],d['joint_finite_vertices'],d['joint_simple_cycles'])==(4004,2,4004,144,68,186)
    assert d['candidate_ground_field']=='Q(sqrt(5))'and d['arbitrary_passive_lift_allowed_in_search']and d['surviving_passive_norm_zero']
    assert not d['all_mixed_face_normals_classified']and not d['arbitrary_completion_excluded']and not d['compact_or_finite_volume_example_found']
    assert d['exported_columns']==[]and not d['old_matrix_changed']and d['optimizer_invocations']==0
    print('PASS_SAVED_RAW23_MIXED_ROOTS',dict(tuples=24000,new_roots=4,changed_active_sets=4004,infinite_cones=2),flush=True)

if __name__=='__main__':main()
