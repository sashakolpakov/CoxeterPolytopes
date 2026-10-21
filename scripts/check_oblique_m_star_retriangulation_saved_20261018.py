#!/usr/bin/env python3
"""Independent saved M-star polar rejection; no project modules.

Inherits the frozen physical J/L/X enclosures and original scaffold.
Reconstructs every continuous-N box and both-sign old52 square independently.
Does not reconstruct the original scaffold, normal banks, or cones.
"""
from fractions import Fraction as F
from hashlib import sha256
from itertools import permutations,combinations
from collections import Counter
from math import isqrt
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_oblique_m_star_retriangulation_interval_20261018'
PRIOR='fresh_oblique_low_s_continuous_n_interval_20261018'
SCALE=10**70


def iv(x):return tuple(map(F,x))if isinstance(x,(list,tuple))else(F(x),F(x))
def rounding(v):
    a,b=v;assert a<=b
    return F(a.numerator*SCALE//a.denominator,SCALE),F(-((-b.numerator*SCALE)//b.denominator),SCALE)
def add(a,b):return rounding((a[0]+b[0],a[1]+b[1]))
def neg(a):return -a[1],-a[0]
def sub(a,b):return add(a,neg(b))
def mul(a,b):
    values=[x*y for x in a for y in b]
    return rounding((min(values),max(values)))
def div(a,b):
    assert b[0]*b[1]>0
    return mul(a,(1/b[1],1/b[0]))
def sqrt(a):
    assert a[0]>=0
    lo=isqrt(a[0].numerator*SCALE*SCALE//a[0].denominator)
    hi=isqrt(a[1].numerator*SCALE*SCALE//a[1].denominator)
    if F(hi,SCALE)**2<a[1]:hi+=1
    return F(lo,SCALE),F(hi,SCALE)
def total(values):
    out=iv(0)
    for v in values:out=add(out,v)
    return out
def dot(x,g,y):return total(mul(mul(a,g[i][j]),b)for i,a in enumerate(x)for j,b in enumerate(y))
def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return{str(k):enc(v)for k,v in x.items()}
    if isinstance(x,(list,tuple)):return[enc(v)for v in x]
    return x
def determinant(g):
    terms=[]
    for p in permutations(range(len(g))):
        v=iv((-1)**sum(p[i]>p[j]for i in range(len(g))for j in range(i+1,len(g))))
        for i in range(len(g)):v=mul(v,g[i][p[i]])
        terms.append(v)
    return total(terms)
def inverse(g):
    n=len(g);d=determinant(g);assert d[0]*d[1]>0
    return[[div(mul(iv((-1)**(i+j)),determinant([[g[a][b]for b in range(n)if b!=i]
        for a in range(n)if a!=j])),d)for j in range(n)]for i in range(n)],d


def restore(g,step):
    w=step['wall'];out={int(j):iv(z)for j,z in step['new_products'].items()}
    assert set(out)==set(g)
    for j,z in out.items():g[j][w]=z
    g[w]={**out,w:iv(step['projected_square'])}


def ngram(g,n,a):
    frame=[10,67,69];inv,det=inverse([[g[i][j]for j in frame]for i in frame]);assert det[1]<0
    cv=[g[i][64]for i in frame];rhs=[iv(-a),iv(0),neg(n)]
    sig=sub(g[64][64],dot(cv,inv,cv));rho=sub(iv(2),dot(rhs,inv,rhs));assert sig[0]>0 and rho[0]>0
    fac=sqrt(div(rho,sig));products={};other=None
    for j in g:
        cj=[g[i][j]for i in frame];linear=dot(rhs,inv,cj);residual=mul(fac,sub(g[64][j],dot(cv,inv,cj)))
        products[j]=rhs[frame.index(j)]if j in frame else add(linear,residual)
        if j==8:other=sub(linear,residual)
    assert other[0]>0
    out={i:{**r,73:products[i]}for i,r in g.items()};out[73]={**products,73:iv(2)}
    return out,dict(frame=frame,inverse=inv,determinant=det,domain=n,N10_magnitude=a,witness=64,
        witness_residual=sig,norm_remainder=rho,opposite_N_sheet_old8=other,new_products=products)


def theta(m):
    if m==6:return sqrt(iv(3))
    assert m in(16,17,19,20)
    def atan(x):
        s=sum((-1)**k*x**(2*k+1)/F(2*k+1)for k in range(100))
        return s,s+x**201/F(201)
    pi=sub(mul(iv(16),atan(F(1,5))),mul(iv(4),atan(F(1,239))))
    x=div(pi,iv(m));xx=mul(x,x);term=value=iv(1)
    for k in range(1,60):
        term=div(mul(neg(xx),term),iv((2*k-1)*(2*k)));value=add(value,term)
    error=xx[1]*max(abs(q)for q in term)/F(119*120)
    return mul(iv(2),(value[0]-error,value[1]+error))


def polar(g,n,a,side):
    gn,np=ngram(g,n,a);frame=[69,73,68]
    gram=[[gn[i][j]for j in frame]for i in frame];inv,det=inverse(gram)
    assert gram[0][0][0]>0 and det[1]<0
    c=[gn[i][52]for i in frame];res=sub(gn[52][52],dot(c,inv,c));square=mul(iv(2),res)
    assert (1<square[0]<=square[1]<2)if side=='below_gap'else(0<square[0]<=square[1]<1)
    return dict(N=np,frame=frame,Gram=gram,inverse=inv,determinant=det,old52_covector=c,
        old52_projected_square=gn[52][52],old52_residual_square=res,
        both_R_signs_old52_product_square=square,rejection=side)


def check_topology(t,old):
    assert(t['V'],t['E'],t['F'],t['euler'])==(9,18,10,1)
    removed={tuple(sorted(x))for x in[[68,69,71],[69,70,71],[70,71,72],[68,71,72]]}
    assert {tuple(sorted(x))for x in t['removed_triangles']}==removed
    expected={tuple(sorted(x))for x in old['triangles']}-removed|{(68,69,70),(68,70,72)}
    tri=t['triangles'];assert len(tri)==len(expected)==10 and {tuple(sorted(x))for x in tri}==expected
    vs=sorted({v for f in tri for v in f});assert vs==t['vertices']and len(vs)==9 and 71 not in vs
    boundary=[8,10,73,70,72,52];assert t['boundary']==boundary
    edge=Counter(tuple(sorted(e))for f in tri for e in combinations(f,2))
    assert len(edge)==18 and set(edge.values())=={1,2}
    assert {e for e,c in edge.items()if c==1}=={tuple(sorted(e))for e in zip(boundary,boundary[1:]+boundary[:1])}
    for v in vs:
        le=[sorted(set(f)-{v})for f in tri if v in f];assert t['links'][str(v)]==le
        adj={i:[]for e in le for i in e}
        for a,b in le:adj[a].append(b);adj[b].append(a)
        seen=set();todo=[next(iter(adj))]
        while todo:
            x=todo.pop()
            if x in seen:continue
            seen.add(x);todo+=adj[x]
        assert seen==set(adj)
        degree=sorted(map(len,adj.values()))
        assert degree==([1,1]+[2]*(len(degree)-2)if v in boundary else[2]*len(degree))


def check():
    d=json.loads((ROOT/f'results/{TAG}.json').read_text());p=json.loads((ROOT/f'results/{PRIOR}.json').read_text())
    for doc in(d,p):
        for path,h in doc['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==h,path
    assert d['status']=='PASS_M_STAR_RETRIANGULATION_NATIVE52_POLAR_PREFLIGHT'and d['interval_scale']==SCALE
    assert d['deleted_original_walls']==p['deleted_walls']and d['exact_deletion_debit']==p['exact_deletion_debit']
    assert d['removed_new_wall']==71 and d['retained_new_walls_before_R']==[67,68,69,73]
    check_topology(d['topology'],p['topology'])
    profiles=[[iv(x)for x in row]for row in p['profiles']]
    assert len(profiles)==7 and all(len(row)==6 for row in profiles)
    assert profiles[0]==[iv(0)]*6
    assert all(row[:4]==[iv(0)]*4 and sum(x!=iv(0)for x in row[4:])==1
        and all(x==iv(0)or x[1]<0 for x in row[4:])for row in profiles[1:])
    base={int(i):{int(j):iv(z)for j,z in r.items()}for i,r in p['refined_primary']['projected_Gram'].items()}
    parents={};domains=[]
    assert len(p['records'])==2 and sum(len(x['origins'])for x in p['records'])==4
    for row in p['records']:
        g={i:r.copy()for i,r in base.items()};steps=row['refined_parent_steps'][:3]
        assert [s['wall']for s in steps]==[67,68,69]
        for step in steps:restore(g,step)
        assert 70 not in g and 71 not in g
        assert steps[0]['profile']==steps[2]['profile']==0 and steps[1]['profile']in(1,4)
        assert steps[2]['frame']==[8,67,68]and steps[2]['labels']==[2,2,5]
        lp=[iv(x)for x in p['profiles'][steps[1]['profile']]]
        assert lp[:4]==[iv(0)]*4 and sorted(lp[4:])==[iv(-1),iv(0)]
        assert g[68][68][0]<=F(4,3)<=g[68][68][1]and g[69][69]==iv(2)
        minor=sub(mul(g[68][68],g[69][69]),mul(g[68][69],g[69][68]));assert minor[0]>0
        duty=dict(L_profile=steps[1]['profile'],L_white_covector=lp,H4_path_labels=[3,3,5],
            projected_LX_minor=minor,forced_R_profile=0,forced_R_frame=[69,73,68],forced_R_products=[0,0,0])
        parents[row['index']]=g
        for a,(k,l)in{0:(19,20),1:(16,17)}.items():
            left,right=theta(k),theta(l);assert left[1]<right[0]<2
            ranges=[('below_gap',(theta(6)[0],left[1])),('above_gap',(right[0],F(2)))]
            domains.append(dict(parent=row['index'],N10_magnitude=a,gap_labels=[k,l],gap_left_theta=left,
                gap_right_theta=right,ranges=ranges,duty=duty))
    assert enc(domains)==d['domains']and len(domains)==4
    assert len(d['records'])==724 and [r['index']for r in d['records']]==list(range(724))
    squares=[]
    for r in d['records']:
        assert r['parent']in parents and r['N10_magnitude']in(0,1)and r['side']in('below_gap','above_gap')
        n=iv(r['domain']);assert n[0]<n[1]
        proof=polar(parents[r['parent']],n,r['N10_magnitude'],r['side'])
        assert enc(proof)==r['proof'];squares.append(proof['both_R_signs_old52_product_square'])
    for r in domains:
        for side,(a,b)in r['ranges']:
            parts=sorted(iv(x['domain'])for x in d['records']if x['parent']==r['parent']and x['N10_magnitude']==r['N10_magnitude']and x['side']==side)
            assert parts[0][0]==a and parts[-1][1]==b and all(x[1]==y[0]for x,y in zip(parts,parts[1:]))
    s=d['summary'];fixed=dict(parents=2,parent_origins=4,N_families=2,range_domains=8,attempted_boxes=1440,
        certified_boxes=724,unresolved=0,accepted_R=0,S_banks=0,new_cones=0)
    assert all(s[k]==v for k,v in fixed.items())
    assert s['minimum_square']==str(min(x[0]for x in squares))and s['maximum_square']==str(max(x[1]for x in squares))
    assert s['minimum_unit_gap']==str(min(x[0]-1 if r['side']=='below_gap'else 1-x[1]for r,x in zip(d['records'],squares)))
    print('PASS_SAVED_M_STAR_RETRIANGULATION_POLAR',len(d['records']),flush=True)


if __name__=='__main__':check()

