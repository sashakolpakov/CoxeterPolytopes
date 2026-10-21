#!/usr/bin/env python3
"""Independent saved polar interval arithmetic; no project modules.

Inherits the frozen physical parent/S enclosures and prior native rejections.
Reconstructs every continuous-N box and both polar arguments independently.
Does not reconstruct the original scaffold, normal banks, or cones.
"""
from fractions import Fraction as F
from hashlib import sha256
from itertools import permutations
from math import isqrt
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_oblique_low_s_high_r_polar_interval_20261018'
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


def rm(g):
    frame=[69,73,72];gram=[[g[i][j]for j in frame]for i in frame];inv,det=inverse(gram)
    assert gram[0][0][0]>0 and det[1]<0
    v=[g[i][71]for i in frame];assert v[0]==v[2]==iv(0) and g[71][71]==iv(2)
    cofactor=sub(mul(gram[0][0],gram[2][2]),mul(gram[0][2],gram[2][0]))
    margin=sub(mul(sub(iv(3),mul(iv(2),g[71][71])),neg(det)),mul(iv(2),mul(mul(v[1],v[1]),cofactor)))
    assert margin[0]>0
    return dict(frame=frame,Gram=gram,inverse=inv,determinant=det,M_covector=v,middle_cofactor=cofactor,
        strict_adjugate_margin=margin,RM_square=mul(iv(2),sub(g[71][71],dot(v,inv,v))))


def rs(g):
    frame=[69,71,73];gram=[[g[i][j]for j in frame]for i in frame];inv,det=inverse(gram)
    assert gram[0][0][0]>0 and det[1]<0 and g[72][72]==iv(2) and g[72][71]==iv(0)
    w=[g[i][8]for i in frame];s=[g[i][72]for i in frame]
    q=sub(g[8][8],dot(w,inv,w));assert q[0]>0
    positive=sqrt(mul(iv(2),q));assert positive[0]>0
    sr=neg(mul(sqrt(div(iv(2),q)),sub(g[8][72],dot(w,inv,s))))
    kind='positive'if sr[0]>0 else'strictly_dotted'
    assert sr[0]>0 or sr[1]<-2
    return dict(frame=frame,Gram=gram,inverse=inv,determinant=det,old8_covector=w,S_covector=s,
        old8_residual_square=q,positive_R_sheet_old8=positive,negative_R_sheet_SR=sr,rejection=kind)


def check():
    d=json.loads((ROOT/f'results/{TAG}.json').read_text());p=json.loads((ROOT/f'results/{PRIOR}.json').read_text())
    for doc in(d,p):
        for path,h in doc['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==h,path
    assert d['status']=='PASS_PHYSICAL_LOW_S_ALL_HIGH_RM_RS_CONTINUOUS_N_POLAR_EXCLUSION'and d['interval_scale']==SCALE
    assert d['parameter_domain']==enc([sqrt(iv(3))[0],F(2)])
    assert d['N_families']==[[2,2,'m>=6'],[3,2,'m>=6']]
    assert d['exact_deletion_debit']==p['exact_deletion_debit']and d['topology']==p['topology']and d['deleted_walls']==p['deleted_walls']
    assert p['summary']['S_recipes']==84 and p['summary']['S_branches']==168
    assert p['summary']['S_strict_rejections']==92 and p['summary']['S_overdomain_survivors']==76
    base={int(i):{int(j):iv(z)for j,z in r.items()}for i,r in p['refined_primary']['projected_Gram'].items()}
    parents={};selected=[]
    assert len(p['records'])==2 and [r['index']for r in p['records']]==[0,1]
    for parent in p['records']:
        g={i:r.copy()for i,r in base.items()}
        assert len(parent['refined_parent_steps'])==5
        for step in parent['refined_parent_steps']:restore(g,step)
        del g[70]
        for row in g.values():del row[70]
        ss=[s for s in parent['S']if s['blocker']is None and s['labels'][2]==2]
        assert len(ss)==12 and sum(s['profile']==0 for s in ss)==7
        parents[parent['index']]=(g,ss)
        selected.append(dict(parent=parent['index'],S=[[s['recipe'],s['branch'],s['profile'],s['labels']]for s in ss]))
    assert selected==d['selected_S']
    assert len(d['records'])==16 and [r['index']for r in d['records']]==list(range(16))
    rms=[];rss=[]
    for row in d['records']:
        g,ss=parents[row['parent']];a=row['N10_magnitude'];assert a in(0,1)
        assert len(row['S'])==12 and iv(row['domain'])[0]<iv(row['domain'])[1]
        for s,proof in zip(ss,row['S']):
            gs={i:r.copy()for i,r in g.items()};restore(gs,s['extension'])
            assert gs[69][71]==gs[72][71]==iv(0)
            assert gs[69][69]==gs[71][71]==iv(2)
            gn,np=ngram(gs,iv(row['domain']),a);rp=rm(gn);sp=rs(gn)if s['profile']==0 else None
            expected=dict(S_recipe=s['recipe'],S_branch=s['branch'],S_profile=s['profile'],S_labels=s['labels'],N=np,high_RM=rp,high_RS=sp)
            assert enc(expected)==proof
            rms.append(rp)
            if sp:rss.append(sp)
    for parent in range(2):
        for a in(0,1):
            parts=sorted(iv(r['domain'])for r in d['records']if r['parent']==parent and r['N10_magnitude']==a)
            assert parts[0][0]==sqrt(iv(3))[0] and parts[-1][1]==2
            assert all(x[1]==y[0]for x,y in zip(parts,parts[1:]))
    summary=d['summary']
    fixed=dict(parents=2,parent_origins=4,low_S_recipes=84,low_S_branches=168,inherited_strict_S_rejections=92,
        inherited_S_overdomain=76,SM_right_survivors=24,SM_right_profile0_survivors=14,
        high_RM_families=48,high_RS_families=28,covering_boxes=16,high_RM_box_proofs=192,high_RS_box_proofs=112,
        unresolved=0,accepted=0,new_cones=0)
    assert all(summary[k]==v for k,v in fixed.items())and len(rms)==192 and len(rss)==112
    bounds=dict(minimum_high_RM_adjugate_margin=min(r['strict_adjugate_margin'][0]for r in rms),
        maximum_high_RM_determinant=max(r['determinant'][1]for r in rms),maximum_high_RS_determinant=max(r['determinant'][1]for r in rss),
        minimum_positive_R_old8=min(r['positive_R_sheet_old8'][0]for r in rss),
        minimum_positive_SR=min(r['negative_R_sheet_SR'][0]for r in rss if r['rejection']=='positive'),
        maximum_dotted_SR=max(r['negative_R_sheet_SR'][1]for r in rss if r['rejection']=='strictly_dotted'))
    assert all(summary[k]==str(v)for k,v in bounds.items())
    print('PASS_SAVED_LOW_S_HIGH_R_POLAR_INTERVAL',len(d['records']),len(rms),len(rss),flush=True)


if __name__=='__main__':check()
