#!/usr/bin/env python3
"""Exact incidence-graph replacement of the central pair in the H22 search.

Integer/Fraction core verification; the optional full-bank check uses checked
int64 products and inherits the previously certified bank's internal pairs.
"""
import argparse
from collections import Counter, deque
from fractions import Fraction as F
import hashlib, json
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = 'data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
ESCAPE = 'results/fresh_h22_new_cusp_completion_20261013_lower_escape.json'
OUT = 'results/fresh_h22_incidence_core_20261013.json'
def read(p): return json.loads((ROOT/p).read_text())
def sha(p): return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def dot(a,b): return sum(x*y for x,y in zip(a,b,strict=True))
def mv(M,v): return [dot(r,v) for r in M]
def pair(Q,a,b): return dot(a,mv(Q,b))
def inertia(M):
    A = [[F(x) for x in r] for r in M]; counts = [0,0,0]
    while A:
        n=len(A); k=next((i for i in range(n) if A[i][i]),None)
        if k is None:
            ij=next(((i,j) for i in range(n) for j in range(i+1,n) if A[i][j]),None)
            if ij is None: counts[2]+=n; break
            i,j=ij
            for a in range(n): A[i][a]+=A[j][a]
            for a in range(n): A[a][i]+=A[a][j]
            k=i
        A[0],A[k]=A[k],A[0]
        for r in A: r[0],r[k]=r[k],r[0]
        p=A[0][0]; counts[0 if p>0 else 1]+=1
        A=[[A[i][j]-A[i][0]*A[0][j]/p for j in range(1,n)] for i in range(1,n)]
    return counts
def classify(p,a,b):
    assert p<=0
    if p==0: return 'right'
    t=F(p*p,a*b)
    if t in (F(1,4),F(1,2),F(3,4)):
        return {F(1,4):'pi/3',F(1,2):'pi/4',F(3,4):'pi/6'}[t]
    assert t>=1, (p,a,b,t)
    return 'parallel' if t==1 else 'dotted'
def core():
    raw=read(SOURCE)['coxeter_input']
    Q=[[-int(x) for x in r] for r in raw['basis_gram']]
    W=[[int(x) for x in r] for r in raw['wall_coordinates_as_columns_transposed']]
    cov=[mv(Q,r) for r in W]; G=[[dot(r,c) for c in cov] for r in W]
    assert inertia(Q)==[21,1,0]
    S=[i for i in range(210) if G[i][i]==2]
    assert len(S)==42 and Counter(G[i][i] for i in range(210))=={2:42,4:168}
    oldcounts=Counter(classify(G[i][j],G[i][i],G[j][j]) for i,j in combinations(range(210),2))
    color={S[0]:0}; todo=deque(color)
    while todo:
        i=todo.popleft()
        for j in S:
            if i==j: continue
            assert G[i][j] in (0,-1)
            if G[i][j]==-1:
                if j not in color: color[j]=1-color[i]; todo.append(j)
                assert color[j]!=color[i]
    assert len(color)==42
    parts=[[i for i in S if color[i]==c] for c in (0,1)]
    assert list(map(len,parts))==[21,21]
    A=[[-G[i][j] for j in parts[1]] for i in parts[0]]
    assert all(sum(r)==5 for r in A)
    assert all(dot(A[i],A[j])==4*(i==j)+1 for i in range(21) for j in range(21))
    a,b=[[sum(W[i][j] for i in part) for j in range(22)] for part in parts]
    ns=[[2*x+5*y for x,y in zip(a,b)],[5*x+2*y for x,y in zip(a,b)]]
    assert all(x%21==0 for r in ns for x in r)
    V=[[x//21 for x in r] for r in ns]; h=[x+y for x,y in zip(*V)]
    assert all(3*h[j]==a[j]+b[j] for j in range(22))
    assert pair(Q,h,h)==-14 and [pair(Q,v,v) for v in V]==[-2,-2]
    assert pair(Q,*V)==-5
    R=[v+[z] for v in V for z in (2,-2)]
    Qn=[r+[0] for r in Q]+[[0]*22+[1]]
    assert inertia(Qn)==[22,1,0]
    assert all(pair(Qn,r,r)==2 for r in R)
    newcounts=Counter(classify(pair(Qn,r,w+[0]),2,G[i][i]) for r in R for i,w in enumerate(W))
    newcounts.update(classify(pair(Qn,r,s),2,2) for r,s in combinations(R,2))
    assert newcounts=={'right':84,'pi/3':86,'dotted':676}
    assert max(pair(Q,h,w) for w in W)<0 and [pair(Q,h,v) for v in V]==[-7,-7]
    # The former m=3 central normals are sums of same-sign new normals.
    assert all(R[0][j]+R[2][j]==(h+[4])[j] for j in range(23))
    assert all(R[1][j]+R[3][j]==(h+[-4])[j] for j in range(23))
    # For every previous center height T>=sqrt(14.5), -7+2T>0.
    assert 4*F(29,2)>49
    s=read(ESCAPE); X=s['numerator']; den=s['denominator']
    xp=[pair(Qn,X,r) for r in R]; assert max(xp)<0
    assert F(pair(Qn,X,X),den**2)==F(s['square'])>0
    assert pair(Qn,X,h+[0])<0 and max(pair(Qn,X,w+[0]) for w in W)<0
    return {'status':'PASS_EXACT_H22_INCIDENCE_CORE_REPLACEMENT',
        'sources':{p:sha(p) for p in (SOURCE,ESCAPE)},
        'incidence_parts':parts,'incidence_identity':'A A^t = 4 I + J; row weight 5',
        'projection_formula':'v=(2 sum(P)+5 sum(L))/21; vprime=(5 sum(P)+2 sum(L))/21',
        'old_projections':V,'new_normals':R,'center':h,'ambient_signature':[22,1,0],
        'old_pair_counts':dict(oldcounts),'new_pair_counts':dict(newcounts),
        'core_facets':214,'actual_facets':'All pairings are nonpositive and the center is strict timelike. Projection of the center onto each wall remains strict on every other wall.',
        'old_m3_center_redundancy':'(h,+/-4) is the sum of the two new same-sign normals.',
        'old_center_incompatibility':'Each old (h,+/-T), sqrt(14.5)<T<=4, has positive pairing -7+2T with a same-sign new wall.',
        'surviving_escape':{'numerator':X,'denominator':den,'norm':s['square'],'new_pair_numerators':xp},
        'scope':'Exact 214-facet Coxeter chamber, infinite volume. This construction is new to this search; no literature novelty or finite-volume claim.'},Q,W,V
def full_bank(Q,V,new_core_facets=214):
    import numpy as np
    from fresh_h22_new_cusp_completion_20261013_verify import setup
    from fresh_h22_new_cusp_completion_20261013_tight import shifted,DEN
    lower=read('results/fresh_h22_new_cusp_completion_20261013_lower.json')
    assert sha(lower['source'])==lower['source_sha256']
    inherited=read(lower['source'])
    assert all(sha(p)==value for p,value in inherited['sources'].items())
    prior,qn,W,h,meta,gens,D,Vt,B,A,ids,sizes=setup()
    cov=np.array(Q,dtype=np.int64)@np.array(V,dtype=np.int64).T
    bounds={}; sources={}
    def products(rows):
        assert 22*int(abs(rows).max())*int(abs(cov).max())<2**63
        return rows@cov
    for name,key,den in [('symmetry_orbit','twice_old_projections',2),('cusp_orbits','old_null_projections',1),('gap_face_orbit','six_old_projections',6)]:
        path='data/enumerations/fresh_h22_'+name+'_20261013.npz'; sources[path]=sha(path)
        z=np.load(ROOT/path); rows=z[key]; p=products(rows)
        if name=='symmetry_orbit':
            m=F(int(p.max()),den); assert m+F(13,4)<-F(3,2)
            # sqrt(10)<13/4 and sqrt(2)<3/2.
            bounds[name]={'max_old_pair':str(m),'signed_new_pair_bound':'max_old_pair+13/4<-3/2<-sqrt(2)','count':2*len(rows)}
        elif name=='gap_face_orbit':
            m=F(int(p.max()),den); assert m+7<-2
            # 2sqrt(11)<7 and sqrt(4/3)<2.
            bounds[name]={'max_old_pair':str(m),'signed_new_pair_bound':'max_old_pair+7<-2<-sqrt(4/3)','count':2*len(rows)}
        else:
            ds={}; counts=Counter()
            for t in np.unique(z['four_heights']):
                height=F(int(t),4); m=F(int(p[z['four_heights']==t].max()))
                vals,amounts=np.unique(p[z['four_heights']==t],return_counts=True)
                for value,amount in zip(vals,amounts):
                    for sign in (-1,1):
                        counts[classify(F(int(value))+sign*2*height,2,height*height)]+=2*int(amount)
                ds[str(height)]=str(m)
            bounds[name]={'height_to_max_old_pair':ds,'pair_counts':dict(counts),'count':2*len(rows)}
    repair={}
    for name,rows,t,ix,index in shifted(D,Vt,B,A,ids,sizes,qn,True):
        p=products(rows)
        norms=np.einsum('ij,ij->i',A,A@qn) if index==0 else np.full(len(rows),[0,279/14,279/14,55/62496,55/62496][index])
        groups={}
        for height in np.unique(t):
            for norm in np.unique(norms[t==height]):
                mask=(t==height)&(norms==norm); m=F(int(p[mask].max()),DEN)
                # New coordinate is t*sqrt(11)/10; cap height is 2.
                # Bound sqrt(11)<10/3, then compare rational squares.
                upper=m+F(2*int(height),3)
                q=F(int(norm)) if index==0 else [None,F(279,14),F(279,14),F(55,62496),F(55,62496)][index]
                assert upper<0 and upper*upper>2*q
                groups[f't={int(height)},norm={q}']={'max_old_pair':str(m),'signed_pair_upper':str(upper)}
        repair[name]={'groups':groups,'count':2*len(rows)}
    for p in ('results/fresh_h22_new_cusp_completion_20261013_lower.json','results/fresh_h22_new_cusp_completion_20261013_lower_escape.json'):
        sources[p]=sha(p)
    total=new_core_facets+sum(x['count'] for x in bounds.values())+sum(x['count'] for x in repair.values())
    assert total==5397336+new_core_facets
    return {'status':'PASS_ALL_NEW_TO_RETAINED_BANK_PAIRS_COXETER','sources':sources,'old_bank_bounds':bounds,'repair_bank_bounds':repair,'facets':total,
        'scope':'Checks the supplied norm-two height-two projections against every retained bank family. Retained-retained compatibility is inherited through the checked manifests. This helper makes no volume assertion and does not check an escape against the supplied projections.'}
if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--verify',action='store_true'); p.add_argument('--full-bank',action='store_true'); a=p.parse_args()
    out,Q,W,V=core()
    if a.full_bank: out['full_bank']=full_bank(Q,V)
    target=OUT.replace('.json','_full.json') if a.full_bank else OUT
    if a.verify: assert out==read(target)
    else: (ROOT/target).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out.get('full_bank',{}).get('facets',out['core_facets']))
