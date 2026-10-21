#!/usr/bin/env python3
"""Independent mixed D4/P2 audit: stdlib Fraction quartic arithmetic and support census.

No producer/project code imported. All cone rays are enumerated by support
nullspaces; the producer primarily uses incremental cone clipping.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import isqrt
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
BASE='fresh_h9_mixed_d4_coupling_20261014'
TAG=BASE+'_review'

class Q:
    __slots__ = ('a', 'b')
    def __init__(self, a=0, b=0):
        if isinstance(a,Q): self.a,self.b=a.a,a.b
        else: self.a,self.b=F(a),F(b)
    def __add__(self,y):
        y=Q(y);return Q(self.a+y.a,self.b+y.b)
    __radd__=__add__
    def __neg__(self):return Q(-self.a,-self.b)
    def __sub__(self,y):return self+-Q(y)
    def __rsub__(self,y):return Q(y)+-self
    def __mul__(self,y):
        y=Q(y);return Q(self.a*y.a+5*self.b*y.b,self.a*y.b+self.b*y.a)
    __rmul__=__mul__
    def __truediv__(self,y):
        y=Q(y);d=y.a*y.a-5*y.b*y.b
        assert d
        return self*Q(y.a/d,-y.b/d)
    def __rtruediv__(self,y):return Q(y)/self
    def __bool__(self):return bool(self.a or self.b)
    def __eq__(self,y):
        y=Q(y);return self.a==y.a and self.b==y.b
    def __hash__(self):return hash((self.a,self.b))
    def sign(self):
        a,b=self.a,self.b
        if not b:return (a>0)-(a<0)
        if not a:return (b>0)-(b<0)
        sa=(a>0)-(a<0);sb=(b>0)-(b<0)
        if sa==sb:return sa
        d=a*a-5*b*b
        return sa*((d>0)-(d<0))
    def enc(self):return [str(self.a),str(self.b)]

R=Q(F(16,11),F(4,11))
class E:
    __slots__=('a','b')
    def __init__(self,a=0,b=0):
        if isinstance(a,E):self.a,self.b=a.a,a.b
        else:self.a,self.b=Q(a),Q(b)
    def __add__(self,y):
        y=E(y);return E(self.a+y.a,self.b+y.b)
    __radd__=__add__
    def __neg__(self):return E(-self.a,-self.b)
    def __sub__(self,y):return self+-E(y)
    def __rsub__(self,y):return E(y)+-self
    def __mul__(self,y):
        y=E(y);return E(self.a*y.a+R*self.b*y.b,self.a*y.b+self.b*y.a)
    __rmul__=__mul__
    def __truediv__(self,y):
        y=E(y);d=y.a*y.a-R*y.b*y.b
        assert d
        return self*E(y.a/d,-y.b/d)
    def __rtruediv__(self,y):return E(y)/self
    def __bool__(self):return bool(self.a or self.b)
    def __eq__(self,y):
        y=E(y);return self.a==y.a and self.b==y.b
    def __hash__(self):return hash((self.a,self.b))
    def sign(self):
        sa,sb=self.a.sign(),self.b.sign()
        if not sa:return sb
        if not sb or sa==sb:return sa
        return sa*(self.a*self.a-R*self.b*self.b).sign()
    def sigma(self):return E(self.a,-self.b)
    def enc(self):return [self.a.enc(),self.b.enc()]

Z=E();ONE=E(1);PHI=E(Q(F(1,2),F(1,2)))
def dec(x):return E(Q(*x[0]),Q(*x[1]))
def dv(v):return [dec(x) for x in v]
def dm(M):return [dv(v) for v in M]
def enc(x):
    if isinstance(x,(Q,E)):return x.enc()
    if isinstance(x,(list,tuple)):return [enc(y) for y in x]
    if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
    return x

def dot(a,b):return sum((x*y for x,y in zip(a,b)),Z)
def mv(M,v):return [dot(row,v) for row in M]
def pair(a,b,M):return dot(a,mv(M,b))
def eye(n):return [[E(i==j) for j in range(n)] for i in range(n)]
def inverse(M):
    n=len(M);A=[row[:]+r for row,r in zip(M,eye(n))]
    for j in range(n):
        p=next(i for i in range(j,n) if A[i][j]);A[j],A[p]=A[p],A[j]
        v=A[j][j];A[j]=[x/v for x in A[j]]
        for i in range(n):
            if i!=j and A[i][j]:
                c=A[i][j];A[i]=[x-c*y for x,y in zip(A[i],A[j])]
    return [r[n:] for r in A]

def inertia(M):
    A=[r[:] for r in M];sg=[]
    for i in range(len(A)):
        d=A[i][i];assert d;sg.append(d.sign())
        for j in range(i+1,len(A)):
            for k in range(i+1,len(A)):A[j][k]=A[j][k]-A[j][i]*A[i][k]/d
    return [sg.count(1),sg.count(-1)]

def kernel_one(M):
    """Forward elimination/back-substitution, independent of producer RREF."""
    A=[r[:] for r in M];n=len(A[0]);piv=[]
    for j in range(n):
        k=len(piv);p=next((i for i in range(k,len(A)) if A[i][j]),None)
        if p is None:continue
        A[k],A[p]=A[p],A[k]
        for i in range(k+1,len(A)):
            if A[i][j]:
                c=A[i][j]/A[k][j]
                for a in range(j,n):A[i][a]=A[i][a]-c*A[k][a]
        piv.append(j)
    if len(piv)!=n-1:return None
    free=next(j for j in range(n) if j not in piv);v=[Z]*n;v[free]=ONE
    for i,j in reversed(list(enumerate(piv))):v[j]=-dot(A[i][j+1:],v[j+1:])/A[i][j]
    assert all(not dot(r,v) for r in M)
    return v

def normkey(v):
    a=next(x for x in v if x);a=a if a.sign()>0 else -a
    return tuple(x/a for x in v)

def census(M,N,expected):
    n=len(M);rows=[mv(M,v) for v in N]
    dual=[sum((v[j] for v in N),Z) for j in range(n)]
    assert pair(dual,dual,M).sign()<0
    found={};rankdrop=dualzero=trials=0
    for ids in combinations(range(len(N)),n-1):
        trials+=1;v=kernel_one([rows[i] for i in ids])
        if v is None:rankdrop+=1;continue
        s=pair(v,dual,M).sign()
        if not s:dualzero+=1;continue
        if s>0:v=[-x for x in v]
        vals=[dot(row,v) for row in rows]
        if any(x.sign()>0 for x in vals):continue
        v=normkey(v);q=pair(v,v,M);active=[i for i,row in enumerate(rows) if not dot(row,v)]
        assert pair(v,dual,M).sign()<0
        found[v]=({-1:'finite',0:'ideal',1:'spacelike'}[q.sign()],active,q)
    old={normkey(dv(row['vector'])):(row['kind'],row['incident'],dec(row['norm'])) for row in expected['rays']}
    assert len(old)==len(expected['rays']) and old==found
    counts=dict(Counter(row[0] for row in found.values()));assert counts==expected['counts']
    return dict(counts=counts,tested_supports=trials,rank_deficient_supports=rankdrop,
        zero_dual_nullspaces=dualzero,dual_sum=dual,dual_sum_norm=pair(dual,dual,M),
        rays=[dict(vector=v,kind=k,incident=a,norm=q) for v,(k,a,q) in sorted(found.items(),key=lambda x:str(enc(x[0])))])

COS={2:Z,3:E(F(1,4)),4:E(F(1,2)),5:E(Q(F(3,8),F(1,8))),6:E(F(3,4)),10:E(Q(F(5,8),F(1,8)))}
def relation(g,q,r):
    if not g:return ('right',2)
    if g.sign()>0:return ('obtuse',None)
    c=g*g/(q*r)
    if (c-1).sign()>=0:return ('parallel' if c==1 else 'dotted',None)
    for m,w in COS.items():
        if c==w:return ('finite',m)
    return ('forbidden_finite',None)

def polar_audit(M,N,h,roundrow):
    rays=roundrow['complete_ray_census']['rays'];bank=roundrow['polar_candidates']
    ids=[r['ray'] for r in bank]
    assert len(ids)==len(set(ids)) and set(ids)=={i for i,r in enumerate(rays) if r['kind']=='spacelike'}
    passes=[];rejected=[]
    for i,row in enumerate(bank):
        v=dv(rays[row['ray']]['vector']);q=pair(v,v,M)
        assert v==dv(row['normal']) and q==dec(row['norm']) and q.sign()>0
        rels=[relation(pair(v,w,M),q,pair(w,w,M)) for w in N]
        bad=[j for j,(t,_) in enumerate(rels) if t in ('obtuse','forbidden_finite')]
        assert bad==row['bad_old_walls'] and (not bad)==row['passes']
        for w,(t,m) in zip(row['relations'],rels):
            assert w['type']==t
            if m is not None:assert w['label']==m
        for j,w in enumerate(row['relations']):
            g=pair(v,N[j],M)
            if 'squared_weight' in w:assert dec(w['squared_weight'])==g*g/(q*pair(N[j],N[j],M))
            if 'product' in w:assert dec(w['product'])==g
        assert pair(v,h,M)==dec(row['center_pairing']) and pair(v,h,M).sign()<0
        if not bad:passes.append(i)
        else:rejected.append(dict(candidate=i,first_bad_wall=bad[0],relation=rels[bad[0]][0]))
    assert passes==roundrow['passing_polar_indices']
    return dict(candidate_count=len(bank),passing_indices=passes,rejections=rejected)

def square_rat(x):
    return x>=0 and isqrt(x.numerator)**2==x.numerator and isqrt(x.denominator)**2==x.denominator

def run():
    data=json.loads((ROOT/f'results/{BASE}.json').read_text())
    for p,h in data['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h
    assert R==Q(*data['field']['R']) and R.sign()>0 and Q(R.a,-R.b).sign()>0
    nr=R.a*R.a-5*R.b*R.b;assert nr==F(16,11)
    assert not square_rat(nr) and not square_rat(nr/5)
    G=eye(5)
    for i in range(5):G[i][i]=E(2)
    for i,x in enumerate([PHI,ONE,ONE,ONE]):G[i][i+1]=G[i+1][i]=-x
    Gi=inverse(G);d=Gi[4][4];D0=2-d
    assert D0==E(Q(F(7,4),F(1,4))) and (4-4*d)/D0==E(R)
    M5=[r+[Z] for r in G]+[[Z]*5+[D0]]
    g=[-Gi[i][4] for i in range(5)]
    W=[v+[Z] for v in eye(5)]+[g+[ONE],[2*x for x in g]+[E(0,-1)]]
    assert M5==dm(data['source_form']) and W==dm(data['source_normals'])
    assert [pair(w,w,M5) for w in W]==[E(2)]*6+[E(4)]
    assert inertia(G)==[4,1] and inertia(M5)==[5,1]
    D=eye(4)
    for i in range(4):D[i][i]=E(2)
    for j in (1,2,3):D[0][j]=D[j][0]=-ONE
    Di=inverse(D);omega=[Di[i][1] for i in range(4)]
    assert pair(omega,omega,D)==1 and D==dm(data['support_gram']) and omega==dv(data['support_weight'])
    M=[r+[Z]*4 for r in M5]+[[Z]*6+r for r in D]
    N=[w+[Z]*4 for w in W];N[6]=[x/2 for x in W[6]]+[-x for x in omega]
    N += [[Z]*6+v for v in eye(4)]
    h5=[-sum(row,Z) for row in Gi]+[Z];h=h5+[-sum(row,Z)/100 for row in Di]
    assert M==dm(data['form']) and N==dm(data['seed_normals']) and h==dv(data['center'])
    assert inertia(M)==[9,1] and all(pair(v,v,M)==2 for v in N)
    # Literal independent defining rows prove pointedness of both ambient cones.
    inverse([mv(M,N[i]) for i in list(range(6))+list(range(7,11))])
    assert pair(h,h,M)==dec(data['center_norm']) and pair(h,h,M).sign()<0
    assert all(pair(h,v,M).sign()<0 for v in N)
    assert pair(h5,h5,M5).sign()<0 and all(pair(h5,w,M5).sign()<0 for w in W)
    source=census(M5,W,data['source_complete_ray_census']);print('source',source['counts'],flush=True)
    seed=census(M,N,data['seed_complete_ray_census']);print('seed',seed['counts'],flush=True)
    rounds=[];assert len(data['polar_rounds'])==2
    for k,row in enumerate(data['polar_rounds']):
        assert row['wall_count']==len(N)
        if k==0:assert row['complete_ray_census']==data['seed_complete_ray_census']
        else:assert row['complete_ray_census']==data['final_complete_ray_census']
        r=polar_audit(M,N,h,row);rounds.append(r)
        for i in r['passing_indices']:N.append(dv(row['polar_candidates'][i]['normal']))
    assert rounds[0]['candidate_count']==20 and len(rounds[0]['passing_indices'])==1
    assert rounds[1]['candidate_count']==38 and rounds[1]['passing_indices']==[]
    assert data['exhausted_iterated_polar_closure'] and len(N)==12 and N==dm(data['final_normals'])
    final=census(M,N,data['final_complete_ray_census']);print('final',final['counts'],flush=True)
    pairrows={tuple(row['pair']):row for row in data['final_pairs']};counts=Counter();facets=[]
    for i,v in enumerate(N):
        q=pair(v,v,M);assert q.sign()>0
        p=[x-pair(h,v,M)*y/q for x,y in zip(h,v)]
        assert p==dv(data['final_actual_facet_points'][i])
        assert pair(p,p,M).sign()<0 and pair(p,h,M).sign()<0
        assert all(pair(p,w,M).sign()==(0 if i==j else -1) for j,w in enumerate(N))
        facets.append(p)
        for j in range(i):
            g0=pair(v,N[j],M);t,m=relation(g0,q,pair(N[j],N[j],M));row=pairrows[j,i]
            assert row['type']==t and t in ('right','finite','parallel','dotted')
            if m is not None:assert row['label']==m
            if 'squared_weight' in row:assert dec(row['squared_weight'])==g0*g0/(q*pair(N[j],N[j],M))
            counts[t+(str(m) if m else '')]+=1
    assert len(pairrows)==66 and dict(counts)==data['final_pair_counts']
    er=data['strict_escape'];v=dv(er['vector']);q=pair(v,v,M);vals=[pair(v,w,M) for w in N]
    assert all(not x.b and not x.a.b for x in v)
    assert q.sign()>0 and pair(v,h,M).sign()<0 and all(x.sign()<0 for x in vals)
    assert q==dec(er['norm']) and vals==dv(er['all_pairings']) and pair(v,h,M)==dec(er['center_pairing'])
    seedN=dm(data['seed_normals']);end=pair(seedN[5],seedN[6],M);cycle=end*end
    ob=data['nonquasiarithmetic_obstruction']
    assert cycle==dec(ob['end_cartan_two_cycle']) and cycle.b
    explicit=E(Q(2,F(1,2)),Q(F(-1,8),F(3,8)));assert cycle==explicit
    assert G==dm(ob['fixed_principal_gram']) and inertia(G)==[4,1]
    assert all(x.sigma()==x for row in G for x in row)
    assert cycle.sigma()!=cycle
    basecycle=pair(seedN[0],seedN[1],M);basecycle=basecycle*basecycle
    assert basecycle==E(Q(F(3,2),F(1,2)))
    files=[f'scripts/{TAG}.py',f'scripts/{BASE}.py',f'results/{BASE}.json']+list(data['sources'])
    return enc(dict(status='PASS_INDEPENDENT_QUARTIC_MIXED_D4_SUPPORT_CENSUS',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(set(files))},
        method='Independent stdlib Fraction Q(sqrt5)(sqrtR), forward nullspace elimination over every support; no producer/project imports',
        field=dict(R=R,norm_R=str(nr),nonnormality_square_tests_pass=True,
            rationale='The nonsquares normR and normR/5 exclude normality; the sole nontrivial abelian subfield is Qsqrt5'),
        source_census=source,seed_census=seed,polar_round_audits=rounds,final_census=final,
        final_pair_counts=dict(counts),all_twelve_actual_facet_points=facets,
        strict_escape=dict(vector=v,norm=q,center_pairing=pair(v,h,M),all_pairings=vals),
        cyclic_field_witness=dict(base_two_cycle=basecycle,end_two_cycle=cycle,
            sigma_cycle=cycle.sigma(),fixed_five_block=G,fixed_five_block_signature=[4,1]),
        scope='Exact complete source/seed/final cone ray and iterative extreme-ray polar census. No arbitrary-cap exhaustion or finite-volume claim. Conditional non-quasi-arithmetic conclusion uses Vinberg conjugate-semidefiniteness as an external theorem.'))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    out=run();path=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],flush=True)
