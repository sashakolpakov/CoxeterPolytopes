#!/usr/bin/env python3
"""Independent Fraction/Qsqrt5 construction and complete double-description audit.

No producer or repository arithmetic is imported.  Incremental halfspace
clipping enumerates rays, independently of the producer's support RREF census.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = 'fresh_h9_d4_prism_coupling_20261014'
TAG = BASE + '_review'

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

Z=Q();ONE=Q(1);PHI=Q(F(1,2),F(1,2))
def dec(x):return Q(*x)
def dv(v):return [dec(x) for x in v]
def dm(M):return [dv(v) for v in M]
def enc(x):
    if isinstance(x,Q):return x.enc()
    if isinstance(x,(tuple,list)):return [enc(y) for y in x]
    if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
    return x

def dot(a,b):return sum((x*y for x,y in zip(a,b)),Z)
def mv(M,v):return [dot(row,v) for row in M]
def pair(a,b,M):return dot(a,mv(M,b))
def eye(n):return [[Q(i==j) for j in range(n)] for i in range(n)]
def rank(M):
    if not M:return 0
    A=[r[:] for r in M];k=0
    for j in range(len(A[0])):
        p=next((i for i in range(k,len(A)) if A[i][j]),None)
        if p is None:continue
        A[k],A[p]=A[p],A[k];v=A[k][j]
        A[k]=[x/v for x in A[k]]
        for i in range(k+1,len(A)):
            if A[i][j]:
                c=A[i][j];A[i]=[x-c*y for x,y in zip(A[i],A[k])]
        k+=1
        if k==len(A):break
    return k

def inverse(M):
    n=len(M);A=[row[:]+r for row,r in zip(M,eye(n))]
    for j in range(n):
        p=next(i for i in range(j,n) if A[i][j]);A[j],A[p]=A[p],A[j]
        v=A[j][j];A[j]=[x/v for x in A[j]]
        for i in range(n):
            if i!=j and A[i][j]:
                c=A[i][j];A[i]=[x-c*y for x,y in zip(A[i],A[j])]
    return [r[n:] for r in A]

def determinant(M):
    A=[r[:] for r in M];out=ONE
    for j in range(len(A)):
        p=next((i for i in range(j,len(A)) if A[i][j]),None)
        if p is None:return Z
        if p!=j:A[j],A[p]=A[p],A[j];out=-out
        v=A[j][j];out=out*v
        for i in range(j+1,len(A)):
            c=A[i][j]/v
            for k in range(j+1,len(A)):A[i][k]=A[i][k]-c*A[j][k]
    return out

def inertia(M):
    # Every leading pivot in these explicit forms is nonzero.
    A=[r[:] for r in M];sg=[]
    for i in range(len(A)):
        d=A[i][i];assert d;sg.append(d.sign())
        for j in range(i+1,len(A)):
            for k in range(i+1,len(A)):A[j][k]=A[j][k]-A[j][i]*A[i][k]/d
    return dict(positive=sg.count(1),negative=sg.count(-1),zero=sg.count(0))

def normalize(v):
    a=next(x for x in v if x);a=a if a.sign()>0 else -a
    return tuple(x/a for x in v)

def clipped_rays(M,N,expected):
    """Pointed full-dimensional cone: start with n independent rows, add all rest.

    At each step an extreme new ray crosses a two-dimensional old face.  Its
    generators have common active row rank n-2; conversely such a crossing is
    extreme.  Thus checking every positive/negative pair is complete, including
    nonsimple rays and nonadjacent pairs (the latter fail the rank test).
    """
    n=len(M);rows=[mv(M,v) for v in N];base=[]
    dual=[sum((v[j] for v in N),Z) for j in range(n)]
    assert pair(dual,dual,M).sign()<0
    for i,row in enumerate(rows):
        if rank([rows[j] for j in base]+[row])>len(base):base.append(i)
        if len(base)==n:break
    assert len(base)==n
    Ainv=inverse([rows[j] for j in base])
    rays={normalize([-Ainv[j][i] for j in range(n)]) for i in range(n)}
    order=base+[i for i in range(len(N)) if i not in base]
    used=base[:];steps=[];cached={}
    for idx in order[n:]:
        active={v:frozenset(i for i in used if not dot(rows[i],v)) for v in rays}
        vals={v:dot(rows[idx],v) for v in rays}
        pos=[v for v in rays if vals[v].sign()>0]
        neg=[v for v in rays if vals[v].sign()<0]
        keep={v for v in rays if vals[v].sign()<=0};made=0
        for p in pos:
            for q in neg:
                common=active[p]&active[q]
                key=tuple(sorted(common))
                if key not in cached:cached[key]=rank([rows[i] for i in key])
                if cached[key]!=n-2:continue
                v=normalize([vals[p]*b-vals[q]*a for a,b in zip(p,q)])
                assert not dot(rows[idx],v)
                assert all(dot(rows[i],v).sign()<=0 for i in used)
                assert rank([rows[i] for i in key]+[rows[idx]])==n-1
                keep.add(v);made+=1
        steps.append(dict(added_wall=idx,old_rays=len(rays),positive=len(pos),negative=len(neg),
                          crossing_pairs=len(pos)*len(neg),adjacent_crossings=made,new_rays=len(keep)))
        rays=keep;used.append(idx)
    found={}
    for v in rays:
        vv=pair(v,v,M);inc=[i for i,row in enumerate(rows) if not dot(row,v)]
        assert rank([rows[i] for i in inc])==n-1
        assert all(dot(row,v).sign()<=0 for row in rows)
        assert pair(v,dual,M).sign()<0
        found[v]=( {-1:'finite',0:'ideal',1:'spacelike'}[vv.sign()],inc,vv)
    archived={normalize(dv(r['vector'])):(r['kind'],r['incident'],dec(r['norm'])) for r in expected['rays']}
    assert len(archived)==len(expected['rays']) and found==archived
    counts=dict(Counter(x[0] for x in found.values()));assert counts==expected['counts']
    return dict(counts=counts,dual_sum=dual,dual_sum_norm=pair(dual,dual,M),initial_independent_walls=base,clipping_steps=steps,
                complete_ray_count=len(found),cached_face_rank_tests=len(cached),
                rays=[dict(vector=v,kind=k,incident=i,norm=q) for v,(k,i,q) in sorted(found.items(),key=lambda x:str(enc(x[0])))])

COS={2:Z,3:Q(F(1,4)),4:Q(F(1,2)),5:Q(F(3,8),F(1,8)),6:Q(F(3,4)),10:Q(F(5,8),F(1,8))}
def relation(g,q,r):
    if not g:return ('right',2)
    if g.sign()>0:return ('obtuse',None)
    c=g*g/(q*r)
    if (c-1).sign()>=0:return ('parallel' if c==1 else 'dotted',None)
    for m,w in COS.items():
        if c==w:return ('finite',m)
    return ('forbidden_finite',None)

def run():
    d=json.loads((ROOT/f'results/{BASE}.json').read_text())
    for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h
    J=eye(6)
    for i in range(5):J[i][i]=Q(2)
    for i,x in enumerate([PHI,ONE,ONE,ONE,ONE]):J[i][i+1]=J[i+1][i]=-x
    Ji=inverse(J);N0=eye(6)+[[-PHI*Ji[i][5] for i in range(6)]]
    D=eye(4)
    for i in range(4):D[i][i]=Q(2)
    for j in (1,2,3):D[0][j]=D[j][0]=-ONE
    Di=inverse(D);omega=[Di[i][1] for i in range(4)]
    M=[r+[Z]*4 for r in J]+[[Z]*6+r for r in D]
    assert J==dm(d['source_form']) and N0==dm(d['source_normals'])
    assert D==dm(d['support_gram']) and omega==dv(d['support_weight']) and M==dm(d['form'])
    assert inertia(J)==dict(positive=5,negative=1,zero=0)
    assert inertia(M)==dict(positive=9,negative=1,zero=0)
    N=[v+[Z]*4 for v in N0];N[5]=N0[5]+[-x for x in omega]
    N += [[Z]*6+v for v in eye(4)]
    h=[-sum(r,Z) for r in Ji]+[-sum(r,Z)/10 for r in Di]
    assert N==dm(d['normals']) and h==dv(d['center'])
    assert pair(h,h,M).sign()<0 and all(pair(h,v,M).sign()<0 for v in N)
    assert all(pair(v,v,M)==2 for v in N)
    G=[[pair(a,b,M) for b in N] for a in N];assert G==dm(d['gram'])
    for item in d['seed_pairs']:
        i,j=item['pair'];assert relation(G[i][j],Q(2),Q(2))==(item['type'],item['label'])
    source=clipped_rays(J,N0,d['source_complete_ray_census']);print('source',source['counts'],flush=True)
    seed=clipped_rays(M,N,d['complete_ray_census']);print('seed',seed['counts'],flush=True)
    accepted=[];polar_checks=[]
    candidate_rays=[r['ray'] for r in d['polar_candidates']]
    assert len(candidate_rays)==len(set(candidate_rays))
    assert set(candidate_rays)=={i for i,r in enumerate(d['complete_ray_census']['rays']) if r['kind']=='spacelike'}
    for k,item in enumerate(d['polar_candidates']):
        r=d['complete_ray_census']['rays'][item['ray']]
        assert r['kind']=='spacelike'
        v=dv(r['vector']);q=pair(v,v,M);assert v==dv(item['normal']) and q==dec(item['norm'])
        rels=[relation(pair(v,w,M),q,Q(2)) for w in N]
        bad=[i for i,(t,_) in enumerate(rels) if t in ('obtuse','forbidden_finite')]
        assert bad==item['bad_old_walls'] and (not bad)==item['passes_all_old_coxeter_pairs']
        for a,(t,m) in zip(item['relations'],rels):
            assert a['type']==t
            if m is not None:assert a['label']==m
        assert pair(v,h,M)==dec(item['center_pairing']) and pair(v,h,M).sign()<0
        polar_checks.append(dict(index=k,bad_old_walls=bad))
        if not bad:accepted.append(k)
    assert accepted==d['passing_polar_indices'] and len(accepted)==4
    extra=[dv(d['polar_candidates'][i]['normal']) for i in accepted]
    for a,b in combinations(extra,2):assert relation(pair(a,b,M),pair(a,a,M),pair(b,b,M))[0]=='dotted'
    E=N+extra;assert E==dm(d['extension_normals'])
    facet_points=[];relations=Counter()
    for i,v in enumerate(E):
        q=pair(v,v,M);assert q.sign()>0
        p=[a-pair(h,v,M)*b/q for a,b in zip(h,v)]
        assert pair(p,p,M).sign()<0 and pair(p,h,M).sign()<0
        assert all(pair(p,w,M).sign()==(0 if i==j else -1) for j,w in enumerate(E))
        facet_points.append(p)
        assert p==dv(d['extension_actual_facet_points'][i])
        for j in range(i):
            t,m=relation(pair(v,E[j],M),q,pair(E[j],E[j],M));assert t in ('right','finite','parallel','dotted')
            relations[t+(f'_{m}' if m else '')]+=1
    extended=clipped_rays(M,E,d['extension_complete_ray_census']);print('extended',extended['counts'],flush=True)
    # Recheck the frozen rational witness with no decimal approximation.
    escape_checks={}
    for key in ('strict_escape','extension_strict_escape'):
        if key not in d:continue
        row=d[key];v=dv(row.get('vector',row.get('point')))
        q=pair(v,v,M);values=[pair(v,w,M) for w in E]
        assert q.sign()>0 and pair(v,h,M).sign()<0 and all(x.sign()<0 for x in values)
        assert all(not x.b for x in v)
        assert q==dec(row['norm']) and pair(v,h,M)==dec(row['center_pairing'])
        assert values==dv(row['all_pairings'])
        escape_checks[key]=dict(vector=v,norm=q,center_pairing=pair(v,h,M),all_pairings=values)
    assert 'extension_strict_escape' in escape_checks
    # Literal module verification, separate from any external classification.
    host=d['literal_arithmetic_host'];B=eye(10)
    for i in range(4):B[6+i][9]=omega[i]
    cols=[list(x) for x in zip(*B)];Bi=inverse(B)
    H=[[pair(a,b,M) for b in cols] for a in cols]
    oi=lambda x:(2*x.b).denominator==1 and (x.a-x.b).denominator==1
    assert all(oi(x) for row in H for x in row)
    det=determinant(H);assert oi(det) and oi(ONE/det)
    HC=[[Q(x.a,-x.b) for x in row] for row in H]
    assert inertia(H)==dict(positive=9,negative=1,zero=0)
    assert inertia(HC)==dict(positive=10,negative=0,zero=0)
    assert H[5][5]==1
    assert cols==dm(host['basis_columns']) and H==dm(host['gram']) and det==dec(host['determinant'])
    matrices=[]
    for n,row in zip(E,host['reflections']):
        c=mv(Bi,n);dual=mv(H,c);q=dot(c,dual)
        T=[[Q(i==j)-2*c[i]*dual[j]/q for j in range(10)] for i in range(10)]
        assert c==dv(row['coordinates']) and T==dm(row['matrix']) and row['integral']
        assert all(oi(x) for line in T for x in line)
        Tc=[list(x) for x in zip(*T)]
        for i in range(10):
            for j in range(10):
                assert dot(T[i],Tc[j])==Q(i==j)
                assert pair(Tc[i],Tc[j],H)==H[i][j]
        matrices.append(T)
    assert len(matrices)==len(host['reflections'])==15
    literal_host=dict(basis_columns=cols,gram=H,determinant=det,
        signature=inertia(H),conjugate_signature=inertia(HC),norm_one_basis_index=5,
        integral_involutive_form_preserving_reflections=len(matrices),
        scope='Literal Z[phi] odd unimodular host and reflection matrices only; no independent classification or global nonreflectivity theorem proved here')
    paths=[f'scripts/{TAG}.py',f'scripts/{BASE}.py',f'results/{BASE}.json']+list(d['sources'])
    out=dict(status='PASS_INDEPENDENT_QSQRT5_COMPLETE_DD_COUPLED_D4_P2',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(set(paths))},
        arithmetic='Independent stdlib Fraction Q(sqrt5), no producer/project imports',
        enumeration='Independent incremental double description; every positive-negative ray pair tested for common-active rank dimension-2',
        signatures=dict(source=inertia(J),ambient=inertia(M)),center=h,center_norm=pair(h,h,M),
        source_census=source,seed_census=seed,polar_checks=polar_checks,accepted_polars=accepted,
        extension_pairs=dict(relations),extension_actual_facet_points=facet_points,
        extension_census=extended,strict_escape_checks=escape_checks,literal_host=literal_host,
        scope='Independent complete linear-cone ray census and actual hyperbolic facets. The strict spacelike witness certifies infinite volume, not compactness or finite volume.')
    return enc(out)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    out=run();path=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],flush=True)
