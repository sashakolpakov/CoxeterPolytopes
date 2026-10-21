#!/usr/bin/env python3
"""Complete finite-contact one-root profiles at literal E7/E8 over Q(sqrt2).

Necessary local arithmetic profiles, not whole polytope constructions.
All comparisons use rational coefficients and exact quadratic-field signs.
"""
import argparse,json
from fractions import Fraction as Q
from hashlib import sha256
from itertools import product,combinations
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_arithmetic_sqrt2_local_20261015'

def sign(x):return (x>0)-(x<0)
def qs(a,b):
    if not b:return sign(a)
    if not a or sign(a)==sign(b):return sign(b)
    return sign(a)*sign(a*a-2*b*b)
def inv(M):
    n=len(M);A=[[Q(v) for v in row]+[Q(i==j) for j in range(n)] for i,row in enumerate(M)]
    for i in range(n):
        pivot=A[i][i];assert pivot>0
        A[i]=[x/pivot for x in A[i]]
        for j in range(n):
            if i!=j:
                c=A[j][i];A[j]=[x-c*y for x,y in zip(A[j],A[i])]
    return [row[n:] for row in A]
def encode(a,b):return [str(a),str(b)]

def check_rank(n):
    edges=[(i,i+1) for i in range(n-2)]+[(2,n-1)]
    G=[[2*int(i==j) for j in range(n)] for i in range(n)]
    for i,j in edges:G[i][j]=G[j][i]=-1
    J=inv(G)
    assert all(sum(Q(G[i][k])*J[k][j] for k in range(n))==int(i==j) for i in range(n) for j in range(n))
    counts={};survivors=[];pair_profiles={pair:[] for pair in combinations(range(n),2)}
    for kind,alphabet in [('mixed34',(0,1,2)),('pure6',(0,1)),('pure8',(0,1))]:
        counts[kind]=dict(total=0,hyperbolic=0,conjugate_positive=0,with_two_finite_replacements=0)
        for row in product(alphabet,repeat=n):
            counts[kind]['total']+=1
            if kind=='mixed34':
                u=[Q(x==1) for x in row];v=[Q(x==2) for x in row]
                a=sum(J[i][j]*(u[i]*u[j]+2*v[i]*v[j]) for i in range(n) for j in range(n))
                b=sum(J[i][j]*(u[i]*v[j]+v[i]*u[j]) for i in range(n) for j in range(n))
            else:
                quadratic=sum(J[i][j]*row[i]*row[j] for i in range(n) for j in range(n))
                a,b=(3*quadratic,Q(0)) if kind=='pure6' else (2*quadratic,quadratic)
            if qs(a-2,b)<=0:continue
            counts[kind]['hyperbolic']+=1
            if qs(2-a,b)<=0:continue
            counts[kind]['conjugate_positive']+=1
            finite=[];replacement_schur={}
            for i in range(n):
                if kind=='mixed34':
                    x=sum(J[i][j]*u[j] for j in range(n));y=sum(J[i][j]*v[j] for j in range(n))
                    da=(x*x+2*y*y)/J[i][i];db=2*x*y/J[i][i]
                else:
                    z=sum(J[i][j]*row[j] for j in range(n))**2/J[i][i]
                    da,db=(3*z,Q(0)) if kind=='pure6' else (2*z,z)
                sa,sb=2-a+da,-b+db
                if qs(sa,sb)>0:finite.append(i)
                replacement_schur[str(i)]=encode(sa,sb)
            rec=dict(kind=kind,labels=[(2 if not x else 3 if x==1 else 4) if kind=='mixed34' else (2 if not x else int(kind[4:])) for x in row],
                     projection_cost=encode(a,b),finite_replacements=finite,replacement_schur=replacement_schur)
            index=len(survivors);survivors.append(rec)
            if len(finite)>=2:counts[kind]['with_two_finite_replacements']+=1
            for pair in combinations(finite,2):pair_profiles[pair].append(index)
    distances=[]
    for i in range(n):
        d={i:0};queue=[i]
        for u in queue:
            for v in range(n):
                if G[u][v]==-1 and v not in d:d[v]=d[u]+1;queue.append(v)
        distances.append(d)
    return dict(type=f'E{n}',gram=G,inverse=[[str(v) for v in row] for row in J],counts=counts,
                surviving_profiles=survivors,marked_pairs=[dict(pair=list(pair),distance=distances[pair[0]][pair[1]],profiles=ix) for pair,ix in pair_profiles.items()],
                pairs_with_triangle=sum(bool(ix) for ix in pair_profiles.values()),
                pairs_without_triangle=[list(pair) for pair,ix in pair_profiles.items() if not ix])

def run():
    records=[check_rank(n) for n in (7,8)]
    files=[f'scripts/{TAG}.py']
    return dict(status='PASS_EXACT_COMPLETE_LITERAL_E7_E8_SQRT2_FINITE_CONTACT_PROFILE_CENSUS',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in files},
                finite_labels=[2,3,4,6,8],records=records,
                scope='Complete necessary one-root finite-contact profiles at these literal spherical components. Conjugate positivity and pair-product cyclefield enforced. Triangle witnesses concern the local cone only; quadrilateral completions and extra spherical components remain unclassified. No global exclusion or new finite-volume example.')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    data=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert data==json.loads(out.read_text())
    else:out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'])
    for rec in data['records']:
        print(rec['type'],rec['counts'],'triangle pairs',rec['pairs_with_triangle'],'missing',rec['pairs_without_triangle'])
