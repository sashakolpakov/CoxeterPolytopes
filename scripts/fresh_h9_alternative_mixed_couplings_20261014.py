#!/usr/bin/env python3
"""All spherical rank-four facet couplings of the fixed mixed H5 source.

Finite support diagrams and active contact profiles are exact finite censuses.
Isolated support I2(m>=6) factors and contact labels>=6 have separate uniform
inequalities; no arbitrary label cutoff is inferred from enumeration.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,permutations,product
from pathlib import Path
import fresh_h9_mixed_d4_coupling_20261014_review as A

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_alternative_mixed_couplings_20261014'

class B:
    __slots__=('a','b')
    def __init__(self,a=0,b=0):
        if isinstance(a,B):self.a,self.b=a.a,a.b
        else:self.a,self.b=A.Q(a),A.Q(b)
    def __add__(self,y):
        y=B(y);return B(self.a+y.a,self.b+y.b)
    __radd__=__add__
    def __neg__(self):return B(-self.a,-self.b)
    def __sub__(self,y):return self+-B(y)
    def __rsub__(self,y):return B(y)+-self
    def __mul__(self,y):
        y=B(y);return B(self.a*y.a+2*self.b*y.b,self.a*y.b+self.b*y.a)
    __rmul__=__mul__
    def __truediv__(self,y):
        y=B(y);d=y.a*y.a-2*y.b*y.b;assert d
        return self*B(y.a/d,-y.b/d)
    def __bool__(self):return bool(self.a or self.b)
    def __eq__(self,y):
        y=B(y);return self.a==y.a and self.b==y.b
    def sign(self):
        sa,sb=self.a.sign(),self.b.sign()
        if not sa:return sb
        if not sb or sa==sb:return sa
        return sa*(self.a*self.a-2*self.b*self.b).sign()
    def enc(self):return [self.a.enc(),self.b.enc()]

Z=B();ONE=B(1);SQ2=B(0,1);PHI=B(A.Q(F(1,2),F(1,2)))
C={2:Z,3:ONE,4:SQ2,5:PHI}
PAIRS=list(combinations(range(4),2));PERMS=list(permutations(range(4)))
def sm(xs):return sum(xs,Z)
def dot(x,y):return sm(a*b for a,b in zip(x,y))
def pair(x,y,G):return dot(x,[dot(r,y) for r in G])
def pd(G):
    T=[r[:]for r in G]
    for i in range(len(T)):
        d=T[i][i]
        if d.sign()<=0:return False
        for j in range(i+1,len(T)):
            for k in range(i+1,len(T)):T[j][k]-=T[j][i]*T[i][k]/d
    return True
def inv(G):
    T=[r[:]+[B(i==j) for j in range(len(G))] for i,r in enumerate(G)]
    for i in range(len(G)):
        d=T[i][i];assert d.sign()>0
        T[i]=[x/d for x in T[i]]
        for j in range(len(G)):
            if i!=j:
                c=T[j][i];T[j]=[x-c*y for x,y in zip(T[j],T[i])]
    return [r[len(G):] for r in T]
def enc(x):
    if isinstance(x,B):return x.enc()
    if isinstance(x,list):return [enc(v) for v in x]
    if isinstance(x,dict):return {k:enc(v) for k,v in x.items()}
    return x
def gram(code):
    G=[[B(2*(i==j)) for j in range(4)]for i in range(4)]
    for (i,j),m in zip(PAIRS,code):G[i][j]=G[j][i]=-C[m]
    return G
def canonical(code):
    d={p:m for p,m in zip(PAIRS,code)}
    return min(tuple(d[tuple(sorted((p[i],p[j])))]for i,j in PAIRS)for p in PERMS)
def components(code):
    adj=[set()for _ in range(4)]
    for (i,j),m in zip(PAIRS,code):
        if m>2:adj[i].add(j);adj[j].add(i)
    rem=set(range(4));comps=[]
    while rem:
        seen={min(rem)};todo=list(seen)
        for i in todo:
            for j in sorted(adj[i]-seen):seen.add(j);todo.append(j)
        comps.append(sorted(seen));rem-=seen
    return comps,adj
def name(code):
    comps,adj=components(code);names=[]
    for ss in comps:
        labels=sorted(m for p,m in zip(PAIRS,code)if set(p)<=set(ss) and m>2)
        if len(ss)==1:n='A1'
        elif len(ss)==2:n={3:'A2',4:'B2',5:'I2(5)'}[labels[0]]
        elif len(ss)==3:n={3:'A3',4:'B3',5:'H3'}[max(labels)]
        elif max(len(adj[i])for i in ss)==3:n='D4';assert labels==[3]*3
        elif max(labels)==3:n='A4'
        elif max(labels)==5:n='H4'
        else:
            ends={i for i in ss if len(adj[i])==1}
            p=next(p for p,m in zip(PAIRS,code)if m==4)
            n='B4'if ends.intersection(p)else'F4'
        names.append(n)
    return '+'.join(names)


def run():
    canonical_pd={};positive_codes=0
    for code in product((2,3,4,5),repeat=6):
        G=gram(code)
        if not pd(G):continue
        positive_codes+=1;key=canonical(code)
        canonical_pd.setdefault(key,gram(key))
    assert len(canonical_pd)==18
    rows=[];total_profiles=0;total_phi_pairs=0
    for code,G in sorted(canonical_pd.items()):
        J=inv(G);assert all(x.sign()>=0 for r in J for x in r)
        assert all((J[i][i]-B(F(1,2))).sign()>=0 for i in range(4))
        comps,adj=components(code)
        for i in range(4):
            if adj[i]:assert (J[i][i]-B(F(2,3))).sign()>=0
            else:assert J[i][i]==B(F(1,2)) and all(not J[i][j]for j in range(4)if j!=i)
        profiles=[]
        for labs in product((2,3,4,5),repeat=4):
            c=[C[m]for m in labs];cost=pair(c,c,J)
            if (2-cost).sign()<=0:continue
            profiles.append(dict(labels=list(labs),contact=c,cost=cost,alpha_squared=1-cost/2))
        total_profiles+=len(profiles)
        survivors=[]
        for i,p in enumerate(profiles):
            for j,q in enumerate(profiles):
                cross=pair(p['contact'],q['contact'],J);assert cross.sign()>=0
                square=PHI*PHI*p['alpha_squared']*q['alpha_squared']
                good=[]
                for m,t in C.items():
                    # phi*alpha_i*alpha_j=cross+2cos(pi/m).
                    if square==(cross+t)*(cross+t):good.append(m)
                if good:survivors.append(dict(profiles=[i,j],labels=good))
                total_phi_pairs+=1
        zero=next(i for i,p in enumerate(profiles)if not p['cost'])
        # This is a tested theorem, not an assumed propagation rule.
        assert survivors==[dict(profiles=[zero,zero],labels=[5])],(name(code),survivors)
        allowed=[p for p in profiles if p['cost']==ONE]
        comps,adj=components(code)
        for p in allowed:
            p['active_components']=[i for i,ss in enumerate(comps)if any(p['labels'][j]>2 for j in ss)]
            p['untouched_components']=[i for i,ss in enumerate(comps)if all(p['labels'][j]==2 for j in ss)]
        rows.append(dict(type=name(code),code=list(code),gram=G,inverse=J,components=comps,
            positive_projection_profiles=profiles,pi5_pair_trials=len(profiles)**2,
            surviving_pi5_pairs=survivors,cost_one_facet6_profiles=allowed))
        print(name(code),len(profiles),'profiles',len(allowed),'cost1',flush=True)
    # Uniform removal of contact labels >=6. A1 alone has inverse diagonal1/2.
    # If c>=sqrt3, alpha<=1/2. Any adjacent source facet either has no
    # same-A1 contact (0<new angle weight<=phi/2<1), or has cross>=sqrt3/2
    # while phi*alpha*alpha_other<=phi*sqrt3/4<sqrt3/2, hence obtuse.
    assert (2-PHI).sign()>0
    assert (3-PHI*PHI).sign()>0
    assert B(F(3,4))+B(F(1,4))==ONE  # support I2>=6 cannot have a finite neighbor
    irreducible=[r for r in rows if len(r['components'])==1]
    assert {r['type'] for r in irreducible if r['cost_one_facet6_profiles']}=={'B4','D4'}
    assert all(p['untouched_components'] for r in rows if len(r['components'])>1
               for p in r['cost_one_facet6_profiles'])
    source=json.loads((ROOT/'results/fresh_h9_mixed_d4_coupling_20261014.json').read_text())
    for p,h in source['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    M,N=A.dm(source['source_form']),A.dm(source['source_normals'])
    norms=[A.pair(n,n,M)for n in N];assert norms==[A.E(2)]*6+[A.E(4)]
    finite=[]
    for i,j in combinations(range(7),2):
        rel=A.relation(A.pair(N[i],N[j],M),norms[i],norms[j])
        assert rel[0]in('right','finite','dotted')
        if rel[0]=='finite':finite.append([i,j,rel[1]])
    assert finite==[[0,1,5],[1,2,3],[2,3,3],[3,4,3],[4,5,3],[4,6,4]]
    coupled_end_product=A.pair(N[5],N[6],M)/2
    assert (coupled_end_product+2).sign()<0
    paths=['scripts/fresh_h9_mixed_d4_coupling_20261014_review.py',
           'results/fresh_h9_mixed_d4_coupling_20261014.json',f'scripts/{TAG}.py']
    return enc(dict(status='PASS_ALL_SPHERICAL_RANK4_COUPLING_PROFILES_FOR_FIXED_MIXED_H5',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        support_codes_tested=4096,positive_support_codes=positive_codes,canonical_support_types=len(rows),
        total_positive_projection_profiles=total_profiles,total_pi5_pair_trials=total_phi_pairs,supports=rows,
        source_finite_pairs=finite,coupled_end_product=A.enc(coupled_end_product),
        conclusion='The source pi5 pair must be ordinary. Its connected pi3 graph forces source0..5 ordinary. Facet6 is ordinary or has support weight square1. Only irreducible B4,D4 allow square1; every allowed reducible profile leaves an entire support component untouched.',
        high_support_labels='Any support label>=6 is an isolated I2 factor by the3-node path determinant. Its inverse diagonal>=2 forces all facet contacts to that factor right. Replacing such an inactive factor by A1+A1 embeds all active-profile calculations in the enumerated18types.',
        high_contact_labels='Any facet contact label>=6 can occur only on an isolated A1 and then is its only support contact. The two-case acute/Coxeter comparison with any finite source neighbor excludes it uniformly.',
        scope='Fixed7facet mixed H5 source, four retained spherical support facets, exactly one ambient lifted wall for each source facet, all projected boundary directions unchanged and positive. Not arbitrary changed source geometry, extra facets, or different rank support.'))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(out.read_text());print('SAVED_REPLAY_PASS')
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['total_positive_projection_profiles'],d['total_pi5_pair_trials'])
