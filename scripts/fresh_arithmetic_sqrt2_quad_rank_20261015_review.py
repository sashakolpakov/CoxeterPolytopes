#!/usr/bin/env python3
"""Independent exact direct-bridge quadrilateral diagnostic and proof audit."""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations,product
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_arithmetic_sqrt2_quad_rank_20261015_review'
PRIMARY='fresh_arithmetic_sqrt2_quad_rank_20261015'


class Q:
    def __init__(self,a=0,b=0):
        if isinstance(a,Q):self.a,self.b=a.a,a.b
        else:self.a,self.b=F(a),F(b)
    def __add__(self,v):v=Q(v);return Q(self.a+v.a,self.b+v.b)
    __radd__=__add__
    def __neg__(self):return Q(-self.a,-self.b)
    def __sub__(self,v):return self+-Q(v)
    def __rsub__(self,v):return Q(v)+-self
    def __mul__(self,v):v=Q(v);return Q(self.a*v.a+2*self.b*v.b,self.a*v.b+self.b*v.a)
    __rmul__=__mul__
    def __truediv__(self,v):
        v=Q(v);d=v.a*v.a-2*v.b*v.b;assert d
        return self*Q(v.a/d,-v.b/d)
    def __eq__(self,v):v=Q(v);return self.a==v.a and self.b==v.b
    def __bool__(self):return bool(self.a or self.b)
    def sign(self):
        sa=(self.a>0)-(self.a<0);sb=(self.b>0)-(self.b<0)
        if not sb:return sa
        if not sa or sa==sb:return sb
        d=self.a*self.a-2*self.b*self.b
        return sa*((d>0)-(d<0))
    def enc(self):return[str(self.a),str(self.b)]


def cycle(n):
    G=[[Q(i==j)for j in range(n)]for i in range(n)]
    for i in range(n):
        j=(i+1)%n;v=Q(0,F(1,2))if i==0 else Q(F(1,2))
        G[i][j]=G[j][i]=-v
    return G


def sub(G,I):return[[G[i][j]for j in I]for i in I]


def determinant(G):
    A=[row[:]for row in G];answer=Q(1)
    for i in range(len(A)):
        p=next((j for j in range(i,len(A))if A[j][i]),None)
        if p is None:return Q()
        if p!=i:A[p],A[i]=A[i],A[p];answer=-answer
        pivot=A[i][i];answer=answer*pivot
        for j in range(i+1,len(A)):
            ratio=A[j][i]/pivot
            for k in range(i+1,len(A)):A[j][k]=A[j][k]-ratio*A[i][k]
    return answer


@lru_cache(None)
def pd(n,entries):
    A=[[Q(*entries[i*n+j])for j in range(n)]for i in range(n)]
    for i in range(n):
        if A[i][i].sign()<=0:return False
        for j in range(i+1,n):
            for k in range(j,n):
                A[k][j]=A[j][k]=A[j][k]-A[j][i]*A[i][k]/A[i][i]
    return True


def ispd(G):return pd(len(G),tuple((v.a,v.b)for row in G for v in row))


def matching_polynomial(n):
    # Rational determinant terms from fixed points and disjoint edge2-cycles.
    # The two oriented whole cycles give the sole irrational term separately.
    edges=[(i,(i+1)%n)for i in range(n)];terms=[]
    for mask in range(1<<n):
        chosen=[i for i in range(n)if mask>>i&1];vertices=[v for i in chosen for v in edges[i]]
        if len(set(vertices))!=len(vertices):continue
        coefficient=F((-1)**len(chosen))
        for i in chosen:coefficient*=F(1,2)if i==0 else F(1,4)
        terms.append(dict(unmatched_diagonal_indices=[i for i in range(n)if i not in vertices],coefficient=str(coefficient)))
    return terms


def run():
    primary=json.loads((ROOT/f'results/{PRIMARY}.json').read_text())
    for p,h in primary.get('sources',{}).items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    cycles={n:cycle(n)for n in(3,4,5)};dc={n:determinant(G)for n,G in cycles.items()}
    rows=[]
    for n,G in cycles.items():
        assert dc[n]==Q(0 if n==3 else-F(1,16),-F(1,2**(n-1)))
        terms=matching_polynomial(n)
        assert sum(F(t['coefficient'])for t in terms)==dc[n].a
        minors=[]
        for i in range(n):
            d=determinant(sub(G,[j for j in range(n)if j!=i]))
            assert d.b==0 and d.a>0;minors.append(d.a)
        for size in range(1,n):
            for support in combinations(range(n),size):assert ispd(sub(G,list(support)))
        rows.append(dict(rank=n,determinant=dc[n].enc(),rational_diagonal_polynomial=terms,
                         proper_cofactors=list(map(str,minors))))
    counts={};survivors=[];tested=0
    for n,m in product((3,4,5),repeat=2):
        if n==m==3:continue
        key=f'{n},{m}';counts[key]=dict(tested=0,all_four_corners_PD=0)
        for black1 in combinations(range(n),2):
            for black2 in combinations(range(m),2):
                for u in set(range(n))-set(black1):
                    for v in set(range(m))-set(black2):
                        for label in(3,4):
                            G=[[Q()for _ in range(n+m)]for _ in range(n+m)]
                            for i in range(n):
                                for j in range(n):G[i][j]=cycles[n][i][j]
                            for i in range(m):
                                for j in range(m):G[n+i][n+j]=cycles[m][i][j]
                            delta=Q(F(1,2))if label==3 else Q(0,F(1,2))
                            G[u][n+v]=G[n+v][u]=-delta
                            tested+=1;counts[key]['tested']+=1
                            if not all(ispd(sub(G,[k for k in range(n+m)if k not in(a,n+b)]))for a,b in product(black1,black2)):continue
                            counts[key]['all_four_corners_PD']+=1
                            dl=determinant(sub(cycles[n],[i for i in range(n)if i!=u]))
                            dr=determinant(sub(cycles[m],[i for i in range(m)if i!=v]))
                            expected=dc[n]*dc[m]-delta*delta*dl*dr
                            actual=determinant(G);assert actual==expected and actual.b>0
                            survivors.append(dict(ranks=[n,m],black_pairs=[black1,black2],bridge=[u,v,label],determinant=actual.enc()))
    assert tested==4032 and len(survivors)==28
    assert {k:v['all_four_corners_PD']for k,v in counts.items()}=={'3,4':6,'3,5':8,'4,3':6,'4,4':0,'4,5':0,'5,3':8,'5,4':0,'5,5':0}
    saved=primary['complete_direct_placement_census']['counts']
    assert saved['total']==tested and saved['four_PD']==len(survivors)
    for k,v in counts.items():
        assert saved['total '+k]==v['tested']
        assert saved.get('four_PD '+k,0)==v['all_four_corners_PD']
    # Literal nonexcluded singular configurations, independently reconstructed.
    c3=[[Q()for _ in range(6)]for _ in range(6)]
    for offset in(0,3):
        for i in range(3):
            for j in range(3):c3[offset+i][offset+j]=cycles[3][i][j]
    c3[2][5]=c3[5][2]=Q(0,-F(1,2))
    assert determinant(c3)==0
    for a,b in product((0,1),(3,4)):assert ispd(sub(c3,[a,2,5,b]))
    previous,current=Q(1),Q(1)
    path_squared=[Q(F(1,2),F(1,4)),Q(F(1,4)),Q(F(1,2)),Q(F(1,4)),Q(F(1,2),F(1,4))]
    for square in path_squared:previous,current=current,current-square*previous
    assert current==0
    paths=[f'scripts/{TAG}.py',f'scripts/{PRIMARY}.py',f'results/{PRIMARY}.json']
    return dict(status='PASS_INDEPENDENT_SCOPED_QSQRT2_QUADRILATERAL_RANK_REVIEW',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
                cycle_algebra=rows,placement_counts=counts,total_placements=tested,survivors=survivors,
                outside_scope_examples=dict(two_C3_bridge4_determinant='0',two_C3_four_corners='PD F4',
                                            path_8_3_4_3_8_determinant='0'),
                pd_cache_entries=pd.cache_info().currsize,
                proof_scope='Two disjoint chordless cycles ranks3,4,5, one label4 each; at least one rank4/5; one direct finite white bridge. Arbitrary pendant white trees at one attachment and orthogonal white components are covered analytically by conjugate-PD diagonal monotonicity. A connector path of length>=2 is ruled out by orthogonal negative subspaces. Shared white vertices, additional cross edges, other Lanner types and two bare rank3 cycles are not claimed exhausted.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['total_placements'],len(out['survivors']))
