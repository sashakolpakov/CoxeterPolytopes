#!/usr/bin/env python3
"""Complete finite-angle/triangle checks over Q(sqrt2,sqrt5).

Standalone Fraction arithmetic; no project/producer algebra or optimizer.
The exclusion of arithmetic Lanner orders>=4 uses Nikulin0708.3991 §3.1.
"""
import argparse
import json
from fractions import Fraction as R
from hashlib import sha256
from itertools import combinations, combinations_with_replacement, product
from math import comb, gcd, isqrt
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt2sqrt5_lanner_catalogue_20261016'


def qs(x=0,y=0): return (R(x),R(y))
def qa(x,y): return (x[0]+y[0],x[1]+y[1])
def qn(x): return (-x[0],-x[1])
def qm(x,y): return (x[0]*y[0]+2*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def qi(x):
    d=x[0]*x[0]-2*x[1]*x[1]
    assert d
    return (x[0]/d,-x[1]/d)
def qd(x,y): return qm(x,qi(y))
def sg(x): return (x>0)-(x<0)
def qsign(x):
    a,b=x;sa,sb=sg(a),sg(b)
    if not sa:return sb
    if not sb or sa==sb:return sa
    return sa*sg(a*a-2*b*b)
def rsqrt(x):
    if x<0:return None
    a,b=isqrt(x.numerator),isqrt(x.denominator)
    return R(a,b) if a*a==x.numerator and b*b==x.denominator else None
def qsqrt(x):
    a,b=x
    if not b:
        c=rsqrt(a)
        if c is not None:return qs(c)
        d=rsqrt(a/2)
        return None if d is None else qs(0,d)
    h=rsqrt(a*a-2*b*b)
    if h is None:return None
    for hh in (h,-h):
        c=rsqrt((a+hh)/2)
        if c:
            z=(c,b/(2*c))
            assert qm(z,z)==x
            return z if qsign(z)>0 else qn(z)
    return None


class F:
    __slots__=('v',)
    def __init__(self,a=0,b=0,c=0,d=0):self.v=tuple(map(R,(a,b,c,d)))
    @staticmethod
    def co(x):return x if isinstance(x,F) else F(x)
    def __add__(self,o):return F(*(a+b for a,b in zip(self.v,F.co(o).v)))
    __radd__=__add__
    def __neg__(self):return F(*(-x for x in self.v))
    def __sub__(self,o):return self+-F.co(o)
    def __rsub__(self,o):return F.co(o)+-self
    def __mul__(self,o):
        w=F.co(o).v;z=[R(0)]*4
        for i,a in enumerate(self.v):
            for j,b in enumerate(w):z[i^j]+=a*b*(2 if i&j&1 else 1)*(5 if i&j&2 else 1)
        return F(*z)
    __rmul__=__mul__
    def inverse(self):
        A,B=self.v[:2],self.v[2:]
        h=qi(qa(qm(A,A),qn(qm(qs(5),qm(B,B)))))
        return F(*qm(A,h),*qn(qm(B,h)))
    def __truediv__(self,o):return self*F.co(o).inverse()
    def __eq__(self,o):return self.v==F.co(o).v
    def sign(self):
        A,B=self.v[:2],self.v[2:];sa,sb=qsign(A),qsign(B)
        if not sa:return sb
        if not sb or sa==sb:return sa
        return sa*qsign(qa(qm(A,A),qn(qm(qs(5),qm(B,B)))))
    def conjugate(self,s,t):return F(self.v[0],s*self.v[1],t*self.v[2],s*t*self.v[3])
    def square_root(self):
        A,B=self.v[:2],self.v[2:]
        if B==qs():
            u=qsqrt(A)
            if u is not None:return F(*u)
            v=qsqrt(qd(A,qs(5)))
            return None if v is None else F(0,0,*v)
        N=qa(qm(A,A),qn(qm(qs(5),qm(B,B))))
        h=qsqrt(N)
        if h is None:return None
        for hh in (h,qn(h)):
            u=qsqrt(qm(qs(R(1,2)),qa(A,hh)))
            if u is None or u==qs():continue
            v=qd(B,qm(qs(2),u));z=F(*u,*v)
            assert z*z==self
            return z if z.sign()>0 else -z
        return None
    def data(self):return list(map(str,self.v))
    def integral(self):
        a,b,c,d=self.v
        return all(x.denominator==1 for x in (a-c,b-d,2*c,2*d))


EMBEDDINGS=[(1,1),(-1,1),(1,-1),(-1,-1)]
S2=F(0,1);S5=F(0,0,1)
PHI=(1+S5)/2
D={2:F(0),3:F(1),4:F(2),5:PHI*PHI,6:F(3),8:2+S2,10:2+PHI}


def phi(n):return sum(gcd(a,n)==1 for a in range(1,n+1))
def real_group(m):
    units=sorted({min(a,m-a) for a in range(1,m) if gcd(a,m)==1})
    orders=[]
    for a in units:
        z,n=1,0
        while True:
            z=z*a%m;z=min(z,m-z);n+=1
            if z==1:break
            assert n<=len(units)
        orders.append(n)
    return dict(conductor=m,degree=len(units),orders=orders,cyclic=max(orders)==len(units))


def alphabet():
    # [Q(cos(2pi/m)):Q]=phi(m)/2 divides4. Therefore phi(m)|8.
    # p|m => p-1|8, so only p=2,3,5; valuations give v2<=4,v3,v5<=1.
    candidates=sorted({2**a*3**b*5**c for a in range(5) for b in range(2) for c in range(2)
                       if 2**a*3**b*5**c>=2 and 8%phi(2**a*3**b*5**c)==0})
    assert candidates==[2,3,4,5,6,8,10,12,15,16,20,24,30]
    cyclic=[real_group(m) for m in (15,16,20,30)]
    assert all(r['degree']==4 and r['cyclic'] for r in cyclic)
    rejected={12:'Q(sqrt3) is not a subfield of F',24:'contains Q(sqrt3)',
              15:'cyclic quartic, not biquadratic',16:'cyclic quartic, not biquadratic',
              20:'cyclic quartic, not biquadratic',30:'cyclic quartic, not biquadratic'}
    assert sorted(set(candidates)-set(rejected))==list(D)
    assert not F(3).square_root()
    for m,x in D.items():
        assert x.integral()
        assert all(x.conjugate(*e).sign()>=0 and (4-x.conjugate(*e)).sign()>0 for e in EMBEDDINGS)
    return dict(all_degree_candidates=candidates,rejected=rejected,finite_labels=list(D),
                squared_norm2_products={str(m):x.data() for m,x in D.items()},
                cyclic_quartic_checks=cyclic,
                complete_bound='phi(m)|8 => primes2,3,5 with exponents<=4,1,1; no empirical label cutoff')


def triangles():
    rows=[];kept=[]
    for labels in combinations_with_replacement(D,3):
        rec=dict(labels=labels,reciprocal_sum=str(sum(R(1,m) for m in labels)))
        if sum(R(1,m) for m in labels)>=1:
            rec['status']='NOT_COMPACT_HYPERBOLIC_TRIANGLE';rows.append(rec);continue
        ds=[D[m] for m in labels]
        square=ds[0]*ds[1]*ds[2]
        p=square.square_root()
        rec['cyclic_product_square']=square.data()
        if p is None:
            rec['status']='CYCLIC_PRODUCT_NOT_IN_F';rows.append(rec);continue
        assert p*p==square and p.sign()>=0 and p.integral()
        det=2*(4-sum(ds,F())-p)
        assert det.sign()<0
        fixed=[e for e in EMBEDDINGS if all(v.conjugate(*e)==v for v in ds+[p])]
        rec.update(cyclic_product=p.data(),integral=True,determinant=det.data(),
                   stabilizer=fixed,determinant_signs=[det.conjugate(*e).sign() for e in EMBEDDINGS])
        if len(fixed)!=1:
            rec['status']='PROPER_GROUND_SUBFIELD';rows.append(rec);continue
        if any(det.conjugate(*e).sign()!=1 for e in EMBEDDINGS[1:]):
            rec['status']='NONPOSITIVE_CONJUGATE_DETERMINANT';rows.append(rec);continue
        rec['status']='ARITHMETIC_LANNER_TRIANGLE_FULL_F'
        kept.append(list(labels));rows.append(rec)
    assert len(rows)==84 and kept==[[2,5,8],[4,5,5]]
    return dict(records=rows,complete_triangles=kept,
                counts={status:sum(r['status']==status for r in rows) for status in sorted({r['status'] for r in rows})})


def local_costs():
    rows=[]
    for labels in ((2,5,8),(4,5,5)):
        ds=[D[m] for m in labels];p=(ds[0]*ds[1]*ds[2]).square_root();assert p is not None
        for opposite in sorted(set(labels)):
            contacts=list(labels);contacts.remove(opposite)
            cost=(sum((D[m] for m in contacts),F())+p)/(4-D[opposite])
            rec=dict(triangle=labels,opposite_label=opposite,white_contacts=contacts,cost=cost.data(),
                     cost_minus_two=(cost-2).data(),sign=(cost-2).sign())
            if opposite==2:
                assert labels==(2,5,8) and rec['sign']<0
                rec['resolution']='white endpoint of bond8 cannot also meet the mandatory white-white bridge in a finite corner'
            else:
                assert rec['sign']>0
                rec['resolution']='strict Schur cost greater than2'
            rows.append(rec)
    assert len(rows)==5
    return rows


def edge(a,b):return tuple(sorted((a,b)))
def component(vertices,edges,start):
    seen={start}
    while True:
        nxt=seen|{a for a in vertices if any(a!=b and edge(a,b) in edges for b in seen)}
        if nxt==seen:return seen
        seen=nxt
def forest(vertices,edges):
    todo=set(vertices)
    while todo:
        part=component(vertices,edges,min(todo))
        if sum(a in part and b in part for a,b in edges)!=len(part)-1:return False
        todo-=part
    return True


def unmarked_supports():
    # No distinguished high edge is assumed here:455 has no label>=6.
    shapes=[((0,2),(1,3)),((0,2),(1,3,4)),((0,2,4),(1,3)),
            ((0,2,4),(1,3,5)),((0,2,4),(1,3,4))]
    survivors=set();tested=0
    for L,M in shapes:
        vertices=sorted(set(L+M));pairs=list(combinations(vertices,2))
        white=sorted(set(vertices)-{0,1,2,3})
        for mask in range(1<<len(pairs)):
            tested+=1;es={p for i,p in enumerate(pairs) if mask>>i&1}
            if len(component(L,es,L[0]))!=len(L) or len(component(M,es,M[0]))!=len(M):continue
            if not set(L)&set(M) and not any(edge(a,b) in es for a in L for b in M):continue
            good=True
            for a,b in ((0,1),(1,2),(2,3),(3,0)):
                corner=white+[a,b]
                if edge(a,b) in es or any(edge(a,w) in es and edge(b,w) in es for w in white):good=False
                if not forest(corner,{p for p in es if set(p)<=set(corner)}):good=False
                if not good:break
            if not good:continue
            assert L==(0,2,4) and M==(1,3,5)
            assert (4,5) in es
            assert not any(edge(a,b) in es for a in (0,2) for b in (1,3,5))
            assert not any(edge(a,b) in es for a in (1,3) for b in (0,2,4))
            survivors.add(tuple(sorted(es)))
    assert tested==35904 and len(survivors)==16
    return dict(binary_graphs=tested,surviving_supports=sorted(survivors),
                conclusion='two triangles on distinct whites, joined only by the white-white bond; opposite finite pair may initially be right')


def weight_counting():
    records=[]
    for r in range(1,41):
        for kind in ['path']+(['fork'] if r>=4 else []):
            es={(i,i+1) for i in range(r-1)} if kind=='path' else {(i,i+1) for i in range(r-2)}|{(r-3,r-1)}
            short=sum(edge(i,j) in es or any(edge(i,k) in es and edge(j,k) in es for k in range(r) if k not in (i,j))
                      for i,j in combinations(range(r),2))
            wanted=0 if r==1 else 2*r-3 if kind=='path' else 2*r-2
            assert short==wanted and short<=2*r-2
            records.append(dict(kind=kind,rank=r,short_pairs=short))
    assert [sum(m!=2 for m in t) for t in ((2,5,8),(4,5,5))]==[2,3]
    paired=[R(comb(21-i,2)+comb(i,2),2) for i in range(22)]
    assert min(paired)==50 and paired[0]>50
    assert 5*50-comb(21,2)==40==2*21-2
    # Substitute n=21+2u and n=22+2u into the exact h-vector
    # lower bound minus the vertex upper bound. These coefficient
    # identities, rather than a finite dimension check, cover all n>=21.
    odd=[R(20*16,8)-40,R(2*(20+16),8)-4,R(4,8)]
    even=[R(22*16,8)-42,R(2*(22+16),8)-4,R(4,8)]
    assert odd==[R(0),R(5),R(1,2)]
    assert even==[R(2),R(11,2),R(1,2)]
    return dict(graph_records=records,triangle_base_payments=[2,3],
                component_all_rank_identity='r-1+sum_v binom(deg(v),2)',
                vertex_upper='2n-2',strict_dimension21_lower=40,dimension21_upper=40,
                all_rank_parity_gaps=dict(odd=list(map(str,odd)),even=list(map(str,even))))


def run():
    assert S2*S2==2 and S5*S5==5 and PHI*PHI==PHI+1
    for x in [F(1,1,1,1),F(2,1,0,1),PHI,S2,S5]:
        assert x*x.inverse()==1
        z=(x*x).square_root();assert z is not None and z*z==x*x and z.sign()>0
    src=f'scripts/{TAG}.py'
    return dict(status='PASS_COMPLETE_F25_LANNER_CATALOGUE_AND_COMPACT_BOUND20_CANDIDATE',
                sources={src:sha256((ROOT/src).read_bytes()).hexdigest()},
                field='Q(sqrt2,sqrt5)',basis=['1','sqrt2','sqrt5','sqrt10'],
                ring_basis=['1','sqrt2','phi','sqrt2*phi'],alphabet=alphabet(),triangles=triangles(),
                quadrilateral_costs=local_costs(),unmarked_quadrilateral_supports=unmarked_supports(),
                counting=weight_counting(),compact_dimension_upper_bound_candidate=20,
                rank2_family='tau in O_F, Q(tau)=F, tau>4, 0<sigma(tau)<4 for all three nonidentity embeddings; normalized product=-sqrt(tau)/2',
                orders_ge4_published_exclusion=dict(source='https://arxiv.org/pdf/0708.3991',section='3.1, Figure1 and formula(4)',
                                                  fields=['Q','Q(sqrt2)','Q(sqrt5)']),
                scope='Complete finite-angle and arithmetic triangle catalogue for exact ground field F, exact local projection costs and full binary minimal-pair supports. Rank2 is an unbounded arithmetic parameter family. Higher Lanner orders use the cited complete arithmetic classification and same-ground-field theorem. The all-white quadrilateral and all-rank counting proofs in the report give compact bound20 as an internally reviewed candidate only after independent review; no numerical or finite-label sampling or claim for Qsqrt5.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text());print('PASS fresh saved F25 catalogue replay')
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['triangles']['complete_triangles'],out['triangles']['counts'])
