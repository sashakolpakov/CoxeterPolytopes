#!/usr/bin/env python3
"""A restricted quadrilateral pattern is impossible over Q(sqrt5).

This checks explicit polynomial identities and arithmetic cases. The all-rank
inequalities and the geometric domain are stated in the accompanying report.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import permutations
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_long_white_two_angle_quad_gate_20261016'

def add(a,b):
    d=a.copy()
    for k,v in b.items():
        d[k]=d.get(k,F(0))+v
        if not d[k]:del d[k]
    return d

def mul(a,b):
    d={}
    for i,x in a.items():
        for j,y in b.items():
            k=tuple(u+v for u,v in zip(i,j));d[k]=d.get(k,F(0))+x*y
    return {k:v for k,v in d.items()if v}

def scale(a,c):return {k:v*c for k,v in a.items()if v*c}
def constant(x):return {(0,0,0):F(x)}if x else{}
def encode(p):return [[list(k),str(v)]for k,v in sorted(p.items())]
def det(M):
    d={}
    for perm in permutations(range(len(M))):
        term=constant((-1)**sum(perm[i]>perm[j]for i in range(len(M))for j in range(i+1,len(M))))
        for i,j in enumerate(perm):term=mul(term,M[i][j])
        d=add(d,term)
    return d

def phi_mul(x,y):
    a,b=x;c,d=y;return(a*c+b*d,a*d+b*c+b*d)
def phi_add(x,y):return tuple(a+b for a,b in zip(x,y))
def phi_inv(x):
    a,b=x;n=a*a+a*b-b*b;assert n
    return((a+b)/n,-b/n)
def phi_scale(x,c):return tuple(c*a for a in x)
def phi_div(x,y):return phi_mul(x,phi_inv(y))
def ring_integral(x):return all(a.denominator==1 for a in x)

def run():
    a={(1,0,0):F(1)};p={(0,1,0):F(1)};q={(0,0,1):F(1)};two=constant(2);minus_one=constant(-1);zero={}
    M=[[a,minus_one,scale(p,-1),minus_one],[minus_one,two,zero,scale(q,-1)],
       [scale(p,-1),zero,two,zero],[minus_one,scale(q,-1),zero,two]]
    determinant=det(M)
    wanted=mul(add(two,q),add(mul(add(scale(a,2),scale(mul(p,p),-1)),add(two,scale(q,-1))),constant(-4)))
    assert determinant==wanted
    # phi(n)<=4 is exhausted by prime divisors2,3,5 and these exponent bounds.
    cyclotomic=[]
    for e2 in range(4):
        for e3 in range(2):
            for e5 in range(2):
                n=2**e2*3**e3*5**e5
                tot=(1 if e2==0 else 2**(e2-1))*(1 if e3==0 else 2)*(1 if e5==0 else 4)
                if tot<=4:cyclotomic.append(n)
    assert sorted(cyclotomic)==[1,2,3,4,5,6,8,10,12]
    label5=[]
    for m in range(5,11):
        # z=1+phi; q=2+4m/[m phi-m-2].
        exact=phi_add((F(2),F(0)),phi_div((F(4*m),F(0)),(F(-m-2),F(m))))
        D=m*m-2*m-4;formula=(F(2)+F(8*m,D),F(4*m*m,D))
        assert exact==formula and not ring_integral(exact)
        label5.append(dict(m=m,q=list(map(str,exact)),integral=False))
    # m>=11: D-8m=(m-11)(m+1)+7>0, so 0<8m/D<1.
    assert [1,-10,-4]==[1,1-11,7-11]
    # For label10, conjugate PSD requires phi>=2-2/m. At m>=6,
    # RHS>=5/3>phi; the last strict comparison uses x²-x-1>0 at x=5/3.
    assert F(5,3)**2-F(5,3)-1==F(1,9)>0
    q10=phi_add((F(2),F(0)),phi_div((F(20),F(0)),(F(-2),F(5))))
    assert q10==(F(2,31),F(100,31))and not ring_integral(q10)
    # Dotted/parallel p>=2: |N(t)|<=8m/(m-1). Odd m>=5 has
    # m²|N(t), while even m>=8 has (m/2)²|N(t); both are impossible.
    assert F(8*5,5-1)==10<25 and F(8*8,8-1)<10<16
    # Only m6 remains: t=3s, N(t) multiple9, |N(t)|<=48/5<18.
    # Thus N(s)=-1, but 3s(z-2)=s+4 gives s=-4 mod3 and N(s)=16 mod3.
    assert 9<=F(48,5)<18 and (-1)%3!=16%3
    paths=[f'scripts/{TAG}.py']
    return dict(status='PASS_LONG_WHITE_TWO_ANGLE_QUADRILATERAL_EXCLUSION',
        sources={x:sha256((ROOT/x).read_bytes()).hexdigest()for x in paths},
        domain=dict(white='A_r, r>=4, attached by one3-edge at its endpoint only to black0',
            black_adjacent_labels={'01':3,'03':3,'12':2,'23':2},black_opposite_products={'02':'-p','13':'-q'},
            spectators='All other common white components orthogonal to this local graph',
            assumptions=['Acute compact Coxeter quadrilateral face','Physical Gram has Lorentz index1','Cyclic products integral in Z[phi]',
                'Conjugate positive semidefinite','Local Gram corank at least1']),
        determinant=dict(variable_order=['a','p','q'],polynomial=encode(determinant),factorization='(2+q)*((2*a-p*p)*(2-q)-4)'),
        schur_diagonal='a=1+1/m, m=r+1>=5',relation='t=q-2=4/(p^2-2-2/m)>0',
        integral_relation='m*t*(p^2-2)=2*t+4*m, hence m divides2*t in Z[phi]',
        conjugate_bound='-4<=t_conjugate<0',
        physical_q_bound='q<2 would make the (1,3) block positive definite, and corank would then make the whole Gram positive semidefinite; q=2 contradicts the determinant. Thus q>2.',
        finite_label_domain=dict(cyclotomic_totient_at_most4=sorted(cyclotomic),
            excluded_quadratic_fields={'8':'Q(sqrt2)','12':'Q(sqrt3)'},remaining_labels=[2,3,4,5,6,10]),
        dotted_branch=dict(norm_upper_bound='8*m/(m-1)',odd_m_at_least5='m^2 divisibility exceeds bound<=10',
            even_m_at_least8='(m/2)^2 divisibility exceeds bound<10',m6='Norm(s)=-1 contradicts s=-4 modulo3'),
        finite_branches=dict(labels_2_3_4='p^2<=2<2+2/m',label6='q rational>2 contradicts conjugate PSD',
            label5_small_m=label5,label5_large_m='m>=11: 0<8m/(m^2-2m-4)<1',
            label10_m5_q=list(map(str,q10)),label10_large_m='m>=6: necessary phi>=2-2/m>=5/3 is false'),
        consequence='This specific two-discounted-corner quadrilateral is impossible; it cannot impose the proposed unconditional adjacent-discount cap1/2.',
        scope='Restricted geometric pattern only. No claim that all quadrilaterals, all white attachments, or a proposed weight rule are covered. No dimension-bound improvement.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(p.read_text())
    else:p.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'])
