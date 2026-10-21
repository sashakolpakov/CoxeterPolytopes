#!/usr/bin/env python3
"""Exact arithmetic supporting the all-label octic-field angle theorem.

The universal subfield argument is in the matching report.  This checker uses
only stdlib rational quadratic arithmetic; it does not certify any new cap.
"""
import argparse,json,math
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_octic_angle_20261015'
SOURCE='results/fresh_h9_after_fourth_rescore_20261015.json'

class F:
    def __init__(self,a=0,b=0):self.a,self.b=Q(a),Q(b)
    def __add__(self,o):
        o=o if isinstance(o,F)else F(o);return F(self.a+o.a,self.b+o.b)
    __radd__=__add__
    def __neg__(self):return F(-self.a,-self.b)
    def __sub__(self,o):return self+-o if isinstance(o,F)else self+F(-o)
    def __mul__(self,o):
        o=o if isinstance(o,F)else F(o)
        return F(self.a*o.a+5*self.b*o.b,self.a*o.b+self.b*o.a)
    __rmul__=__mul__
    def conj(self):return F(self.a,-self.b)
    def norm(self):return self.a*self.a-5*self.b*self.b
    def __truediv__(self,o):
        o=o if isinstance(o,F)else F(o);z=self*o.conj();n=o.norm()
        return F(z.a/n,z.b/n)
    def __pow__(self,n):
        assert n>=0;r=F(1)
        for _ in range(n):r=r*self
        return r
    def __eq__(self,o):
        o=o if isinstance(o,F)else F(o);return (self.a,self.b)==(o.a,o.b)
    def sign(self):
        if not self.b:return (self.a>0)-(self.a<0)
        if not self.a:return (self.b>0)-(self.b<0)
        if self.a*self.b>0:return (self.a>0)-(self.a<0)
        z=self.a*self.a-5*self.b*self.b
        return ((z>0)-(z<0))*((self.a>0)-(self.a<0))
    def out(self):return [str(self.a),str(self.b)]

def square(q):
    if q<0:return False
    a,b=math.isqrt(q.numerator),math.isqrt(q.denominator)
    return a*a==q.numerator and b*b==q.denominator
def phi(n):return sum(math.gcd(i,n)==1 for i in range(1,n+1))
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()

def run():
    src=json.loads((ROOT/SOURCE).read_text())
    for p,h in src['sources'].items():assert sha(p)==h,p
    D=F(6,Q(16,5));d3=F(50,22);d5=F(Q(245,4),Q(111,4));gold=F(Q(1,2),Q(1,2))
    R20=F(10,2);R15=F(30,-6);R5=D*d5
    assert d3==gold**4*R20
    assert R15/R20==F(3)/(gold**2)
    assert R5==F(Q(1623,2),Q(725,2))and R5.norm()==1501
    assert [x.sign()for x in(D,d3,d5)]==[1,1,1]
    assert [x.conj().sign()for x in(D,d3,d5)]==[-1,1,-1]
    rad=[D,d3,d5];independence=[]
    for mask in range(1,8):
        x=F(1)
        for j in range(3):
            if mask>>j&1:x=x*rad[j]
        n=x.norm();assert not square(n)
        independence.append(dict(mask=mask,radicand=x.out(),norm=str(n),norm_is_rational_square=False))
    # Every rational square in F is either a rational square or five times one.
    nongalois=[]
    for name,x in [('R5',R5),('delta3_R5',d3*R5)]:
        n=x.norm();assert not square(n)and not square(n/5)
        nongalois.append(dict(radicand_name=name,norm=str(n),norm_over_5=str(n/5),
            norm_is_square_in_F=False))
    assert not square(Q(3))and not square(Q(3,5))
    # phi(m)^2 >= m/2 bounds all possible cyclotomic degrees <=8 by m<=512.
    degree_candidates={str(d):[]for d in [1,2,4,8]}
    for m in range(2,513):
        t=phi(m);assert 2*t*t>=m
        deg=1 if m==2 else t//2
        if deg in [1,2,4,8]:degree_candidates[str(deg)].append(m)
    assert degree_candidates=={'1':[2,3,4,6],'2':[5,8,10,12],
        '4':[15,16,20,24,30],'8':[17,32,34,40,48,60]}
    records=src['exact_nine_contact_recipes']
    assert [r['delta']for r in records]==[[d3.out(),['0','0']],[d5.out(),['0','0']]]
    common={2:F(0),3:F(1),4:F(2),5:(F(3,1)/2),6:F(3),10:F(5,1)/2}
    for m,a in common.items():
        assert a.sign()>=0 and (F(4)-a).sign()>0
    return dict(status='PASS_EXACT_ALL_LABEL_FIELD_GATE_NO_NEW_CAP_CERTIFICATION',
        sources={p:sha(p)for p in [SOURCE,f'scripts/{TAG}.py']},
        field_encoding='[a,b] means a+b*sqrt5',D=D.out(),delta3=d3.out(),delta5=d5.out(),
        physical_radicand_signs=[1,1,1],conjugate_radicand_signs=[-1,1,-1],
        seven_squareclass_independence_checks=independence,
        quartic_nongalois_checks=nongalois,
        identities=dict(delta3='phi^4*(10+2sqrt5)',R5=R5.out(),
            conductor15_over_conductor20=(R15/R20).out()),
        maximum_totally_real_fields=dict(L3='F(sqrt(delta3)) = Q(zeta20+zeta20^-1)',
            L5='F(sqrt(D*delta5))',compositum='F(sqrt(delta3),sqrt(D*delta5))'),
        finite_label_candidates_by_degree=degree_candidates,
        complete_finite_alphabets=dict(old=[2,3,4,5,6,10],delta3=[2,3,4,5,6,10,20],
            delta5=[2,3,4,5,6,10],mutual_compositum=[2,3,4,5,6,10,20]),
        cartan_two_cycle_formula='A=4*(r.s)^2/(r^2*s^2)=4*cos(pi/m)^2; also require r.s<=0',
        common_cycle_values={str(m):a.out()for m,a in common.items()},
        label20_cycle='2 + sqrt(delta3)/(2*phi^2)',parallel_cycle=4,
        scope='Complete necessary field alphabet, with every listed value present. Actual pair compatibility, signs, facet incidence and finite volume are separate geometric checks. No arithmetic-lattice hypothesis.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    data=run();out=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert data==json.loads(out.read_text());print('PASS saved exact all-label octic field gate')
    else:out.write_text(json.dumps(data,indent=2)+'\n');print(data['status'])
