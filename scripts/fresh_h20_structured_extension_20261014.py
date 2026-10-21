#!/usr/bin/env python3
"""All-real escape for four frozen paired curves and their height flips.

No producer imports, numerical acceptance or finite angle-label cutoff.
Rational interval bounds certify a finite set of analytic monotonicity gates.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from math import isqrt
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_structured_extension_20261014'
SCALE=1<<160


class I:
    def __init__(self,a,b=None):
        self.a=F(a);self.b=F(a if b is None else b)
        assert self.a<=self.b
    def __add__(self,v):
        v=as_i(v);return I(self.a+v.a,self.b+v.b)
    __radd__=__add__
    def __neg__(self):return I(-self.b,-self.a)
    def __sub__(self,v):return self+-as_i(v)
    def __rsub__(self,v):return as_i(v)+-self
    def __mul__(self,v):
        v=as_i(v);z=[self.a*v.a,self.a*v.b,self.b*v.a,self.b*v.b]
        return I(min(z),max(z))
    __rmul__=__mul__
    def __truediv__(self,v):
        v=as_i(v);assert v.a>0 or v.b<0
        return self*I(1/v.b,1/v.a)
    def abs(self):
        return I(0 if self.a<=0<=self.b else min(abs(self.a),abs(self.b)),max(abs(self.a),abs(self.b)))
    def enc(self):return [str(self.a),str(self.b)]


def as_i(x):return x if isinstance(x,I) else I(x)
def sqrtq(x):
    x=F(x);assert x>0
    k=isqrt(x.numerator*SCALE*SCALE//x.denominator)
    assert F(k*k,SCALE*SCALE)<=x<F((k+1)**2,SCALE*SCALE)
    return I(F(k,SCALE),F(k+1,SCALE))
def mv(A,x):return[sum(a*b for a,b in zip(row,x)) for row in A]
def dot(x,y):return sum(a*b for a,b in zip(x,y))


def run():
    paths=['results/fresh_h20_coupled_lift_20261014.json',
           'results/fresh_dimension_construction_routes_20261014.json',
           f'scripts/{TAG}.py']
    d=json.loads((ROOT/paths[0]).read_text());old=json.loads((ROOT/paths[1]).read_text())
    for src in (d,old):
        for p,h in src['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    Q=[[F(x)for x in row]for row in d['form']]
    ip=lambda x,y:dot(x,mv(Q,y))
    H=list(map(F,d['strict_center']));X=list(map(F,d['strict_escape']['point']))
    v=[F(z['a'])for z in old['strict_escape']['boundary_point']]
    assert all(F(z['b'])==0 for z in old['strict_escape']['boundary_point'])
    assert ip(v,v)==25
    E=[x/5 for x in v];assert ip(E,E)==1 and ip(E,H)<0 and ip(H,H)==-4296
    assert ip(X,X)==F(d['strict_escape']['norm'])>0 and ip(X,H)<0
    assert all(F(r['d'])==2 for r in d['normals'][:60])
    roots=[[I(F(a))+sqrtq(r['d'])*F(b)for a,b in zip(r['a'],r['b'])]for r in d['normals'][:60]]
    assert all(ip(r,X).b<0 for r in roots)
    pre=[7,12,9,4];main=[29,1,29,20];neg=[1,1,20,20]
    lows=[F(3689,2000),F(11157,2500),F(3689,2000),F(11157,2500)]
    final_lows=[F(5),lows[1],F(5),lows[3]]
    frames=[];base_frames=[];range_cert=[]
    for i,rec in enumerate(d['base_recipes']):
        S=rec['support'];h=list(map(F,rec['h0']));w=list(map(F,rec['w']))
        z=[a-ip(h,E)*b for a,b in zip(h,E)]
        qz=-ip(z,z);qw=ip(w,w)
        assert qz>0 and qw>0 and ip(E,z)==ip(E,w)==ip(z,w)==0 and ip(z,H)<0 and ip(w,H)==0
        for j in S:
            r=d['normals'][j];assert all(F(x)==0 for x in r['b'])
            a=list(map(F,r['a']));assert ip(a,E)==ip(a,z)==ip(a,w)==0
        T=[I(a)/sqrtq(qz)for a in z];W=[I(a)/sqrtq(qw)for a in w]
        base_frames.append((E,T,W))
        # Both sheets of the unit spacelike hyperbola are explicitly covered.
        # On the negative-E sheet even the pair's mean has a positive old pairing.
        aa=ip([-a+b for a,b in zip(E,T)],roots[neg[i]])
        cc=ip([-a-b for a,b in zip(E,T)],roots[neg[i]])
        assert aa.a>0 and cc.a>0
        pol=[]
        for j,L0,L1 in ((pre[i],F(0),F(1)),(main[i],F(1),lows[i])):
            a=ip([x+y for x,y in zip(E,T)],roots[j])
            b=2*ip(W,roots[j]).abs()
            c=ip([x-y for x,y in zip(E,T)],roots[j])
            left=a*L0*L0+b*L0+c;right=a*L1*L1+b*L1+c
            assert a.b<0 and left.a>0 and right.a>0
            pol.append(dict(wall=j,interval=[str(L0),str(L1)],a=a.enc(),b=b.enc(),c=c.enc(),left=left.enc(),right=right.enc()))
        gap=None
        if i in (0,2):
            j=1 if i==0 else 20
            assert F(d['norms'][j])==2
            e=ip(E,roots[j]);t=ip(T,roots[j]);w0=ip(W,roots[j])
            assert e.b<0 and t.b<0
            def value(p):return e*((p+1/p)/2)+t*((p-1/p)/2)-w0
            low=value(lows[i]);high=value(F(5))
            assert low.b<0 and high.a>-1
            gap=dict(wall=j,old_norm='2',at_acute_lower=low.enc(),at_five=high.enc(),E_pair=e.enc(),T_pair=t.enc(),
                     conclusion='For the minus cap and every parameter between the acute lower bound and5, pairing lies strictly in(-1,0), so normalized cosine lies strictly in(0,1/2).')
        L=final_lows[i];c=(L+1/L)/2;s=(L-1/L)/2
        M=[c*a+s*b for a,b in zip(E,T)];U=[s*a+c*b for a,b in zip(E,T)]
        frames.append((M,U,W))
        range_cert.append(dict(support=i,negative_sheet_wall=neg[i],negative_sheet_coefficients=[aa.enc(),cc.enc()],
                               concave_positive_polynomials=pol,cosine_gap=gap,necessary_parameter_lower=str(L),
                               exact_frame=dict(E_norm='1',T_norm='-1',W_norm='1',all_cross='0',T_raw_norm=str(-qz),W_raw_norm=str(qw))))

    transport=[];family=[]
    for signs in product((1,-1),repeat=2):
        D=[1]*19+list(signs)
        assert [x*y for x,y in zip(H,D)]==H
        assert all(D[i]*D[j]*Q[i][j]==Q[i][j]for i in range(21)for j in range(21))
        flip=lambda a:[x*y for x,y in zip(a,D)]
        oldkeys={(tuple(r['a']),tuple(r['b']))for r in d['normals'][:60]}
        for r in d['normals'][:60]:
            aa=[str(F(x)*y)for x,y in zip(r['a'],D)];bb=[str(F(x)*y)for x,y in zip(r['b'],D)]
            assert (tuple(aa),tuple(bb))in oldkeys
        for i,(M,U,W)in enumerate(frames):
            M,U,W=map(flip,(M,U,W));m=ip(M,X);u=ip(U,X);w=ip(W,X)
            # cosh(v)>=1, sinh(v)>=0: these inequalities control the full
            # unbounded v=log(parameter/lower)>=0 interval, not samples.
            assert (m+w.abs()).b<0 and u.b<0
            transport.append(dict(signs=list(signs),support=i,mean_pairing=m.enc(),derivative_pairing=u.enc(),
                                  spatial_pairing=w.enc(),uniform_upper=(m+w.abs()).enc()))
            family.append((M,U,W))
    # A supplementary monotonicity proof: any different named family has
    # strictly dotted cross pairs throughout these necessary parameter tails.
    joins=[]
    for i,j in combinations(range(16),2):
        M,U,W=family[i];N,V,Z=family[j]
        tests=[ip(M,N)+ip(M,Z).abs(),ip(U,N)+ip(U,Z).abs(),ip(M,V),ip(U,V),
               ip(N,M)+ip(N,W).abs(),ip(V,M)+ip(V,W).abs(),ip(N,U),ip(V,U)]
        assert all(x.b<0 for x in tests),(i,j)
        pairs=[]
        for a,b in product((1,-1),repeat=2):
            p=ip([x+a*y for x,y in zip(M,W)],[x+b*y for x,y in zip(N,Z)])
            assert p.b < -2,(i,j,a,b)
            pairs.append(dict(signs=[a,b],pair=p.enc()))
        joins.append(dict(families=[i,j],derivative_gates=[x.enc()for x in tests],corner_pairs=pairs))
    # A concrete acute but non-Coxeter cutter shows why the cosine gap is
    # essential on the reflected complements: mere acuteness does not suffice.
    E0,T0,W0=base_frames[0];p=F(37,20);c=(p+1/p)/2;s=(p-1/p)/2
    D=[1]*19+[-1,1];perm=d['height_sign_permutations'][2]
    inc={perm[j]for j in d['base_recipes'][0]['support']}
    acute_examples=[]
    for eps in (1,-1):
        r=[d*(c*e+s*t+eps*w)for e,t,w,d in zip(E0,T0,W0,D)]
        pairs=[ip(r,oldr)for oldr in roots]
        assert all(pairs[j].b<0 for j in range(60)if j not in inc)
        # Incidences follow the exact orthogonal frame identities checked above.
        cut=ip(r,X)
        if cut.a>0:
            bad=pairs[perm[1]];assert bad.a>-1 and bad.b<0
            acute_examples.append(dict(support=0,height_signs=[-1,1],parameter=str(p),cap_sign=eps,
                                       escape_pairing=cut.enc(),forbidden_old_wall=perm[1],forbidden_pairing=bad.enc(),
                                       strict_nonincident_old_pairings=[dict(wall=j,bound=pairs[j].enc())for j in range(60)if j not in inc]))
    assert acute_examples
    return dict(status='PASS_ALL_REAL_PAIRED_CURVES_SHARE_STRICT_ESCAPE',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
                interval_scale=SCALE,range_certificates=range_cert,transported_families=transport,
                cross_family_dotted=joins,acute_noncoxeter_cutters=acute_examples,old60_escape_strict=True,
                common_strict_center=X_to_strings(H),center_norm=str(ip(H,H)),
                escape=X_to_strings(X),escape_norm=str(ip(X,X)),
                scope='Four fixed-W normalized paired curves determined by the frozen rational-H recipes; all real hyperbola parameters, both sheets and four height flips. Each selected curve parameter contributes BOTH orthogonal cap normals and retains old60 with complete nonpositive Coxeter pair gates. Arbitrary finite joint selections/replacements within this family cannot have finite volume. Single unpaired caps, rotations of W, different supports, or replacement of old60 are not excluded.')


def X_to_strings(X):return list(map(str,X))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run()
    dest=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(dest.read_text())
    else:dest.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],len(d['transported_families']),'families',len(d['cross_family_dotted']),'cross-family gates')
