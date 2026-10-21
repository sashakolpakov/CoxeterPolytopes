#!/usr/bin/env python3
"""Independent literal witness cross-check for the reviewed two-block theorem.

Standard-library only. Reconciles every primary symbolic polynomial with the
separate Fraction determinant output, recomputes Bernstein coefficients, and
certifies their signs using integer square-root enclosures. The full category
census is freshly replayed by the separately authored block_independent.py;
this file does not claim a third census or determinant implementation.
"""
import argparse
import ast
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, gcd, isqrt
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_bugaenko_remaining_20261014_block_gate_review'
PRIMARY='results/fresh_h9_bugaenko_remaining_20261014_block_gate.json'
SEPARATE='results/fresh_h9_bugaenko_remaining_20261014_block_independent.json'

# Polynomials in x, with coefficients in the squarefree-radical basis.
def clean(a): return {k:v for k,v in a.items() if v}
def const(a): return {(0,1):F(a)} if a else {}
def add(a,b):
    out=dict(a)
    for k,v in b.items(): out[k]=out.get(k,F(0))+v
    return clean(out)
def scale(a,c): return clean({k:v*c for k,v in a.items()})
def mul(a,b):
    out={}
    for (i,r),v in a.items():
        for (j,s),w in b.items():
            g=gcd(r,s); k=(i+j,r*s//(g*g))
            out[k]=out.get(k,F(0))+v*w*g
    return clean(out)
def parse(text):
    def read(n):
        if isinstance(n,ast.Constant):
            assert type(n.value) is int
            return const(n.value)
        if isinstance(n,ast.Name):
            assert n.id=='x'; return {(1,1):F(1)}
        if isinstance(n,ast.UnaryOp):
            assert isinstance(n.op,(ast.USub,ast.UAdd))
            return scale(read(n.operand),-1 if isinstance(n.op,ast.USub) else 1)
        if isinstance(n,ast.Call):
            assert isinstance(n.func,ast.Name) and n.func.id=='sqrt' and len(n.args)==1
            r=n.args[0]; assert isinstance(r,ast.Constant) and type(r.value) is int
            assert r.value in [2,3,5,6,10,15,30]
            return {(0,r.value):F(1)}
        assert isinstance(n,ast.BinOp)
        a,b=read(n.left),read(n.right)
        if isinstance(n.op,ast.Add): return add(a,b)
        if isinstance(n.op,ast.Sub): return add(a,scale(b,-1))
        if isinstance(n.op,ast.Mult): return mul(a,b)
        if isinstance(n.op,ast.Div):
            assert set(b)=={(0,1)}; return scale(a,1/b[0,1])
        assert isinstance(n.op,ast.Pow) and set(b)=={(0,1)}
        exponent=b[0,1]; assert exponent.denominator==1 and 0<=exponent<=20
        out=const(1)
        for _ in range(int(exponent)): out=mul(out,a)
        return out
    return read(ast.parse(text,mode='eval').body)
def decode_coeffs(rows):
    return clean({(i,int(r)):F(c) for i,row in enumerate(rows) for r,c in row.items()})
def negative(a):
    assert all(i==0 for i,r in a)
    lo=hi=F(0); D=10**80
    for (i,r),c in a.items():
        k=isqrt(r*D*D); l=F(k,D); h=F(k+(k*k<r*D*D),D)
        lo+=c*(l if c>=0 else h); hi+=c*(h if c>=0 else l)
    assert hi<0, (a,lo,hi)
    return str(hi)
def bernstein(a):
    degree=max((i for i,r in a),default=0)
    power=[{} for _ in range(degree+1)]
    for (k,r),c in a.items():
        for j in range(k+1):
            term={(0,r):c*comb(k,j)*F(9,10)**(k-j)/10**j}
            power[j]=add(power[j],term)
    return [sum_polys(scale(power[j],F(comb(i,j),comb(degree,j)))
                      for j in range(i+1)) for i in range(degree+1)]
def sum_polys(rows):
    out={}
    for row in rows: out=add(out,row)
    return out
def key(orders,a,b,values): return tuple(map(tuple,[orders,a,b,values]))
def run():
    primary=json.loads((ROOT/PRIMARY).read_text())
    separate=json.loads((ROOT/SEPARATE).read_text())
    dependencies={PRIMARY,SEPARATE}
    for data in [primary,separate]:
        for p,h in data['sources'].items():
            assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
            dependencies.add(p)
    independent={key(r['orders'],r['first_code'],r['second_code'],r['cross_values']):r
                 for r in separate['schur_proofs']}
    assert len(independent)==42
    records=[];seen=set();parametric=0;coefficients=0
    for table in primary['tables']:
        for row in table['internal_pairs_allowing_nonright']:
            for c in row['completions']:
                k=key(table['orders'],row['first_code'],row['second_code'],c['cross_values'])
                assert k not in seen; seen.add(k)
                other=independent[k]
                assert max(c['cross_values'])<=6
                assert sum(v==7 for v in row['first_code']+row['second_code'])<=1
                p=parse(c['schur']['polynomial'])
                assert p==decode_coeffs(other['polynomial_coefficients'])
                degree=max(i for i,r in p)
                assert degree==c['schur']['degree']
                if degree:
                    assert c['schur']['interval']==['9/10','1']; parametric+=1
                bb=bernstein(p)
                assert bb==[parse(v) for v in c['schur']['negative_bernstein_coefficients']]
                assert bb==[decode_coeffs([v]) for v in other['negative_bernstein_coefficients']]
                upper=[negative(v) for v in bb]; coefficients+=len(bb)
                records.append(dict(orders=table['orders'],first_code=row['first_code'],
                    second_code=row['second_code'],cross_values=c['cross_values'],
                    degree=degree,strict_negative_upper_bounds=upper))
    assert seen==set(independent)
    # T7(x)+1=(x+1)*(8x^3-4x^2-4x+1)^2; the largest root is cos(pi/7).
    cubic=parse('8*x**3-4*x**2-4*x+1')
    assert mul(parse('x+1'),mul(cubic,cubic))==parse('64*x**7-112*x**5+56*x**3-7*x+1')
    x=F(9,10); low=8*x**3-4*x*x-4*x+1
    assert low==F(-1,125) and 24*x*x-8*x-4==F(206,25)
    assert 48*x-8>0  # derivative is positive and increasing on[.9,1].
    assert sum(t['internal_type_pairs'] for t in primary['tables'])==767
    assert primary['excluded_catalogue_types']==separate['applications_checked']==2365
    assert len({r['case'] for r in primary['applications']})==2365
    dependencies.add(str(Path(__file__).resolve().relative_to(ROOT)))
    return dict(status='PASS_SOURCE_SCOPE_AND_LITERAL_TWO_BLOCK_WITNESS_REVIEW',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(dependencies)},
        internal_type_pairs=767,nonright_completions=42,parametric_completions=parametric,
        bernstein_coefficients_checked=coefficients,catalogue_applications=2365,
        cosine_seven_lower_bound=dict(polynomial_at_nine_tenths=str(low),
            derivative_lower_bound='206/25',chebyshev_identity=True),
        witnesses=records,
        scope='Source/proof review and fresh primary plus separately authored Fraction replays. '
              'This additional stdlib checker independently parses/reconciles all42 determinant '
              'polynomials, recomputes Bernstein coefficients and encloses their signs. '
              'It is not a third category census/determinant engine or a review of the remaining588 searches.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=run();out=ROOT/'results'/f'{TAG}.json'
    if args.verify: assert data==json.loads(out.read_text())
    else: out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['nonright_completions'],data['bernstein_coefficients_checked'])
