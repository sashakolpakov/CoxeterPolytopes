#!/usr/bin/env python3
"""Exact determinant triage and a bound-independent end-weight obstruction."""
import argparse
from itertools import product
import json
from fresh_h9_changed_form_20261013_independent import f,Census,Z,O,T,F,ROOT,sq,AI
OUTPUT=ROOT/'results/fresh_h9_changed_form_20261013_structure.json'

def valuation(n,p):
    n=abs(n);assert n
    k=0
    while n%p==0:n//=p;k+=1
    return k

def run(verify=False):
    B=f.BETA;D=f.D
    assert f.field_norm(T)==1 and abs(f.field_norm(D))==1 and f.field_norm(B)==3
    assert f.mul(B,sq(B))==f.scale(sq(f.divide(D,T)),3)
    # Determinant-only triage: equality of squareclasses is not an isometry.
    a5_ratio=f.divide(f.scale(O,3),f.mul(B,sq(B)))
    assert a5_ratio==sq(f.divide(T,D))
    d5_ratio=f.divide(f.scale(O,4),sq(B));assert d5_ratio==sq(f.divide(f.TWO,B))
    alpha=f.scale(T,2);det=f.scale(f.mul(alpha,B),-14);norm=f.field_norm(det)
    assert norm.denominator==1 and valuation(norm.numerator,7)==3
    comparisons=[]
    for a in(T,f.scale(T,2),f.mul(T,B),f.scale(f.mul(T,B),2)):
        old=f.neg(f.mul(a,sq(sq(B))))
        comparisons.append(('B^4+A1',a,old))
    for a in(T,f.scale(T,2)):
        old=f.scale(f.mul(a,sq(B)),-10)
        comparisons.append(('B^2+A4+A1',a,old))
    rows=[]
    for host,a,old in comparisons:
        old_norm=f.field_norm(old);ratio=f.field_norm(f.divide(det,old))
        assert valuation(old_norm.numerator,7)==0 and old_norm.denominator==1
        assert abs(valuation(ratio.numerator,7)-valuation(ratio.denominator,7))%2==1
        rows.append({'host':host,'alpha':f.encoded(a),'determinant':f.encoded(old),
                     'determinant_norm':str(old_norm),'new_to_old_norm_ratio':str(ratio)})
    mod7=[(x**3-3*x-1)%7 for x in range(7)];assert all(mod7)
    mod7field=lambda x:tuple(F(int(t)%7)for t in x)
    def power7(x,n):
        out=O
        for _ in range(n):out=mod7field(f.mul(out,x))
        return out
    local7=[]
    for q in(O,f.TWO,f.FOUR,B,f.scale(B,2),f.scale(B,4)):
        symbol=power7(f.neg(q),171)
        expected=O if q in(B,f.scale(B,2),f.scale(B,4))else f.scale(O,6)
        assert symbol==expected
        local7.append({'q_divided_by7':f.encoded(q),'minus_q_euler_symbol':f.encoded(symbol),
                       'locally_possible':symbol==O})
    # Unit signatures span all eight signs, so TP units are squares.
    signatures=[tuple(f.sign(u,e)for e in range(3))for u in(T,D,f.neg(O))]
    span=set()
    for powers in product((0,1),repeat=3):
        span.add(tuple(__import__('math').prod(signatures[i][j]**powers[i]for i in range(3))for j in range(3)))
    assert len(span)==8
    # The H9 escape lies in the saturated binary span(time,7omega6).
    source=ROOT/'results/fresh_h9_changed_form_20261013_h9_2theta_h100.json'
    data=json.loads(source.read_text());c=Census(data)
    inv=tuple(tuple(f.neg(f.inverse(alpha))if i==j==0 else
        f.divide(f.TWO if i==j else T,B)if 1<=i<=2 and 1<=j<=2 else
        AI[i-3][j-3]if 3<=i<=8 and 3<=j<=8 else f.scale(O,F(1,2))if i==j==9 else Z
        for j in range(10))for i in range(10))
    assert all(f.dot(c.G[i],tuple(inv[k][j]for k in range(10)))==(O if i==j else Z)for i in range(10)for j in range(10))
    exponent=f.scale(B,14)
    assert all(f.integral(f.mul(exponent,x))for row in inv for x in row)
    for prime in(f.TWO,B,f.scale(O,7)):
        shortened=f.divide(exponent,prime)
        assert any(not f.integral(f.mul(shortened,x))for row in inv for x in row)
    v=f.vecdecode(data['geometry']['hyperideal_witness']['vector'])
    direction=(Z,Z,Z,*(f.scale(O,j)for j in range(1,7)),Z)
    assert c.norm(direction)==f.scale(O,42)
    assert all(v[i]==Z for i in(1,2,9))and all(v[3+j]==f.scale(v[3],j+1)for j in range(6))
    assert f.sign(c.norm(v))>0 and all(f.sign(f.dot(v,c.cov(f.vecdecode(r))))<=0 for r in data['all_walls'])
    # At beta, T maps to1 in F3. For q2 andq14 the divided norm equation
    # is21t^2-Tx^2=1 or7, giving-x^2=1 inF3: impossible.
    residue3=lambda a:sum(int(t)for t in a)%3
    assert residue3(B)==0 and residue3(T)==1
    odd_norm_table=[(-x*x)%3 for x in range(3)];assert 1 not in odd_norm_table
    # For q4 orq28, reflection integrality forces t even. Reducing the
    # divided norm equation modulo2 forces Tx^2=0, hence x even.
    dyadic=[]
    for coords in product((0,1),repeat=3):
        x=tuple(map(F,coords));value=tuple(int(t)%2 for t in f.mul(T,sq(x)))
        assert (not any(value))==(not any(coords))
        dyadic.append({'x':list(coords),'T_x_squared':list(value)})
    result={'status':'PASS_EXACT_STRUCTURE_AND_BINARY_NO_ROOT_CERTIFICATE',
        'new_host':'<-2T>+B+A6+A1','determinant':f.encoded(det),'determinant_norm':str(norm),
        'comparison_scope':rows,'prime7_polynomial_residues':mod7,
        'prime7_reflection_norm_test':local7,
        'prime7_reflection_norm_argument':'Reflection integrality makes time, B and scalar coordinates divisible by7 and the A6 coordinates congruent to t*(1,2,3,4,5,6). Norm/7 is therefore -t^2 modulo7. The residues -1,-2,-4 are nonsquares inF343, excluding norm7/14/28 classes without beta at all heights.',
        'A5_determinant_ratio_to_old_2Tbeta':f.encoded(a5_ratio),'A5_square_witness':f.encoded(f.divide(T,D)),
        'D5_determinant_ratio_to_old_T':f.encoded(d5_ratio),'D5_square_witness':f.encoded(f.divide(f.TWO,B)),
        'triage_scope':'A5 andD5 determinant tests do not establish a new quadratic space. These identities alone assert neither field isometry nor equality of full reflection groups.',
        'unit_signatures':[list(x)for x in signatures],'signature_group_size':len(span),
        'dual_gram':[[f.encoded(x)for x in row]for row in inv],
        'exact_dual_exponent':f.encoded(exponent),'primitive_norm_divisor_bound':f.encoded(f.scale(exponent,2)),
        'source':str(source.relative_to(ROOT)),'end_weight_direction':f.vecencode(direction),
        'end_weight_norm':f.encoded(c.norm(direction)),'saved_escape':f.vecencode(v),'escape_norm':f.encoded(c.norm(v)),
        'saturated_slice_coordinate_minor':[[1,0],[0,1]],'slice_primitive_norm_pool':[2,4,14,28],
        'mod_beta_negative_squares':odd_norm_table,'dyadic_nonzero_square_table':dyadic,
        'binary_no_root_scope':'No primitive integral reflective normal lies in span(time,7omega6(A6)), at any height. Beta norms first force all supported coordinates divisible by beta; q2/q14 fail modulo beta, while q4/q28 force both slice coordinates even. The A6 reversal has the same conclusion. No full-host nonreflectivity or obstruction to cuts outside this binary span is claimed.'}
    if verify:assert json.loads(OUTPUT.read_text())==result
    else:OUTPUT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':result['status'],'compared_hosts':len(rows),'binary_slice':'NO_PRIMITIVE_ROOTS'}));return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();run(a.verify)
