#!/usr/bin/env python3
"""Exact sign ledger distinguishing Lorentz-inertia and dimension-nine rank gates."""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_lorentz_inertia_20261015'

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def sign5(a,b):
    if not a:return(b>0)-(b<0)
    if not b or(a>0)==(b>0):return(a>0)-(a<0)
    z=a*a-5*b*b;return((a>0)-(a<0))*((z>0)-(z<0))
def sign4(v):
    a,b,c,d=map(F,v);x,y=sign5(a,c),sign5(b,d)
    if not x:return y
    if not y or x==y:return x
    return x*sign5(a*a+5*c*c-2*b*b-10*d*d,2*a*c-4*b*d)
def anp(v):
    p,q,r,s=[F(0)]*(4-len(v))+list(map(F,v))
    # (sqrt2+sqrt5)^2=7+2sqrt10; its cube=17sqrt2+11sqrt5.
    return(7*q+s,17*p+r,11*p+r,2*q)

def run():
    sources={};seen={}
    def read(tag,review=True):
        p=f'results/{tag}.json';d=json.loads((ROOT/p).read_text());sources[p]=sha(p)
        assert str(d['status']).startswith('PASS')
        for q,h in d.get('sources',{}).items():assert sha(q)==h,q
        if review:
            rp=f'results/{tag}_review.json';r=json.loads((ROOT/rp).read_text());sources[rp]=sha(rp)
            assert str(r['status']).startswith('PASS')
            assert r['sources'].get(p)==sha(p)
            for q,h in r.get('sources',{}).items():assert sha(q)==h,q
        seen[tag]=d;return d
    ledger=[]
    d=read('fresh_h9_star_zero_20261014_own_three')
    assert d['field']=='Q(sqrt2,sqrt5); ANP descending rational coefficients at primitive generator sqrt2+sqrt5'
    signs=Counter();anchors=0
    for block in d['records']:
        assert all(x in(4,5)for x in block['block_labels'])
        assert all(sign4(anp(x))<0 for x in block['anchor_determinants']);anchors+=len(block['anchor_determinants'])
        for row in block['cases']:signs[sign4(anp(row['determinant']))]+=1
    assert signs=={1:20824,-1:1128}and anchors==784 and sum(signs.values())==21952
    ledger.append(dict(source='fresh_h9_star_zero_20261014_own_three',positive_regions=20824,negative_regions=1128,
        negative_anchor_minors=784,negative_diagonal_reason='own pi3/pi3 cost=1/[2(1-c)]>1 for c>=cos(pi/4)',
        universal_inertia=20824,rank_only=1128,negative_region_signature=[10,1],
        scope='The sixteen finite internal4/5 banks only; their separate high-block reductions are not recounted.'))
    d=read('fresh_h9_star_shared_high_20261014');count=Counter()
    for r in d['cases']:
        count[r['status']]+=1
        if r['status']=='STRICT_POSITIVE_FOR_ALL_D_GREATER_THAN_FOUR':
            assert sign4(anp(r['slope']))>0 and sign4(anp(r['boundary_D4']))>0
    assert sum(count.values())==7840 and d['positive_cases']==7392 and d['actual_vertex_rejections']==448
    assert count['STRICT_POSITIVE_FOR_ALL_D_GREATER_THAN_FOUR']==7392
    ledger.append(dict(source='fresh_h9_star_shared_high_20261014',universal_inertia=7392,
        universal_actual_face_failures=448,rank_only=0,negative_diagonal_reason='1-D/4-F<0, D>4,F>=0',
        scope='The frozen first-three internal4/5 coefficient cases, before separately proved high-block reductions.'))
    d=read('fresh_h9_star_opposite_high_20261014');count=Counter()
    for r in d['cases']:
        count[r['status']]+=1
        if r['status']=='STRICT_INTERIOR_BERNSTEIN_6_7_TO_ONE':
            v=[sign4(anp(x))for x in r['bernstein_6_7']];assert min(v)>=0 and max(v)>0
        elif r['status']=='STRICT_BERNSTEIN_9_10_TO_ONE':
            v=[sign4(anp(x))for x in r['bernstein_9_10']];assert min(v)>0
    assert count=={'STRICT_INTERIOR_BERNSTEIN_6_7_TO_ONE':43412,'NONPOSITIVE_ACTUAL_VERTEX_AT_RELAXED_LAST_INTERNAL_FOUR':3612,'STRICT_BERNSTEIN_9_10_TO_ONE':16}
    ledger.append(dict(source='fresh_h9_star_opposite_high_20261014',universal_inertia=43428,
        universal_actual_face_failures=3612,rank_only=0,negative_diagonal_reason='1-D/4-F<0, D>4,F>=0',
        scope='The frozen first-three internal4/5 coefficient cases, before separately proved high-block reductions.'))
    d=read('fresh_h9_star_third34_20261014');n=0
    for row in d['continuous_cases']:
        for term in row['shifted_terms']:
            a,b=map(F,term['coefficient']);assert sign4((a,b,0,0))>0;n+=1
    assert len(d['continuous_cases'])==6 and n==66
    ledger.append(dict(source='fresh_h9_star_third34_20261014',universal_negative_planes=6,
        positive_coefficients=66,rank_only=0,negative_diagonal_reason='1-X<0 and 2-2mu-Y_last-F-F_middle<0 with X>1,Y_last>7/2'))
    d=read('fresh_h9_star_single_third_20261014');n=0
    for row in d['continuous_cases']:
        for term in row['shifted_terms']:assert sign4(term['coefficient'])>0;n+=1
    assert len(d['continuous_cases'])==11 and n==78
    ledger.append(dict(source='fresh_h9_star_single_third_20261014',universal_inertia=11,positive_coefficients=78,
        rank_only=0,negative_diagonal_reason='1-D/4<0 with D>100/19'))
    d=read('fresh_h9_star_single_third_opposite_20261015');count=Counter();co=Counter()
    for row in d['continuous_cases']:
        count[row['target_sign']]+=1
        for x in row['coefficients']:assert sign4(x['coefficient'])==row['target_sign'];co[row['target_sign']]+=1
    assert count=={1:15,-1:1}and co=={1:165,-1:3}
    ledger.append(dict(source='fresh_h9_star_single_third_opposite_20261015',universal_inertia=15,
        rank_only=1,negative_diagonal_reason='1-D/4<0 with D>100/19',negative_region_signature=[10,1],
        negative_exception='Owner4,third5,last7; foreign at third-right owner root; only the other third-anchor contact pi3.'))
    # The fresh primary and its separate Fraction review are now frozen.
    d=read('fresh_h9_star_third_four_20261015');count=Counter();co=Counter()
    for row in d['foreign_cases']:
        count[row['target_sign']]+=1
        for x in row['bernstein_coefficients']:assert sign4(x)==row['target_sign'];co[row['target_sign']]+=1
    assert count=={1:8,-1:1}and co=={1:24,-1:3}
    ledger.append(dict(source='fresh_h9_star_third_four_20261015',universal_no_foreign='All patterns and every owner>=5',
        universal_foreign_inertia=8,rank_only=1,negative_diagonal_reason='1-D/4<0 with D>100/19',
        negative_region_signature=[10,1],independent_review='generation and saved replay PASS',
        negative_exception='Owner5,third4,last7; opposite anchors, foreign at third-right root; only other third-anchor contact pi3.'))
    # In both negative exceptions F lies in[1/2,1), mu=0. The anchor minor is
    # (1-F)+D*(F/4-7/16) <= 1/2-3D/16 < 0.
    bound=F(1,2)-F(3,16)*F(100,19);assert bound==F(-37,76)<0
    # F=1/2 for owner4; F=(5+sqrt5)/10 for owner5.
    # This simplification follows from 1/[4(1-cos^2(pi/5))].
    F5=(F(1,2),0,F(1,10),0)
    assert sign4((F5[0]-F(1,2),0,F5[2],0))>0
    assert sign4((1-F5[0],0,-F5[2],0))>0
    sources[f'scripts/{TAG}.py']=sha(f'scripts/{TAG}.py')
    return dict(status='PASS_LORENTZ_INERTIA_VERSUS_RANK_LEDGER',sources=sources,
        lemma=dict(positive_support_rank=8,schur_residual_size=3,
            positive_determinant_negative_diagonal='Exactly two negative residual eigenvalues; impossible in every Lorentz space.',
            negative_determinant_indefinite_two_minor='Exactly one negative and two positive residual eigenvalues; full Gram signature(10,1), requiring ambient hyperbolic dimension>=10.',
            two_by_two='Two negative diagonals and positive determinant give a negative-definite plane in every dimension.'),
        ledger=ledger,negative_exception_anchor_minor=dict(formula='1-F+D*(F/4-7/16)',F_interval=['1/2','1'],
            D_lower='100/19',strict_upper=str(bound)),
        arithmetic_scope='Separate Fraction signs of saved exact coefficients and determinant values; inherited producer/review bindings. No recensus of the original matrices or literal face guards.',
        global_scope='Conditional local forbidden configurations only. No proof that dimension or face averaging forces these configurations. No new upper dimension bound, compact polytope or finite-volume construction.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],len(d['ledger']))
