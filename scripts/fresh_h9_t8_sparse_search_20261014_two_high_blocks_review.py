#!/usr/bin/env python3
"""Separate symbolic/literal-incidence review of the two-high-block theorem."""
import argparse,json
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
BASE='fresh_h9_t8_sparse_search_20261014_two_high_blocks'
TAG=BASE+'_review'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def scalar(a):return sum(s.Rational(v)*r for v,r in zip(a,[1,s.sqrt(2),s.sqrt(3),s.sqrt(6)]))
def run():
    primary=f'results/{BASE}.json';d=json.loads((ROOT/primary).read_text())
    for p,h in d['sources'].items():assert sha(p)==h
    geom=json.loads((ROOT/GEOM).read_text())['representatives'][8]['geometry']
    vertices=[set(v['incident']) for v in geom['vertices']]
    missing={frozenset(x) for x in geom['minimal_nonfaces']}
    S=[1,2,4,5,7,8,12,13];blocks=[[1,2],[4,5],[7,8],[12,13]]
    assert d['positive_support']==S and set(S)<=vertices[d['positive_support_actual_vertex']]
    assert d['fixed_labels']==[]
    assert d['finite_label_lower_bounds']==[[1,2,7],[4,5,7],[7,8,4],[12,13,4]]
    core=sorted(S+[11,14]);pairs=set()
    for row in d['all45_core_meeting_pairs']:
        p=tuple(row['pair']); assert set(p)<=vertices[row['actual_vertex']];pairs.add(p)
    assert pairs==set(combinations(core,2))
    expected={}
    for i,j in combinations(range(4),2):
        for a,b in product(blocks[i],blocks[j]):
            if i<2:
                support=[next(x for x in blocks[i] if x!=a),a,b]
                cos=[s.sqrt(3)/2,s.Rational(1,2)];w=s.Matrix([s.sqrt(3)/2,1,s.Rational(1,2)])
            else:
                support=[next(x for x in blocks[i] if x!=a),a,b,next(x for x in blocks[j] if x!=b)]
                cos=[s.sqrt(2)/2,s.Rational(1,2),s.sqrt(2)/2];w=s.Matrix([1,s.sqrt(2),s.sqrt(2),1])
            M=s.eye(len(support))
            for k,c in enumerate(cos):M[k,k+1]=M[k+1,k]=-c
            assert all(s.simplify(x)==0 for x in M*w) and all(x.is_positive for x in w)
            expected[a,b]=(support,M,w)
    assert len(expected)==24
    for row in d['cross_right_guards']:
        support,M,w=expected.pop(tuple(row['pair']))
        assert support==row['support'] and set(support)<=vertices[row['actual_vertex']]
        assert s.Matrix([[scalar(x) for x in rr] for rr in row['gram']])==M
        assert s.Matrix([scalar(x) for x in row['positive_null_vector']])==w
    assert not expected
    needed={frozenset(x) for x in [[1,2,14],[4,5,14],[11,12,13]]}
    assert needed<=missing and needed=={frozenset(x) for x in d['missing_triples']}
    for row in d['anchor_high_right_budgets']:
        own=row['active_own_root'];assert own in [12,13]
        support=(set(S)-({12,13}-{own}))|{11}
        assert support==set(row['support']) and support<=vertices[row['actual_vertex']]
    assert {x['active_own_root'] for x in d['anchor_high_right_budgets']}=={12,13}
    seen=set()
    for row in d['third_other_blocks_right_budgets']:
        a,b=row['active_high_roots'];assert a in [1,2] and b in [4,5];seen.add((a,b))
        support=({a,b,14}|set(blocks[2])|set(blocks[3]))
        assert support==set(row['support']) and support<=vertices[row['actual_vertex']]
    assert seen==set(product([1,2],[4,5]))
    stars=[]
    for row in d['third_anchor_upper_stars']:
        a,b=row['support'][2:]; assert row['support']==[11,14,a,b]
        assert (a,b) in seen and set(row['support'])<=vertices[row['actual_vertex']]
        M=s.eye(4);w=s.Matrix([s.sqrt(2)/2,1,s.Rational(1,2),s.Rational(1,2)])
        for i,c in [(0,s.sqrt(2)/2),(2,s.Rational(1,2)),(3,s.Rational(1,2))]:M[1,i]=M[i,1]=-c
        assert all(s.simplify(x)==0 for x in M*w)
        assert M==s.Matrix([[scalar(x) for x in rr] for rr in row['gram']])
        assert w==s.Matrix([scalar(x) for x in row['positive_null_vector']])
        stars.append(row['support'])
    assert len(stars)==4 and {tuple(x[2:]) for x in stars}==seen
    # Derive U6 independently from the trigonometric multiple-angle identity.
    x=s.symbols('x');U=s.chebyshevu(6,x)
    assert s.Poly(U,x).all_coeffs()==[64,0,-80,0,24,0,-1]
    t=s.symbols('t');f=64*t**3-80*t**2+24*t-1
    assert f.subs(t,s.Rational(4,5))==-s.Rational(29,125)
    assert s.diff(f,t).subs(t,s.Rational(3,4))==12
    assert s.diff(f,t,2).subs(t,s.Rational(3,4))==128
    # f'>0 on[3/4,1]; c=cos(pi/7)>cos(pi/6), and f(c^2)=0,
    # so c^2>4/5. Monotonicity extends this to every integer m>=7.
    c=s.symbols('c',positive=True)
    Q=s.Matrix([[1,-c],[-c,1]]).inv()
    v=s.Matrix([s.Rational(1,2)]*2)
    assert s.simplify((v.T*Q*v)[0]-1/(2*(1-c)))==0
    low=[1+s.sqrt(2)/2,(3+s.sqrt(5))/4,1+s.sqrt(5)/5,s.Integer(2)]
    assert all((z-s.Rational(5,4)).is_positive for z in low)
    assert s.simplify(((1+s.sqrt(5))/4)**2/(1-s.Rational(1,2))-low[1])==0
    assert s.simplify(s.Rational(1,2)/(1-((1+s.sqrt(5))/4)**2)-low[2])==0
    assert (1-(5+s.sqrt(5))/10).is_positive
    assert s.Rational(1,4)/(1-s.Rational(3,4))==1
    assert s.Rational(1,4)+s.Rational(1,4)+s.Rational(1,2)==1
    assert s.Rational(3,2)*s.Rational(1,4)-s.Rational(1,2)**2==s.Rational(1,8)
    paths=[primary,f'scripts/{BASE}.py',GEOM,f'scripts/{TAG}.py']
    return dict(status='PASS_INDEPENDENT_TWO_HIGH_BLOCKS_DOMAIN_AND_INDEX_REVIEW',sources={p:sha(p) for p in paths},
        implementation='No producer/project imports. SymPy exact symbolic matrices, independent support construction and literal source-vertex checks.',
        exact_null_guards=28,cross_guards=24,actual_third_budgets=4,actual_anchor_budgets=2,meeting_pairs=45,
        missing_triples=[sorted(x) for x in sorted(needed,key=lambda x:sorted(x))],
        residual_stronger_bounds=dict(third_square='less than -3/2',anchor_square='less than -1/4',off_diagonal_square='at most1/4',determinant='greater than1/8'),
        scope='Four internal finite-label lower bounds only. No outside contact labels, fixed third contacts, finite label cutoff, polynomial bank, or diagram automorphism assumed. Geometry realizability inherited from its frozen certificate.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'])
