#!/usr/bin/env python3
"""Independent SymPy reconstruction of the universal ideal-corner scalars.

No producer or other project arithmetic code is imported. All profiles and
the formal Schur matrix are rebuilt literally; geometry is reviewed in MD.
"""
import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import sympy as S

ROOT=Path(__file__).resolve().parents[1]
PRIMARY='fresh_distance8_ideal_corner_replacement_20261017'
TAG=PRIMARY+'_review'
q=S.sqrt(2)

def wire(x):
    x=S.expand(S.radsimp(S.sympify(x)));a=x.coeff(q,0);b=x.coeff(q)
    assert a.is_Rational and b.is_Rational and S.expand(x-a-b*q)==0
    return [str(a),str(b)]

def nonnegative(x):
    a,b=map(S.Rational,wire(x))
    if b==0:return a>=0
    if a==0:return b>=0
    if a>0:return b>=0 or a*a>=2*b*b
    return b>0 and 2*b*b>=a*a

def arithmetic():
    saved=json.loads((ROOT/f'results/{PRIMARY}.json').read_text())
    ps=[dict(name='isolated',c0=S.Integer(2),moment=[S.Integer(0)]*2,position=[0,0],minimum_height=S.Integer(1))]
    for name,moment,c0,pos,height in(('end3',1,1,1,1),('end4',q,0,1,q),('next3',2,0,2,1)):
        for j in range(2):
            ps.append(dict(name=f'{name}_{j}',c0=S.Integer(c0),moment=[moment if i==j else S.Integer(0)for i in range(2)],position=[pos if i==j else 0 for i in range(2)],minimum_height=S.sympify(height)))
    ps.append(dict(name='both_ends3',c0=S.Integer(0),moment=[S.Integer(1)]*2,position=[1,1],minimum_height=S.Integer(2)))
    profiles=[]
    for p in ps:
        profiles.append(dict(name=p['name'],c0=wire(p['c0']),moment=[wire(x)for x in p['moment']],position=p['position'],minimum_height=wire(p['minimum_height']),umax=wire(1/p['minimum_height'])))
    assert profiles==saved['cycle_profiles']
    rows=[]
    for c,d in product(ps,repeat=2):
        bounds=[max(9,7+c['position'][j]+d['position'][j])if c['position'][j]*d['position'][j]else 9 for j in range(2)]
        a,b=map(lambda k:S.Rational(1,k),bounds)
        for ui,vi in product((0,1),repeat=2):
            u=ui/c['minimum_height'];v=vi/d['minimum_height']
            L=[u*c['moment'][j]+v*d['moment'][j]for j in range(2)]
            value=S.expand(1-2*(a*L[0]+b*L[1])-(a+b)*(c['c0']*u*u+d['c0']*v*v)-a*b*(L[0]-L[1])**2)
            assert nonnegative(value-S.Rational(5,81))
            rows.append(dict(c=c['name'],d=d['name'],arc_lower_bounds=bounds,upper_height_reciprocals=[ui,vi],normalized_determinant=wire(value)))
    assert rows==saved['cycle_rectangle_checks']and len(rows)==256
    assert any(r['normalized_determinant']==['5/81','0']for r in rows)
    tp=[('isolated',0,2,1),('outer_end3',0,1,1),('outer_end4',0,0,q),('outer_next3',0,0,1),('inner_end3',1,1,1),('inner_next3',2,0,1),('both_ends3',1,0,2),('inner_end4',q,0,q)]
    tree=[];costs={}
    for name,l,c,h in tp:
        l,c,h=map(S.sympify,(l,c,h));g=S.expand(2*l/h+c/h**2)
        assert nonnegative(g)and nonnegative(4-g);costs[name]=g
        tree.append(dict(name=name,moment=wire(l),c0=wire(c),minimum_height=wire(h),maximum_single_root_cost=wire(g)))
    assert tree==saved['tree_profiles']
    pairs=[]
    for c,d in product(costs,repeat=2):
        gap=1-(costs[c]+costs[d])/9;assert nonnegative(gap-S.Rational(1,9))
        pairs.append(dict(c=c,d=d,P_lower_bound=9,gap=wire(gap)))
    assert pairs==saved['tree_pair_checks']and len(pairs)==64
    assert any(r['gap']==['1/9','0']for r in pairs)
    # A genuinely independent symbolic determinant calculation.
    al,A,B,C,D,E,h,k=S.symbols('alpha A B C D E h k')
    G=S.Matrix([[al,-al,B-h,-A],[-al,al,-B,A-k],[B-h,-B,C,-E],[-A,A-k,-E,D]])
    expected=(h*k-k*B-h*A)**2-al*(C*k*k+D*h*h+2*E*h*k)
    determinant=S.Poly(S.expand(G.det(method='domain-ge')),al,A,B,C,D,E,h,k)
    assert S.expand(determinant.as_expr()-expected)==0
    terms=[dict(exponents=list(e),coefficient=str(c))for e,c in sorted(determinant.terms())]
    assert terms==saved['formal_schur_determinant']['terms']
    # Spectator correction, derived directly rather than assumed by name.
    u,v,p,t,r=S.symbols('u v p t r')
    corrected=expected.subs({C:C-p,D:D-r,E:E+t},simultaneous=True)
    correction=S.expand((corrected-expected).subs({h:1/u,k:1/v})*u*u*v*v)
    assert S.expand(correction-al*(p*u*u-2*t*u*v+r*v*v))==0
    # Verify endpoint and next-node identities symbolically for arbitrary P.
    P=S.symbols('P',positive=True)
    endpoint=(P-1)/P;second=2*(P-2)/P
    assert all(S.cancel(z)==0 for z in(2-endpoint-1-1/P,2-2*endpoint-2/P,2-second-4/P,2-1-endpoint-1/P))
    return dict(cycle_profiles=profiles,cycle_rectangle_checks=rows,tree_profiles=tree,
        tree_pair_checks=pairs,formal_schur_terms=terms,
        symbolic_spectator_correction='alpha*(p*u^2-2*t*u*v+r*v^2)',
        symbolic_path_inverse_identities=4,cycle_uniform_gap=['5/81','0'],tree_uniform_gap=['1/9','0'])

def run():
    data=arithmetic()
    paths=[Path(__file__)]+[ROOT/f'{folder}/{PRIMARY}.{ext}'for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    for stem in('fresh_cusp_factor_payment_framework_20261017','fresh_affine_distance8_feature_payment_20261017'):
        paths += [ROOT/f'results/{stem}.md',ROOT/f'results/{stem}.json']
    for p in paths:assert p.is_file(),p
    return dict(status='PASS_INDEPENDENT_SYMPY_IDEAL_DISTANCE9_CORNER_REPLACEMENT_REVIEW',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        counts=dict(cycle_rectangle_values=256,tree_pair_values=64,tree_single_root_costs=8,
            symbolic_path_inverse_identities=4,formal_determinant_identities=1,spectator_identities=1),
        reconstruction=data,scope='Independent exact scalar reconstruction and written mathematical/source review. The finite/affine classification and actual cusp-factor geometry remain explicit inherited inputs. No new ambient bound or formalization.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--arithmetic-only',action='store_true');args=ap.parse_args()
    if args.arithmetic_only:
        d=arithmetic();print('PASS_INDEPENDENT_IDEAL_REPLACEMENT_ARITHMETIC',len(d['cycle_rectangle_checks']),len(d['tree_pair_checks']));raise SystemExit(0)
    d=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,sort_keys=True,indent=2)+'\n')
    print(d['status'],d['counts'],flush=True)
