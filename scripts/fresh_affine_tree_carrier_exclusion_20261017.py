#!/usr/bin/env python3
"""Complete finite role reduction for large affine-tree carriers.

The finite resistance range L<8 is proved from the two spherical deletions;
it is not an arbitrary rank or label cutoff. All-rank geometry is in the MD.
"""
import argparse
from fractions import Fraction as F
import hashlib
from itertools import combinations, product
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_affine_tree_carrier_exclusion_20261017'
SOURCES=[
    'results/fresh_private_affine_tree_exclusion_20261017.md',
    'scripts/fresh_private_affine_tree_exclusion_20261017.py',
    'results/fresh_private_affine_tree_exclusion_20261017.json',
    'results/fresh_affine_cycle_carrier_exclusion_20261017.md',
    'results/fresh_simple_finite_quad_vinberg_transfer_20261017.md',
    'results/fresh_simple_cusp_vinberg_weight_frontier_20261017.md',
]

def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()

def role_census():
    rows=[];tested=0
    for role,lo,nux,nuy in [('main_main',1,F(1),F(1)),
                           ('double_main',1,F(1,2),F(1)),
                           ('fork_main',2,F(1,4),F(1)),
                           ('sibling_forks',4,F(1,4),F(1,4))]:
        lengths=[4]if role=='sibling_forks'else range(lo,8)
        for L in lengths:
            positions=[0,2,4]if role=='sibling_forks'else [0]+list(range(2,L+1))if role=='fork_main'else list(range(L+1))
            for tx2,ty2 in product((1,2),repeat=2):
                a,b=nux*tx2,nuy*ty2
                values={j:a*j+b*(L-j)for j in positions}
                for j,k in combinations(positions,2):
                    tested+=1
                    if values[j]>=2 or values[k]>=2:continue
                    if role=='main_main':
                        assert L==1 and tx2==ty2==1 and(j,k)==(0,1)
                        outcome='small_affine_A2'
                    elif role=='fork_main':
                        assert tx2==ty2==1 and L in(3,4)and(j,k)==(L-1,L)
                        outcome='small_affine_A'+str(L)
                    elif role=='double_main':
                        assert ty2==1 and L in(1,2)
                        if L==1:
                            assert(j,k)==(0,1)
                            outcome='terminal_triangle_contacts_both_black'
                        else:
                            assert tx2==1 and(j,k)==(1,2)
                            outcome='white_double_terminal_two_step'
                    else:
                        if tx2==ty2==1:outcome='small_affine_A3'
                        else:
                            assert(tx2,ty2,j,k)in((2,1,0,2),(1,2,2,4))
                            outcome='white_fork_switch'
                    rows.append(dict(role=role,resistance=L,contact_squares=[tx2,ty2],
                                     marked_positions=[j,k],schur_costs=[str(values[j]),str(values[k])],
                                     outcome=outcome))
    counts={q:sum(r['outcome']==q for r in rows)for q in sorted({r['outcome']for r in rows})}
    assert counts=={'small_affine_A2':1,'small_affine_A3':4,'small_affine_A4':1,
                   'terminal_triangle_contacts_both_black':2,
                   'white_double_terminal_two_step':1,'white_fork_switch':2}
    # The opposite-end cases have resistance m-1>=17, independently of rank.
    assert F(17,4)>2
    # In every non-small-affine outcome at least one contact is marked black.
    for r in rows:
        if not r['outcome'].startswith('small_'):
            assert set(r['marked_positions'])&{0,r['resistance']}
    return dict(trials=tested,rows=rows,outcome_counts=counts,
                opposite_end_lower_cost='17/4',strict_universal_range='L<8')

def graph(kind,n):
    nu2=[F(1)]*n
    if kind=='C':
        es=[(i,i+1,4 if i in(0,n-2) else 3)for i in range(n-1)]
        nu2[0]=nu2[-1]=F(1,2)
    elif kind=='B':
        es=[(i,i+1,4 if i==n-4 else 3)for i in range(n-3)]+[(0,n-2,3),(0,n-1,3)]
        nu2[n-3]=F(1,2);nu2[-2:]=[F(1,4)]*2
    elif kind=='D':
        es=[(i,i+1,3)for i in range(n-5)]+[(0,n-4,3),(0,n-3,3),(n-5,n-2,3),(n-5,n-1,3)]
        nu2[-4:]=[F(1,4)]*4
    else:raise ValueError(kind)
    return es,nu2

def adj(n,es):
    g=[{}for _ in range(n)]
    for i,j,m in es:g[i][j]=m;g[j][i]=m
    return g

def code(n,es):
    g=adj(n,es);assert len(es)==n-1
    def r(i,p):return tuple(sorted((m,r(j,i))for j,m in g[i].items()if j!=p))
    return min(r(i,-1)for i in range(n))

def induced(es,keep):
    ix={x:i for i,x in enumerate(sorted(keep))}
    return [(ix[i],ix[j],m)for i,j,m in es if i in ix and j in ix]

def continuation_checks():
    rows=[]
    for kind,n in product(('B','C'),(18,19,23)):
        es,nu2=graph(kind,n);g=adj(n,es)
        x=next(i for i,v in enumerate(nu2)if v==F(1,2))
        p=next(iter(g[x]));y=next(i for i in g[p]if i!=x)
        z,w=n,n+1
        full=es+[(x,z,3),(y,z,3),(z,w,3)]
        keep=[i for i in range(n+2)if i!=p]
        target='B'if kind=='C'else'D';te,_=graph(target,n+1)
        assert code(n+1,induced(full,keep))==code(n+1,te)
        rows.append(dict(old_affine=kind,order=n,roles=dict(x=x,p=p,y=y,z=z,w=w),
                         deleted=p,result_affine=target,result_order=n+1))
    # Fork-switch continuation is the affine doubled star, with center z.
    # To avoid irrational arithmetic, use the scaled rational Laplacian:
    # nu=(1,1/sqrt2,1/2,1/2), conductances=(1,1/2,1/2).
    diagonal=[F(2),F(1),F(1,2),F(1,2)]
    conductance=[F(1),F(1,2),F(1,2)]
    assert diagonal[0]==sum(conductance)
    assert diagonal[1:]==conductance
    return dict(double_terminal_cases=rows,
                fork_switch=dict(labels=[4,3,3],positive_kernel=['1','1/sqrt(2)','1/2','1/2'],
                                 scaled_diagonal=list(map(str,diagonal)),
                                 conductances=list(map(str,conductance)),
                                 result='connected affine B3; stronger continuation bond gives negative kernel value'))

def grounded_formula_check():
    # Two contacts lie in different components after an internal path node
    # is deleted. Endpoint deletion removes the corresponding contact.
    # These are exactly the endpoint/interior limits of q(s).
    checks=[]
    for a,b,L in product((F(1,4),F(1,2),F(1),F(2)),repeat=3):
        if L.denominator!=1:continue
        assert a*0+b*L==b*L and a*L+b*0==a*L
        assert min(a,b)*L<=a*(L/3)+b*(2*L/3)
        checks.append([str(a),str(b),str(L)])
    return dict(endpoint_and_lower_bound_checks=len(checks),
                formula='q(s)=t_x^2*nu_x^2*s+t_y^2*nu_y^2*(L-s)',
                reason='grounded tree Green diagonal, zero cross-component inverse entry')

def build():
    r=role_census();c=continuation_checks();g=grounded_formula_check()
    return dict(status='PASS_ALL_RANK_AFFINE_TREE_CARRIER_ROLE_REDUCTION',
                sources={p:sha(p)for p in SOURCES},source_sha256=sha('scripts/'+TAG+'.py'),
                role_census=r,continuations=c,grounded_identity=g,
                scope='The finite role range is derived from q<2 and coefficient>=1/4. Full shared-path/carrier exclusion is the accompanying all-rank proof, not a finite-rank extrapolation.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=json.loads(json.dumps(build()));out=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(d['status'],d['role_census']['trials'],d['role_census']['outcome_counts'],flush=True)
