#!/usr/bin/env python3
"""Exact algebra for nonsimple-cusp factor payments and spectator corrections.

All-rank geometric arguments are in the report; finite checks are corroborative.
No search, optimizer, or claim of an unrestricted dimension improvement.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

import sympy as s

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_cusp_factor_payment_framework_20261017'
INPUTS = ['results/fresh_simple_cusp_affine_exclusion_20261017.json',
          'results/fresh_simple_cusp_affine_exclusion_20261017.md',
          'results/fresh_simple_cusp_distance_bound_20261017.md',
          'results/fresh_finite_volume_simplex_truncation_20261017.md',
          'results/fresh_cusp_finite_port_triangle_capacity_20261017.md',
          'results/fresh_simple_cusp_affine_a_quadrilateral_20261017.json',
          'results/fresh_simple_cusp_affine_tree_completion_20261017.md']


def beta(n):
    return F(n * (n - 6), 8) if n % 2 == 0 else F((n - 1) * (n - 5), 8)


def run():
    old = json.loads((ROOT / INPUTS[0]).read_text())
    for p, h in old['sources'].items():
        assert sha256((ROOT / p).read_bytes()).hexdigest() == h
    assert old['minimal_indefinite_core_order_bound'] == 5
    graphs = old['inherited_classical_distance7_checks']
    assert len(graphs) == 245
    assert all(row['cost'] <= 7 * row['order'] for row in graphs)

    n, t, count = s.symbols('n t count', integer=True)
    N = 56 * count + 6
    even = s.expand((N + 2*t) * (N + 2*t - 6) - 56*count*(N + 2*t))
    odd = s.expand((N + 1 + 2*t - 1)*(N + 1 + 2*t - 5) - 56*count*(N + 1 + 2*t))
    assert s.expand(even - (4*t*t + 2*N*t)) == 0
    assert s.expand(odd - (4*t*t + (2*N+4)*t + N+6)) == 0
    thresholds = []
    for k in range(1, 9):
        bound = 56*k + 6
        assert beta(bound) == 7*k*bound
        assert beta(bound+1) > 7*k*(bound+1)
        thresholds.append(dict(max_factors=k, bound=bound-1, excluded_even_endpoint=bound,
            endpoint_finite_vertex_slack=str(beta(bound)-7*bound)))
    incidence_checks = []
    for dim in range(10, 181):
        eta = F(4*(dim-1), dim-2) if dim % 2 == 0 else F(4*dim, dim-1)
        alpha = 5/eta-1
        assert alpha * F(dim*(dim-1), 2) == beta(dim)
        assert alpha * F(dim-1, 2) == beta(dim)/dim
        incidence_checks.append([dim, str(alpha), str(beta(dim)/(7*dim))])

    # Four-root determinant, with a PSD Schur correction only on c,d.
    a,A,B,C,D,E,h,k = s.symbols('a A B C D E h k')
    p,q,r = s.symbols('p q r')
    M = s.Matrix([[a,-a,B-h,-A],[-a,a,-B,A-k],
                  [B-h,-B,C,-E],[-A,A-k,-E,D]])
    det = s.expand(M.det())
    desired = (h*k-k*B-h*A)**2-a*(C*k*k+D*h*h+2*E*h*k)
    assert s.expand(det-desired) == 0
    corrected = M.copy()
    corrected[2,2] -= p
    corrected[3,3] -= r
    corrected[2,3] -= q
    corrected[3,2] -= q
    change = s.expand(corrected.det()-det)
    assert s.expand(change - a*(p*k*k-2*q*h*k+r*h*h)) == 0

    # A~1 x A~_(m-1): rational local cusp models of all ranks, not global
    # finite-volume chambers. The positive kernel sums share the same null nu.
    models = []
    for m in range(3, 25):
        Ar = s.Matrix(m-1, m-1, lambda i,j: 2 if i==j else -1 if abs(i-j)==1 else 0)
        Q = s.diag(-s.Integer(0), -s.Integer(0), s.Integer(1), Ar)
        Q[0,1] = Q[1,0] = -1
        roots = []
        e = s.zeros(m+2,1); e[2]=1
        nu = s.zeros(m+2,1); nu[0]=1
        roots += [e,nu-e]
        for i in range(m-1):
            v=s.zeros(m+2,1);v[3+i]=1;roots.append(v)
        roots.append(nu-sum(roots[2:],s.zeros(m+2,1)))
        R=s.Matrix.hstack(*roots)
        G=R.T*Q*R
        assert sum(roots[:2],s.zeros(m+2,1)) == nu
        assert sum(roots[2:],s.zeros(m+2,1)) == nu
        assert G.rank() == m and R.rank() == m+1
        assert all(G[i,j] == 0 for i in range(2) for j in range(2,m+2))
        H=s.zeros(m+2,1);H[0]=m;H[1]=1;H[2]=-s.Rational(1,2)
        for i in range(1,m):H[2+i]=-s.Rational(i*(m-i),2*m)
        pairings=list(R.T*Q*H)
        assert pairings[:2] == [-s.Rational(1,2)]*2
        assert pairings[2:] == [-s.Rational(1,m)]*m
        norm=(H.T*Q*H)[0]
        expected=-2*m+s.Rational(1,4)+s.Rational(m*m-1,12*m)
        assert norm == expected and norm < 0
        models.append(dict(ambient_dimension=m+1, factor_ranks=[1,m-1],
            gram_rank=m, normal_span_rank=m+1, strict_center_norm=str(norm),
            strict_center_pairings=list(map(str,pairings))))
    paths=[Path(__file__)]+[ROOT/p for p in INPUTS]
    return dict(status='PASS_EXACT_CUSP_FACTOR_TRADEOFF_AND_PSD_SPECTATOR_IDENTITY',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        hypothesis='Every irreducible affine cusp factor has rank at least9.',
        weighted_inequality='beta_n/(7*n) < (n*v_f+sum(s_v*V_v))/(n*v_f+sum(V_v))',
        two_factor_dimension_bound=117, strictness_source="Khovanskii Theorems8-10: local cut ratios are bounded by eta, while a generic global minimum has ratio2<eta.",
        classical_graph_checks=len(graphs), thresholds=thresholds,
        threshold_polynomials_times8=dict(even=str(even),odd=str(odd)),
        incidence_identity_checks=incidence_checks,
        spectator_determinant_change=str(change),
        spectator_normalized_change='alpha*(p*u^2-2*q*u*v+r*v^2)>=0 for [[p,q],[q,r]] PSD',
        local_cusp_models=models,
        local_model_scope='All-rank rational two-factor cusp neighborhoods with a rank1 factor; no claim of finite-volume completion in arbitrary dimensions.',
        scope='Restricted geometric theorem and exact determinant lemma; no unrestricted163 improvement, full3sV face-payment claim, novelty claim, or optimizer.')


if __name__ == '__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();path=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert result==json.loads(path.read_text())
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],len(result['local_cusp_models']),'local models;',
          result['classical_graph_checks'],'inherited graph rows;',
          'one exact PSD spectator identity',flush=True)
