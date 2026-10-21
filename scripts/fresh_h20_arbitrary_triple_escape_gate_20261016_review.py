#!/usr/bin/env python3
"""Independent SymPy reconstruction of all arbitrary-chain triple identities."""
import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import sympy as S

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_arbitrary_triple_escape_gate_20261016_review'
PRIMARY='results/fresh_h20_arbitrary_triple_escape_gate_20261016.json'
BASE='results/fresh_h20_s61_multicap_20261014_exact.json'
UPGRADE_REVIEW='results/fresh_h20_asymmetric_arc_escape_gate_20261016_review.json'


def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):
    d=json.loads((ROOT/p).read_text())
    for q,h in d.get('sources',{}).items():assert sha(q)==h,(p,q)
    return d
def rat(v):return S.Rational(str(v))


def enclosure(poly,lo,hi):
    # Mean-value theorem on0<C<1; unlike producer's interval Horner evaluation.
    assert 0<lo<hi<1
    midpoint=(lo+hi)/2
    radius=(hi-lo)/2
    derivative_bound=sum(k*abs(poly.nth(k)) for k in range(1,poly.degree()+1)) if not poly.is_zero else 0
    value=poly.eval(midpoint)
    return value-radius*derivative_bound,value+radius*derivative_bound


def run():
    primary,base,upgrade=map(read,(PRIMARY,BASE,UPGRADE_REVIEW))
    assert upgrade['degree_over_K_Qab']==50
    assert upgrade['weak_escape_norm']=='190011/110000'
    x=S.symbols('x')
    gp={}
    for k in range(1,26):
        num=S.Poly(2*(S.chebyshevt(k,x)-x),x,domain=S.ZZ)
        quotient,remainder=S.div(num,S.Poly(1-x,x))
        assert remainder.is_zero
        gp[k]=quotient
    assert gp[1].is_zero and gp[2].as_expr()==-4*x-2
    alo,ahi=map(rat,base['parameter']['isolating_interval'])
    clo,chi=alo*alo/(alo*alo+2),ahi*ahi/(ahi*ahi+2)
    lo,hi=map(rat,primary['C_interval'])
    assert 0<lo<clo<chi<hi<1
    gap_records=[]
    for k in range(1,26):
        saved=primary['chain_pair_polynomials'][k-1]
        coeff=[int(gp[k].nth(i)) for i in range(max(0,gp[k].degree())+1)]
        assert coeff==saved['polynomial']
        if k==1:continue
        bounds=enclosure(gp[k],lo,hi)
        difference=enclosure(gp[k]-gp[k-1],lo,hi)
        assert bounds[1]<-2 and difference[1]<0
        gap_records.append(dict(gap=k,bounds=list(map(str,bounds)),difference_bounds=list(map(str,difference))))

    saved_cases={tuple(r['gaps']):r for r in primary['cases']}
    cases=[];counts=Counter();translations=0
    for q in range(24,0,-1):
        for p in range(25-q,0,-1):
            a,b,c=gp[p],gp[q],gp[p+q]
            # Reconstruct the adjugate by2x2 minors, then the determinant by
            # first-row expansion. Poly arithmetic is provided by SymPy.
            adj=[[4-b*b,b*c-2*a,a*b-2*c],
                 [b*c-2*a,4-c*c,a*c-2*b],
                 [a*b-2*c,a*c-2*b,4-a*a]]
            det=2*adj[0][0]+a*adj[0][1]+c*adj[0][2]
            columns=[det,-adj[0][0],-adj[1][1],-adj[2][2],
                     -2*adj[0][1],-2*adj[0][2],-2*adj[1][2]]
            degree=max(poly.degree() for poly in columns if not poly.is_zero)
            assert degree<=2*(p+q)-2<=48
            literal=[[int(poly.nth(i)) for i in range(max(0,poly.degree())+1)] for poly in columns]
            row=saved_cases[p,q]
            assert literal==row['coefficient_columns']
            M=S.Matrix([[poly.nth(i) for poly in columns] for i in range(degree+1)])
            reduced,pivots=M.rref()
            null=M.nullspace()
            assert len(null)+len(pivots)==7
            assert list(pivots)==row['pivot_columns']
            saved_null=[S.Matrix(list(map(rat,v))) for v in row['kernel_basis']]
            assert null==saved_null
            bounds=enclosure(det,lo,hi)
            assert bounds[1]<0
            # Verify the nonlinear square/product reduction, rather than just
            # recognizing the producer's literal kernel class.
            if not null:
                forced='all coefficient variables vanish'
                counts['zero_kernel']+=1
            else:
                square_pair=(2,3) if p==1 and q>1 else (1,2) if q==1 and p>1 else (1,3)
                constraint=S.zeros(7,1)
                for j in square_pair:constraint[j]=1
                # Orthogonality to the full kernel proves this sum of actual
                # nonnegative squares is zero. Both contacts vanish.
                assert all((constraint.T*v)[0]==0 for v in null)
                remaining=next(j for j in (1,2,3) if j not in square_pair)
                reduced_actual=M[:,[0,remaining]]
                if (p,q)==(1,1):
                    assert remaining==2 and reduced_actual.rank()==1
                    assert reduced_actual*S.Matrix([S.Rational(1,2),1])==S.zeros(M.rows,1)
                    forced='outer contacts zero; q_y=middle_contact_squared/2'
                    counts['consecutive_exception']+=1
                else:
                    assert reduced_actual.rank()==2
                    forced='two contacts forced zero by squares; all cross products vanish; last square and q_y forced zero'
                    counts['nonlinear_zero']+=1
            translations+=26-p-q
            cases.append(dict(gaps=[p,q],degree=degree,kernel_dimension=len(null),
                              determinant_bounds=list(map(str,bounds)),actual_contact_reduction=forced))
    assert len(cases)==300 and translations==2600
    assert counts==dict(zero_kernel=242,nonlinear_zero=57,consecutive_exception=1)
    # Independent2x2 inverse and monotonicity proof for any outer pair.
    g,f,h=S.symbols('g f h',real=True)
    pair=S.Matrix([[2,g],[g,2]])
    vec=S.Matrix([-f,-h])
    bound=S.factor(-(vec.T*pair.inv()*vec)[0])
    target=(2*f*f+2*h*h-2*g*f*h)/(g*g-4)
    assert S.simplify(bound-target)==0
    assert S.simplify(S.diff(bound,f)-(4*f-2*g*h)/(g*g-4))==0
    assert S.simplify(S.diff(bound,h)-(4*h-2*g*f)/(g*g-4))==0
    assert S.simplify(bound.subs({f:2,h:2})-8/(-g-2))==0
    assert S.simplify((8/(-g-2)).subs(g,gp[2].as_expr())-2/x)==0
    # The full support-cost proof is already independently checked. Bind
    # every recomputed profile and endpoint in this new primary to it.
    budget=primary['independently_recomputed_budget']
    assert budget['total_cost_upper']==upgrade['support_budget']=='281/70'
    for row in budget['A10_profiles']:
        expected=next(r for r in upgrade['A10_profiles'] if r['labels']==row['labels'])
        assert rat(row['cost'])==rat(expected['cost'])
        for key,old in [('pay','pairing'),('margin','margin')]:
            value=rat(row[key][0])+rat(row[key][1])*S.sqrt(2)
            assert S.simplify(value-S.sympify(expected[old]))==0
    assert len(budget['A10_profiles'])==8 and len(budget['excluded_profiles'])==28
    for row in budget['endpoint_cases']:
        expected=next(r for r in upgrade['continuous_endpoint_checks'] if r['wall']==row['wall'])
        assert rat(row['lead'])==rat(expected['linear_coefficient'])
        assert rat(row['margin'])==rat(expected['margin'])
    paths=[f'scripts/{TAG}.py',PRIMARY,BASE,UPGRADE_REVIEW]
    return dict(status='PASS_INDEPENDENT_ALL_2600_CHAIN_TRIPLES_AND_ESCAPE_CUTTER_GATE',
                sources={p:sha(p) for p in paths},gap_cases=len(cases),actual_triples=translations,
                maximum_degree=max(r['degree'] for r in cases),class_counts=dict(counts),
                chain_gap_bounds=gap_records,cases=sorted(cases,key=lambda r:r['gaps']),
                independent_outer_pair_projection=str(bound),
                minimum_dotted_contacts_for_E0_cutter=24,
                scope='Independent symbolic Chebyshev/minor construction, all300 coefficient matrices and kernels, actual-square/product reductions, derivative-based rational gap/determinant enclosures and symbolic outer-pair norm bound. Inherits independently reviewed degree50 over KQab, physical retained134 walls, weak E0 feasibility and complete support-budget proof. Cutter-only condition; no dimension bound or finite-volume example.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['gap_cases'],result['actual_triples'])
