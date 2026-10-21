#!/usr/bin/env python3
"""Exact finite equality cases for the spectator difference energy through3/4."""
import argparse,json
from hashlib import sha256
from pathlib import Path
import sympy as S

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_small_spectator_difference_gap_20261017'

def pd(g):
    pivots=[];a=g.copy()
    while a.rows:
        p=S.simplify(a[0,0]);assert p.is_positive;pivots.append(str(p))
        a=S.simplify(a[1:,1:]-a[1:,0:1]*a[0:1,1:]/p)
    return pivots

def run():
    rows=[]
    specs=[('A1_single',1,[],[1],[0],S.Rational(1,2)),
        ('A2_single',2,[1],[1,0],[0,0],S.Rational(2,3)),
        ('A3_endpoint_single',3,[1,1],[1,0,0],[0,0,0],S.Rational(3,4)),
        ('A2_adjacent',2,[1],[1,0],[0,1],S.Rational(2,3)),
        ('A3_adjacent',3,[1,1],[1,0,0],[0,1,0],S.Rational(3,4)),
        ('I2_4_adjacent',2,[S.sqrt(2)],[1,0],[0,1],2-S.sqrt(2)),
        ('I2_5_adjacent',2,[(1+S.sqrt(5))/2],[1,0],[0,1],1-S.sqrt(5)/5)]
    for name,n,weights,x,y,want in specs:
        g=2*S.eye(n)
        for i,w in enumerate(weights):g[i,i+1]=g[i+1,i]=-w
        x=S.Matrix(x);y=S.Matrix(y);z=x-y;J=g.inv();energy=S.simplify((z.T*J*z)[0])
        assert S.simplify(energy-want)==0
        assert energy.is_positive and(S.Rational(3,4)-energy).is_nonnegative
        ep=[]
        for contact in(x,y):
            h=g.row_join(-contact).col_join((-contact.T).row_join(S.Matrix([[2]])))
            ep.append(pd(h))
        rows.append(dict(name=name,gram=[[str(v)for v in r]for r in g.tolist()],
            x=list(map(str,x)),y=list(map(str,y)),inverse=[[str(v)for v in r]for r in J.tolist()],
            energy=str(energy),positive_pivots=pd(g),separate_extension_pivots=ep,
            inverse_difference=list(map(str,J*z))))
    # For a third node of attachment magnitude u appended at an A2 end,
    # Schur complementation makes the difference energy monotone in u^2.
    v=S.symbols('v',positive=True)
    extension=S.Rational(2,3)+v/(6*(3-v))
    assert S.simplify(extension-S.Rational(3,4)-(v-1)/(4*(3-v)))==0
    assert S.simplify(S.diff(extension,v)-1/(2*(v-3)**2))==0
    a3=next(r for r in rows if r['name']=='A3_adjacent')
    assert a3['inverse_difference']==['1/4','-1/2','-1/4']
    # Exact rational bounds: sqrt2>7/5, sqrt3<7/4 force a nonunit
    # contact pair to have energy>96/125>3/4.
    assert S.Rational(7,5)**2<2 and S.Rational(7,4)**2>3
    assert (1+S.Rational(7,5))**2/(4+2*S.Rational(7,4))==S.Rational(96,125)>S.Rational(3,4)
    values=[S.Integer(0),S.Rational(1,2),S.Rational(2,3),S.Rational(3,4),2-S.sqrt(2),1-S.sqrt(5)/5]
    assert all(S.simplify(S.Rational(4,7)-q)!=0 for q in values)
    # The exceptional noncrystallographic/doubled values straddle4/7 exactly.
    assert S.simplify(S.Rational(4,7)-(1-S.sqrt(5)/5)).is_positive
    assert S.simplify((2-S.sqrt(2))-S.Rational(4,7)).is_positive
    stems=['fresh_selective_distance9_dimension108_20261017','fresh_ideal_distance6_mixed_replacement_20261017']
    paths=[Path(__file__)]+[ROOT/f'{folder}/{stem}.{ext}'for stem in stems
        for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_DISCRETE_SPECTATOR_ENERGY_AT_MOST_THREE_QUARTERS',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        hypotheses='Z and each separate norm2 extension Z+c,Z+d positive definite; native acute Coxeter contacts; nonnegative contact vectors x,y have disjoint supports.',
        possible_energies=list(map(str,values)),literal_patterns=rows,
        A2_third_node_energy=str(extension),A2_excess_over_three_quarters=str(S.factor(extension-S.Rational(3,4))),
        nonunit_pair_lower_bound='strictly greater than96/125>3/4',
        excluded_cycle14_energy='4/7',
        equality_three_quarters='Exactly one unit contact at an endpoint of a whole simple A3 component, or one unit contact from each side at adjacent vertices of a whole simple A3 component. All other spectator contacts vanish.',
        scope='Universal small spectator-energy classification under the explicit strict edge and disjoint-support hypotheses. Completeness is proved analytically in the report; finite cases check its exact values and sharpness. No complete face-payment or dimension theorem asserted.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],len(d['literal_patterns']),'literal equality patterns',flush=True)
