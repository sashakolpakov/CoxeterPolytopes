#!/usr/bin/env python3
"""Exact finite-contact E7/E8 triangle screen over Q(sqrt5).

This is a local diagnostic, not a dimension exclusion or a global chamber.
No existing producer is imported or modified.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import product,combinations
import json
from pathlib import Path
import sympy as s
from sympy.polys.polyerrors import IsomorphismFailed

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_arithmetic_quadratic_bound_20261015'


def add(x,y):return(x[0]+y[0],x[1]+y[1])
def mul(x,y):return(x[0]*y[0]+x[1]*y[1],x[0]*y[1]+x[1]*y[0]+x[1]*y[1])
def scale(x,k):return(x[0]*k,x[1]*k)
def neg(x):return(-x[0],-x[1])
def sign(x):
    # a+b phi=(a+b/2)+(b/2)sqrt5.
    a=x[0]+x[1]/2;b=x[1]/2
    sa=(a>0)-(a<0);sb=(b>0)-(b<0)
    if not sb:return sa
    if not sa or sa==sb:return sb
    d=a*a-5*b*b
    return sa*((d>0)-(d<0))
def conj(x):return(x[0]+x[1],-x[1])
def enc(x):return[str(a)for a in x]


def gram(n):
    G=s.eye(n)*2
    for i,j in [(i,i+1)for i in range(n-2)]+[(2,n-1)]:G[i,j]=G[j,i]=-1
    return G


def old_quadrilaterals():
    """Reject only these literal saved recipes, not all quad completions."""
    old=json.loads((ROOT/'results/fresh_compact29_small_faces_20261015.json').read_text())
    def decode(v):
        return sum(s.Rational(a)*b for a,b in zip(v,[1,s.sqrt(2),s.sqrt(5),s.sqrt(10)]))
    phi=(1+s.sqrt(5))/2
    mag={2:s.Integer(0),3:s.Integer(1),4:s.sqrt(2),5:phi}
    rows=[]
    for row in old['cases']:
        rec=row['recipe'];proof=row['proof']
        if not rec['type'].startswith('E')or rec['kind']!='quadrilateral':continue
        a,b=map(decode,proof['dotted_magnitude'])
        y=a+b*s.sqrt(decode(proof['positive_radicand']))
        if rec['cd_label']==2:
            expression=s.expand(y*y);kind='dotted two-cycle'
        else:
            # A nonzero c-to-old contact plus the unique old E-tree path
            # closes a literal cycle. Every old path edge has magnitude1.
            contact_label=rec['ac_label']
            if contact_label==2:
                assert rec['c_contacts'];contact_label=rec['c_contacts'][0][1]
            assert contact_label!=2
            expression=s.expand(y*mag[rec['cd_label']]*mag[contact_label])
            kind='fundamental cycle through dotted edge'
        try:
            exact=s.to_number_field(expression,s.sqrt(5)).as_expr()
            polynomial=s.Poly(s.minpoly(exact))
            assert polynomial.LC()!=1,(rec,polynomial)
            reason='cyclic product in Q(sqrt5) but not an algebraic integer'
            detail=str(polynomial.as_expr())
        except IsomorphismFailed:
            reason='cyclic product outside Q(sqrt5)';detail=None
        rows.append(dict(type=rec['type'],marked=rec['marked'],cycle_kind=kind,
                         obstruction=reason,primitive_integer_minimal_polynomial=detail))
    assert len(rows)==23
    assert sum('outside'in row['obstruction']for row in rows)==20
    return rows


def run():
    records=[]
    for n in(7,8):
        G=gram(n);J=[[F(v)for v in row]for row in G.inv().tolist()]
        assert all(v>0 for row in J for v in row)
        minimum=min(J[i][i]for i in range(n))
        assert minimum==(F(3,2)if n==7 else F(2))
        # The other three uniform squareclasses (labels4,6,10) have
        # conjugate magnitudes squared2,3,3-phi respectively.
        pure_class_margins=[add(scale(q,minimum),(-F(2),F(0)))for q in[(F(2),F(0)),(F(3),F(0)),(F(3),-F(1))]]
        assert all(sign(v)>0 for v in pure_class_margins)
        qmask=[]
        for mask in range(1<<n):
            ids=[i for i in range(n)if mask>>i&1]
            qmask.append(sum((J[i][j]for i in ids for j in ids),F(0)))
        accepted=[];subcritical=[];counts={'total':0,'conjugate_not_positive':0,'physical_not_indefinite':0,'accepted':0}
        mark_witnesses={pair:[]for pair in combinations(range(n),2)}
        for labels in product((0,1,2),repeat=n):
            counts['total']+=1
            ones=sum(1<<i for i,v in enumerate(labels)if v==1)
            phis=sum(1<<i for i,v in enumerate(labels)if v==2)
            cross=(qmask[ones|phis]-qmask[ones]-qmask[phis])/2
            cost=(qmask[ones]+qmask[phis],qmask[phis]+2*cross)
            if sign(add((F(2),F(0)),neg(conj(cost))))<=0:
                counts['conjugate_not_positive']+=1;continue
            if sign(add(cost,(-F(2),F(0))))<=0:
                counts['physical_not_indefinite']+=1
                subcritical.append(dict(contacts=[[i,3 if v==1 else 5]for i,v in enumerate(labels)if v],cost=enc(cost)))
                continue
            counts['accepted']+=1
            b=[(F(0),F(0))if v==0 else(F(1),F(0))if v==1 else(F(0),F(1))for v in labels]
            dual=[tuple(sum(J[i][j]*b[j][k]for j in range(n))for k in range(2))for i in range(n)]
            deletion_costs=[add(cost,neg(scale(mul(dual[i],dual[i]),1/J[i][i])))for i in range(n)]
            good=[i for i in range(n)if sign(add((F(2),F(0)),neg(deletion_costs[i])))>0]
            pairs=list(combinations(good,2))
            index=len(accepted)
            for pair in pairs:mark_witnesses[pair].append(index)
            accepted.append(dict(contacts=[[i,3 if v==1 else 5]for i,v in enumerate(labels)if v],
                                 cost=enc(cost),conjugate_cost=enc(conj(cost)),
                                 finite_deletions=good,triangle_marked_pairs=pairs,
                                 deletion_costs=list(map(enc,deletion_costs))))
        records.append(dict(type=f'E{n}',inverse=[[str(v)for v in row]for row in J],
                            pure_class_exclusion_margins=list(map(enc,pure_class_margins)),counts=counts,
                            profiles=accepted,subcritical_profiles=subcritical,
                            marked_pairs=[dict(pair=pair,profiles=ix)for pair,ix in mark_witnesses.items()],
                            triangle_covered_pairs=sum(bool(ix)for ix in mark_witnesses.values())))
    assert records[0]['subcritical_profiles']==[{'contacts':[],'cost':['0','0']},{'contacts':[[5,3]],'cost':['3/2','0']}]
    assert records[1]['subcritical_profiles']==[{'contacts':[],'cost':['0','0']}]
    # Ambient spectator check: component of leaf5 after each E7 deletion.
    edges=[(i,i+1)for i in range(5)]+[(2,6)]
    component_sizes=[]
    for omitted in range(7):
        reach=set()if omitted==5 else{5}
        while True:
            nxt=reach|{v for u,v in edges if u in reach and v!=omitted}|{u for u,v in edges if v in reach and u!=omitted}
            if nxt==reach:break
            reach=nxt
        component_sizes.append(len(reach))
    assert component_sizes==[6,5,3,2,1,0,6]
    forced=[i for i,z in enumerate(component_sizes)if z<=2]
    assert forced==[3,4,5]
    covered={tuple(r['pair'])for r in records[0]['marked_pairs']if r['profiles']}
    assert set(combinations(forced,2))<=covered
    # Explicit genuine subcritical rescue by an extra A1 with a pi5 contact.
    extra=(F(1,2),F(1,2)) # phi^2/2
    total=add((F(3,2),F(0)),extra)
    assert sign(add(total,(-F(2),F(0))))>0 and sign(add((F(2),F(0)),neg(conj(total))))>0
    J=[[F(a)for a in row]for row in records[0]['inverse']]
    rescued_deletions=[]
    for a in range(7):
        cost=add((J[5][5]-J[5][a]**2/J[a][a],F(0)),extra)
        if sign(add((F(2),F(0)),neg(cost)))>0:rescued_deletions.append(a)
    assert rescued_deletions==forced
    source_paths=[f'scripts/{TAG}.py','scripts/fresh_compact29_small_faces_20261015.py',
                  'results/fresh_compact29_small_faces_20261015.json']
    return dict(status='PASS_EXACT_QSQRT5_FINITE_CONTACT_TRIANGLE_SCREEN',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in source_paths},
                records=records,literal_prior_quad_obstructions=old_quadrilaterals(),
                spectator_lemma=dict(E7_leaf5_component_sizes=component_sizes,possible_marked_nodes=forced,
                                     explicit_A1_rescue_total_cost=enc(total),rescue_conjugate_cost=enc(conj(total)),
                                     rescued_finite_deletions=rescued_deletions,
                                     proof='A nonrational spectator cost forces a label5 in the connected finite corner containing leaf5, the new root and a spectator; rank<=4 forces the retained E7 leaf component to have size<=2.'),
                scope='Complete finite-contact triangle marked-pair necessities at literal E7/E8 in arithmetic ground field Q(sqrt5), including arbitrary orthogonal spherical spectator components by the analytic lemma. The 23 saved quad recipes are excluded individually. No exhaustive dotted opposite-side profiles, simultaneous stars, global chamber or dimension bound is certified.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'])
    for r in out['records']:print(r['type'],r['counts'],'triangle pairs',r['triangle_covered_pairs'])
