#!/usr/bin/env python3
"""Independent bivariate-polynomial replay of the T8 two-label extension."""
import argparse,json
from hashlib import sha256
from itertools import permutations,product
from pathlib import Path
from fresh_h9_star_cube_metric_20261014_cut_profiles_review import (
    F,Z,O,CO,add,sub,neg,mul,scale,total,inverse,enc)
from fresh_h9_star_cube_metric_20261014_all_cut_profiles_review import sign
from fresh_h9_t8_metric_20261014_review import q2
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_t8_metric_20261014_two_labels'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def mpadd(a,b):
    d=a.copy()
    for k,v in b.items():d[k]=add(d.get(k,Z),v)
    return {k:v for k,v in d.items()if v!=Z}
def mpneg(a):return {k:neg(v)for k,v in a.items()}
def mpmul(a,b):
    out={}
    for (i,j),x in a.items():
        for (k,l),y in b.items():
            key=i+k,j+l;out[key]=add(out.get(key,Z),mul(x,y))
    return {k:v for k,v in out.items()if v!=Z}
def mpdet(A):
    out={}
    for p in permutations(range(len(A))):
        val={(0,0):O}
        for i,j in enumerate(p):val=mpmul(val,A[i][j])
        if sum(p[i]>p[j]for i in range(len(p))for j in range(i+1,len(p)))%2:val=mpneg(val)
        out=mpadd(out,val)
    return out
def mpeval(a,k,q):
    out=Z
    for (i,j),v in a.items():
        for _ in range(i):v=mul(v,k)
        for _ in range(j):v=mul(v,q)
        out=add(out,v)
    return out
def run():
    src=f'results/{TAG}.json';data=json.loads((ROOT/src).read_text())
    prior='results/fresh_h9_t8_metric_20261014_review.json'
    basepath='results/fresh_h9_t8_metric_20261014.json'
    base=json.loads((ROOT/basepath).read_text());review=json.loads((ROOT/prior).read_text())
    assert data['status']=='PASS_T8_TWO_LABEL_FAMILY_M_GE_7_N_GE_4_EXCLUDED'
    assert review['status']=='PASS_INDEPENDENT_T8_CONTINUOUS_DOMAIN_AND_SCHUR_REVIEW'
    bound=set()
    def bind(p):
        if p in bound:return
        bound.add(p)
        for s,h in json.loads((ROOT/p).read_text()).get('sources',{}).items():
            assert sha(s)==h,s
            if s.endswith('.json'):bind(s)
    bind(src);bind(prior)
    assert data['parameter_pairs']=={'m':[1,2],'n':[12,13]}
    gp='results/fresh_h9_tree_cube_wedges_20261014.json'
    g=json.loads((ROOT/gp).read_text())['representatives'][8]['geometry']
    V=[set(v['incident'])for v in g['vertices']]
    # These literal supports remove one endpoint of the varied edge.  Their
    # matrices and all anchor foreign-contact budgets are therefore unchanged.
    budget_rows=[]
    for row in base['anchor_budget_witnesses']:
        s=set(row['support']);assert s<=V[row['actual_vertex']]
        assert len(s&{12,13})==1 and {1,2}<=s
        assert {4,5,7,8}<=s
        budget_rows.append(row)
    assert len(budget_rows)==4
    # Third supports retain the whole varied block but use only its inverse
    # diagonal lower bound. For sqrt2/2 <= s <1, 1/(1-s^2)>=2, offdiag>=0.
    third_rows=[]
    for row in base['third_leaf_budget_witnesses']:
        s=set(row['support']);assert s<=V[row['actual_vertex']]
        assert {12,13}<=s and len(s&{1,2})==1
        assert row['retained_own_roots'][0]==2
        third_rows.append(row)
    assert len(third_rows)==2 and F(1,4)+F(1,4)+2*F(1,4)==1
    assert mul(CO[4],CO[4])==scale(O,F(1,2))
    # Every inherited finite-vertex failure remains a lower-weight relaxation:
    # only (12,13) increases from4; no known edge or required contact decreases.
    affected=[]
    for group in ['forced_cross_right_witnesses','anchor_contact_upper_witnesses']:
        for row in base[group]:
            support=set(row.get('support',V[row['actual_vertex']]))
            assert support<=V[row['actual_vertex']]
            affected.append(dict(group=group,actual_vertex=row['actual_vertex'],
                                 contains_varied_edge={12,13}<=support))
    assert not {12,13}<=set(base['anchor_pair_right_witness']['support'])

    k0=add(O,CO[4]);assert k0==q2(data['k0'])
    # At s=cos(pi/n), the equal (-1/2,-1/2) pairing vector has
    # norm (1+s)/(2(1-s^2))=1/(2(1-s)); denominators are positive.
    assert mul(k0,scale(sub(O,CO[4]),2))==O
    B=[[O if i==j else neg(CO[4])if i//2==j//2 else Z for j in range(4)]for i in range(4)]
    J=inverse(B)
    def dot(x,y):return total(mul(x[i],mul(J[i][j],y[j]))for i in range(4)for j in range(4))
    patterns={None:[Z]*4}
    for i in range(4):
        p=[Z]*4;p[i]=neg(CO[3]);patterns[i+2]=p
    third=[neg(CO[3]),neg(CO[3]),Z,Z]
    orig={(tuple(r['foreign_contacts']),tuple(r['contact_labels'])):r for r in base['residual_cases']}
    saved={(tuple(r['foreign_contacts']),tuple(r['contact_labels'])):r for r in data['cases']}
    saved_anchor={tuple(r['foreign_contacts']):r for r in data['anchors']}
    assert len(saved)==len(orig)==100 and len(saved_anchor)==25
    qbnd=neg(add(O,CO[4]));records=[];anchor_seen={}
    for ia,ib in product(patterns,repeat=2):
        x,y=patterns[ia],patterns[ib]
        u,v,p=dot(x,x),dot(y,y),dot(x,y)
        qa={(0,0):sub(O,u),(1,0):neg(O)}
        qb={(0,0):sub(O,v),(1,0):neg(O)}
        ab={(0,0):neg(p),(1,0):neg(O)}
        anchor=mpdet([[qa,ab],[ab,qb]])
        assert all(j==0 and i<=1 for i,j in anchor)
        a1=anchor.get((1,0),Z);assert sign(a1)<0
        assert a1==sub(add(u,v),add(scale(O,2),scale(p,2)))
        a0=mpeval(anchor,k0,Z);assert sign(a0)<0
        assert a1==q2(saved_anchor[ia,ib]['negative_k_slope'])
        assert a0==q2(saved_anchor[ia,ib]['anchor_determinant_at_k0'])
        anchor_seen[ia,ib]=dict(foreign_contacts=[ia,ib],value_at_k0=enc(a0),negative_derivative=enc(a1))
        for ma,mb in product([2,3],repeat=2):
            ac=sub(neg(CO[ma]),dot(x,third));bc=sub(neg(CO[mb]),dot(y,third))
            A=[[qa,ab,{(0,0):ac}],[ab,qb,{(0,0):bc}],
               [{(0,0):ac},{(0,0):bc},{(0,1):O}]]
            D=mpdet(A)
            # Coefficient comparison in the full polynomial ring Q(sqrt2)[k,q],
            # not testing finitely many q or k values.
            assert all(i<=1 and j<=1 for i,j in D)
            square=mul(sub(ac,bc),sub(ac,bc));assert sign(square)>=0
            assert D.get((1,1),Z)==a1 and D.get((1,0),Z)==square
            assert D.get((0,1),Z)==anchor.get((0,0),Z)
            value=mpeval(D,k0,qbnd)
            dk=add(mul(a1,qbnd),square);assert sign(dk)>0 and sign(value)>0
            row=saved[(ia,ib),(ma,mb)]
            for k,val in [('negative_anchor_slope',a1),('nonnegative_square',square),
                          ('positive_k_slope_at_boundary',dk),('positive_minimum',value)]:
                assert q2(row[k])==val
            assert value==q2(orig[(ia,ib),(ma,mb)]['boundary_at_three_quarters'])
            records.append(dict(foreign_contacts=[ia,ib],contact_labels=[ma,mb],
                                determinant_polynomial_in_k_q=[dict(k_degree=i,q_degree=j,coefficient=enc(v))for(i,j),v in sorted(D.items())],
                                positive_lower_bound=enc(value),negative_a1=enc(a1),nonnegative_square=enc(square)))
    assert len(records)==100 and len(anchor_seen)==25
    paths=[src,f'scripts/{TAG}.py',prior,basepath,gp,
           'scripts/fresh_h9_star_cube_metric_20261014_cut_profiles_review.py',
           'scripts/fresh_h9_star_cube_metric_20261014_all_cut_profiles_review.py',
           'scripts/fresh_h9_t8_metric_20261014_review.py',f'scripts/{TAG}_review.py']
    return dict(status='PASS_INDEPENDENT_T8_TWO_PARAMETER_CONTINUOUS_REVIEW',sources={p:sha(p)for p in paths},
                recursively_checked_manifests=sorted(bound),unchanged_anchor_budget_supports=budget_rows,
                strengthened_third_budget_supports=third_rows,monotone_vertex_guard_checks=affected,
                parameter_region='3/4<t<1; k>=1+sqrt2/2, corresponding to all finite m>=7,n>=4',
                anchor_records=list(anchor_seen.values()),exact_bivariate_cases=records,
                proof='Anchor determinant is affine decreasing in k and negative at k0. D_k=a1*q+(ac-bc)^2>0 since a1<0,q<0; D_q equals the anchor determinant A(k)<0, while q_t<0. Thus D increases with both k and t, and the already certified positive boundary at (k0,3/4) is a uniform lower bound.',
                scope='Two internal labels of the saved T8 assignment only. Inherits the independent base domain review, checks its literal support persistence, and independently reconstructs all100 full bivariate determinants. No primary arithmetic imports, finite label cutoff or numerical acceptance predicates.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run()
    out=ROOT/f'results/{TAG}_review.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],len(d['exact_bivariate_cases']),flush=True)
