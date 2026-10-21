#!/usr/bin/env python3
"""Independent literal-face and polynomial audit of two local H9 gates.

Uses only the standard library, and imports no producer arithmetic. The
complete original incidence census is inherited; saved vertices, all gate
placements and the decisive determinant identities are checked directly.
"""
import argparse,json
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,permutations

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_dotted_central_20261014_extra_gates_review'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json'
BASE='results/fresh_h9_dotted_central_20261014_gate.json'
STRONG='results/fresh_h9_dotted_central_20261014_strong_contact_gate.json'
DENSE='results/fresh_h9_dotted_central_20261014_dense_pair_gate.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def rank(M):
    a=[list(map(F,r))for r in M];i=0
    for j in range(len(a[0])):
        k=next((k for k in range(i,len(a))if a[k][j]),None)
        if k is None:continue
        a[i],a[k]=a[k],a[i];q=a[i][j];a[i]=[x/q for x in a[i]]
        for k in range(i+1,len(a)):
            q=a[k][j]
            if q:a[k]=[x-q*y for x,y in zip(a[k],a[i])]
        i+=1
        if i==len(a):break
    return i

# An independent small exact polynomial ring for determinant identities.
NV=7
def poly(x):return{(0,)*NV:F(x)}if x else{}
def var(j):e=[0]*NV;e[j]=1;return{tuple(e):F(1)}
def add(a,b):
    r=a.copy()
    for m,c in b.items():r[m]=r.get(m,F(0))+c
    return{m:c for m,c in r.items()if c}
def neg(a):return{m:-c for m,c in a.items()}
def sub(a,b):return add(a,neg(b))
def mul(a,b):
    r={}
    for u,c in a.items():
        for v,d in b.items():
            m=tuple(x+y for x,y in zip(u,v));r[m]=r.get(m,F(0))+c*d
    return{m:c for m,c in r.items()if c}
def square(a):return mul(a,a)
def det(M):
    n=len(M);ans={}
    for p in permutations(range(n)):
        s=(-1)**sum(p[i]>p[j]for i in range(n)for j in range(i+1,n))
        term=poly(s)
        for i in range(n):term=mul(term,M[i][p[i]])
        ans=add(ans,term)
    return ans
def encoded(p):return[dict(exponents=list(m),coefficient=str(c))for m,c in sorted(p.items())]

def run():
    geo,base,strong,dense=[json.loads((ROOT/p).read_text())for p in[GEOM,BASE,STRONG,DENSE]]
    for data in[geo,base,strong,dense]:
        for p,h in data['sources'].items():assert sha(p)==h,p
    g=geo['representatives'][5]['geometry'];assert g['dimension']==9 and g['facet_count']==18
    normals=[list(map(F,r))for r in g['normals']];offsets=list(map(F,g['offsets']))
    V=[frozenset(v['incident'])for v in g['vertices']]
    assert len(V)==444 and len(set(V))==444
    for row,v in zip(g['vertices'],V):
        p=list(map(F,row['point']));assert len(p)==9 and len(v)==9
        slack=[b-sum(a*x for a,x in zip(n,p))for n,b in zip(normals,offsets)]
        assert all(x>=0 for x in slack)
        assert frozenset(i for i,x in enumerate(slack)if x==0)==v
        assert rank([normals[i]for i in v])==9
    def containing(S):return next((i for i,v in enumerate(V)if frozenset(S)<=v),None)
    def missing_triple(S):
        return containing(S)is None and all(containing(t)is not None for t in combinations(S,2))
    blocks=[tuple(b)for b in g['blocks']if len(b)==2]
    outside=sorted(set(range(18))-set(sum((list(b)for b in blocks),[])))
    placements=[]
    for A,Bb in permutations(blocks,2):
        for x in outside:
            if not missing_triple(A+(x,)):continue
            for y in outside:
                if y==x or not missing_triple(Bb+(y,)):continue
                tests=[A+Bb,Bb+(x,),A+(y,)]+[(a,x,y,b)for a in A for b in Bb]
                if all(containing(s)is not None for s in tests):placements.append((A,Bb,x,y))
    assert len(placements)==44
    def key(r):return(tuple(r['dense_block']),tuple(r['other_block']),r['dense_leaf'],r['other_leaf'])
    assert {key(r)for r in base['maps']}==set(placements)
    strong_expected=set();dense_expected=set()
    for A,Bb,x,y in placements:
        for b in Bb:
            assert containing(A+(y,b))is not None
            strong_expected.add((A,Bb,x,y,b))
        assert all(containing((x,y,b))is not None for b in Bb)
        dense_expected.add((A,Bb,x,y))
    assert len(strong_expected)==88 and len(dense_expected)==44
    assert {(key(r)+(r['strong_contact_wall'],))for r in strong['maps']}==strong_expected
    assert {key(r)for r in dense['maps']}==dense_expected
    used=set();support_occurrences=0
    for data,kind in[(strong,'strong'),(dense,'dense')]:
        expected_conditions=set()
        for r in data['maps']:
            A,Bb,x,y=key(r)
            for s in r['actual_supports']:
                i=s['actual_vertex'];assert frozenset(s['support'])<=V[i]
                used.add(i);support_occurrences+=1
            if kind=='strong':
                b=r['strong_contact_wall']
                cond=[sorted(A)+[4],sorted(Bb)+[6],sorted((x,A[0]))+[3],sorted((x,A[1]))+[3],sorted((y,b))+[4]]
            else:
                cond=[sorted(A)+[6],sorted(Bb)+[4]]+[sorted((x,a))+[3]for a in A]+[sorted((y,b))+[3]for b in Bb]
            assert r['bounds']==sorted(cond)
            expected_conditions.add(tuple(map(tuple,sorted(cond))))
        assert {tuple(map(tuple,r['bounds']))for r in data['families']}==expected_conditions
        assert len(data['families'])==len(expected_conditions)==len(data['maps'])
        for f in data['families']:
            assert f['fixed']==[] and f['bounds']==data['maps'][f['witness_map']]['bounds']
    one=poly(1);zero=poly(0)
    alpha,beta,delta,H,Y,FF,K=[var(i)for i in range(NV)]
    path=[[one,neg(alpha),zero,zero],[neg(alpha),one,neg(delta),zero],
          [zero,neg(delta),one,neg(beta)],[zero,zero,neg(beta),one]]
    path_det=sub(mul(sub(one,square(alpha)),sub(one,square(beta))),square(delta))
    assert det(path)==path_det
    triangle=[[one,neg(alpha),neg(beta)],[neg(alpha),one,neg(delta)],[neg(beta),neg(delta),one]]
    triangle_det=sub(sub(sub(sub(one,square(alpha)),square(beta)),square(delta)),mul(poly(2),mul(alpha,mul(beta,delta))))
    assert det(triangle)==triangle_det
    dense_min_numerator=mul(poly(F(1,2)),add(one,alpha))
    assert mul(dense_min_numerator,mul(poly(2),sub(one,alpha)))==sub(one,square(alpha))
    residual=[[sub(sub(one,H),FF),neg(K)],[neg(K),sub(one,Y)]]
    simple_lower=sub(mul(sub(H,one),sub(Y,one)),FF)
    assert sub(det(residual),simple_lower)==sub(mul(FF,Y),square(K))
    # All thresholds are strict/exact. No label cutoff is substituted.
    assert 1-F(3,4)-F(1,4)==0 # A/B label >=6 cannot join a finite third root.
    assert F(1,4)/(1-F(1,2))+F(1,2)==1 # Strong contact rules out y-A.
    assert (1-F(1,4))*(1-F(1,2))==F(3,8)<F(1,2) # Strong path: xy<=3.
    assert F(1,2)/(1-F(3,4))==2 # Strong y projection lower bound.
    assert 2*F(1,2)**2-F(1,4)**2==F(7,16)>0 # sqrt2/2-1/4>0.
    assert F(3,4)>F(25,36) and F(1,2)>F(4,9) # H>3,Y>3/2.
    assert 4*5>9 # cos^2(pi/5)>9/16: no-foreign xy<=4.
    assert 1-F(1,2)-F(1,2)==0 # Foreign label >=4 impossible.
    assert 5<F(49,9) # cos^2(pi/5)<2/3: F<3/4.
    assert (3-1)*(F(3,2)-1)-F(1,2)==F(1,2)
    assert (3-1)*(F(3,2)-1)-F(3,4)==F(1,4)
    paths=[GEOM,BASE,STRONG,DENSE,
        'scripts/fresh_h9_dotted_central_20261014_strong_contact_gate.py',
        'scripts/fresh_h9_dotted_central_20261014_dense_pair_gate.py',f'scripts/{TAG}.py']
    return dict(status='PASS_INDEPENDENT_LITERAL_FACE_AND_SCHUR_AUDIT',sources={p:sha(p)for p in paths},
        saved_vertices_checked=444,placements_reconstructed=44,strong_placements=88,dense_placements=44,
        support_occurrences=support_occurrences,distinct_containing_vertices=len(used),
        polynomial_variables=['alpha','beta','delta','H','Y','F','K'],
        identities=dict(four_path=encoded(path_det),three_node=encoded(triangle_det),
            foreign_residual_minus_lower=encoded(sub(mul(FF,Y),square(K)))),
        strict_determinant_bounds=dict(strong_Qsqrt2=['-1/4','1/2'],dense_no_foreign='1/2',dense_foreign='1/4'),
        scope='Independent literal coordinates, placements and exact determinant/threshold arithmetic; complete original incidence census inherited. All-label sufficient local gates only, not whole case5 exclusion or a rerun of the search.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert data==json.loads(p.read_text())
    else:p.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['strong_placements'],data['dense_placements'],flush=True)
