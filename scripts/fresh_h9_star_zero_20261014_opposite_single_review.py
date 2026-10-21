#!/usr/bin/env python3
"""Independent opposite-single branch audit via tensor interpolation.

No producer imports or symbolic-algebra dependency.  Arithmetic uses the
ordered basis 1,sqrt(2),sqrt(5),sqrt(10), and determinants use permutations.
"""
import argparse
import json
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import permutations, product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_star_zero_20261014_opposite_single'
Z = (F(0),) * 4
O = (F(1), F(0), F(0), F(0))
S = [3, 4, 6, 7, 9, 10, 12, 13]
P = [2, 5, 8, 11]
C = [14, 15, 16, 17]

def add(x, y): return tuple(a+b for a, b in zip(x, y))
def neg(x): return tuple(-a for a in x)
def sub(x, y): return add(x, neg(y))
def scale(x, a): return tuple(a*b for b in x)
def total(xs):
    v = Z
    for x in xs: v = add(v, x)
    return v

@lru_cache(None)
def mul(x, y):
    out = [F(0)]*4
    for i, a in enumerate(x):
        if not a: continue
        for j, b in enumerate(y):
            if b: out[i ^ j] += a*b*(2 if i & j & 1 else 1)*(5 if i & j & 2 else 1)
    return tuple(out)

def sign_quad(a, b, d):
    sg = lambda q: (q > 0)-(q < 0)
    if not b: return sg(a)
    if not a: return sg(b)
    if sg(a) == sg(b): return sg(a)
    return sg(a)*sg(a*a-d*b*b)

@lru_cache(None)
def sign(x):
    a, b, c, d = x
    u, v = sign_quad(a, b, 2), sign_quad(c, d, 2)
    if not v: return u
    if not u: return v
    if u == v: return u
    return u*sign_quad(a*a+2*b*b-5*c*c-10*d*d, 2*a*b-10*c*d, 2)

@lru_cache(None)
def inverse(x):
    conjugates = [(x[0], s*x[1], t*x[2], s*t*x[3])
                  for s, t in [(-1, 1), (1, -1), (-1, -1)]]
    n = O
    for y in conjugates: n = mul(n, y)
    den = mul(x, n)
    assert den[0] and den[1:] == (0, 0, 0)
    return scale(n, 1/den[0])

def dec_anp(row):
    p = [F(v) for v in row]
    p = [F(0)]*(4-len(p))+p
    a3, a2, a1, a0 = p
    return a0+7*a2, a1+17*a3, a1+11*a3, 2*a2

COS = {2: Z, 3: (F(1, 2), F(0), F(0), F(0)),
       4: (F(0), F(1, 2), F(0), F(0)),
       5: (F(1, 4), F(0), F(1, 4), F(0))}
PERMS = [(p, (-1)**sum(p[i] > p[j] for i in range(3) for j in range(i+1, 3)))
         for p in permutations(range(3))]
def determinant3(M):
    return total(scale(mul(mul(M[0][p[0]], M[1][p[1]]), M[2][p[2]]), sg)
                 for p, sg in PERMS)
def sha(p): return sha256((ROOT/p).read_bytes()).hexdigest()

def power_interpolate(values):
    """Degree<=2 polynomial at0,1,2 to its coefficients."""
    f0,f1,f2=values
    return [f0,sub(scale(f1,2),scale(add(scale(f0,3),f2),F(1,2))),
            scale(add(sub(f2,scale(f1,2)),f0),F(1,2))]

def bernstein_interpolate(values):
    """Degree2 Bernstein coefficients from evaluations at0,1/2,1."""
    f0,f1,f2=values
    return [f0,sub(scale(f1,2),scale(add(f0,f2),F(1,2))),f2]

def tensor_interpolate(table,transform):
    first=[transform(row)for row in table]
    columns=[transform([first[i][j]for i in range(3)])for j in range(3)]
    return {(i,j):columns[j][i]for i in range(3)for j in range(3)}

def vertex_failure(labels,record,V):
    nodes=sorted(V[record['actual_vertex']])
    M=[[O if a==b else neg(COS[labels.get(tuple(sorted((a,b))),2)])for b in nodes]for a in nodes]
    for k in range(9):
        pivot=M[k][k]
        if sign(pivot)<=0:
            assert k+1==record['prefix_order'] and pivot==dec_anp(record['last_pivot'])
            return k+1
        for i in range(k+1,9):
            for j in range(k+1,9):M[i][j]=sub(M[i][j],mul(mul(M[i][k],M[k][j]),inverse(pivot)))
    raise AssertionError('claimed nonpositive vertex is positive definite')

def run():
    src=f'results/{TAG}.json';data=json.loads((ROOT/src).read_text())
    pending=[src];seen=set();bound={}
    while pending:
        p=pending.pop()
        if p in seen:continue
        seen.add(p);d=json.loads((ROOT/p).read_text())
        for q,hh in d.get('sources',{}).items():
            assert sha(q)==hh,q;bound[q]=hh
            if q.endswith('.json'):pending.append(q)
    geo=json.loads((ROOT/'results/fresh_h9_tree_cube_wedges_20261014.json').read_text())['representatives'][0]['geometry']
    V=[set(r['incident'])for r in geo['vertices']];missing={frozenset(m)for m in geo['minimal_nonfaces']}
    assert len(V)==378
    perm=data['literal_twin_swap'];expected=list(range(18));expected[12],expected[13]=13,12
    assert perm==expected and {frozenset(perm[i]for i in v)for v in V}==set(map(frozenset,V))
    assert {frozenset(perm[i]for i in m)for m in missing}==missing
    own={tuple(sorted((leaf,a))):3 for i in range(3)for leaf in(P[i],C[i])for a in S[2*i:2*i+2]}
    own[11,12]=own[13,17]=2
    assert data['fixed_other_own_contacts']==sorted([list(p)+[m]for p,m in own.items()])
    assert data['nonright_own_contacts']==[[11,13,4],[12,17,4]]
    assert data['internal_label_lower_bounds']==[S[2*i:2*i+2]+[4]for i in range(4)]
    for i in range(4):
        for leaf in(P[i],C[i]):assert frozenset(S[2*i:2*i+2]+[leaf])in missing
    gates=json.loads((ROOT/'results/fresh_h9_star_zero_20261014_gates.json').read_text())
    assert len(data['arbitrary_high_block_reduction'])==4
    for i,r in enumerate(data['arbitrary_high_block_reduction']):
        assert r['high_block']==S[2*i:2*i+2]
        lab=own|{(11,13):4,(12,17):4}|{tuple(S[2*j:2*j+2]):6 if j==i else 4 for j in range(4)}
        g=next(g for g in gates['families']if g['id']==r['reviewed_gate'])
        assert g['origin']in('two_dense','strong_contact')
        assert all(lab.get((a,b),0)>=m for a,b,m in g['bounds'])
    cross={tuple(sorted((a,b)))for i in range(4)for j in range(i+1,4)for a in S[2*i:2*i+2]for b in S[2*j:2*j+2]}
    assert {tuple(r['pair'])for r in data['cross_right_guards']}==cross
    sqrt2=scale(COS[4],2);null=[O,sqrt2,sqrt2,O]
    path=[[O if i==j else neg(COS[[4,3,4][min(i,j)]])if abs(i-j)==1 else Z for j in range(4)]for i in range(4)]
    assert all(total(mul(path[i][j],null[j])for j in range(4))==Z for i in range(4))
    for row in data['cross_right_guards']:
        a,b=row['pair'];i=S.index(a)//2;j=S.index(b)//2
        assert set(row['support'])==set(S[2*i:2*i+2]+S[2*j:2*j+2])
        assert set(row['support'])<=V[row['actual_vertex']] and row['relaxed_path']==[4,3,4]
    wanted=[(11,13,12),(17,12,13),(14,4,3)]
    assert len(data['projection_budgets'])==3
    for row,(leaf,keep,drop)in zip(data['projection_budgets'],wanted):
        assert row['leaf']==leaf and row['own_retained']==keep
        assert set(row['support'])==(set(S)-{drop})|{leaf}
        assert set(row['support'])<=V[row['actual_vertex']]
    # Anchors have retained own cost>=1/2, plus foreign minimum2*(1/2)^2.
    # Third has own cost1/4: two foreign contacts or one>=4 both violatePD.
    assert F(1,2)+2*F(1,4)==1
    assert F(1,4)+2*2*F(1,4)>1 and F(1,4)+2*F(1,2)>1
    mr=data['mutual_upper_guard'];assert mr['support']==[11,13,17]
    assert set(mr['support'])<=V[mr['actual_vertex']] and mr['test_label']==4
    M=[[O,neg(COS[4]),neg(COS[4])],[neg(COS[4]),O,Z],[neg(COS[4]),Z,O]]
    assert determinant3(M)==Z
    upper_records=[]
    assert len(data['third_contact_upper_guards'])==4
    for row in data['third_contact_upper_guards']:
        last=row['last_internal'];low=5 if last==4 else 4
        assert last in(4,5)and low==row['minimum_own_label']
        assert tuple(row['pair'])in((11,14),(14,17))and row['test_label']==4
        lab=own|{(11,13):low,(12,17):low}|{tuple(S[2*j:2*j+2]):last if j==3 else 4 for j in range(4)}
        lab[tuple(row['pair'])]=4
        upper_records.append(vertex_failure(lab,row,V))
    # The missing last triples supply the strict lower own-label restrictions.
    assert sub(sub(O,mul(COS[4],COS[4])),mul(COS[4],COS[4]))==Z
    assert sign(sub(sub(O,mul(COS[5],COS[5])),mul(COS[3],COS[3])))>0
    seen_cases=set();positive=rejected=coeff_count=0;minimum=None;records=[]
    for row in data['cases']:
        blocks=tuple(row['block_labels']);extra=row['foreign_third'];mutual=row['mutual_label'];ma,mb=row['third_labels']
        key=(blocks,extra,mutual,ma,mb);assert key not in seen_cases;seen_cases.add(key)
        assert blocks in set(product((4,5),repeat=4))and extra in(None,2,3,4,5,6,7)
        assert mutual in(2,3)and ma in(2,3)and mb in(2,3)
        low=5 if blocks[-1]==4 else 4;lo=COS[low]
        assert row['lower_own_label']==low and dec_anp(row['interval_lower'])==lo and dec_anp(row['interval_upper'])==O
        J=[[Z]*8 for _ in range(8)]
        for i,m in enumerate(blocks):
            c=COS[m];v=inverse(sub(O,mul(c,c)));assert sign(sub(v,scale(O,2)))>=0
            J[2*i][2*i]=J[2*i+1][2*i+1]=v
            J[2*i][2*i+1]=J[2*i+1][2*i]=mul(c,v)
        third=[Z]*8;third[0]=third[1]=neg(COS[3])
        if extra is not None:third[extra]=neg(COS[3])
        def inner(v,w):return total(mul(v[i],total(mul(J[i][j],w[j])for j in range(8)))for i in range(8))
        def evaluate(alpha,beta):
            a=[Z]*8;b=[Z]*8;a[7]=neg(alpha);b[6]=neg(beta)
            P=[a,b,third]
            raw=[[O,neg(COS[mutual]),neg(COS[ma])],[neg(COS[mutual]),O,neg(COS[mb])],[neg(COS[ma]),neg(COS[mb]),O]]
            M=[[sub(raw[i][j],inner(P[i],P[j]))for j in range(3)]for i in range(3)]
            return determinant3(M)
        # Each variable belongs to just one added root; a determinant has
        # degree at most2 in it. Nine evaluations therefore determine it.
        grid=[Z,O,scale(O,2)]
        poly=tensor_interpolate([[evaluate(a,b)for b in grid]for a in grid],power_interpolate)
        saved={tuple(r['exponents']):dec_anp(r['coefficient'])for r in row['polynomial']}
        assert {k:v for k,v in poly.items()if v!=Z}==saved
        if row['status']=='STRICT_POSITIVE_BERNSTEIN':
            grid=[lo,scale(add(lo,O),F(1,2)),O]
            bern=tensor_interpolate([[evaluate(a,b)for b in grid]for a in grid],bernstein_interpolate)
            saved={tuple(r['index']):dec_anp(r['coefficient'])for r in row['bernstein']}
            assert bern==saved and all(sign(v)>0 for v in bern.values())
            for v in bern.values():
                if minimum is None or sign(sub(v,minimum))<0:minimum=v
            positive+=1;coeff_count+=len(bern)
        else:
            assert row['status']=='NONPOSITIVE_ACTUAL_VERTEX_AT_MINIMUM_OWN_LABELS'
            lab=own|{tuple(S[2*i:2*i+2]):blocks[i]for i in range(4)}|{(11,13):low,(12,17):low,(11,17):mutual,(11,14):ma,(14,17):mb}
            if extra is not None:lab[tuple(sorted((14,S[extra])))]=3
            vertex_failure(lab,row['vertex_failure'],V);rejected+=1
        records.append(dict(blocks=blocks,foreign=extra,contacts=[mutual,ma,mb],status=row['status']))
    assert seen_cases==set((b,e,x,y,z)for b in product((4,5),repeat=4)for e in(None,2,3,4,5,6,7)for x,y,z in product((2,3),repeat=3))
    assert positive==data['positive_polynomial_cases']==832 and rejected==data['actual_vertex_rejections']==64
    assert coeff_count==data['positive_bernstein_coefficients']==7488
    assert minimum==dec_anp(data['minimum_bernstein_coefficient'])
    paths=[src,f'scripts/{TAG}.py',f'scripts/{TAG}_review.py']
    return dict(status='PASS_INDEPENDENT_OPPOSITE_SINGLE_RECTANGLE_AND_DOMAIN_REVIEW',
        sources={p:sha(p)for p in paths},method='Standalone Fraction basis; full Schur matrices; permutation determinants; polynomial and Bernstein coefficients reconstructed by tensor interpolation, not producer expansion.',
        recursive_source_manifests=len(seen),actual_vertices=len(V),cross_right_null_witnesses=24,
        foreign_budget_actual_supports=3,high_label_inherited_gate_placements=4,mutual_upper_actual_supports=1,
        fresh_third_upper_vertex_failures=upper_records,polynomials_reconstructed=896,
        positive_Bernstein_cases=positive,positive_Bernstein_coefficients=coeff_count,
        actual_vertex_rejections=rejected,literal_twin_swap_checked=True,records=records,
        scope='Only last internal>=4 and opposite single last-leaf contacts>=4, with first3 blocks own3/3. Last internal3 is excluded from this review. Checks all new literal placements, algebra and certificates; inherits the previously reviewed strong/two-dense gate theorems through verified placements. No whole case0 claim.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}_review.json'
    if a.verify:assert d==json.loads(out.read_text());print('SAVED_REPLAY_PASS')
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['polynomials_reconstructed'],d['positive_Bernstein_coefficients'])

