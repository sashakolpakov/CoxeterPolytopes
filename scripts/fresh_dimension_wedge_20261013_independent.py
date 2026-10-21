#!/usr/bin/env python3
"""Stdlib exact replay of every recorded local wedge obstruction."""
import argparse,hashlib,json
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'results/fresh_dimension_wedge_20261013.json'
PRODUCER=ROOT/'scripts/fresh_dimension_wedge_20261013.py'
OUT=ROOT/'results/fresh_dimension_wedge_20261013_independent.json'
ZERO=(Q(0),Q(0));ONE=(Q(1),Q(0))

def mul(x,y):
    # phi²=phi+1; no symbolic or numerical algebra dependency.
    a,b=x;c,d=y
    return a*c+b*d,a*d+b*c+b*d

def scale(x,c):return tuple(c*v for v in x)

def sign(x):
    a,b=x;u=a+b/2;v=b/2
    if not v:return (u>0)-(u<0)
    if not u:return (v>0)-(v<0)
    if (u>0)==(v>0):return (u>0)-(u<0)
    delta=u*u-5*v*v
    assert delta
    return ((u>0)-(u<0))*((delta>0)-(delta<0))

def run(verify=False):
    d=json.loads(SOURCE.read_text())
    for p,h in d['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
    b=json.loads((ROOT/'data/known_examples/bugaenko_compact_h8.json').read_text())
    g8=[[(Q(x['a']),Q(x['b']))for x in row]for row in b['wall_gram']]
    assert len(g8)==11 and all(sign(g8[i][i])>0 for i in range(11))
    assert all(g8[i][j]==g8[j][i]and sign(g8[i][j])<=0 for i in range(11)for j in range(i))
    c8=d['H8_to_H9']['certificates'];assert sorted(r['wall']for r in c8)==list(range(11))
    for row in c8:
        i,j,m=row['wall'],row['neighbor'],row['forbidden_old_label'];assert i!=j and m in(3,5)
        assert b['coxeter_relations'][i][j]==m and sign(g8[i][j])<0
        ratio=scale(ONE,Q(1,4))if m==3 else(Q(1,4),Q(1,4))
        assert mul(g8[i][j],g8[i][j])==mul(ratio,mul(g8[i][i],g8[j][j]))
    v=json.loads((ROOT/'data/enumerations/vinberg_x3_h19.json').read_text())['wall_gram']
    assert all(type(x)is int for row in v for x in row)
    g19=[[-x for x in row]for row in v]
    b=json.loads((ROOT/'data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json').read_text())['coxeter_input']
    W=[[int(x)for x in row]for row in b['wall_coordinates_as_columns_transposed']]
    F=[[-int(x)for x in row]for row in b['basis_gram']]
    cov=[[sum(a*b for a,b in zip(row,col,strict=True))for col in zip(*F)]for row in W]
    g21=[[sum(a*b for a,b in zip(row,col,strict=True))for col in W]for row in cov]
    results={}
    for key,g,n in(('H19_to_H20',g19,36),('H21_to_H22',g21,210)):
        assert len(g)==n and all(len(row)==n for row in g)
        assert all(g[i][i]>0 for i in range(n))and all(g[i][j]==g[j][i]and g[i][j]<=0 for i in range(n)for j in range(i))
        rows=d[key]['certificates'];assert sorted(r['wall']for r in rows)==list(range(n))
        for row in rows:
            i=row['wall']
            if row['obstruction']=='old finite label3':
                j=row['neighbor'];assert i!=j and g[i][j]<0 and 4*g[i][j]**2==g[i][i]*g[j][j]
            else:
                j,k=row['finite6_neighbor'],row['parallel_neighbor'];assert len({i,j,k})==3
                assert g[i][j]<0 and 4*g[i][j]**2==3*g[i][i]*g[j][j]
                assert g[i][k]<0 and g[i][k]**2==g[i][i]*g[k][k]
                assert row['forced_forbidden_squared_cosine']=='1/3'and Q(1,4)<Q(1,3)<Q(1,2)
        counts=dict(Counter(r['obstruction']for r in rows));assert counts==d[key]['obstruction_counts']
        results[key]={'walls':n,'obstruction_counts':counts}
    # The three unrestricted finite-label possibilities have squared
    # cosine sums 1/2, 3/4, (2+phi)/4, i.e. old labels 4,6,10.
    assert Q(1,4)+Q(1,4)==Q(1,2)and Q(1,4)+Q(1,2)==Q(3,4)
    assert tuple(x+y for x,y in zip(scale(ONE,Q(1,4)),(Q(1,4),Q(1,4))))==(Q(1,2),Q(1,4))
    assert d['general_necessary_old_finite_labels']==[4,6,10]
    out={'status':'PASS_INDEPENDENT_EXACT_WEDGE_SOURCE_REPLAY','sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(SOURCE,PRODUCER)},
      'H8_to_H9':{'walls':11,'forbidden_old_labels':dict(Counter(str(r['forbidden_old_label'])for r in c8))},**results,
      'total_excluded_facets':257,'arithmetic':'Independent stdlib integers/Fractions and two-coefficient golden-field operations; no NumPy or SymPy.',
      'mathematical_review':'The full-rank positive-definite triple forces the replacement pair to be right. The resulting cosine-square sum has exactly the possibilities (3,3),(3,4),(3,5), without a label cutoff. The parallel-neighbor obstruction for an old label six is valid. The corrected deletion and cusp side observations give genuine spacelike escapes within their explicitly stated scope.',
      'scope':'All 257 recorded local obstructions independently replayed. This rules out the specified single-facet two-slope lift, not arbitrary dimension raising; no new root search or polytope classification is performed.'}
    if verify:assert json.loads(OUT.read_text())==out
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['total_excluded_facets'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');run(p.parse_args().verify)
