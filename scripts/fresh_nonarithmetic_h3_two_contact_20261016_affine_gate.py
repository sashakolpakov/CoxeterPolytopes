#!/usr/bin/env python3
"""Exact permanent-cusp exclusion of the pilot's null-line recipes."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
BASE='fresh_nonarithmetic_h3_two_contact_20261016';TAG=BASE+'_affine_gate'

def run():
    data=json.loads((ROOT/f'results/{BASE}.json').read_text())
    old=json.loads((ROOT/'results/fresh_h9_persistent_cusp_screen_20261016.json').read_text())
    cap=json.loads((ROOT/'results/fresh_nonarithmetic_escape52_20261016_cap.json').read_text())
    matrix=old['pair_matrix'];relations={r['wall']:r for r in cap['relations']}
    assert matrix[0][1]==5
    def label(i,j):
        if i==j:return 1
        if max(i,j)<49:return matrix[i][j]
        assert max(i,j)==52
        r=relations[min(i,j)];assert r['kind']in('right','finite');return r['label']
    cos={2:s.Integer(0),3:s.Integer(1),4:s.sqrt(2),6:s.sqrt(3)}
    rows=[];families=[]
    recipes=sorted({(tuple(f['anchors']),tuple(f['labels']))for f in data['families']if f['kind']=='line'})
    for anchors,labs in recipes:
        i,j=anchors;m=label(i,j)
        assert all(label(k,h)==2 for k in anchors for h in(0,1))
        assert sum((F(1,l)for l in(m,*labs)),F(0))==1
        G=s.Matrix([[2,-cos[m],-cos[labs[0]]],[-cos[m],2,-cos[labs[1]]],[-cos[labs[0]],-cos[labs[1]],2]])
        K=G.nullspace();assert len(K)==1 and G.rank()==2
        v=K[0]
        if v[0].is_negative:v=-v
        assert all(x.is_positive for x in v)and G*v==s.zeros(3,1)
        assert all(G.extract(list(p),list(p)).det().is_positive for p in((0,1),(0,2),(1,2)))
        ids=[f['index']for f in data['families']if tuple(f['anchors'])==anchors and tuple(f['labels'])==labs and f['kind']=='line']
        families.extend(ids);rows.append(dict(anchors=anchors,new_labels=labs,old_label=m,
                                              norm_two_gram=[[str(x)for x in row]for row in G.tolist()],
                                              positive_null_kernel=[str(x)for x in v],family_indices=ids))
    assert len(rows)==23 and len(families)==25
    paths=[f'scripts/{TAG}.py',f'results/{BASE}.json','results/fresh_h9_persistent_cusp_screen_20261016.json',
           'results/fresh_nonarithmetic_escape52_20261016_cap.json','results/fresh_nonarithmetic_escape52_20261016_cap_review.json',
           'results/fresh_h9_fixed_core_cusp_obstruction_20261013.md',
           'results/fresh_h20_weighted_binary_20261014_dotted_completion_review.md']
    return dict(status='PASS_EXACT_ALL_PARAMETER_NULL_LINE_CUSP_OBSTRUCTION',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
                forbidden_retained_edge=[0,1],forbidden_label=5,recipes=rows,excluded_family_indices=sorted(families),
                proof='For every real norm-two new normal right toS6 with any listed two finite contacts, the displayed three-root Gram has a positive null kernel. A common strict future interior makes its positive null vector nonzero and future. All acute retained/additional walls have nonpositive product with it, so it persists in every acute completion. Both oldpi5 roots0,1 are orthogonal to it. A finite-volume Euclidean cusp would have a full translation lattice, but the product of these two reflections has irrational linear trace d-2+2cos(2pi/5), impossible for an integral lattice automorphism.',
                scope='Exact conditional family obstruction for23 recipes, covering all25 numerical null-line branches and all their parameters. Does not classify the remaining1119 conics or exclude arbitrary rightS6 extensions.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(p.read_text());print('PASS saved exact affine-recipe gate')
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],len(d['recipes']),len(d['excluded_family_indices']))
