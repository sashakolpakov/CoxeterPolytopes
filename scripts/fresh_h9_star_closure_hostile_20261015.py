#!/usr/bin/env python3
"""Small independent literal/control audit of the star0 closure.

Does not rerun the163894-node DFS or import any project arithmetic.
"""
import argparse
import ast
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_star_closure_hostile_20261015'
BASE='fresh_h9_star_large_block_resume_20261015'


def run():
    paths=[f'scripts/{BASE}_{s}.py' for s in ('gates','search','independent')]+[
        f'results/{BASE}_gates.json',f'results/{BASE}_search_batch_016.json',
        f'results/{BASE}_independent.json',
        'results/fresh_h9_tree_cube_wedges_20261014.json',
        'scripts/fresh_h9_tree_cube_wedges_20261014.py',
        'results/fresh_h9_tree_cube_wedges_20261014_summary.json',
        'results/fresh_h9_tree_cube_wedges_20261014_label_replay.json',
        'results/fresh_h9_dotted_central_20261014_next_targets.json',
        'scripts/fresh_h9_dotted_central_20261014_next_targets.py',
    ]
    seen={}
    def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
    def read(p):
        if p in seen:return seen[p]
        d=json.loads((ROOT/p).read_text());seen[p]=d
        for q,h in d.get('sources',{}).items():
            assert sha(q)==h,q
            if q.endswith('.json'):read(q)
        return d
    for p in paths:
        if p.endswith('.json'):read(p)
    g=read('results/fresh_h9_tree_cube_wedges_20261014.json')
    rep=g['representatives'][0]['geometry']
    missing={frozenset(m) for m in rep['minimal_nonfaces']}
    vertices={tuple(v['incident']) for v in rep['vertices']}
    reconstructed={v for v in combinations(range(18),9) if not any(m<=set(v) for m in missing)}
    assert reconstructed==vertices and len(vertices)==378
    triangles=sorted(tuple(sorted(m)) for m in missing if len(m)==3)
    assert len(triangles)==8
    categories=[t for t in product(range(2,8),repeat=3) if sum((F(1,m) for m in t),F(0))<1]
    assert len(categories)==175
    # Every category containing7 has constant spherical/affine/hyperbolic
    # status throughout all represented finite labels. The sole spherical
    # case is(2,2,m), and every other such triple is strictly hyperbolic.
    assert F(1,2)+F(1,3)+F(1,7)<1
    for row in product(range(2,8),repeat=3):
        if 7 not in row:continue
        expected=sorted(row)!=[2,2,7]
        assert (row in categories)==expected
    twins=[tuple(b) for b in rep['blocks'] if len(b)==2]
    assert len(twins)==4
    all_twins=set(sum((list(b) for b in twins),[]))
    for a,b in twins:
        t=next(t for t in triangles if a in t and b in t)
        assert next(v for v in t if v not in (a,b)) not in all_twins
    symmetry_count=0
    for choices in product((0,1),repeat=4):
        p=list(range(18))
        for (a,b),take in zip(twins,choices):
            if take:p[a],p[b]=b,a
        assert {tuple(sorted(p[i] for i in v)) for v in vertices}==vertices
        assert {frozenset(p[i] for i in m) for m in missing}==missing
        symmetry_count+=1
    clauses=[]
    for a,b in combinations(missing,2):
        if a&b or any(c<=a|b and c not in (a,b) for c in missing):continue
        cross={tuple(sorted((i,j))) for i in a for j in b}
        assert all(any(set(p)<=set(v) for v in vertices) for p in cross)
        clauses.append(sorted(cross))
    assert len(clauses)==38
    binder=read(f'results/{BASE}_gates.json')
    old=read('results/fresh_h9_star_pi4_resume_20261015_gates.json')
    new=read('results/fresh_h9_star_three_large_blocks_20261015.json')
    assert binder['families']==old['families']+new['families']
    assert len(old['families'])==493 and len(new['families'])==4
    families=binder['families']
    assert len({(tuple(map(tuple,f['fixed'])),tuple(map(tuple,f['bounds']))) for f in families})==497
    fixed_labels=sorted({m for f in families for a,b,m in f['fixed']})
    lower_labels=sorted({m for f in families for a,b,m in f['bounds']})
    assert fixed_labels==[2,3,4,5] and lower_labels==[3,4,5,6,7]
    assert all(any({a,b}<=set(t) for t in triangles)
               for f in families for a,b,m in f['fixed']+f['bounds'])
    blocks=[(3,4),(6,7),(9,10),(12,13)]
    assert {tuple(map(tuple,f['bounds'])) for f in new['families']}=={
        tuple(b+(6,) for j,b in enumerate(blocks) if j!=i) for i in range(4)}
    assert all(not f['fixed'] for f in new['families'])
    # Review link binds the exact new primary; no assertion by status alone.
    nr=read('results/fresh_h9_star_three_large_blocks_20261015_review.json')
    assert nr['sources']['results/fresh_h9_star_three_large_blocks_20261015.json']==sha('results/fresh_h9_star_three_large_blocks_20261015.json')
    independent=read(f'results/{BASE}_independent.json')
    primary=read(f'results/{BASE}_search_batch_016.json')
    assert independent['status']=='PASS_INDEPENDENT_COMPLETE_STAR0_EXCLUSION'
    assert independent['complete'] and independent['cutoff'] is None
    assert independent['survivor'] is None and independent['unexplored_frontier'] is None
    assert independent['nodes']==163894 and independent['analytic_conditions']==497
    assert primary['status']=='EXACT_ALL_MANDATORY_LABELS_AND_BRIDGES_EMPTY'
    assert primary['solution'] is None and primary['pending_frames']==[]
    assert primary['new_nodes']==46153 and primary['cumulative_additional_nodes']==166030
    assert primary['new_nodes']<primary['settings']['node_budget']
    syntax=ast.parse((ROOT/f'scripts/{BASE}_independent.py').read_text())
    imported=[]
    for node in ast.walk(syntax):
        if isinstance(node,ast.Import):imported.extend(a.name for a in node.names)
        if isinstance(node,ast.ImportFrom):imported.append(node.module)
    assert set(imported)=={'argparse','json','collections','fractions','functools','hashlib','itertools','pathlib'}
    # Literal aggregate replay: all24 star cut orders and5 wedge choices.
    # Derive the stellar vertex replacement directly, then wedge its facets.
    star_rows=[r for r in g['recipes'] if r['tree']=='star']
    expected={(tuple(order),u) for order in permutations(((0,1),(0,2),(0,3),(0,4))) for u in range(5)}
    assert {(tuple(map(tuple,r['cut_order'])),r['unwedged_coordinate']) for r in star_rows}==expected
    mapped_vertices=0
    for recipe in star_rows:
        v={frozenset(2*i+s for i,s in enumerate(signs)) for signs in product((0,1),repeat=5)}
        for step,(i,j) in enumerate(recipe['cut_order']):
            a,b=2*i,2*j; out=set()
            for face in v:
                if {a,b}<=face:
                    out.add((face-{a})|{10+step});out.add((face-{b})|{10+step})
                else:out.add(face)
            v=out
        # The cut base is flag. Verify this from its complete vertex set.
        miss2=[frozenset(p) for p in combinations(range(14),2) if not any(set(p)<=vv for vv in v)]
        assert {frozenset(w) for w in combinations(range(14),5) if not any(m<=set(w) for m in miss2)}==v
        chosen={2*i+1 for i in range(5) if i!=recipe['unwedged_coordinate']}
        blocks=[];next_id=0
        for i in range(14):
            size=2 if i in chosen else 1
            blocks.append(tuple(range(next_id,next_id+size)));next_id+=size
        w=set()
        for vv in v:
            whole=set(j for i in vv for j in blocks[i])
            for choices in product(*(blocks[i] for i in sorted(chosen-vv))):
                w.add(tuple(sorted(whole|set(choices))))
        clutter={frozenset(j for i in m for j in blocks[i]) for m in miss2}
        p=recipe['facet_isomorphism']; assert sorted(p)==list(range(18))
        target=g['representatives'][recipe['representative']]['geometry']
        assert {tuple(sorted(p[i] for i in vv)) for vv in w}=={tuple(vv['incident']) for vv in target['vertices']}
        assert {frozenset(p[i] for i in m) for m in clutter}=={frozenset(m) for m in target['minimal_nonfaces']}
        mapped_vertices+=len(w)
    multiplicities=Counter(r['representative'] for r in star_rows)
    assert multiplicities=={i:24 for i in range(5)}
    summary=read('results/fresh_h9_tree_cube_wedges_20261014_summary.json')
    assert summary['excluded_cases']==[1,2,3,4]
    replay=read('results/fresh_h9_tree_cube_wedges_20261014_label_replay.json')
    for i in range(1,5):
        row=summary['cases'][i]
        assert row['status']=='EXACT_ALL_MANDATORY_LABELS_AND_BRIDGES_EMPTY'
        assert row['nodes']<row['node_limit'] and row['partial_survivor'] is None
    assert {1,2,3,4}<={r['case'] for r in replay['completed_cases']}
    targets=read('results/fresh_h9_dotted_central_20261014_next_targets.json')
    remaining=set(range(33))-set(targets['closed_cases'])-{0}
    assert len(remaining)==22 and all(g['representatives'][i]['tree']=='T' for i in remaining)
    paths.append(f'scripts/{TAG}.py')
    return dict(status='PASS_SOURCE_AND_LITERAL_AUDIT_STAR_CLOSURE',
                sources={p:sha(p) for p in paths},
                recursively_checked_json_manifests=len(seen),
                independent_traversal_reexecuted=False,
                inherited_complete_empty_run_nodes=independent['nodes'],
                vertices=378,mandatory_triangles=8,triangle_categories=175,
                twin_symmetries=symmetry_count,bridge_clauses=len(clauses),
                predicates=497,predicate_fixed_labels=fixed_labels,predicate_lower_bounds=lower_labels,
                independently_reconstructed_star_recipes=len(star_rows),
                literal_vertex_occurrences=mapped_vertices,
                star_representative_multiplicities={str(k):v for k,v in sorted(multiplicities.items())},
                remaining_T_representatives=sorted(remaining),
                scope='Source-level DFS/binder proof review plus independent small combinatorial replay. Prior analytic predicate proofs, rational incidence realization, and saved complete traversals are inherited. No new expensive DFS or metric search; all120-star conclusion is limited to this240-recipe catalogue.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=run();dest=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert data==json.loads(dest.read_text())
    else:dest.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['independently_reconstructed_star_recipes'])
