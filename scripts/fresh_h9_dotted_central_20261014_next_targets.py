#!/usr/bin/env python3
"""Literal eligibility table for three proved gates; no label search."""
import argparse
import json
from hashlib import sha256
from itertools import permutations
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_dotted_central_20261014_next_targets'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json'
SUMMARY='results/fresh_h9_tree_cube_wedges_20261014_summary.json'
GATES={
    'dense_high':'results/fresh_h9_dotted_central_20261014_gate.json',
    'strong_contact':'results/fresh_h9_dotted_central_20261014_strong_contact_gate.json',
    'two_dense':'results/fresh_h9_dotted_central_20261014_dense_pair_gate.json'}
CLOSED={1,2,3,4,5,8,17,22,25,26}
CLOSURES=['results/fresh_h9_dotted_central_20261014_independent.json',
    'results/fresh_h9_t8_sparse_search_20261014_independent_dfs.json']+[
    f'results/fresh_h9_other_trees_20261014_independent_dfs_case_{i:02d}.json' for i in (17,22,25,26)]
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def mask(v):return sum(1<<x for x in v)


def placements(g):
    vertices=[set(v['incident']) for v in g['vertices']]
    vm=[mask(v) for v in vertices]
    missing={frozenset(s) for s in g['minimal_nonfaces']}
    blocks=[b for b in g['blocks'] if len(b)==2]
    assert len(blocks)==4
    outside=sorted(set(range(g['facet_count']))-set(sum(blocks,[])))
    cache={}
    def actual(kind,s):
        s=sorted(s);assert len(set(s))==len(s)
        key=tuple(s)
        if key not in cache:
            m=mask(s);cache[key]=next((i for i,v in enumerate(vm) if m&v==m),None)
        i=cache[key]
        if i is None:return None
        assert set(s)<=vertices[i]
        return dict(kind=kind,support=s,actual_vertex=i)
    maps={k:[] for k in GATES}
    for A,B in permutations(blocks,2):
        for x in outside:
            if frozenset(A+[x]) not in missing:continue
            for y in outside:
                if y==x or frozenset(B+[y]) not in missing:continue
                specs=[('positive_four_block',A+B),('x_foreign_right',B+[x]),('y_foreign_right',A+[y])]
                specs += [('middle_label_upper',[x,y,a,b]) for a in A for b in B]
                supports=[actual(k,s) for k,s in specs]
                if any(s is None for s in supports):continue
                def record(name,bounds,extras):
                    rows=supports+extras
                    if any(s is None for s in rows):return
                    bounds=sorted(bounds)
                    # All thresholded pairs really meet, independent of labels.
                    pairs=[]
                    for a,b,m in bounds:
                        assert a<b and frozenset((a,b)) not in missing
                        w=actual('threshold_pair',[a,b]);assert w is not None
                        pairs.append(w)
                    maps[name].append(dict(A=A,B=B,x=x,y=y,
                        missing_own_triples=[sorted(A+[x]),sorted(B+[y])],
                        supports=rows,threshold_pair_witnesses=pairs,bounds=bounds))
                record('dense_high',[sorted(A)+[6],sorted(B)+[7]]+
                    [sorted((x,a))+[3] for a in A],[])
                for b in B:
                    record('strong_contact',[sorted(A)+[4],sorted(B)+[6]]+
                        [sorted((x,a))+[3] for a in A]+[sorted((y,b))+[4]],
                        [actual('strong_contact_forces_foreign_right',A+[y,b])])
                record('two_dense',[sorted(A)+[6],sorted(B)+[4]]+
                    [sorted((x,a))+[3] for a in A]+[sorted((y,b))+[3] for b in B],
                    [actual('foreign_contact_forces_xy_right',[x,y,b]) for b in B])
    out={}
    for name,rr in maps.items():
        unique={}
        for i,r in enumerate(rr):unique.setdefault(tuple(map(tuple,r['bounds'])),i)
        out[name]=dict(literal_placements=len(rr),unique_conditions=len(unique),maps=rr,
            conditions=[dict(bounds=[list(x) for x in key],witness_map=i) for key,i in sorted(unique.items())])
    return out


def run():
    paths=[GEOM,SUMMARY,*GATES.values(),*CLOSURES,
        'results/fresh_h9_dotted_central_20261014_gate_review.json',
        'results/fresh_h9_dotted_central_20261014_extra_gates_review.json']
    source={}
    for p in paths:
        d=json.loads((ROOT/p).read_text());source[p]=digest(p)
        for q,h in d.get('sources',{}).items():assert digest(q)==h,q
    geo=json.loads((ROOT/GEOM).read_text())
    summary=json.loads((ROOT/SUMMARY).read_text())
    for path,case in zip(CLOSURES,(5,8,17,22,25,26)):
        proof=json.loads((ROOT/path).read_text())
        assert proof['complete'] and proof['cutoff'] is None
        assert proof['survivor'] is None and proof['unexplored_frontier'] is None
        assert proof.get('case',case)==case
        assert proof['vertex_count']==len(geo['representatives'][case]['geometry']['vertices'])
    old={r['case']:r for r in summary['cases']}
    assert len(geo['representatives'])==33 and set(summary['excluded_cases'])=={1,2,3,4}
    rows=[]
    for rep in geo['representatives']:
        case=rep['id'];g=rep['geometry'];gate=placements(g)
        row=dict(case=case,tree=rep['tree'],cut_order=rep['cut_order'],
            unwedged_coordinate=rep['unwedged_coordinate'],vertices=len(g['vertices']),
            missing_size_counts=g['missing_size_counts'],closed=case in CLOSED,
            earlier_status=old[case]['status'],earlier_nodes=old[case]['nodes'],
            gates=gate,total_unique_conditions=sum(v['unique_conditions'] for v in gate.values()))
        rows.append(row)
        print(case,'CLOSED'if row['closed']else'ELIGIBLE',row['tree'],
            {k:(v['literal_placements'],v['unique_conditions']) for k,v in gate.items()},flush=True)
    # The formerly proved case5 counts are an exact template regression.
    for name,expected in [('dense_high',(44,24)),('strong_contact',(88,88)),('two_dense',(44,44))]:
        g=rows[5]['gates'][name]
        assert (g['literal_placements'],g['unique_conditions'])==expected
        oldgate=json.loads((ROOT/GATES[name]).read_text())
        assert {tuple(map(tuple,r['bounds'])) for r in g['maps']}=={tuple(map(tuple,r['bounds'])) for r in oldgate['maps']}
    eligible=[r for r in rows if not r['closed']]
    assert len(eligible)==23
    ranking=sorted(eligible,key=lambda r:(-r['total_unique_conditions'],r['vertices'],r['case']))
    source[f'scripts/{TAG}.py']=digest(f'scripts/{TAG}.py')
    return dict(status='PASS_LITERAL_GATE_ELIGIBILITY_ONLY_NO_SEARCH',sources=source,
        closed_cases=sorted(CLOSED),eligible_cases=[r['case'] for r in eligible],
        ranked_by_gate_count_then_vertex_count=[r['case'] for r in ranking],cases=rows,
        template_scope='Ordered pairs of the four declared two-facet wedge blocks, with distinct outside leaves and their missing own triples. Every required finite support has a literal containing saved vertex. No global incidence symmetry or extension to arbitrary non-block pairs is assumed.',
        theorem_scope='Counts certify where the three previously proved all-label gates may be applied IF their threshold hypotheses hold. They do not assert those thresholds, exclude a case, establish a surviving label assignment, or provide a metric. No new label search was run.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(out.read_text());print('SAVED_REPLAY_PASS')
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['ranked_by_gate_count_then_vertex_count'])
