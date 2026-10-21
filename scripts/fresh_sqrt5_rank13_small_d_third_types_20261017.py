#!/usr/bin/env python3
"""Literal non-ADE third-vertex costs incident to three high-budget D forests."""
import argparse,hashlib,json
from collections import Counter,defaultdict
from fractions import Fraction as F
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_rank13_small_d_third_types_20261017'
INPUT='results/fresh_sqrt5_rank13_triangle_third_costs_20261017.json'
ENDPOINT='results/fresh_sqrt5_rank13_high_ade_endpoints_20261017.json'
SELECTED=[['D4','D4','D5'],['A3','D5','D5'],['A1','D4','D4','D4']]

def run():
    data=json.loads((ROOT/INPUT).read_text());end=json.loads((ROOT/ENDPOINT).read_text())
    for doc in(data,end):
        for p,h in doc['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    selected={tuple(f)for f in SELECTED};lookup={tuple(t['forest']):t for t in end['targets']}
    all_thirds={tuple(t['forest']):t for t in data['third_forests']}
    rows=[];targets=[]
    for i,p in enumerate(data['triangle_patterns']):
        hits=[j for j,c in enumerate(p['selected_corners'])if tuple(c['forest'])in selected]
        if not hits:continue
        third=p['third_corner'];forest=tuple(third['forest'])
        exceptional=[name for name in forest if name[0]not in 'ADE']
        assert len(exceptional)==1 and exceptional[0]in ['H3','H4','I2(5)','I2(10)']
        corner=F(third['corner_baseline']);vertex=F(third['vertex_baseline'])
        assert vertex==F(all_thirds[forest]['baseline'])
        required=F(p['required_selected_charge']);assert required==2-corner
        assert (corner,required)==((F(2),F(0))if exceptional==['H4']else(F(1),F(1)))
        rows.append(dict(pattern_index=i,selected_occurrences=hits,
            selected_corners=p['selected_corners'],third_corner=third,exceptional_component=exceptional[0],
            third_vertex_margin_over12=str(vertex-12),required_selected_charge=str(required),
            occurrences=p['occurrences'],witness=p['witness']))
    for forest in SELECTED:
        matching=[r for r in rows if any(c['forest']==forest for c in r['selected_corners'])]
        thirdmap=defaultdict(list)
        for r in matching:thirdmap[tuple(r['third_corner']['forest'])].append(r)
        thirds=[]
        for f,patterns in sorted(thirdmap.items()):
            base=F(all_thirds[f]['baseline'])
            thirds.append(dict(forest=list(f),baseline=str(base),margin_over12=str(base-12),
                pattern_indices=[r['pattern_index']for r in patterns],
                exceptional_component=patterns[0]['exceptional_component'],
                required_selected_charges=sorted({r['required_selected_charge']for r in patterns})))
        targets.append(dict(forest=forest,index=lookup[tuple(forest)]['index'],
            baseline=lookup[tuple(forest)]['baseline'],pattern_count=len(matching),third_forest_count=len(thirds),
            exceptional_pattern_counts=dict(sorted(Counter(r['exceptional_component']for r in matching).items())),
            pattern_indices=[r['pattern_index']for r in matching],third_forests=thirds))
    assert [(t['pattern_count'],t['third_forest_count'])for t in targets]==[(19,11),(20,12),(11,6)]
    union={tuple(r['third_corner']['forest'])for r in rows}
    paths=[f'scripts/{TAG}.py',INPUT,ENDPOINT]
    return dict(status='PASS_EXACT_SMALL_D_THIRD_VERTEX_COST_EXTRACTION',
        sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},targets=targets,
        matching_pattern_count=len(rows),distinct_third_forest_count=len(union),records=rows,
        scope='Complete extraction from the frozen marked triangle catalogue for any occurrence of the three named initial types. The other ADE endpoint still ranges over the entire215-type set. These are literal necessary local triangle patterns and their third-vertex baseline costs, not new geometric realizations, a quadrilateral classification, or a dimension bound.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run()
    out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],[(t['forest'],t['pattern_count'],t['third_forest_count'])for t in d['targets']],flush=True)
