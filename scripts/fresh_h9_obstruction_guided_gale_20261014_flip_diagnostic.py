#!/usr/bin/env python3
"""Exact restricted 5/5 flip diagnostic; no realization claim for new spheres."""
import argparse,json
from collections import Counter
from itertools import combinations
from fresh_h9_rank5_gale_20261014_core import ROOT,sha,simplex_link_factors
from fresh_h9_bugaenko_surgery_20261014_incidence_independent import clutter
from fresh_h9_obstruction_guided_gale_20261014_supports import SupportGates
TAG='fresh_h9_obstruction_guided_gale_20261014_flip_diagnostic'
SOURCE='results/fresh_h9_obstruction_guided_gale_20261014.json'


def run():
    d=json.loads((ROOT/SOURCE).read_text())
    for p,h in d['sources'].items():assert sha(p)==h,p
    selected=sorted(d['cases'],key=lambda r:(r['recognized_support_count'],r['case']))[:20]
    gates=SupportGates();counts=Counter();rows=[]
    for ref in selected:
        assert sha(ref['path'])==ref['sha256']
        parent=json.loads((ROOT/ref['path']).read_text());g=parent['geometry']
        M=g['minimal_nonfaces'];V=set(g['facets'])
        carriers=sorted({a|b for a,b in combinations(V,2) if (a|b).bit_count()==10})
        for X in carriers:
            P=sum(1<<i for i in range(15) if X>>i&1 and X^(1<<i) in V);N=X^P
            if P.bit_count()!=5 or N.bit_count()!=5 or P not in M:continue
            old={v for v in V if v&N==N}
            expected={X^(1<<i) for i in range(15) if P>>i&1}
            if old!=expected:continue
            new=sorted((V-old)|{X^(1<<i) for i in range(15) if N>>i&1})
            mm=list(clutter(new,15));counts['exact_bistellar_recipes']+=1
            row=dict(seed=ref['path'],carrier=X,old_missing_face=P,replaced_face=N,
                     old_vertex_count=len(V),new_vertex_count=len(new))
            if max(map(int.bit_count,mm))>5:
                counts['large_minimal_nonface']+=1
                row.update(verdict='LARGE_MINIMAL_NONFACE',witness=next(m for m in mm if m.bit_count()>5))
            elif simplex_link_factors(new):
                counts['simplex_link_factor']+=1;row['verdict']='TRUNCATION_LINK_FILTER'
            else:
                obs=gates.all(mm);counts['passes_incidence_filters']+=1
                row.update(minimal_nonfaces=mm,recognized_support_count=len(obs))
                assert obs  # This completed bounded diagnostic found none free.
                counts['exact_literal_obstruction']+=1
                row.update(verdict='EXACT_INHERITED_OBSTRUCTION',gate=obs[0])
            rows.append(row)
    assert dict(counts)==dict(exact_bistellar_recipes=625,large_minimal_nonface=543,
                             passes_incidence_filters=82,exact_literal_obstruction=82)
    paths=[SOURCE,f'scripts/{TAG}.py','scripts/fresh_h9_obstruction_guided_gale_20261014_supports.py',
           'scripts/fresh_h9_rank5_gale_20261014_core.py','scripts/fresh_h9_bugaenko_surgery_20261014_incidence_independent.py']+gates.paths
    return dict(status='PASS_RESTRICTED_COMBINATORIAL_FLIP_DIAGNOSTIC',sources={p:sha(p) for p in paths},
                selected_seed_paths=[r['path'] for r in selected],counts=dict(counts),records=rows,
                scope='All combinatorially valid 5/5 bistellar flips of the selected20 saved nerves are inspected, not all flip sizes or all parents. New complexes are simplicial spheres by the exact star replacement, but no Euclidean polytopality or Gale realization is asserted. These625 recipes are not added to the74 newly realized Gale nerves.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
    path=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'])
