#!/usr/bin/env python3
"""Exact 120-case certificate for the written maximum-energy carrier defect."""
import argparse
from fractions import Fraction as Q
from itertools import combinations
from math import comb
import json
from pathlib import Path
from check_right_angled_mixed_small_three_face_dimension_research_20261018 import sign_sqrt2

ROOT=Path(__file__).resolve().parents[1]
CERT=ROOT/'results/right_angled_mixed_max_energy_defect_research_20261018.json'


def generate():
    limits=[8,7,5,4,2,2];forbidden=[9,8,6,5,3,3]
    # r/9+g*(sqrt2-1)^2-1 >=0 rules out the indicated clique.
    for g,r in enumerate(forbidden):
        assert sign_sqrt2(Q(r,9)+3*g-1,-2*g)>=0
    assert sign_sqrt2(17,-12)>0  # 6*(sqrt2-1)^2-1
    dimensions={}
    for n in (11,12):
        rows=[]
        for k in range(n+1):
            m=n-k
            for g in range(min(m,5)+1):
                r=limits[g];small,remainder=divmod(k,r)
                parts=[small+1]*remainder+[small]*(r-remainder)
                edges=sum(a*b for a,b in combinations(parts,2))
                triangles=sum(a*b*c for a,b,c in combinations(parts,3))
                lower=comb(k,3)+m*comb(k,2)+comb(m,3)-comb(g,3)-triangles-g*edges
                rows.append(dict(k=k,m=m,g=g,clique_limit=r,
                    balanced_parts=[part for part in parts if part],maximum_edges=edges,
                    maximum_triangles=triangles,large_carrier_lower=lower))
        bound=min(row['large_carrier_lower'] for row in rows)
        assert bound=={11:27,12:40}[n]
        dimensions[str(n)]=dict(cases=len(rows),local_large_carrier_lower=bound,
            minimizers=[{key:row[key] for key in ('k','m','g')} for row in rows if row['large_carrier_lower']==bound],rows=rows)
    assert dimensions['11']['cases']==57 and dimensions['12']['cases']==63
    return dict(status='EXACT_MAXIMUM_ENERGY_LARGE_CARRIER_DEFECT',
        clique_limits_by_good_finite_edge_count=limits,forbidden_cliques=forbidden,
        dimensions=dimensions,
        scope=dict(finite_vertex_required=True,universal_geometric_proof_written=True,
            global_positive_excess_fraction_proved=False,new_example_claim=False,
            dimension_exclusion=False,optimizer_used=False,lean_compilation_performed=False))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write',action='store_true')
    args=parser.parse_args();record=generate()
    if args.write:CERT.write_text(json.dumps(record,indent=2)+'\n')
    else:
        with CERT.open() as stream:saved=json.load(stream)
        assert saved==record
    print('PASS 120 exact local cases: max-energy vertex has >=27/40 large carriers in n11/12; no dimension exclusion')


if __name__=='__main__':main()
