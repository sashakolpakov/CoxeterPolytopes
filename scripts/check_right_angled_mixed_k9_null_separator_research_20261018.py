#!/usr/bin/env python3
"""Tiny exact checks supporting the written K9 null-separator argument."""
import argparse
from fractions import Fraction as Q
from itertools import combinations
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
CERT=ROOT/'results/right_angled_mixed_k9_null_separator_research_20261018.json'


def generate():
    # Coordinates are three times the Euclidean h-perp/h root vectors.
    differences={tuple(3*(int(k==i)-int(k==j)) for k in range(9))
                 for i in range(9) for j in range(9) if i!=j}
    triples={tuple(1-3*int(k in chosen) for k in range(9))
             for chosen in combinations(range(9),3)}
    roots=differences|triples|{tuple(-x for x in root) for root in triples}
    assert len(differences)==72 and len(triples)==84 and len(roots)==240
    assert all(sum(root)==0 and sum(x*x for x in root)==18 for root in roots)
    products=set()
    for root in roots:
        for other in roots:
            product=sum(x*y for x,y in zip(root,other));products.add(product)
            assert product%9==0
            assert tuple(x-(product//9)*y for x,y in zip(other,root)) in roots
    assert products=={-18,-9,0,9,18}
    assert all(sum(root[k] for root in triples)==0 for k in range(9))
    assert all(Q(m,9-m)<4 for m in range(1,8))
    assert Q(15,14)==Q(1,2)+Q(4,7) and 14<16
    # A permitted positive-height normal: a=sqrt35, all nine c_i=2.
    assert 9*2**2-35==1 and 35>5**2 and 18**2>9*35
    return dict(status='EXACT_K9_NULL_SEPARATOR_SUPPORTING_CHECKS',
        positive_support_maximum_if_at_most_seven='7/2 < 4',
        orthogonal_pair_time_square_upper='14',
        scaled_root_configuration=dict(roots=240,differences=72,triple_directions=84,
            squared_norm=18,reflection_pairs_checked=240**2,possible_products=sorted(products)),
        algebraically_permitted_cutoff_normal=dict(time_coefficient='sqrt(35)',
            nine_spatial_coefficients='2',h_product='18-3sqrt(35)>0',
            products_with_old_normals='-2',products_with_bad_triangle_normals='4-sqrt(35)<-1'),
        scope=dict(universal_separator_proof_written=True,finite_volume_completion_constructed=False,
            complete_k9_corner_excluded=False,dimension_ceiling_improved=False,
            optimizer_used=False,lean_compilation_performed=False))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write',action='store_true')
    args=parser.parse_args();record=generate()
    if args.write:CERT.write_text(json.dumps(record,indent=2)+'\n')
    else:
        with CERT.open() as stream:saved=json.load(stream)
        assert saved==record
    print('PASS null-separator arithmetic and 57600 finite root reflections; no K9 or dimension exclusion')


if __name__=='__main__':main()
