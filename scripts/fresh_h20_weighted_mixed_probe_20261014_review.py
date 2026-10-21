#!/usr/bin/env python3
"""Independent small arithmetic review; no producer or project imports.

The prior 88-cutter classification/adjacency theorem is a named dependency,
not repeated here.  This checks the new Hamming covering counts, integer
graph bounds, strictness constants, and all-label coordinate extension.
"""
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from math import comb
from pathlib import Path
import argparse
import json

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_weighted_mixed_probe_20261014'


def sign_phi(a, b):
    """Sign of a+b*phi, using only rational operations."""
    x, y = F(a) + F(b, 2), F(b, 2)
    if not y:
        return (x > 0) - (x < 0)
    if not x:
        return (y > 0) - (y < 0)
    if x*y > 0:
        return (x > 0) - (x < 0)
    norm = x*x - 5*y*y
    assert norm != 0
    return ((x > 0) - (x < 0))*((norm > 0) - (norm < 0))


def run():
    primary = json.loads((ROOT/f'results/{TAG}.json').read_text())
    assert primary['status'] == 'PASS_ALL_LABEL_MIXED_R2_AND_ALL_VERTEX_RADII_OBSTRUCTION'
    for path, digest in primary['sources'].items():
        assert sha256((ROOT/path).read_bytes()).hexdigest() == digest
    prior_path = ROOT/'results/fresh_h20_seed_deformation_20261014.json'
    prior = json.loads(prior_path.read_text())
    for path, digest in prior['sources'].items():
        assert sha256((ROOT/path).read_bytes()).hexdigest() == digest
    assert prior['R4_vertex_family']['maximum_compatible_caps'] == 512
    assert prior['R2_nonvertex_family']['cutters_per_half_probe'] == 88
    assert prior['R2_nonvertex_family']['adjacent_golden_probe_pair_counts'] == {
        'forbidden': 6484, 'right': 1260}

    # All-label coordinate extension: the smaller opposite gap is <=b<2,
    # hence its finite label is <=3,4,5 for b=1,sqrt2,phi respectively.
    # Only the following two complementary gaps need exclusion beyond the
    # old finite alphabet. Both fall strictly between labels7 and8.
    sqrt2_lo, sqrt2_hi = F(1414213, 10**6), F(1414214, 10**6)
    phi_lo, phi_hi = F(1618033, 10**6), F(1618034, 10**6)
    assert sqrt2_lo**2 < 2 < sqrt2_hi**2
    assert phi_lo**2-phi_lo-1 < 0 < phi_hi**2-phi_hi-1
    forbidden_gaps = [(2*sqrt2_lo-1, 2*sqrt2_hi-1),
                      (2*phi_lo-sqrt2_hi, 2*phi_hi-sqrt2_lo)]
    assert all(F(91, 50) < lo < hi < F(46, 25)
               for lo, hi in forbidden_gaps)
    x = F(3, 7)
    cos7_upper = 2*(1-x*x/2+x**4/24)
    assert cos7_upper < F(91, 50)
    assert 2+sqrt2_lo > F(46, 25)**2

    # |vertex-halfvertex|² = N/4+2H, N/4=5+2phi.
    # Radius-four covering is H <= 11/2-phi.
    rows = []
    by_golden = [0]*9
    boundary = 0
    for h1, h2, h5 in product(range(11), range(2), range(9)):
        sign = sign_phi(F(h1+2*h2+h5)-F(11, 2), h5+1)
        boundary += sign == 0
        if sign <= 0:
            multiplicity = comb(10, h1)*comb(8, h5)
            by_golden[h5] += multiplicity
            rows.append([h1, h2, h5, multiplicity])
    assert boundary == 0
    assert by_golden == [187, 88, 0, 0, 0, 0, 0, 0, 0]
    assert sum(by_golden) == 275
    assert primary['R4_half_probe_cover']['cover_rows'] == rows
    assert primary['R4_half_probe_cover']['counts_by_golden_difference'] == by_golden
    assert primary['all_label_R2_gate']['smaller_gap_max_labels'] == [3, 4, 5]

    # For A in Q8: e(A)<=|A| and e(A)>=8(|A|-128).
    max_uncovered = (8*128)//7
    assert max_uncovered == 146
    assert 8*(max_uncovered-128) <= max_uncovered
    assert 8*(max_uncovered+1-128) > max_uncovered+1
    fibers = 2**11
    required = fibers*(256-max_uncovered)
    possible = 512*sum(by_golden)
    assert required == 225280 and possible == 140800
    assert required > possible
    assert primary['mixed_union']['minimum_uncovered_half_probes'] == required-possible
    assert primary['mixed_union']['maximum_total_covered'] == 2**19-(required-possible)

    # Every non-R4 seed-compatible vertex radius has square<=8.
    # The closest half-probe distance is N/4=5+2phi>8.
    assert sign_phi(-3, 2) > 0
    # Half probes are strictly outside the central radius-two ball.
    assert sign_phi(1, 2) > 0
    # Every covering R2 recipe is a strict cutter, not a boundary case.
    assert sign_phi(8, -5) < 0

    paths = [f'scripts/{TAG}.py', f'results/{TAG}.json',
             f'results/{TAG}.md', 'scripts/fresh_h20_seed_deformation_20261014.py',
             'results/fresh_h20_seed_deformation_20261014.json',
             f'scripts/{TAG}_review.py']
    return dict(
        status='PASS_INDEPENDENT_MIXED_PROBE_COUNTS_AND_PROOF_REVIEW',
        sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        hamming_trials=198, covering_rows=rows,
        half_probes_per_R4_by_golden_difference=by_golden,
        boundary_probe_occurrences=boundary,
        max_uncovered_per_golden_fiber=max_uncovered,
        fibers=fibers, required_R4_covered=required,
        maximum_R4_covered=possible, contradiction_margin=required-possible,
        all_label_extension=dict(
            smaller_gap_max_labels=[3, 4, 5],
            complementary_gap_bounds=[[str(a), str(b)] for a, b in forbidden_gaps],
            twice_cos_pi7_upper=str(cos7_upper),
            twice_cos_pi8_lower='46/25'),
        vertex_radius_gate='R=4 or R^2<=8; all latter miss every half-probe strictly',
        scope='Independent Fraction count/strictness checks and mathematical source review. '
              'The prior R2 coordinate/cutter classification and adjacency gates are inherited; '
              'the two complementary-gap exclusions extend it to every finite label. '
              'No second full cutter enumeration is claimed. The family allows arbitrary '
              'seed-compatible vertex-centered radii and all-label R2 nonvertex caps, but '
              'not arbitrary nonvertex radii, changed seed/central core, or arbitrary caps.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    data = run()
    path = ROOT/f'results/{TAG}_review.json'
    if args.verify:
        assert json.loads(path.read_text()) == data
    else:
        path.write_text(json.dumps(data, indent=2)+'\n')
    print(data['status'])
