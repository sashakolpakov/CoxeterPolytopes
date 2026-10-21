#!/usr/bin/env python3
"""Exact four-profile quadratic fork census; prepare now, root guards execution.

The complete actual-field and weak-prefix lemmas and a separately pinned proof
of the four anchor domains are mathematical inputs. This script reuses ONLY
the audited arithmetic and prefix-list functions of cac3bc1d; its old isolated
fork hypothesis, anchor generator, census and main are NOT called.

Four finite domains (all endpoints in the integer loops are inclusive):
 FAR3: T=2..10, N=-5..3, 2<t<10,
       H=(19*t^2+4*t-4)/4;                           81 rows.
 NEAR3: T=1..9, N=-13..-1, 2<t<10,
       H=(19*t^2+34*t+11)/4;                       117 rows.
 FAR_DPOS3: T=2..9, N=-4..-1, 2<t<10,
       H=t*(19*t+8)/4;                             32 rows.
 FAR4: T=3..14, N=-5..-1, t>0 and 8<t^2<200,
       H=t*(19*t+8)/4;                             60 rows.
Here t is the larger root of X^2-T*X+N. Keep positive NONSQUARE discriminants,
positive physical H and strictly negative conjugate H. This visits all290
integer rows, including N=0 in FAR3 (reducible, therefore filtered).

Normalizations: t=2gamma in the three native3-helper families, raw norm2;
t=2sqrt2*gamma in FAR4, raw norm4. H is the corresponding RAW squared time,
not its unit-wall value. The proof input must justify integrality, actual
K=Q(t), these trace/norm/physical ranges and the complete remaining profiles.

The257 prefix rows come from the unchanged complete d8d8 inventory. Canonical
squarefree d and its positive radical identify the physical quadratic fields;
test U/H using the full exact K-square test, including both rational norm
signs and the rational/pure-radical branches. Every positive result has an
exact multiplied square-root witness. Loop termination and square-test
correctness are the audited cac3bc1d algorithm, not a numerical approximation.

Output groups ALL necessary matches by anchor and weak index. No empty outcome
or geometric realization is assumed. This is no whole-D19/dimension theorem.
No large JSON input, scientific package, shell job, old main or data worker is
imported. Run ONLY via the existing guarded runner with --parent-lock. The
base source is verified BEFORE its normal Python import and AGAIN afterward;
all mathematical inputs and this worker are pinned before/after the census.
The required JSON destination must be new; stdout contains bounded counters.
"""

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE = 'scripts/check_d19_isolated_quadratic_fork_time_compatibility_20261021.py'
PROOF = 'results/arithmetic_d19_remaining_dotted_fork_physical_bounds_and_census_20261021.md'
PINS = {
    BASE: 'cac3bc1d4b9c0a1c98667917f1603f81f76f411e9c378b0bd3880b87145ebc40',
    'results/arithmetic_d19_actual_field_is_quadratic_20261021.md':
        '790b70eba68daf53c286ec5bf19b0b911b4089431638267dd9c8ebd8e4404646',
    'results/arithmetic_d19_quadratic_weak_prefix_inventory_20261021.md':
        'd8d80322d9d35d4df65ca72d897b6d60958648ab267e2562a9c803b085d5f8db',
    PROOF: '80d67e3b48b33534f7e716e2011323123d4b05a4f1aa6fe9fd770919bb9b930f',
}
FAMILIES = (
    ('FAR3', range(2, 11), range(-5, 4)),
    ('NEAR3', range(1, 10), range(-13, 0)),
    ('FAR_DPOS3', range(2, 10), range(-4, 0)),
    ('FAR4', range(3, 15), range(-5, 0)),
)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(65536), b''):
            h.update(block)
    return h.hexdigest()


def anchor_candidates(core):
    rows, domain_counts, candidate_counts = [], {}, {}
    for family, traces, norms in FAMILIES:
        domain_counts[family] = len(traces)*len(norms)
        candidate_counts[family] = 0
        for trace in traces:
            for norm in norms:
                t = core.positive_root(trace, norm)
                if t is None:
                    continue
                if family == 'FAR4':
                    if not (t.sign() > 0 and (t*t-8).sign() > 0
                            and (200-t*t).sign() > 0):
                        continue
                elif not ((t-2).sign() > 0 and (10-t).sign() > 0):
                    continue
                if family == 'FAR3':
                    time = (19*t*t+4*t-4)/4
                elif family == 'NEAR3':
                    time = (19*t*t+34*t+11)/4
                else:
                    time = t*(19*t+8)/4
                if not (time.sign() > 0 and time.conjugate().sign() < 0):
                    continue
                rows.append({'family': family, 'trace': trace, 'norm': norm,
                             't': t, 'time': time})
                candidate_counts[family] += 1
    assert sum(domain_counts.values()) == 290
    return rows, domain_counts, candidate_counts


def census(core):
    anchors, domains, candidates = anchor_candidates(core)
    prefixes = core.prefix_candidates()
    assert len(prefixes) == 257
    survivors, zero_match, result_counts = [], [], {}
    prefix_counts, survivor_counts, weak_union = {}, {}, set()
    same_field_pairs = 0
    for prefix in prefixes:
        family = prefix['family']
        prefix_counts[family] = prefix_counts.get(family, 0)+1
    for anchor in anchors:
        matches = []
        for prefix in prefixes:
            if prefix['time'].d != anchor['time'].d:
                continue
            same_field_pairs += 1
            ratio = prefix['time']/anchor['time']
            assert ratio.sign() > 0 and ratio.conjugate().sign() > 0
            witness, reason = core.square_in_field(ratio)
            result_counts[reason] = result_counts.get(reason, 0)+1
            if witness is not None:
                assert witness*witness == ratio
                matches.append({'prefix_family': prefix['family'],
                    'index': prefix['index'], 'prefix_parameter': prefix['parameter'],
                    'prefix_time': prefix['time'].encode(), 'ratio': ratio.encode(),
                    'verified_square_root': witness.encode(), 'square_branch': reason})
        base = {'anchor_family': anchor['family'], 'trace': anchor['trace'],
                'norm': anchor['norm'], 't': anchor['t'].encode(),
                'anchor_raw_time': anchor['time'].encode()}
        if matches:
            indices = sorted({match['index'] for match in matches})
            survivors.append(dict(base, possible_weak_indices=indices,
                                  compatible_weak_candidates=matches))
            weak_union.update(indices)
            family = anchor['family']
            survivor_counts[family] = survivor_counts.get(family, 0)+1
        else:
            zero_match.append(base)
    return {
        'status': 'COMPLETE_EXACT_FOUR_PROFILE_NECESSARY_COMPATIBILITY_CENSUS',
        'scope': 'Quadratic local-data compatibility only; no actual-wall realization or dimension conclusion.',
        'integer_domain_rows': sum(domains.values()),
        'integer_domain_rows_by_family': domains,
        'anchor_models': len(anchors), 'anchor_models_by_family': candidates,
        'prefix_candidate_rows': len(prefixes), 'prefix_rows_by_family': prefix_counts,
        'same_field_pairs': same_field_pairs, 'exact_field_square_tests': same_field_pairs,
        'field_square_result_counts': result_counts,
        'no_compatible_pair': not survivors,
        'remaining_anchor_models': len(survivors),
        'remaining_anchor_models_by_family': survivor_counts,
        'remaining_compatible_pairs': sum(len(row['compatible_weak_candidates']) for row in survivors),
        'union_of_possible_weak_indices': sorted(weak_union),
        'remaining_necessary_models_not_geometric_realizations': survivors,
        'zero_match_anchor_models': zero_match,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--proof', required=True,
                        help='Complete new four-profile domain proof directly under results/')
    parser.add_argument('--proof-sha', required=True)
    parser.add_argument('--out', required=True, help='New JSON file directly under results/')
    args = parser.parse_args()
    proof = (ROOT/args.proof).resolve()
    assert proof.parent == ROOT/'results' and proof.suffix == '.md'
    assert args.proof == str(proof.relative_to(ROOT))
    assert args.proof == PROOF and args.proof_sha == PINS[PROOF]
    destination = (ROOT/args.out).resolve()
    assert destination.parent == ROOT/'results' and destination.suffix == '.json'
    assert not destination.exists(), 'Output must be new'
    pins = dict(PINS)
    for name, pin in pins.items():
        assert digest(ROOT/name) == pin, 'Input changed: '+name
    original = digest(Path(__file__))

    import check_d19_isolated_quadratic_fork_time_compatibility_20261021 as core
    assert Path(core.__file__).resolve() == ROOT/BASE
    assert digest(ROOT/BASE) == PINS[BASE]
    output = census(core)
    output['checker_sha256'] = original
    output['mathematical_and_arithmetic_source_pins'] = pins
    assert digest(Path(__file__)) == original
    for name, pin in pins.items():
        assert digest(ROOT/name) == pin, 'Input changed: '+name
    with destination.open('x') as stream:
        json.dump(output, stream, indent=2, sort_keys=True)
        stream.write('\n')
    summary = {key: value for key, value in output.items()
               if key not in ('remaining_necessary_models_not_geometric_realizations',
                              'zero_match_anchor_models',
                              'mathematical_and_arithmetic_source_pins')}
    print(json.dumps(summary, sort_keys=True))


if __name__ == '__main__':
    main()
