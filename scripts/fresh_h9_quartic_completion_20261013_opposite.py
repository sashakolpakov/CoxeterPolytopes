#!/usr/bin/env python3
"""Exact opposite-end cap, followed by a 46-wall volume test.

--verify imports no numerical libraries or search producers.
"""
import argparse
from fractions import Fraction as Q
import json
from scripts.fresh_h9_quartic_completion_20261013_polar import (
    ROOT, read, K, Quadratic, Z, O, P, parse, enc, add, sub, mul, div, scale,
    dot, matvec, pairing, lift_K, attach, build as polar_build,
    mirror_check, check_escape, propose_escape)

SOURCE = 'data/enumerations/fresh_h9_quartic_completion_20261013_polar.json'
OUTPUT = ROOT/'data/enumerations/fresh_h9_quartic_completion_20261013_opposite.json'
NEXT = ROOT/'data/enumerations/fresh_h9_quartic_completion_20261013_next_targets.json'


def build():
    polar, ctx = polar_build()
    saved = read(SOURCE)
    assert all(saved[k] == v for k, v in polar.items())
    G, I, center, old, case, Fpolar, papolar, cpolar, prior = ctx
    assert I[6][9] == I[9][6] == Z
    d9 = div(sub(scale(O, 2), scale(I[6][6], 2)), I[9][9])
    assert K.sign(d9) > 0
    E = Quadratic(K, d9)
    F = Quadratic(E, E.rational(2))
    w = [F.zero]*10
    w[6] = E.zero, E.one
    w[9] = (Z, O), E.zero
    pa = tuple(F.neg(p) for p in w)
    c = matvec(F, tuple(tuple(lift_K(F, t) for t in row) for row in I), pa)
    norm = dot(F, c, pa)
    assert F.sign(F.sub(norm, F.rational(2))) == 0
    cp = dot(F, c, tuple(lift_K(F, x) for x in center))
    assert F.sign(cp) > 0
    assert all(F.sign(p) <= 0 for p in pa)
    assert F.sign(F.sub(F.mul(pa[6], pa[6]), F.rational(2))) == 0
    assert F.sign(F.sub(F.mul(pa[9], pa[9]), F.rational(4))) > 0
    oldchecks = []
    for i in case['retained_old_caps']:
        n = old[i]
        T = Quadratic(F, lift_K(F, n['d']))
        p = (dot(F, pa, tuple(lift_K(F, x) for x in n['a'])),
             dot(F, pa, tuple(lift_K(F, x) for x in n['b'])))
        assert T.sign(p) < 0
        assert T.sign(T.sub(T.mul(p, p), lift_K(T, scale(n['q'], 2)))) > 0
        oldchecks.append({'index': i, 'kind': 'dotted', 'pairing': T.encode(p)})
    extra_prior = prior+[
        ({'source': SOURCE, 'kind': 'polar_pi3'}, Fpolar, cpolar, papolar)]
    newchecks = []
    for ref, f, cs, ps in extra_prior:
        T, csT = attach(F, f.base.d, f.d, cs)
        p = dot(T, tuple(T.lift(T.base.lift(x)) for x in pa), csT)
        assert T.sign(p) < 0
        assert T.sign(T.sub(T.mul(p, p), T.rational(4))) > 0
        newchecks.append({'reference': ref, 'kind': 'dotted', 'pairing': T.encode(p)})
    y = [Z]*10
    y[6], y[9] = sub(scale(O, 2), P), O
    ynorm, yc = pairing(K, y, I, y), pairing(K, y, I, center)
    assert ynorm == sub(O, scale(P, Q(1, 2))) and K.sign(ynorm) > 0
    assert yc == scale(O, Q(-7, 2))
    oldsigns = [Quadratic(K, old[i]['d']).sign(
        (dot(K, old[i]['a'], y), dot(K, old[i]['b'], y)))
        for i in case['retained_old_caps']]
    assert all(s >= 0 for s in oldsigns)
    assert [i for i, s in zip(case['retained_old_caps'], oldsigns, strict=True) if not s] == [67, 310]
    newsigns = [f.sign(dot(f, cs, tuple(lift_K(f, t) for t in y)))
                for ref, f, cs, ps in extra_prior]
    assert all(s > 0 for s in newsigns)
    target_pair = dot(F, c, tuple(lift_K(F, t) for t in y))
    assert F.sign(target_pair) < 0
    target_data = read(polar['target_source'])
    targets = [tuple(map(parse, r['certificate']['y'])) for r in target_data['targets']]
    target_signs = [F.sign(dot(F, c, tuple(lift_K(F, t) for t in v))) for v in targets]
    assert target_signs[1] < 0 and target_signs[30] < 0
    context = G, I, center, old, case, F, pa, c, extra_prior
    mirrors = mirror_check(context)
    assert [i for i, p in enumerate(papolar) if Fpolar.sign(p)] == [4, 8]
    mirrors['preceding_entries'] += 1
    record = {
        'status': 'EXACT_OPPOSITE_END_CAP', 'polar_source': SOURCE,
        'D9': enc(d9), 'field_degree_upper_bound': 8,
        'seed_pairings': list(map(F.encode, pa)),
        'seed_coefficients': list(map(F.encode, c)),
        'norm': F.encode(norm), 'center_inequality': F.encode(cp),
        'seed_angles': [2, 2, 2, 2, 2, 2, 4, 2, 2, 'dotted'],
        'old_checks': oldchecks, 'new_checks': newchecks,
        'exact_opposite_ray': {'y': list(map(enc, y)), 'norm': enc(ynorm),
            'center_pairing': enc(yc), 'old_zero_caps': [67, 310],
            'old_signs': oldsigns, 'new_signs': newsigns,
            'new_cap_inequality': F.encode(target_pair)},
        'target_source': polar['target_source'], 'target_signs': target_signs,
        'target_cuts': [i for i, s in enumerate(target_signs) if s < 0],
        'mirror_check': mirrors, 'facet_count': 46,
        'scope': 'Exact seed6 pi/4 cap, dotted to every retained nonseed facet; '
                 'cuts the explicit opposite gap while preserving all 45 previous facets.'}
    return record, context


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--verify', action='store_true')
    args = p.parse_args()
    record, context = build()
    if args.verify:
        saved = json.loads(OUTPUT.read_text())
        volume = saved.pop('volume')
        assert saved == record
        assert volume['status'] == 'EXACT_STRICT_SPACELIKE_ESCAPE'
        check_escape(volume['certificate'], context)
        record['volume'] = volume
    else:
        OUTPUT.write_text(json.dumps(record, indent=2)+'\n')
        print({'status': record['status'], 'target_cuts': record['target_cuts'],
               'mirror_check': record['mirror_check']}, flush=True)
        record['volume'] = propose_escape(context)
        OUTPUT.write_text(json.dumps(record, indent=2)+'\n')
    print({'status': 'PASS_EXACT_OPPOSITE_END_COUPLING', 'facets': record['facet_count'],
           'volume_status': record['volume']['status']}, flush=True)
    if record['volume']['status'] == 'EXACT_STRICT_SPACELIKE_ESCAPE':
        cert = record['volume']['certificate']
        assert cert['projectively_equal_saved_target_indices'] == []
        earlier = read(record['target_source'])['targets']
        assert len(earlier) == 36
        next_targets = {
            'source_target_archive': record['target_source'],
            'source_configuration': str(OUTPUT.relative_to(ROOT)),
            'appended_target_index': 36,
            'suggested_dominant_dual_coordinates': [3, 8, 9],
            'targets': earlier+[{'certificate': cert, 'name': '46_facet_escape'}],
            'scope': '37 saved target entries, with repetitions inherited from the first36. '
                     'The new index36 is exactly projectively distinct from every earlier '
                     'entry. Coordinates3,8,9 dominate numerically; the strict witness '
                     'itself has all10 coordinates positive. An independent cusp obstruction '
                     'rules out finite-volume acute wall-addition completion of the retained '
                     '46-wall core; this target is preserved for continuations that replace '
                     'obstructing retained walls.'}
        if args.verify:
            assert json.loads(NEXT.read_text()) == next_targets
        else:
            NEXT.write_text(json.dumps(next_targets, indent=2)+'\n')


if __name__ == '__main__':
    main()
