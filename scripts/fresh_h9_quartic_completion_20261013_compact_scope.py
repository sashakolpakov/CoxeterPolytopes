#!/usr/bin/env python3
"""Exact inherited cusp obstructions and a 43-wall deletion checkpoint.

No numerical modules or search producers are imported.
"""
import argparse
import json
from itertools import combinations
from scripts.fresh_h9_quartic_completion_20261013_opposite import (
    build as opposite_build, OUTPUT as OPPOSITE)
from scripts.fresh_h9_quartic_completion_20261013_polar import (
    ROOT, SOURCE as SOURCE44, read, K, Quadratic, Z, O, parse, enc, scale,
    dot, lift_K, check_escape)
from scripts.fresh_h9_quartic_completion_20261013_geometry_verify import (
    check as geometry_check, inertia_K)

OUTPUT = ROOT/'data/enumerations/fresh_h9_quartic_completion_20261013_compact_scope.json'


def embed_tower(target, source):
    """Return a common tower and an exact embedding of the source."""
    if source is K:
        return target, lambda x: lift_K(target, x)
    middle, embed = embed_tower(target, source.base)
    common = Quadratic(middle, embed(source.d))
    return common, lambda x: (embed(x[0]), embed(x[1]))


def lift_existing(target, source, x):
    if target is source:
        return x
    return target.lift(lift_existing(target.base, source, x))


def run():
    geometry = geometry_check(SOURCE44)
    assert geometry['facet_counts'] == [44]
    record, context = opposite_build()
    saved = json.loads(OPPOSITE.read_text())
    assert all(saved[k] == v for k, v in record.items())
    cert = saved['volume']['certificate']
    check_escape(cert, context)
    G, I, center, old, case, F, pa, c, prior = context
    removed = [31, 51, 197]
    profiles = {}
    for i in range(10):
        profiles[f'seed{i}'] = (K, tuple(O if j == i else Z for j in range(10)),
                                tuple(G[j][i] for j in range(10)), scale(O, 2))
    for i in case['retained_old_caps']:
        n = old[i]
        f = Quadratic(K, n['d'])
        profiles[f'old{i}'] = (f, tuple(zip(n['a'], n['b'])),
                              tuple(zip(n['pa'], n['pb'])), n['q'])
    references = {}
    for j, (ref, f, cs, ps) in enumerate(prior):
        profiles[f'new{j}'] = f, cs, ps, scale(O, 2)
        references[f'new{j}'] = ref
    profiles['opposite'] = F, c, pa, scale(O, 2)
    assert len(profiles) == 46
    inherited_parallel = []
    for index in case['retained_old_caps']:
        f, cs, ps, q = profiles[f'old{index}']
        for i, p in enumerate(ps):
            if f.sign(f.sub(f.mul(p, p), lift_K(f, scale(q, 2)))) == 0:
                assert f.sign(p) < 0
                inherited_parallel.append([f'old{index}', f'seed{i}'])
    assert inherited_parallel == [
        ['old31', 'seed8'], ['old51', 'seed8'], ['old51', 'seed9'], ['old197', 'seed1']]
    retained = [name for name in profiles if name not in {f'old{i}' for i in removed}]
    assert len(retained) == 43
    # Every other pair has already been exactly checked by geometry_check,
    # opposite_build and polar_build; their gates reject new parallel pairs.
    assert all(a not in retained or b not in retained for a, b in inherited_parallel)
    y = tuple(map(parse, cert['y']))
    assert all(K.sign(x) > 0 for x in y)
    for name in retained:
        f, cs, ps, q = profiles[name]
        assert f.sign(dot(f, cs, tuple(lift_K(f, x) for x in y))) > 0
    # Any PSD seed principal matrix must omit8 or9, since that two-wall
    # determinant is negative. Each such subset lies in one PD maximal.
    seed_checks = []
    for omitted in (8, 9):
        ids = [i for i in range(10) if i != omitted]
        gram = [[G[i][j] for j in ids] for i in ids]
        assert inertia_K(gram) == (9, 0, 0)
        seed_checks.append({'omitted': omitted, 'indices': ids, 'inertia': [9, 0, 0]})
    bad2 = K.sub(K.mul(G[8][8], G[9][9]), K.mul(G[8][9], G[8][9]))
    assert K.sign(bad2) < 0
    circuits = []
    proposed = [
        ['seed0', 'seed1', 'seed2', 'new4'],
        ['seed2', 'seed3', 'seed4', 'old1028'],
        ['seed5', 'seed6', 'old310', 'old837'],
        ['seed4', 'seed5', 'old496'],
        ['seed3', 'seed4', 'old496']]
    for names in proposed:
        assert all(name in retained for name in names)
        gram = [[2 if i == j else 0 for j in range(len(names))]
                for i in range(len(names))]
        for name in names:
            f, cs, ps, q = profiles[name]
            assert q == scale(O, 2)
            assert f.sign(dot(f, cs, tuple(lift_K(f, x) for x in center))) > 0
        for i, j in combinations(range(len(names)), 2):
            fi, ci, pi, qi = profiles[names[i]]
            fj, cj, pj, qj = profiles[names[j]]
            common, embed_j = embed_tower(fi, fj)
            value = dot(common, tuple(lift_existing(common, fi, x) for x in pi),
                        tuple(embed_j(x) for x in cj))
            if common.sign(value) == 0:
                entry = 0
            else:
                assert common.sign(common.sub(value, common.rational(-1))) == 0
                entry = -1
            gram[i][j] = gram[j][i] = entry
        assert all(sum(row) == 0 for row in gram)
        assert inertia_K([[scale(O, x) for x in row] for row in gram]) == (len(names)-1, 0, 1)
        circuits.append({'walls': names, 'gram': gram,
                         'positive_kernel': [1]*len(names),
                         'type': 'affine_A'+str(len(names)-1)})
    nonseed = sorted({name for circuit in circuits for name in circuit['walls']
                      if not name.startswith('seed')})
    hitting = []
    for size in range(len(nonseed)+1):
        hitting = [list(ids) for ids in combinations(nonseed, size)
                   if all(set(ids).intersection(circuit['walls']) for circuit in circuits)]
        if hitting:
            break
    assert size == 4 and len(hitting) == 2
    return {
        'status': 'PASS_EXACT_COMPACTNESS_SCOPE',
        'source_46': str(OPPOSITE.relative_to(ROOT)),
        'source_44': SOURCE44,
        'inherited_parallel_pairs': inherited_parallel,
        'removed_old_caps': removed, 'retained_wall_names': retained,
        'references_for_new_names': references, 'facet_count': 43,
        'transferred_strict_escape': cert,
        'remaining_parallel_pairs': [],
        'seed_no_affine_proof': {'positive_maximals': seed_checks,
                                'pair_8_9_determinant': enc(bad2)},
        'remaining_positive_affine_circuits': circuits,
        'minimum_nonseed_deletions_for_displayed_circuits': size,
        'minimum_nonseed_hitting_sets_for_displayed_circuits': hitting,
        'scope': 'The specified43-wall deletion chamber is an acute Coxeter chamber '
                 'of infinite volume, with no parallel wall pairs. The five displayed '
                 'positive affine circuits still forbid compactification by adding '
                 'acute walls while keeping this core. This is not an exhaustive '
                 'affine-circuit census and does not obstruct finite-volume completion. '
                 'The ten-wall seed alone contains no affine principal subdiagram.'}


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--verify', action='store_true')
    args = p.parse_args()
    result = run()
    if args.verify:
        assert json.loads(OUTPUT.read_text()) == result
    else:
        OUTPUT.write_text(json.dumps(result, indent=2)+'\n')
    print({k: result[k] for k in ('status', 'facet_count', 'remaining_parallel_pairs')})
    print({'positive_affine_circuits': len(result['remaining_positive_affine_circuits'])})
