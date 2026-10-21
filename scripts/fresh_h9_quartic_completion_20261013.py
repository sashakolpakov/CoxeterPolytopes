#!/usr/bin/env python3
"""Full quartic anchor equations, multi-gap cuts, and controlled replacements."""
import argparse
from fractions import Fraction as Q
import heapq
from itertools import combinations, product
import json
import numpy as np
import sympy as sp

from scripts.fresh_h9_joint_completion_20261013 import setup, ROOT, read
from scripts.fresh_h8_alternate_caps_20261013 import real
from scripts.fresh_h8_alternate_core_angles_independent_20261013 import totient
from scripts.fresh_h9_quartic_completion_20261013_arithmetic import (
    K, Quadratic, Z, O, P, add, scale, parse, enc, dot, pairing, matvec,
    old_pair_gate, new_pair_gate, finite_label)

OUTPUT = ROOT / 'data/enumerations/fresh_h9_quartic_completion_20261013.json'
ANGLE_OUTPUT = ROOT / 'data/enumerations/fresh_h9_quartic_completion_20261013_angles.json'


def angle_polynomials():
    # phi(m)^2/m >=1/2 implies m<=2048 when phi(m)<=32.
    labels = [m for m in range(2, 2049)
              if m == 2 or totient(m) in (2, 4, 8, 16, 32)]
    z, v = sp.symbols('z v')
    out = []
    for m in labels:
        polynomial = (sp.Poly(v, v) if m == 2 else
                      sp.Poly(sp.resultant(sp.cyclotomic_poly(m, z),
                                          z*z-(v-2)*z+1, z), v).sqf_part().monic())
        intervals = polynomial.intervals(eps=sp.Rational(1, 10**25))
        lower, upper = intervals[-1][0]
        assert polynomial.degree() == (1 if m == 2 else totient(m)//2)
        assert sum(multiplicity for _, multiplicity in intervals) == polynomial.degree()
        out.append((m, tuple(Q(str(c)) for c in polynomial.all_coeffs()),
                    Q(str(lower)), Q(str(upper))))
    ANGLE_OUTPUT.write_text(json.dumps([
        {'label': m, 'coefficients_descending': list(map(str, p)),
         'largest_root_interval': [str(lo), str(hi)]}
        for m, p, lo, hi in out], indent=2) + '\n')
    return out


def run(min_retained=16, keep=2, limit=160, seed_finite_free=False,
        targeted_anchor_distances=False, min_cuts=2, focus=(5, 22),
        target_source='data/enumerations/fresh_h9_joint_completion_20261013_targets.json',
        required_source=None, required_indices=(), include_seed_anchors=False):
    G, I, center, old, _, cn, In, _, original = setup()
    anchor = 1039
    chosen = sorted([i for i in original if i != 1038] + [anchor])
    rational_anchors = [i for i in chosen if old[i]['d'] == Z]
    if include_seed_anchors:
        rational_anchors += list(range(-10, 0))
    targets = read(target_source)
    ys = [tuple(map(parse, r['certificate']['y'])) for r in targets['targets']]
    Y = np.array([[real(x) for x in y] for y in ys])
    IY, IC = In @ Y.T, In @ np.array([real(x) for x in center])
    mags = (Z, O, P, add(O, P), scale(P, 2),
            add(scale(O, 2), P), add(O, scale(P, 2)))
    anchor_targets = []
    if targeted_anchor_distances:
        # A saved spacelike witness itself gives a norm-two acute normal.
        # Its pairing with the retained quartic anchor estimates a useful
        # *dotted* relation; do not constrain that relation to small angles.
        distances = []
        for target in focus:
            y = Y[target]
            magnitude = float(y @ cn[anchor] * np.sqrt(2/(y @ In @ y)))
            integers = sorted({int(np.floor(magnitude)), int(np.ceil(magnitude))})
            integers = [n for n in integers if n > 2]
            if not integers:
                continue
            distances.extend(integers)
            anchor_targets.append({'target': target, 'numerical_pairing_magnitude': magnitude,
                                   'rational_proposals': integers})
        mags = (Z, O, P) + tuple(scale(O, j) for j in sorted(set(distances)))
    cases = [(j, m, n) for j in rational_anchors
             for m, n in product(range(len(mags)), range(3))]
    aa = np.broadcast_to(cn[anchor], (len(cases), 10))
    bb = np.array([cn[j] if j >= 0 else np.eye(10)[-1-j] for j, _, _ in cases])
    rhs = np.array([[real(mags[m]), real(mags[n])] for _, m, n in cases])
    oldcn, oldnorms = cn[chosen], np.array([real(old[i]['q']) for i in chosen])
    required = []
    if required_source:
        from scripts.fresh_h9_quartic_completion_20261013_geometry import numeric
        import mpmath as mp
        mp.mp.dps = 90
        extra = read(required_source)
        assert parse(extra['anchor_radicand']) == old[anchor]['d']
        E0 = Quadratic(K, old[anchor]['d'])
        for i in required_indices:
            c = extra['caps'][i]
            f = Quadratic(E0, E0.parse(c['discriminant_over_E']))
            cc = tuple(map(f.parse, c['seed_coefficients']))
            required.append((i, f.d, cc, f))
        required_cn = np.array([[float(numeric(f, x)) for x in cc] for _, _, cc, f in required])
    labels = [m for m in range(2, 2049)
              if m == 2 or totient(m) in (2, 4, 8, 16, 32)]
    allowed = np.array([np.cos(np.pi/m)**2 for m in labels])
    patterns = [(supp, labs) for k in range(4)
                for supp in combinations(range(7), k)
                for labs in product((0, 1), repeat=k)]
    stats = {'free_triples': 0, 'center_seed_candidates': 0,
             'targeted_candidates': 0, 'retention_candidates': 0,
             'best_retained': 0, 'best_targets_cut': 0}
    heaps, serial = {}, 0
    for free in combinations(range(10), 3):
        stats['free_triples'] += 1
        remaining = [i for i in range(10) if i not in free]
        assigned = np.zeros(len(cases), dtype=bool)
        for u, v, w in (free, (free[1], free[0], free[2]),
                        (free[2], free[0], free[1])):
            determinant = aa[:, v]*bb[:, w]-aa[:, w]*bb[:, v]
            active = (abs(determinant) > 1e-9) & ~assigned
            assigned |= active
            if not active.any():
                continue
            ids = np.flatnonzero(active)
            a, b, det, m = aa[ids], bb[ids], determinant[ids], rhs[ids]
            w1 = np.zeros((len(ids), 10))
            w1[:, u] = 1
            w1[:, v] = (-a[:, u]*b[:, w]+a[:, w]*b[:, u])/det
            w1[:, w] = (-a[:, v]*b[:, u]+a[:, u]*b[:, v])/det
            A = np.einsum('bi,ij,bj->b', w1, In, w1)
            for positions, fixed_labels in patterns:
                fixed = np.zeros(10)
                for i, label in zip(positions, fixed_labels):
                    fixed[remaining[i]] = (1., real(P))[label]
                r0, r1 = m[:, 0]-a @ fixed, m[:, 1]-b @ fixed
                w0 = np.broadcast_to(fixed, (len(ids), 10)).copy()
                w0[:, v] = (r0*b[:, w]-a[:, w]*r1)/det
                w0[:, w] = (a[:, v]*r1-r0*b[:, v])/det
                B = 2*np.einsum('bi,ij,bj->b', w1, In, w0)
                C = np.einsum('bi,ij,bj->b', w0, In, w0)-2
                D = B*B-4*A*C
                good = (D > 1e-9) & (abs(A) > 1e-9)
                for branch in (-1, 1):
                    x = np.full(len(ids), np.nan)
                    x[good] = (-B[good]+branch*np.sqrt(D[good]))/(2*A[good])
                    ws = w0+x[:, None]*w1
                    cuts = ws @ IY
                    counts = (cuts > 1e-7).sum(axis=1)
                    freevalues = ws[:, list(free)]
                    if seed_finite_free:
                        nearseed = np.min(abs(freevalues[:, :, None]**2/4-allowed), axis=2) < 1e-6
                        seedgood = np.all((freevalues >= -1e-7) &
                                          ((freevalues > 2+1e-7) | nearseed), axis=1)
                    else:
                        seedgood = np.all(freevalues > 2+1e-7, axis=1)
                    ok = good & seedgood
                    ok &= ws @ IC < -1e-7
                    stats['center_seed_candidates'] += int(ok.sum())
                    focused = np.any(cuts[:, list(focus)] > 1e-7, axis=1)
                    ok &= focused & (counts >= min_cuts)
                    if required:
                        ep = -ws @ required_cn.T
                        er = ep*ep/4
                        enear = np.min(abs(er[:, :, None]-allowed), axis=2) < 1e-6
                        ok &= np.all((ep <= 1e-7) & ((er > 1+1e-7) | enear), axis=1)
                    stats['targeted_candidates'] += int(ok.sum())
                    indices = np.flatnonzero(ok)
                    if not len(indices):
                        continue
                    pairings = -ws[indices] @ oldcn.T
                    ratios = pairings*pairings/(2*oldnorms)
                    near = np.min(abs(ratios[:, :, None]-allowed), axis=2) < 1e-6
                    compatible = (pairings <= 1e-7) & ((ratios > 1+1e-7) | near)
                    retained = compatible.sum(axis=1)
                    stats['best_retained'] = max(stats['best_retained'], int(retained.max()))
                    stats['best_targets_cut'] = max(stats['best_targets_cut'], int(counts[indices].max()))
                    for jj in np.flatnonzero(retained >= min_retained):
                        j = int(indices[jj])
                        if not compatible[jj, chosen.index(anchor)]:
                            continue
                        stats['retention_candidates'] += 1
                        other, mi, ni = cases[int(ids[j])]
                        mask = tuple(chosen[k] for k in range(len(chosen)) if compatible[jj, k])
                        focus_mask = tuple(t for t in focus if cuts[j, t] > 1e-7)
                        score = (len(focus_mask), int(counts[j]), int(retained[jj]),
                                 float(max(cuts[j, t] for t in focus)/np.linalg.norm(ws[j] @ In)))
                        spec = (u, v, w, tuple(remaining[i] for i in positions),
                                fixed_labels, other, mi, ni, branch)
                        serial += 1
                        item = score, serial, spec
                        heap = heaps.setdefault((other, mask, focus_mask), [])
                        if len(heap) < keep:
                            heapq.heappush(heap, item)
                        elif score > heap[0][0]:
                            heapq.heapreplace(heap, item)
        if stats['free_triples'] % 20 == 0:
            print(stats, flush=True)
    proposals = sorted([item for heap in heaps.values() for item in heap], reverse=True)
    print({'retained_proposals': len(proposals), 'exact_limit': limit}, flush=True)
    polynomials = angle_polynomials()
    E = Quadratic(K, old[anchor]['d'])
    IE = tuple(tuple(E.lift(x) for x in row) for row in I)
    ar = tuple(zip(old[anchor]['a'], old[anchor]['b']))
    caps, rejected, seen = [], [], set()

    def save():
        result = {'schema': 'full-quartic-anchor-H9-v1',
                  'base_source': 'data/enumerations/fresh_h8_alternate_3_6_20261013.json',
                  'target_source': target_source, 'focus_target_indices': list(focus),
                  'required_cap_source': required_source,
                  'required_cap_indices': list(required_indices),
                  'include_seed_anchors': include_seed_anchors,
                  'angle_source': str(ANGLE_OUTPUT.relative_to(ROOT)),
                  'anchor': anchor, 'anchor_radicand': enc(E.d),
                  'selected_configuration': chosen, 'minimum_retained': min_retained,
                  'anchor_pairing_magnitudes': list(map(enc, mags)),
                  'target_guided_anchor_distances': anchor_targets,
                  'minimum_targets_cut': min_cuts,
                  'solved_seed_pairings_may_be_finite': seed_finite_free,
                  'proposal_statistics': stats, 'scored_proposals': len(proposals),
                  'exact_limit': limit, 'caps': caps, 'rejected_proposals': rejected,
                  'scope': 'Full quartic anchor equation, without coefficient splitting. '
                  'Numerical proposals only; exact accepted roots, seed relations, '
                  'retained old-cap relations, center orientation and target cuts. '
                  'No exhaustive-family or volume conclusion.'}
        OUTPUT.write_text(json.dumps(result, indent=2)+'\n')

    for pi, (score, _, spec) in enumerate(proposals[:limit]):
        u, v, w, support, fixed_labels, other, mi, ni, branch = spec
        br = (tuple(E.lift(x) for x in old[other]['a']) if other >= 0 else
              tuple(E.one if i == -1-other else E.zero for i in range(10)))
        fixed = [E.zero]*10
        for i, label in zip(support, fixed_labels):
            fixed[i] = E.lift((O, P)[label])
        M, N = E.lift(mags[mi]), E.lift(mags[ni])
        determinant = E.sub(E.mul(ar[v], br[w]), E.mul(ar[w], br[v]))
        r0, r1 = E.sub(M, dot(E, ar, fixed)), E.sub(N, dot(E, br, fixed))
        w0, w1 = list(fixed), [E.zero]*10
        w0[v] = E.div(E.sub(E.mul(r0, br[w]), E.mul(ar[w], r1)), determinant)
        w0[w] = E.div(E.sub(E.mul(ar[v], r1), E.mul(r0, br[v])), determinant)
        w1[u] = E.one
        w1[v] = E.div(E.sub(E.mul(ar[w], br[u]), E.mul(ar[u], br[w])), determinant)
        w1[w] = E.div(E.sub(E.mul(ar[u], br[v]), E.mul(ar[v], br[u])), determinant)
        assert E.sign(E.sub(dot(E, ar, w0), M)) == E.sign(dot(E, ar, w1)) == 0
        assert E.sign(E.sub(dot(E, br, w0), N)) == E.sign(dot(E, br, w1)) == 0
        A = pairing(E, w1, IE, w1)
        B = E.mul(E.rational(2), pairing(E, w0, IE, w1))
        C = E.sub(pairing(E, w0, IE, w0), E.rational(2))
        D = E.sub(E.mul(B, B), E.mul(E.rational(4), E.mul(A, C)))
        assert E.sign(D) > 0
        field = Quadratic(E, D)
        x = (E.div(E.neg(B), E.mul(E.rational(2), A)),
             E.div(E.rational(branch), E.mul(E.rational(2), A)))
        pa = tuple(field.neg(field.add(field.lift(a), field.mul(x, field.lift(b))))
                   for a, b in zip(w0, w1))
        IF = tuple(tuple(field.lift(x) for x in row) for row in IE)
        ca = matvec(field, IF, pa)
        key = D, ca
        if key in seen:
            continue
        seen.add(key)
        assert field.sign(field.sub(dot(field, pa, ca), field.rational(2))) == 0
        seed_relations = []
        seed_passed = True
        for j in range(10):
            assert field.sign(pa[j]) <= 0
            square = field.mul(pa[j], pa[j])
            if j in (u, v, w):
                comparison = field.sign(field.sub(square, field.rational(4)))
                label = None if comparison >= 0 else finite_label(field, square, polynomials)
                kind = 'dotted' if comparison > 0 else 'parallel' if comparison == 0 else 'finite' if label is not None else 'noncoxeter'
                seed_relations.append({'seed': j, 'kind': kind, 'label': label})
                if kind not in ('dotted', 'finite'):
                    seed_passed = False
            else:
                assert any(field.sign(field.sub(square, field.lift(E.lift(t)))) == 0
                           for t in (Z, O, add(O, P)))
                label = finite_label(field, square, polynomials)
                assert label in (2, 3, 5)
                seed_relations.append({'seed': j, 'kind': 'finite', 'label': label})
        if not seed_passed:
            rejected.append({'proposal_index': pi, 'seed_relations': seed_relations})
            continue
        assert field.sign(dot(field, ca, tuple(field.lift(E.lift(t)) for t in center))) > 0
        assert field.sign(field.add(dot(field, pa, tuple(field.lift(t) for t in ar)),
                                    field.lift(M))) == 0
        checks = [{'old_node': i, **old_pair_gate(field, pa, old[i], polynomials)}
                  for i in chosen]
        required_checks = [{'cap_index': i, **new_pair_gate(field, pa, dd, cc, polynomials)}
                           for i, dd, cc, _ in required]
        kept = [r['old_node'] for r in checks if r['compatible']]
        dropped = [r['old_node'] for r in checks if not r['compatible']]
        cuts = [i for i, y in enumerate(ys)
                if field.sign(dot(field, ca, tuple(field.lift(E.lift(t)) for t in y))) < 0]
        if (len(kept) < min_retained or anchor not in kept or (other >= 0 and other not in kept) or
                any(not r['compatible'] for r in required_checks)):
            rejected.append({'proposal_index': pi, 'retained': kept, 'dropped': dropped})
            continue
        assert set(cuts) & set(focus) and len(cuts) >= min_cuts
        caps.append({'free_seed_indices': [u, v, w], 'rational_anchor': other,
                     'anchor_pairing_magnitudes': [enc(mags[mi]), enc(mags[ni])],
                     'line_constant': list(map(E.encode, w0)),
                     'line_direction': list(map(E.encode, w1)),
                     'norm_polynomial_descending': list(map(E.encode, (A, B, C))),
                     'discriminant_over_E': E.encode(D), 'parameter': field.encode(x),
                     'seed_pairings': list(map(field.encode, pa)),
                     'seed_coefficients': list(map(field.encode, ca)),
                     'seed_relations': seed_relations,
                     'retained_nodes': kept, 'dropped_nodes': dropped,
                     'pair_checks': checks, 'cut_target_indices': cuts,
                     'required_cap_pair_checks': required_checks,
                     'numeric_score': list(score)})
        save()
        print({'proposal': pi, 'accepted': len(caps), 'retained': len(kept),
               'cuts': cuts}, flush=True)
    save()
    print({'accepted': len(caps), 'rejected': len(rejected)}, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--min-retained', type=int, default=16)
    parser.add_argument('--keep', type=int, default=2)
    parser.add_argument('--limit', type=int, default=160)
    parser.add_argument('--output')
    parser.add_argument('--seed-finite-free', action='store_true')
    parser.add_argument('--targeted-anchor-distances', action='store_true')
    parser.add_argument('--min-cuts', type=int, default=2)
    parser.add_argument('--focus', default='5,22')
    parser.add_argument('--target-source', default='data/enumerations/fresh_h9_joint_completion_20261013_targets.json')
    parser.add_argument('--required-source')
    parser.add_argument('--required-caps', default='')
    parser.add_argument('--include-seed-anchors', action='store_true')
    args = parser.parse_args()
    if args.output:
        OUTPUT = ROOT / args.output
    run(args.min_retained, args.keep, args.limit, args.seed_finite_free,
        args.targeted_anchor_distances, args.min_cuts,
        tuple(map(int, args.focus.split(','))), args.target_source,
        args.required_source,
        tuple(map(int, args.required_caps.split(','))) if args.required_caps else (),
        args.include_seed_anchors)
