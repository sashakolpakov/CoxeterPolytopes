#!/usr/bin/env python3
"""Exact streamed scout for fixed-D16 cusp profiles; no geometry inferred.

Domain: ranks sum to56; retain one varying rank16; inert ranks sum to40,
at most one inert rank below8, and at least one inert rank>=11.  The
integer partition generator has finite decreasing remainder and sorted
tails, hence no omitted/duplicate profiles in precisely this domain.
The mathematical extension from this domain requires the inherited
low8 concentration theorem with the retained factors held fixed.

No repository modules or large result files are loaded.  Run only under
run_bounded_research_check_20261018.py --parent-lock.  This worker writes
nothing; the wrapper preserves its bounded stdout in a fresh ledger.
"""

from fractions import Fraction as F
import json


def parts(total, low=8):
    if total == 0:
        yield ()
    for first in range(low, total + 1):
        for tail in parts(total - first, first):
            yield (first,) + tail


def envelope(rank):
    if rank <= 8:
        return F(rank, 2)
    table = {
        9: F(71, 16), 10: F(53, 11), 11: F(31, 6),
        12: F(11, 2), 13: F(40, 7), 14: F(367, 60),
        15: F(411, 64), 16: F(227, 34), 17: F(119, 18),
        18: F(128, 19), 19: F(137, 20), 20: F(146, 21),
        21: F(155, 22), 22: F(163, 23), 23: F(85, 12),
    }
    return table[rank] if rank < 24 else F(5 * rank + 56, rank + 1)


def evaluate(profile):
    h = [1]
    vertices = 1
    for rank in profile:
        out = [0] * (len(h) + rank)
        for j, value in enumerate(h):
            for shift in range(rank + 1):
                out[j + shift] += value
        h = out
        vertices *= rank + 1
    assert len(h) == 57 and sum(h) == vertices
    g = sum(value * min(j, 56 - j) for j, value in enumerate(h))
    curvature = F(5 * g - 56 * vertices, 2)
    k = profile.count(16)  # Conservative: every rank16 may be D16.
    margin = curvature / vertices - sum(map(envelope, profile), F()) - F(3 * k, 136)
    return margin, margin * vertices


def main():
    count = 0
    min_margin = min_absolute = None
    margin_profile = absolute_profile = None
    by_k = {}
    for small in range(8):
        for tail in parts(40 - small):
            if not any(rank >= 11 for rank in tail):
                continue
            profile = tuple(sorted((16,) + ((small,) if small else ()) + tail))
            margin, absolute = evaluate(profile)
            count += 1
            if min_margin is None or margin < min_margin:
                min_margin, margin_profile = margin, profile
            if min_absolute is None or absolute < min_absolute:
                min_absolute, absolute_profile = absolute, profile
            k = profile.count(16)
            if k not in by_k or margin < by_k[k][0]:
                by_k[k] = margin, profile
    assert count and min_margin is not None
    print(json.dumps({
        'status': 'EXACT_SCOUT_OUTPUT_NOT_A_GEOMETRIC_THEOREM',
        'domain': 'fixed16; inert total40; at most one inert rank<8; some inert rank>=11',
        'profiles': count,
        'minimum_normalized_margin': str(min_margin),
        'normalized_minimizer': margin_profile,
        'minimum_absolute_surplus': str(min_absolute),
        'absolute_minimizer': absolute_profile,
        'positive_on_this_domain': min_margin > 0,
        'absolute_surplus_above_43_over_128': min_absolute > F(43, 128),
        'by_rank16_count': {k: {'margin': str(v[0]), 'profile': v[1]} for k, v in by_k.items()},
        'ideal_rule': 'manuscript envelope; effective extra charge3/136 per rank16',
        'fresh_geometric_completeness_proved': False,
    }, sort_keys=True))


if __name__ == '__main__':
    main()
